import collections
import zipfile
import pathlib
import os.path
import tempfile
import shutil

import openpathsampling as paths
# from openpathsampling.experimental.storage import Storage

# situations we need to track:
# 1. We need some way to disambiguate context for a mover; a given mover
#    might show up in more than one context.
# 2. We need to keep the context to a minimal lengthls
# this should handle checkpoint context writing
class EmptyCheckpointer:
    def __enter__(self):
        # this way, submovers get None
        return self

    def __exit__(self, exc, val, tb):
        pass

    def load_checkpoint(self):
        return None, None

    def save_checkpoint(self, data, additional_files=None):
        pass

    def next_context(self, context="", sequence=0, *, parent=True):
        return EmptyCheckpointer()

    def delete_checkpoint(self):
        pass

    @property
    def engine_extras(self):
        return []

def default_context_label(checkpointer, context_label, sequence):
    return str(len(checkpointer.children))

def uuid_labeler(checkpointer, context_label, sequence):
    if hasattr(context, '__uuid__'):
        context = self.encoding.encode(context.__uuid__)
    elif isinstance(context, int):
        context = str(context)
    elif isinstance(context, str):
        pass
    else:
        context = default_context_label(checkpointer, context_label)
    return context


class Checkpointer:
    """

    Parameters
    ----------
    storage_handler: :class:`.StorageHandler`
        Object to handle the actual backend where checkpoints are stored
        (filesystem, database, or other).
    context: os.Pathlike
        the path to the checkpoint
    frequency: Union[Dict[Engine, int], int]
        Dictionary mapping engine instance to the checkpoint frequency (as
        number of frames between checkpoints). An int can also be given, in
        which case that value is used for all engines.
    engine_extras : Optional[list[StorableObject]]
        extra things to store when an engine is checkpointed (typically
        collective variables, which enables storing of cached results)
    """
    def __init__(self, storage_handler, context="", frequency=1,
                 labeler=None, *, engine_extras=None, parent=None):
        if isinstance(frequency, int):
            frequency = collections.defaultdict(lambda: frequency)

        if labeler is None:
            labeler = default_context_label

        if engine_extras is None:
            engine_extras = []

        self.storage_handler = storage_handler
        self.context = pathlib.Path(context)
        self.frequency = frequency
        self.labeler = labeler
        self.engine_extras = engine_extras
        self.parent = parent

        self.children = []
        self._tempdir_manager = None

    def __enter__(self):
        return self

    def __exit__(self, exc, value, tb):
        if self._tempdir_manager:
            self._teardown_tempdir(exc, value, tb)
        if self.parent is None and exc is None:
            # we only delete the checkpoint on a successful run
            self.delete_checkpoint()

    def _get_storage_class(self):
        from openpathsampling.experimental.storage import (
            Storage, monkey_patch_all
        )
        monkey_patch_all(paths)
        return Storage

    def _setup_tempdir(self):
        """Create the temporary directory this checkpointer will use.
        """
        if self._tempdir_manager is not None:
            raise RuntimeError(f"The checkpoint at label '{self.context}' "
                               "is already in use.")
        self._tempdir_manager = tempfile.TemporaryDirectory()
        tempdir = pathlib.Path(self._tempdir_manager.__enter__())
        return pathlib.Path(tempdir.__enter__())

    def _teardown_tempdir(self, exc, value, tb):
        self._tempdir_manager.__exit__(exc, value, tb)
        self._tempdir_manager = None

    def next_context(self, context=None, sequence=0, *, parent=True):
        """Generate the next context within this checkpoint.

        This should always be used as a context manager.
        """
        parent = self if parent else None
        label = self.labeler(self, context, sequence)
        child = self.__class__(self.storage_handler,
                               self.context / label,
                               frequency=self.frequency,
                               labeler=self.labeler,
                               engine_extras=self.engine_extras,
                               parent=parent)
        self.children.append(child)
        return child

    def load_checkpoint(self):
        """Load data for this checkpoint.

        Returns
        -------
        data, files : Tuple[Dict | None, Dict | None]
            The first element of this 2-tuple is the ``data`` dictionary
            given when saving the checkpoint. The second element is a
            variant of the ``files`` dictionary. If no checkpoint has been
            saved for this context, this returns ``(None, None)``. NOTE: the
            files that are returned are in a temporary directory that will
            be deleted when this checkpoint's context ends. You should move
            them back to the correct location after calling this function.
        """
        if not self.context / "check.zip" in self.storage_handler:
            return None, None

        Storage = self._get_storage_class()

        tempdir = self._setup_tempdir()
        checkpoint = self.storage_handler.load(self.context / "check.zip",
                                               tempdir / "check.zip")
        with zipfile.ZipFile(tempdir / "check.zip", mode="r") as zipf:
            zipf.extractall(tempdir)

        st = Storage(tempdir / "check.db")
        data = st.tags["checkpoint"]
        file_labels = st.tags["files"]
        st.close()
        full_paths = {
            key: tempdir / archive_loc
            for key, archive_loc in file_labels.items()
        }

        for path in full_paths.values():
            assert path.exists()

        return data, full_paths

    def save_checkpoint(self, data, additional_files=None):
        """Save checkpoint data.

        This should only be called once per checkpoint context.

        Parameters
        ----------
        data : dict[str, Any]
            dictionary mapping variable names to values; values can be
            anything storable by OpenPathSampling
        additional_files: dict[str, os.Pathlike]
            dictionary mapping file label to Pathlike for the file. Note
            that the labels here are arbitrary, but will be used on reload
            to identify the same files, despite different paths
        """
        # annoyance for now (until SimStore is official storage)
        Storage = self._get_storage_class()

        if additional_files is None:
            additional_files = {}
        # complications here are based on the concern that the files given
        # could be arbitrary files, and we want them to extract into our
        # temp directory and ask the checkpoint caller to move them to where
        # they need to be.
        if not additional_files:
            arcnames = {}
        elif len(additional_files) == 1:
            arcnames = {key: str(pathlib.Path("files") / value.name)
                        for key, value in additional_files.items()}
        elif len(additional_files) > 1:
            common = os.path.commonpath(additional_files.values())
            arcnames = {
                key: str("files" / pathlib.Path(value).relative_to(common))
                for key, value in additional_files.items()
            }
        else:
            raise RuntimeError("This should never happen")

        tempdir = self._setup_tempdir()
        # zip everything into a checkpoint file
        st = Storage(tempdir / "check.db", mode='w')
        st.tags['checkpoint'] = data
        st.tags['files'] = arcnames
        st.close()
        with zipfile.ZipFile(tempdir / "check.zip", mode='x') as zipf:
            zipf.write(tempdir / "check.db", arcname="check.db")
            for key, file in additional_files.items():
                zipf.write(file, arcname=arcnames[key])

        self.storage_handler.store(self.context / "check.zip",
                                   tempdir / "check.zip")
        self._teardown_tempdir(None, None, None)

    def delete_checkpoint(self):
        """Delete this checkpoint and all its children.
        """
        for child in self.children:
            child.delete_checkpoint()

        contents = self.storage_handler.list_directory(self.context)
        if contents and contents != [self.context / 'check.zip']:
            raise RuntimeError("This should not happen")

        if contents:
            self.storage_handler.delete(self.context / "check.zip")


# TODO: Future PR
# class EnginesOnlyCheckpointer(Checkpointer):
#     """This is a checkpointer that only checkpoints engine movers.

#     The basic implementation is that the checkpoint data is stored
#     internally until it registers that it is in an EngineMover; then it
#     dumps all the stored state.

#     The idea is that engine movers are the only things that actually need
#     checkpointing; writing files for other movers is kind of overkill.
#     """
