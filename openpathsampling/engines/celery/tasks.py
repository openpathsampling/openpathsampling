from celery import Celery
from celery.signals import worker_process_init
import logging

logger = logging.getLogger(__name__)


app = Celery('tasks', backend='redis://localhost:6379/0', broker='redis://localhost:6379/0')
app.config_from_object('celeryconfig')


# This will force each child process to have a seperate INSTANCE_UUID
# only necessary if you fork your processes in which case python
# will reuse exising variables that _seem_ untouched and the ops package
# import will not be executed again
# on remote machine this should actually not be necessary

def _redo_uuid(**kwargs):
    import openpathsampling.netcdfplus as npl
    # logger.info('OLD UUID `%s`' % npl.StorableObject.INSTANCE_UUID)
    npl.StorableObject.initialize_uuid()
    logger.info('NEW UUID `%s`' % npl.StorableObject.get_instance_uuid())

# need to use registered instance for sender argument.
worker_process_init.connect(_redo_uuid)


@app.task(name='openpathsampling.engine.celery.tasks.generate')
def generate(engine, template, ensemble, init_args=None, init_kwargs=None):
    if init_args is None:
        init_args = []
    if init_kwargs is None:
        init_kwargs = {}

    engine.initialize(*init_args, **init_kwargs)
    traj = engine.generate(template, ensemble.can_append)
    return traj


@app.task(name='openpathsampling.engine.celery.tasks.list_cache')
def list_cache():
    """
    Return a list of UUIDs containing all objects present in the cache

    This is useful for checking which objects should actually be
    transferred. This is not save yet, becuase in the meantime the
    cache could loose some objects in which case you will not transmit these.
    This needs a function to lock the cache in that time and release after
    transmission


    Returns
    -------
    list of UUID
        the list of UUIDs of the objects that do not beed to be transmitted
    """
    return []

@app.task
def uuid(obj):
    return obj.__uuid__

@app.task
def identity(obj):
    return obj


@app.task
def attr(obj, attr):
    return getattr(obj, attr)


@app.task
def ls(obj):
    return dir(obj)


@app.task
def temp(obj):
    return obj._lazy.values()[0].__dict__


@app.task
def exev(template, s1, s2):
    exec(s1)
    return eval(s2)

