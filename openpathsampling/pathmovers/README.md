# Adding new path movers

Path movers are among the more complicated things to add to OpenPathSampling.
This is because path movers may involve multiple classes. Typically, a new path
mover includes at least one `PathMover` subclass. But most also include a
`MoveStrategy` subclass, and some include a specialized `MoveScheme` subclass.

We recommend putting all these class definitions in the same file. The easiest
way to learn how to implement something is to have a complete example of
something similar, and by keeping all the relevant code together, your
contributions can be good examples for future contributors.

Then you can import the relevant classes into the correct subpackage handlers.
Your subclass of `PathMover` (and any associated helper classes) should be
imported into `pathmovers/__init__.py`.  If you have a `MoveStrategy` subclass
that should be publicly exposed, you should import it into
`pathmovers/move_strategies.py`. If you have a `MoveScheme` subclass to expose,
import it into `pathmovers/move_schemes.py`.  The classes that you expose in
this way will be automatically exposed in the root `openpathsampling` namespace
(or `openpathsampling.strategies`, in the case of `MoveStrategy` subclasses).
