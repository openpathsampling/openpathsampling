So you want to add a new dynamics engine to OpenPathSampling? Great! We've
tried to make it as easy as possible, because supporting more engines is a
quick way to make OpenPathSampling useful to more people.

There are two main approaches to controlling a dynamics engine from
OpenPathSampling: direct control and indirect control. The first thing you
need to decide is which approach is best for your engine.

Direct control is when you have an API that directly controls the dynamics:
for example, you can ask the dynamics to generate 10 frames, and it gives
you those 10 frames and then just waits for the next request.  Direct
control works best when there's already a Python API for your engine, or
when start-up time is nearly nothing compared to obtaining a single MD step.

With indirect control, we take a different approach: we launch a trajectory
as a separate process, and check in on the trajectory file it outputs until
we hit the stopping condition. Then we kill the subprocess and clean up.

---

For both APIs, the main ideas are the same: your subclass needs to translate
information back and forth between our `Snapshot` class and your engine's
internal data storage, and your subclass needs to be able to tell your
engine to run its dynamics. The main difference is that the indirect
approach uses the file system and external processes as intermediates for
this translation.

## API elements for both control models

Most parts of creating a new engine are independent of the approach to
controlling the dynamics. 

**TODO: list this stuff**

## Direct control API

If you intend to use direct control, your engine class should inherit from
`paths.engines.DynamicsEngine`. You will need to implement these methods:

**TODO: more on the direct API**

## Indirect control API

If you intend to use indirect control, your engine class should inherit from
`paths.engines.ExternalEngine`, and you should consider overriding the
following methods:

* `read_frame_from_file(filename, frame_num)`: reads the frame from the
  external engine's file format
* `write_frame_to_file(filename, snapshot, mode="a")`: writes the frame to
  the external engine's format (used to initiate trajectories)
* `engine_command()`: returns a string of the command to be called by the
  operating system
* `set_filenames(number)`: sets the file names for step `number`
* `cleanup()`: does any clean-up tasks needed after killing the subprocess
  (e.g., removing temporary files)

Additionally, you may wish to override the following options:

* `killsig` (class variable): the signal sent to terminate the process
  (default is `signal.SIGTERM`).
* `default_sleep_ms` (set in `options`): time the engines sleeps before
  checking again whether a new frame has been written. Note that
  `ExternalEngine`s will automatically optimize the sleep time until you set
  the option `auto_optimize_sleep` to `False`. 

## Testing your new engine

We strongly recommend developing thorough unit tests for your engine, and we
require unit tests for any engines added to the OpenPathSampling core. In
particular, there are some special cases for each of direct and indirect
control that your tests should include.

### Extra tests for direct control

### Extra tests for indirect control

Indirect control requires a few extra tests, to ensure that your file
reading is working correctly.

* Write a test which writes half a frame to a file, then tries to
  `get_next_frame` (this part should return `None`), then finish writing the
  frame to the file, and try to `get_next_frame` (now it should work). This
  simulates the case that `get_next_frame` is called while the frame is
  being written. (You might also check for various EOL options, in case
  the engine uses platform-dependent EOL choices.)
