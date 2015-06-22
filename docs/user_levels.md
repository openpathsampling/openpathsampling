# Levels of Users

OpenPathSampling aims to be easy for a beginning user to learn to use, but
also powerful enough to meet the needs of the most advanced users. To
accomplish that, we have tried to split the documentation into several
levels of experience. Our hope is that this will allow us to provide the
detailed information desired by power users without overwhelming novices.


## Novice user

The novice user may not be fully familiar with the methods involved in path
sampling. With OpenPathSampling, a novice user can still set up basic path
sampling calculations, 

## Normal user

The normal user is familiar with path sampling, and may need to do analysis
beyond our default abilities. As such, the normal user will need to learn
not just how to run a simulation in OpenPathSampling, but also some of the
underlying 

## Power user

A power user may be developing new methods which can be framed as
combinations of other methods. The power user will benefit from the
incredibly powerful and general parts that the tools in OpenPathSampling are
built from: in particular, the general treatment of `PathMover`s and of
`Ensemble`s are likely to be things a power user will want to learn about.

## Contributor

A contributor is willing to dive into the code to develop new methods beyond
what even the general tools provided to power users include. In particular,
this might include new simulation schemes or support for new molecular
dynamics engines. In addition to having a deep understanding of the APIs
defined within OpenPathSampling, the contributor should be very comfortable
working in Python.

---

We have tried to develop a framework which is useful for users at all these
levels, and which makes transitioning from one level to the next as easy as
possible when a user needs to. 
