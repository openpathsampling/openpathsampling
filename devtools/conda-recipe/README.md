This is a recipe for building the current development package into a conda
binary.

The installation on travis-ci is done by building the conda package, installing
it, running the tests, and then if successful pushing the package to anaconda
(and the docs to AWS S3). The anaconda auth token is an encrypted environment
variable generated using:

```bash
anaconda auth -n [token_name] -o omnia --max-age 22896000 -c --scopes api
```
and then saved in the environment variable `ANACONDA_TOKEN`. You can limit the rights of the token by using `api:write`. The rights of an anaconda token can also be changed in your anaconda profile at [www.anaconda.org]().

You can set up travis to store an encrypted token via

```bash
gem install travis
travis encrypt ANACONDA_TOKEN=xx
```

where xx is the token output by binstar.  The final command should print a line (containing 'secure') for inclusion in your .travis.yml file.
