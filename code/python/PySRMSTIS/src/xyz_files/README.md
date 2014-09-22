### Quick translation of XYZ file to netCDF

Includes tests suitable for nose testing framework. Currently only handles a
single trajectory at a time.

Converting .xyz to .nc is easy:

    python xyztranslator.py -o output.nc input.xyz

However, that won't work if output.nc already exists.

Converting .nc to .xyz requires an additional .xyz file to provide the
topology (in this case, just atom names, and taken from the first frame of
the .xyz trjectory):

    python xyztranslator.py -o output.xyz -t topol.xyz input.nc

The file `geisslerene000000.xyz` is provided as an example in the .xyz
format.
