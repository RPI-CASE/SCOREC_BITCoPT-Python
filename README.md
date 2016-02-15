# BITCoPT
This is the second version of the BITCoPT code developed by Daniel W. Zaide.

## Documentation 
General documentation is in doc/guide.tex. There is other documentation in the doc directory.

## Dependencies
There are three main dependencies that may need to be installed: pytz for timezones, solar for solar calculations, and shapely for geometries.

Dependencies can be obtained using

python setup.py develop

Or manually installed with

pip install pytz

pip install solar
* Requires pysolar 0.6 for Python 2.X - https://github.com/pingswept/pysolar/releases/tag/0.6

pip install shapely


## Tests
Tests can be run with python test.py.

## Examples

* example.py
* diffusion2D.py
* poisson2D.py

