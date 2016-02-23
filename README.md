# BITCoPT
This is the second version of the BITCoPT code developed by Daniel W. Zaide.

## Documentation 
General documentation is in doc/guide.tex. There is other documentation in the doc directory.

## Dependencies
There are three main dependencies that may need to be installed: pytz for timezones, solar for solar calculations, and shapely for geometries.

### Automatic Install
Update distribute libraries with (may require admin):
```
pip install --upgrade pip
```
```
pip install --upgrade distribute
```
```
pip install --upgrade setuptools
```
Dependencies can be obtained using:
```
  python setup.py develop
```

### Manual Install
This method of installation has been more successful.

Install pytz (http://pytz.sourceforge.net/) using:
```
pip install pytz
```

Install Pysolar (https://github.com/pingswept/pysolar), legacy 0.6 for Python 2.x.
* Download source code from https://github.com/pingswept/pysolar/releases/tag/0.6
* Unzip and move folder to accessible location
* Within CMD or Command Line type:
```
python setup.py install
```

Install Shapely (https://pypi.python.org/pypi/Shapely):
```
pip install shapely
```

## Tests
Tests can be run with python test.py.

## Examples

* example.py
* diffusion2D.py
* poisson2D.py

