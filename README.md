# BITCoPT
This is the second version of the BITCoPT code developed by Daniel W. Zaide, Justin Shultz, and Assad Oberai, through a collaboration between the Scientific Computation Research Center (SCOREC) and Center for Architecture Science and Ecology (CASE) at the Rensselaer Polytechnic Institute.

## Documentation 
General documentation is in doc/guide.tex. There is other documentation in the doc directory. Comments have been added throughout the code to describe functions, inputs, and output. 

## Dependencies
There are three main dependencies that may need to be installed: pytz for timezones, solar for solar calculations, and shapely for geometries. This code works best with Python 2.x.

### Manual Install
This method of installation has been more successful.

#### Install pytz (http://pytz.sourceforge.net/) using:
```
pip install pytz
```
If installed already, you can move to the next step.


#### Install Pysolar (https://github.com/pingswept/pysolar), legacy 0.6 for Python 2.x:
* Download source code from https://github.com/pingswept/pysolar/releases/tag/0.6
* Unzip and move folder to accessible location
* Open the pysolar-0.6 folder
* If on Windows: Click on the navigation bar and type CMD
* If on Mac/Linux: Navigate CMD/Console into the folder
* Within CMD/Console type:
```
python setup.py install
```

#### Install Shapely (https://pypi.python.org/pypi/Shapely):

For Mac OSX and Linux:
```
pip install shapely
```

For Windows:
* Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely
* Search for shapely 
* Download the correct cp27 *.whl file that matches your python installation, 32 or 64 bit
* Navigate to download folder and type "cmd" in the path/url bar at the top OR open CMD and navigate to downloads folder
* Type into CMD:
```
pip install [name of shapely file].whl
```

#### Install joblib (https://pypi.python.org/pypi/joblib):
```
pip install joblib
```

#### Install numdifftools (https://pypi.python.org/pypi/Numdifftools):
```
pip install numdifftools
```

#### Install eppy (http://pythonhosted.org//eppy/):
```
pip install eppy
```


### Automatic Install - not currently working
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

Dependencies can be installed using:
```
  python setup.py develop
```

## Tests
Tests can be run with python test.py.

## Examples

* example.py
* diffusion2D.py
* poisson2D.py

