from setuptools import setup

setup(name='BITCoPT',
      version='2.0',
      description='A simple python code for modelling BITCoPT and dynamic facades.',
      url='https://github.com/RPI-CASE',
      author='Daniel W. Zaide',
      author_email='dan.zaide@gmail.com',
      license='MIT',
      install_requires=[
          'pytz',
          'solar',
          'shapely'
      ],
      zip_safe=False)