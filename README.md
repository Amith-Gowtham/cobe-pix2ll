# ***Amith's PIX2LL*** 
> Convert COBE's Quadrilateralised Spherical Cube (CSC) pixel numbers to Spherical Coordinate Values

Scripts in Python and IDL/GDL to perform CSC number->Lat-Lon of any specified coordinate system (Ecliptic, Galactic or Equatorial)

Adapted from a similar programme written in FORTRAN found in COBE-DIRBE Explanatory Supplement (https://lambda.gsfc.nasa.gov/product/cobe/dirbe_exsup.cfm)

### _How to run the IDL scripts_
Post-cloning, open the working dir in IDL/GDL cmd prompt of choice and type in:

`result=pix2ll(*pixel_number, *resolution, *CoordSys)`

> _**Note**: The IDL programme were written in GDL 0.9.7._

### _How to run the Python scripts_
Make similar necessary changes to the Driver Code at the end of the python script and run

> _**Note**: Requires Python version 9.10 and above_




