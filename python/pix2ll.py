#=====================================================================
# PIX2LL (w/ dependancies)
#      PURPOSE:
#          This routine returns longitude and latitude pointing to center 
#          of pixel from COBE Quadrilaterised Spherical Cube given pixel 
#          number, cube resolution and end coordinate system. 
#          Additionally the spherical quadrilateralized cube face, 
#          x, y positions and face number can be returned if necessary.
#          
#      CALLING SEQUENCE: 
#          result=pix2ll(pixel, res, coord, nerd)    
#     
#      INPUT:
#          pixel - Int. Pixel number
#          resolution - INt. Quad-cube resolution (Note: 6 for DMR/FIRAS 
#                       and 9 for DIRBE)
#          coord - String. output coordinate system
#                  'E' or 'e': Geocentric True Ecliptic (default) 
#                  'Q' or 'q': Equatorial (ICRS) 
#                  'G' or 'g': Galactic. J2000. 
#          nerd - Int. When set to 1, returns face number and pixels 
#                 (x, y) coords. Takes 0 or 1. Default 0        
#     
#      OUTPUT:
#          result - Longitude and Latitude of center of pixel in 
#                      the specified coordinate system. When nerd=1, 
#                      keyword is chosen, an array of 5 elements
#                      is returned: [lon, lat, x, y, face_no]
#     
#      LIBRARIES USED:
#          NUMPY, ASTROPY 
#     
#      REVISION HISTORY:
#          Adapted from the FORTRAN script in COBE-DIRBE Explanatory 
#          Supplement.
#          Written by Gowtham, A. S.; December, 2022              
#=====================================================================

from numpy import *
import astropy.units as u
from astropy.coordinates import SkyCoord

def fc(x, y):
# ============================================================
#      INPUT: X,Y IN RANGE -1 TO +1 ARE DATABASE CO-ORDINATES
#      OUTPUT: [XI, ETA] EACH IN RANGE -1 TO +1 ARE 
#              TANGENT PLANE CO-ORDINATES
#      BASED ON SUBROUTINE FORWARD_CUBE in COBE-DIRBE ES 
#      (A POLYNOMIAL FIT ADAMPTED FROM FCFIT.FOR)
# ============================================================
    XX, YY=x*x, y*y
    P=array([-0.27292696,    -0.07629969,    -0.02819452,    -0.22797056,
            -0.01471565,    0.27058160,     0.54852384,     0.48051509,
            -0.56800938,    -0.60441560,    -0.62930065,    -1.74114454,
            0.30803317,     1.50880086,     0.93412077,     0.25795794,
            1.71547508,     0.98938102,     -0.93678576,    -1.41601920,
            -0.63915306,    0.02584375,     -0.53022337,    -0.83180469,
            0.08693841,     0.33887446,     0.52032238,     0.14381585])

    xi=x*(1.0 + (1.0 - XX)*( 
        P[0]+XX*(P[1]+XX*(P[3]+XX*(P[6]+XX*(P[10]+XX*(P[15]+XX*P[21]))))) +   
        YY*( P[2]+XX*(P[4]+XX*(P[7]+XX*(P[11]+XX*(P[16]+XX*P[22])))) + 
        YY*( P[5]+XX*(P[8]+XX*(P[12]+XX*(P[17]+XX*P[23]))) + 
        YY*( P[9]+XX*(P[13]+XX*(P[18]+XX*P[24])) + 
        YY*( P[14]+XX*(P[19]+XX*P[25]) + 
        YY*( P[20]+XX*P[26] + YY*P[27])))))))
    eta=y*( 1.0 + (1.0 - YY)*( 
        P[0]+YY*(P[1]+YY*(P[3]+YY*(P[6]+YY*(P[10]+YY*(P[15]+YY*P[21]))))) +     
        XX*( P[2]+YY*(P[4]+YY*(P[7]+YY*(P[11]+YY*(P[16]+YY*P[22])))) + 
        XX*( P[5]+YY*(P[8]+YY*(P[12]+YY*(P[17]+YY*P[23]))) + 
        XX*( P[9]+YY*(P[13]+YY*(P[18]+YY*P[24])) + 
        XX*( P[14]+YY*(P[19]+YY*P[25]) + 
        XX*( P[20]+YY*P[26] + XX*P[27])))))))

    return xi, eta   

def xyaxis(face, xi, eta):
# ==============================================================
#       CONVERTS FACE NUMBER (0-5) AND XI, ETA (b/w -1. and +1.) 
#       INTO A UNIT VECTOR 'OUT'
# ==============================================================

# To preserve symmetry, the normalization sum must always 
# have the same ordering (i.e. largest to smallest).
    xi1, eta1 = max([abs(xi), abs(eta)]), min([abs(xi), abs(eta)])
    norm = 1.0/sqrt(1.0 + xi1**2.0 + eta1**2.0)

    out = [0.0, 0.0, 0.0]
    match face:
        case 1:
            out = [1.0, xi, eta]
        case 2:
            out = [-xi, 1.0, eta]
        case 0:
            out = [-eta, xi, 1.0]            
        case 3:
            out = [-1.0, -xi, eta]
        case 4:
            out = [xi, -1.0, eta]
        case 5:
            out = [eta, xi, -1.0]
        case _:
            print("ERR: XYAXIS: Incorrect FACE number")   

    result = array(out)*norm
    
    return result

def pixel_vec(pixel, res):
# ============================================================
#      Routine to return unit vector pointing to center of 
#      pixel given pixel number and resolution of the cube.
# ============================================================
    one, two = int(1), int(2) 
    scale = 2.0**(res-1.0)/2.0
    ppf = two**(two*(res-one)) # Pixels per face
    
    face = int(pixel/ppf)
    fpix = int(pixel - face*ppf)

    bin = binary_repr(fpix, width=16)
    bin = array([bit for bit in reversed(bin)])
    ix, iy, xc, yc = 0.0, 0.0, 0.0, 0.0

    for i in range(0, 15):
        if i % 2 == 0:
            ix = ix+(int(bin[i])*(2.0**xc))
            xc = xc+1.0
        else:
            iy = iy+(int(bin[i])*(2.0**yc))
            yc = yc+1.0     
    ix, iy = int(ix), int(iy)
    x, y = (ix - scale + 0.5)/scale, (iy - scale + 0.5)/scale
    xi, eta=fc(x, y)
    vec=xyaxis(face, xi, eta)

    return ix, iy, face, vec

def pix2ll(pixel, res, coord, nerd):
# ============================================================
#      Main Routine
# ============================================================

# Finding vectors in CSC coordinate system
    x, y, face, vec = pixel_vec(pixel, res)

# Converting CSC unit vectors to Ecliptic Longitudes and Latitudes 
    norm_inv = sqrt(vec[0]**2.0 + vec[1]**2.0 + vec[2]**2.0)
    if norm_inv != 0.0:
        vec=vec/norm_inv
    else:
        print("Err: PIX2LL: Unit vector normalisation is 0")
# Finding Longitude
    if vec[0] == 0 and vec[1] == 0:
        lon = 0.0
    else:
        lon = double(180.0/pi)*arctan2(vec[1], vec[0])
    if lon < 0.0: lon = lon+360.0    

# Finding Latitude  
    lat = arctan2(vec[2], sqrt(vec[0]**2.0 + vec[1]**2.0))
    lat = lat*double(180.0/pi)  

    lonlat = SkyCoord(lon=lon*u.degree, lat=lat*u.degree, frame='geocentrictrueecliptic')

# Converting to Equatorial or Galactic if necessary      
    match coord:
        case 'E' | 'e':
            a, b = lon, lat
        case 'Q' | 'q':
            radec = lonlat.transform_to('icrs')
            a, b = float(radec.ra/u.degree), float(radec.dec/u.degree)
        case 'G' | 'g':
            lb = lonlat.transform_to('galactic')
            a, b = float(lb.l/u.degree), float(lb.b/u.degree)       
        case _:
            print("ERR: PIX2LL: Incorrect CoordSys Choice")  
            print("ERR: PIX2LL: Choose 'G' or'g' for Galactic System")    
            print("ERR: PIX2LL: Choose 'E' or'E' for Ecliptic System")  
            print("ERR: PIX2LL: Choose 'Q' or'q' for Equatorial System")      

    if (nerd == 1): result=[a, b, x, y, face]
    else: result=[a, b]

    return result

#Driver code
print(pix2ll(3, 6, 'e', 1))
