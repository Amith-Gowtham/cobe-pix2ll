function pix2ll, pixel ,res, coord, result, all=all
    ; ============================================================
    ; PURPOSE:
    ;     Routine returns longitude and latitude pointing to center 
    ;     of pixel given pixel number and resolution of the cube and  
    ;     coordinate system. Note, all coordintes are at epoch=J2000.
    ;     Additionally the spherical quadrilateralized cube face, 
    ;     x, y positions can be returned if necessary.
    ;     
    ; CALLING SEQUENCE: 
    ;     result=pix2ll(pixel, res, coord)    
    ;
    ; INPUT:
    ;     pixel - Pixel number
    ;     resolution - Quad-cube resolution
    ;     coord - output coordinate system
    ;
    ; OUTPUT:
    ;     result - Longitude and Latitude of center of pixel in 
    ;                 the specified coordinate system. If /all 
    ;                 keyword is chosen, an array of 5 elements
    ;                 is returned: [lon, lat, x, y, face_no]
    ;
    ; SUBROUTINES CALLED:
    ;     PIXEL_VECTOR, FORWARD_CUBE, XYAXIS, UV2ll, 
    ;     CONV_E2G, CONV_E2Q, EULER     
    ;
    ; REVISION HISTORY:
    ;     Adapted from the FORTRAN script in COBE-DIRBE Explanatory 
    ;     Supplement.
    ;     Written by Gowtham, A. S.; November, 2021              
    ; ============================================================


    out=pixel_vector(pixel, res)

    ; Finding Ecliptic spherical coords
    lonlat=uv2ll(out.vector)

    ; Converting Ecliptic lonlat to either Galactic  
    ; or Equatorial if necessary
    if ((coord eq 'E') or (coord eq 'e')) then begin
        result=lonlat
    endif
    if ((coord eq 'G') or (coord eq 'g')) then begin
        euler, lonlat[0], lonlat[1], gl, gb, 5 
        result=[gl, gb]
    endif
    if ((coord eq 'Q') or (coord eq 'q')) then begin
        euler, lonlat[0], lonlat[1], ra, dec, 4
        result=[ra, dec]
    endif

    if (keyword_set(all) eq 1) then result=[result, out.x, out.y, out.face]

    return, result
end
