function pixel_vector, pixel, res, out
    ; ============================================================
    ; Routine to return unit vector pointing to center of 
    ; pixel given pixel number and resolution of the cube.
    ; ============================================================
    
    one=1l & two=2l
    scale=long(2^(res-2)) 
    ppf=long(two^(two*(res-one))) ; Pixels per face

    face=fix(pixel/ppf) ; Note: Integer division truncates
    fpix=fix(pixel - face*ppf)

    bin=reverse(binary(fpix))
    ix=0.0 & iy=0.0
    jx=[] & jy=[]

    ; Break pixel number down into x and y bits:
    for i=0, 15 do if (i mod 2 eq 0) then jx=[jx, string(bin[i])] $
        else jy=[jy, string(bin[i])]
    for i=0, 7 do begin
        ix=ix+(jx[i]*(2^i)) 
        iy=iy+(jy[i]*(2^i))
    endfor    

    x=(fix(ix) - scale + 0.5)/scale
    y=(fix(iy) - scale + 0.5)/scale
    xieta=forward_cube(x, y)
    vec=xyaxis(face, xieta[0], xieta[1])

    out={face:0.0, x:0, y:0, vector:dblarr(3)}
    out.face=face & out.x=fix(ix) & out.y=fix(iy) & out.vector=vec
    return, out
end
