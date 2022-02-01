function uv2ll, in, out

    norm_inv=sqrt(in[0]^2.0 + in[1]^2.0 + in[2]^2.0)
    if (norm_inv ne 0.0) then norm_inv=1.0/norm_inv $
    else stop, 'Err: UV2LL: Unit vector normalisation is 0'

    v=dblarr(3)
    v[0]=in[0] * norm_inv
    v[1]=in[1] * norm_inv
    v[2]=in[2] * norm_inv

    ; Finding Ecliptic Longitude
    if ((v[1] eq 0.0) and (v[0] eq 0.00)) then lon=0.0 $
    else lon=double(180.0/!PI)*atan(v[1], v[0])
    if (lon lt 0.0) then lon=lon+360.0

    ; Finding Ecliptic Latitude
    lat=double(180.0/!PI)*atan(v[2], sqrt(v[0]^2.0 + v[1]^2.0))

    out=float([lon, lat])
    return, out
end
