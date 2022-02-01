function forward_cube, X, Y, xieta
    ; ============================================================
    ; INPUT: X,Y IN RANGE -1 TO +1 ARE DATABASE CO-ORDINATES
    ; OUTPUT: [XI, ETA] EACH IN RANGE -1 TO +1 ARE 
    ;         TANGENT PLANE CO-ORDINATES
    ; BASED ON SUBROUTINE FORWARD_CUBE in COBE-DIRBE ES (A POLYNO-
    ; MIAL FIT FOUND USING FCFIT.FOR)
    ; ============================================================
    
    P=double([$  
        -0.27292696,    -0.07629969,    -0.02819452,    -0.22797056, $
        -0.01471565,    0.27058160,     0.54852384,     0.48051509, $
        -0.56800938,    -0.60441560,    -0.62930065,    -1.74114454, $
        0.30803317,     1.50880086,     0.93412077,     0.25795794, $
        1.71547508,     0.98938102,     -0.93678576,    -1.41601920, $
        -0.63915306,    0.02584375,     -0.53022337,    -0.83180469, $
        0.08693841,     0.33887446,     0.52032238,     0.14381585])

    xieta=dblarr(2) & XX=X*X & YY=Y*Y

    ; Calculating Xi and Eta
    xieta[0]=X*(1.0 + (1.0 - XX)*( $
        P[0]+XX*(P[1]+XX*(P[3]+XX*(P[6]+XX*(P[10]+XX*(P[15]+XX*P[21]))))) + $    
        YY*( P[2]+XX*(P[4]+XX*(P[7]+XX*(P[11]+XX*(P[16]+XX*P[22])))) + $
        YY*( P[5]+XX*(P[8]+XX*(P[12]+XX*(P[17]+XX*P[23]))) + $
        YY*( P[9]+XX*(P[13]+XX*(P[18]+XX*P[24])) + $
        YY*( P[14]+XX*(P[19]+XX*P[25]) + $
        YY*( P[20]+XX*P[26] + YY*P[27])))))))

    xieta[1]=Y*( 1.0 + (1.0 - YY)*( $
        P[0]+YY*(P[1]+YY*(P[3]+YY*(P[6]+YY*(P[10]+YY*(P[15]+YY*P[21]))))) + $    
        XX*( P[2]+YY*(P[4]+YY*(P[7]+YY*(P[11]+YY*(P[16]+YY*P[22])))) + $
        XX*( P[5]+YY*(P[8]+YY*(P[12]+YY*(P[17]+YY*P[23]))) + $
        XX*( P[9]+YY*(P[13]+YY*(P[18]+YY*P[24])) + $
        XX*( P[14]+YY*(P[19]+YY*P[25]) + $
        XX*( P[20]+YY*P[26] + XX*P[27])))))))

    return, xieta
end        
