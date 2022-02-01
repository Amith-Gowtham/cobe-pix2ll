function xyaxis, face, xi, eta, out
    ; ============================================================
    ; CONVERTS FACE NUMBER NFACE (0-5) AND XI, ETA (-1. - +1.) 
    ; INTO A UNIT VECTOR C
    ; ============================================================

    ; To preserve symmetry, the normalization sum must always 
    ; have the same ordering (i.e. largest to smallest).
    xi1=max([abs(xi), abs(eta)])
    eta1=min([abs(xi), abs(eta)])
    norm=1.0/sqrt(1.0 + xi1^2.0 + eta1^2.0)

    out=fltarr(3)

    case face of
        0: begin
            out[2]= norm
            out[0]= -eta*norm 
            out[1]= xi*norm
        end
        1: begin
            out[0]= norm
            out[2]= eta*norm 
            out[1]= xi*norm
        end
        2: begin
            out[1]= norm
            out[2]= eta*norm 
            out[0]= -xi*norm
        end
        3: begin
            out[0]= -norm
            out[2]= eta*norm 
            out[1]= -xi*norm
        end
        4: begin
            out[1]= -norm
            out[2]= eta*norm 
            out[0]= xi*norm
        end
        5: begin
            out[2]= -norm
            out[0]= eta*norm 
            out[1]= xi*norm
        end
        else: stop, 'ERR: FORWARD_CUBE: Incorrect FACE number. 0-5.'
    endcase
    return, out
end
