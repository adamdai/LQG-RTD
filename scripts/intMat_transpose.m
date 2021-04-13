% returns the transpose of a matrix zonotope

function At = intervalM_transpose(A)
    Ct = center(A)';
    Rt = rad(A)';
    At = intervalMatrix(Ct,Rt);
end