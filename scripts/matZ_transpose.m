% returns the transpose of a matrix zonotope

function At = matZ_transpose(A)
    Ct = A.center';
    if isempty(A.generator)
        At = matZonotope(Ct);
    else
        Gt = cellfun(@transpose,A.generator,'UniformOutput',false);
        At = matZonotope(Ct,Gt);
    end
end