function L = laplacian(f)
%LAPLACIAN   Laplacian of a SURFACEFUNV.
%   LAPLACIAN(F) returns a SURFACEFUNV whose components are the Laplacians
%   of the components of F.
%
%   See also LAP.

if ( isempty(f) )
    L = surfacefunv;
    return
end

comps = cellfun(@laplacian, f.components, 'UniformOutput', false);
L = surfacefunv(comps{:});

end
