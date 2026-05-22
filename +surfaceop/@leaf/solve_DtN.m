function u = solve_DtN(P, bc)
%SOLVE   Solve a leaf patch.
%   U = SOLVE(P, BC) returns a cell array containing the solution values U.

% Extract the domain from the patch:
id = P.id;
n = size(P.domain.x{id}, 1);
nb = size(P.xyz, 1);
ncomp = size(P.S, 2)/nb;

if ( ~isnumeric(bc) )
    % Evaluate the RHS if given a function handle:
    bc = feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3));
    if ( ncomp > 1 && ismatrix(bc) && size(bc, 2) == ncomp )
        bc = reshape(bc, [], 1);
    elseif ( ncomp > 1 && isvector(bc) && numel(bc) == nb )
        bc = repmat(bc(:), ncomp, 1);
    end
elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = repmat(bc, size(P.S, 2), 1);
elseif ( ncomp > 1 && ismatrix(bc) && size(bc, 1) == nb && size(bc, 2) == ncomp )
    bc = reshape(bc, [], 1);
elseif ( ncomp > 1 && isvector(bc) && numel(bc) == nb )
    bc = repmat(bc(:), ncomp, 1);
end

% Evaluate the solution operator for the patch:
u = P.S * bc + P.u_part;

% Return cell output for consistency with PARENT/SOLVE():
if ( ncomp == 1 )
    U = cell(1, size(u, 2));
    for k = 1:size(u, 2)
        U{k} = reshape(u(:,k), n, n);
    end
else
    U = cell(1, ncomp);
    for k = 1:ncomp
        rows = (k-1)*n^2 + (1:n^2);
        U{k} = reshape(u(rows,1), n, n);
    end
end
u = U;

end
