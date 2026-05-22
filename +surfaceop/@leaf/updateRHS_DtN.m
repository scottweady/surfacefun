function P = updateRHS_DtN(P, rhs)
%UPDATERHS   Update RHS of an SURFACEOP.LEAF object.
%   P = UPDATERHS(P, F) replaces the existing RHS of an initialized
%   SURFACEOP.LEAF object P with that given in F, which must be a matrix
%   (or cell array) of tensor-product Chebyshev values. Since only one
%   subproblem needs to be solved on each patch in the update process
%   (rather than O(n) in the original initialization) this can lead to
%   considerable performance gains when solving for multiple RHSs.
%
%   See also SURFACEOP.LEAF.INITIALIZE.

% Developer note: At user levels this is typically called with 
%  >> P.rhs = F

n = P.n;
dom = P.domain;
id = P.id;
ncomp = size(P.u_part, 1)/n^2;

nrhs = 1;
ii = false(n);
ii(2:n-1,2:n-1) = true;

if ( iscell(rhs) && isa(rhs{1}, 'function_handle') )
    rhs = rhs{1};
end

if ( iscell(rhs) )
    if ( numel(rhs) == 1 && isnumeric(rhs{1}) && isvector(rhs{1}) && ...
            numel(rhs{1}) ~= n^2 )
        rhs = rhs{1};
    else
        nrhs = size(rhs, 2);
        rhs = reshape([rhs{1,:}], n^2, nrhs);
    end
end

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    % Constant RHS.
    rhs = repmat(rhs, n^2, 1);
elseif ( isnumeric(rhs) && isvector(rhs) && numel(rhs) ~= n^2 )
    % Constant vector RHS.
    nrhs = numel(rhs);
    rhs = repmat(reshape(rhs, 1, nrhs), n^2, 1);
elseif ( isnumeric(rhs) && ~isscalar(rhs) )
    % We already have the values of the RHS.
elseif ( isa(rhs, 'function_handle') )
    rhs = feval(rhs, dom.x{id}, dom.y{id}, dom.z{id});
    if ( isvector(rhs) && numel(rhs) ~= n^2 )
        nrhs = numel(rhs);
        rhs = repmat(reshape(rhs, 1, nrhs), n^2, 1);
    else
        rhs = reshape(rhs, n^2, []);
    end
end

if ( ncomp > 1 )
    if ( size(rhs, 2) == 1 )
        rhs = repmat(rhs, 1, ncomp);
    elseif ( size(rhs, 2) ~= ncomp )
        error('SURFACEOP:LEAF:updateRHS:rhsDimensions', ...
            'Coupled vector PDE righthand sides must have %d components.', ncomp);
    end
    nrhs = 1;
    rhs = reshape(rhs(ii,:), [], 1);
else
    nrhs = size(rhs, 2);
    rhs = rhs(ii,:);
end

% Scale by the Jacobian if the element is singular
if ( dom.singular(id) )
    rhsScale = repmat(dom.J{id}(ii).^3, ncomp, 1);
    rhs = rhsScale .* rhs;
end

% Apply A^-1
if ( isempty(P.Ainv) )
    error('SURFACEOP:LEAF:updateRHS:operatorNotStored', ...
        'Discretized operator A was not stored. Cannot update RHS.');
end
S = P.Ainv(rhs);

% Append boundary points to solution operator:
tmpS = zeros(ncomp*n^2, nrhs);
if ( ncomp == 1 )
    iiSol = ii;
else
    iiSol = blockIdx(find(ii), ncomp, n^2);
end
tmpS(iiSol,:) = S;
S = tmpS;
P.u_part = S;

% Normal derivative:
P.du_part = P.normal_d * S;

end

function idx = blockIdx(idx, ncomp, nbase)

idx = idx(:);
idx = idx + (0:ncomp-1)*nbase;
idx = idx(:);

end
