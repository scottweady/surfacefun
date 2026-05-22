function u = solve(S, bc)
%SOLVE   Perform a surface solve.
%   U = SOLVE(S, BC) returns a cell array U of function values representing
%   the solution to the PDE, subject to the Dirichlet boundary conditions
%   specified by the function handle BC, contained in the SURFACEOP object
%   S. If S has not yet been built (see BUILD()) then SOLVE() will build
%   it. If S has not yet been initialized (see INITIALIZE()) then an error
%   is thrown.
%
%   The full sequence for solving a problem using a SURFACEOP object S is:
%
%      initialize(S, OP, RHS)
%      build(S)
%      u = S\bc % or u = solve(S, bc)
%
%   See also BUILD, INITIALIZE.
%
%   If S was initialized with a SURFACEFUNV righthand side, or BC is a
%   vector-valued function handle/constant vector, U is returned as a
%   SURFACEFUNV.

% Build if required:
if ( ~isBuilt(S) )
    build(S);
end

if ( ~isInitialized(S) )
    error('The object has not been initialized.');
end

if ( nargin == 1 )
    bc = [];
end

[bc, bcOutputType] = parseBoundaryData(bc, S.patches{1}.xyz, S.numComponents);

switch S.method
    case 'DtN'
        solve = @solve_DtN;
    case 'ItI'
        solve = @solve_ItI;
end

% Solve the patch object:
u = solve(S.patches{1}, bc);

if ( strcmp(S.outputType, 'surfacefunv') || strcmp(bcOutputType, 'surfacefunv') )
    ncomp = size(u, 2);
    if ( ncomp < 2 )
        error('SURFACEOP:solve:vectorDimensions', ...
            'Vector PDE solves require one column per component.');
    end
    components = arrayfun(@(k) surfacefun(u(:,k), S.domain), 1:ncomp, ...
                          'UniformOutput', false);
    u = surfacefunv(components{:});
else
    % Package into a surfacefun:
    u = surfacefun(u, S.domain);
end

end

function [bc, outputType] = parseBoundaryData(bc, xyz, ncomp)
%PARSEBOUNDARYDATA   Normalize vector-valued Dirichlet data at the top level.

outputType = 'surfacefun';
fromFunction = false;

if ( isempty(bc) )
    return
elseif ( isa(bc, 'surfacefunv') )
    error('SURFACEOP:solve:surfacefunvBC', ...
        ['SURFACEFUNV boundary data cannot be evaluated on the solver ' ...
         'skeleton. Use a function handle or a constant vector instead.']);
elseif ( isa(bc, 'function_handle') )
    bc = feval(bc, xyz(:,1), xyz(:,2), xyz(:,3));
    fromFunction = true;
end

if ( isnumeric(bc) && isvector(bc) && numel(bc) == ncomp && ncomp >= 2 )
    bc = repmat(reshape(bc, 1, ncomp), size(xyz, 1), 1);
    outputType = 'surfacefunv';
elseif ( fromFunction && isnumeric(bc) && ismatrix(bc) && size(bc, 2) == ncomp && ncomp >= 2 )
    outputType = 'surfacefunv';
end

if ( ncomp > 1 && strcmp(outputType, 'surfacefunv') )
    bc = reshape(bc, [], 1);
elseif ( ncomp > 1 && isnumeric(bc) && isvector(bc) && numel(bc) == size(xyz, 1) )
    bc = repmat(bc(:), ncomp, 1);
end

end
