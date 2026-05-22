function h = times(f, g)
%.*   Pointwise multiplication for SURFACEFUNV.
%   F.*G multiplies F and G pointwise. Either F or G may be a scalar,
%   SURFACEFUN, or numeric vector.
%
%   See also MTIMES, RDIVIDE.

if ( isnumeric(f) || isa(f, 'surfacefun') )
    h = times(g, f);
    return
end

% f is surfacefunv
ncomp = numel(f.components);
h = f;

if ( isscalar(g) )
    for k = 1:ncomp
        h.components{k} = f.components{k} .* g;
    end
elseif ( isnumeric(g) && numel(g) == ncomp )
    for k = 1:ncomp
        h.components{k} = f.components{k} .* g(k);
    end
elseif ( isa(g, 'surfacefun') )
    for k = 1:ncomp
        h.components{k} = f.components{k} .* g;
    end
elseif ( isa(g, 'surfacefunv') )
    for k = 1:ncomp
        h.components{k} = f.components{k} .* g.components{k};
    end
else
    error('SURFACEFUNV:times:invalid', ...
        'Unsupported operand types for .*');
end

end
