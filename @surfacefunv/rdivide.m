function h = rdivide(f, g)
%./   Pointwise right divide for SURFACEFUNV.
%   F./G divides F by G, where F and G may be SURFACEFUNV objects or scalars.
%
%   See also LDIVIDE, COMPOSE.

h = surfacefunv;
if ( isempty(f) || isempty(g) )
    return
end

if ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    ncomp = numel(f.components);
    h.components = cell(1, ncomp);
    for k = 1:ncomp
        h.components{k} = f.components{k} ./ g.components{k};
    end
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefun') )
    ncomp = numel(f.components);
    h.components = cell(1, ncomp);
    for k = 1:ncomp
        h.components{k} = f.components{k} ./ g;
    end
elseif ( isa(f, 'surfacefun') && isa(g, 'surfacefunv') )
    ncomp = numel(g.components);
    h.components = cell(1, ncomp);
    for k = 1:ncomp
        h.components{k} = f ./ g.components{k};
    end
elseif ( isa(f, 'surfacefunv') && isnumeric(g) )
    ncomp = numel(f.components);
    h.components = cell(1, ncomp);
    if ( isscalar(g) )
        for k = 1:ncomp
            h.components{k} = f.components{k} ./ g;
        end
    elseif ( numel(g) == ncomp )
        for k = 1:ncomp
            h.components{k} = f.components{k} ./ g(k);
        end
    else
        error('SURFACEFUNV:rdivide:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
elseif ( isnumeric(f) && isa(g, 'surfacefunv') )
    ncomp = numel(g.components);
    h.components = cell(1, ncomp);
    if ( isscalar(f) )
        for k = 1:ncomp
            h.components{k} = f ./ g.components{k};
        end
    elseif ( numel(f) == ncomp )
        for k = 1:ncomp
            h.components{k} = f(k) ./ g.components{k};
        end
    else
        error('SURFACEFUNV:rdivide:invalid', ...
            'F and G must be surfacefunv objects, scalars, or constant vectors.');
    end
else
    error('SURFACEFUNV:rdivide:invalid', ...
        'F and G must be surfacefunv objects, scalars, or constant vectors.');
end

end
