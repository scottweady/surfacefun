function h = plus(f, g)
%+   Plus for SURFACEFUNV.
%   F + G adds the SURFACEFUNV F and G. F and G must have the same domains
%   and discretization sizes. F and G may also be scalars or constant vectors.
%
%   See also MINUS.

if ( isnumeric(f) )
    h = plus(g, f);
    return
end

h = f;
ncomp = numel(f.components);

if ( isnumeric(g) )
    if ( isscalar(g) )
        for k = 1:ncomp
            h.components{k} = f.components{k} + g;
        end
    elseif ( numel(g) == ncomp )
        for k = 1:ncomp
            h.components{k} = f.components{k} + g(k);
        end
    else
        error('SURFACEFUNV:plus:invalid', ...
            'F and G must be surfacefunv objects, constant vectors, or scalars.');
    end
elseif ( isa(f, 'surfacefunv') && isa(g, 'surfacefunv') )
    for k = 1:ncomp
        h.components{k} = f.components{k} + g.components{k};
    end
else
    error('SURFACEFUNV:plus:invalid', ...
        'F and G must be surfacefunv objects, constant vectors, or scalars.');
end

end
