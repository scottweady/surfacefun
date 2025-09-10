function dom = image(f, dom)
%IMAGE   Map a planar grid to a surface using a parametric map.
% f: function handle of the form f(u,v) = [x y z], where (u,v) are
%   coordinates in the plane and (x,y,z) are coordinates on the surface.
% dom: reference domain (planar)
%

    if nargin < 2
        error('Not enough input arguments.');
    end
    
    if ~isa(dom, 'surfacemesh')
        error('DOM must be a SURFACEMESH object.');
    end

    % Get reference coordinates
    u = dom.x;
    v = dom.y;

    % Initialize output
    x = cell(size(u));
    y = cell(size(u));
    z = cell(size(u));

    % Get reference grid
    for elem = 1 : length(u)
        shp = size(u{elem});
        xyz = f(u{elem}(:), v{elem}(:));
        x{elem} = reshape(xyz(:,1), shp);
        y{elem} = reshape(xyz(:,2), shp);
        z{elem} = reshape(xyz(:,3), shp);
    end
    
    % Create surfacemesh object
    dom = surfacemesh(x, y, z);

end