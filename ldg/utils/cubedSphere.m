function [x, y, z] = cubedSphere(p, Nrefine)
% Constructs a mesh on the unit sphere.
%   p - mesh order (p + 1) points per patch
%   Nrefine - number of subdivisions of each squared panel

  % Coordinates of patch corners
  I = linspace(-1, 1, Nrefine + 1);
  J = linspace(-1, 1, Nrefine + 1);

  % Patch counter
  n = 1;

  % Initialize
  x1 = cell(1);
  y1 = cell(1);
  z1 = cell(1);

  % Create top patch centered at X = [0,0,1]
  for i = 1 : Nrefine
    for j = 1 : Nrefine
    
      % Grid
      u = chebpts(p + 1, I(i:i+1));
      v = chebpts(p + 1, J(j:j+1));
      w = 1;

      [u, v, w] = meshgrid(u, v, w);
  
      % Inflate points to unit sphere
      uijk = u ./ sqrt(u.^2 + v.^2 + w.^2);
      vijk = v ./ sqrt(u.^2 + v.^2 + w.^2);
      wijk = w ./ sqrt(u.^2 + v.^2 + w.^2);
  
      [x1{n}, y1{n}, z1{n}] = spherePatch(uijk, vijk, wijk);
      n = n + 1;

    end
  end

  % Rotate to cover other faces
  [x2, y2, z2] = rotate(x1, y1, z1, [0 0 1], [1 0 0]);
  [x3, y3, z3] = rotate(x1, y1, z1, [0 0 1], [0 1 0]);
  [x4, y4, z4] = rotate(x1, y1, z1, [0 0 1], [-1 0 0]);
  [x5, y5, z5] = rotate(x1, y1, z1, [0 0 1], [0 -1 0]);
  [x6, y6, z6] = rotate(x1, y1, z1, [0 0 1], [0 0 -1]);

  x = [x1 x2 x3 x4 x5 x6];
  y = [y1 y2 y3 y4 y5 y6];
  z = [z1 z2 z3 z4 z5 z6];

end