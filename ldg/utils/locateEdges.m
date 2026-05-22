
function [xedge, yedge, zedge] = locateEdges(v0, X, x, y, z)
% Locates edges that intersect the inclusion patch 
%   v0 - angle of inclusion
%   X - coordinate of inclusion
%   x, y, z - coordinates within patch

  % Initialize
  xedge = [];
  yedge = [];
  zedge = [];

  % Distance threshold
  tol = 3 * (1 - cos(v0));

  % Squared distance to center
  xsq = (x - X(1)).^2 + (y - X(2)).^2 + (z - X(3)).^2;

  % Bottom edge
  if all(xsq(:, 1) > tol)
    xedge = [xedge; x(:, 1)'];
    yedge = [yedge; y(:, 1)'];
    zedge = [zedge; z(:, 1)'];
  end

  % Top edge
  if all(xsq(:, end) > tol)
    xedge = [xedge; x(:, end)'];
    yedge = [yedge; y(:, end)'];
    zedge = [zedge; z(:, end)'];
  end
  
  % Left edge
  if all(xsq(1, :) > tol)
    xedge = [xedge; x(1, :)];
    yedge = [yedge; y(1, :)];
    zedge = [zedge; z(1, :)];
  end

  % Right edge
  if all(xsq(end, :) > tol)
    xedge = [xedge; x(end, :)];
    yedge = [yedge; y(end, :)];
    zedge = [zedge; z(end, :)];
  end
   
end
