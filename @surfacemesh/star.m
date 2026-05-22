function dom = star(n, a, m, nref)
%DISK   Create a disk mesh.
% n - number of grid points
% nref - number of refinements (default 1)
% r - extent of inner square [-r, r]^2 (default sqrt(2)/4)

  if nargin < 2
    a = 1;
    m = 5;
    nref = 1;
  end

  r = (1 - a) * sqrt(2) / 4;

  % Default number of refinements
  if nargin < 4
    nref = 1;
  end

  % Number of boundary panels
  Nu = 4 * nref;
  Nv = nref;

  % Initialize
  x = cell(1);
  y = cell(1);
  z = cell(1);

  % Counter
  elem = 1;

  % Create boundary panels

  v = chebpts(n, [0 1 / Nv])'; %normal coordinate

  for nv = 1 : Nv

    u = chebpts(n, [0 2 * pi / Nu]); %tangential coordinate

    for nu = 1 : Nu

      % Edges
      [x1, y1] = boundary(u, a, m);
      [x2, y2] = square(u);

      % Interpolate
      x{elem} = (1 - v) .* x1 + r * v .* x2;
      y{elem} = (1 - v) .* y1 + r * v .* y2;

      % Update counter
      elem = elem + 1;

      % Update tangential coordinate
      u = u + 2 * pi / Nu;

    end
 
    % Update normal coordinate
    v = v + 1 / Nv;

  end

  % Horizontal coordinate
  v = chebpts(n, [0 2 * r / nref]) - r;

  for nv = 1 : nref

    % Vertical coordinate
    u = chebpts(n, [0 2 * r / nref]) - r;

    for nu = 1 : nref
  
      % Create grid
      a = cos(3 * pi / 4);
      b = sin(3 * pi / 4);
      x{elem} = a * v - b * u';
      y{elem} = b * v + a * u';

      % Update counter
      elem = elem + 1;

      % Update vertical coordinate
      u = u + 2 * r / nref;

    end

    % Update horizontal coordinate
    v = v + 2 * r / nref;

  end

  % Vertical coordinate (flat)
  for elem = 1 : length(x)
    z{elem} = zeros(n);
  end

  % Create mesh
  dom = surfacemesh(x, y, z);

end

% Square boundary
function [tx, ty] = square(s)

    s = mod(s, 2*pi);

    % Preallocate
    x = zeros(size(s));
    y = zeros(size(s));

    % Constants
    k = 4 / pi;

    % Quadrant masks
    q1 = s >= 0       & s < pi/2;
    q2 = s >= pi/2    & s < pi;
    q3 = s >= pi      & s < 3*pi/2;
    q4 = s >= 3*pi/2  & s < 2*pi;

    % Assign x and y for each quadrant
    x(q1) = k * s(q1) - 1;              y(q1) = -1;
    x(q2) = 1;                          y(q2) = k * (s(q2) - pi/2) - 1;
    x(q3) = 1 - k * (s(q3) - pi);       y(q3) = 1;
    x(q4) = -1;                         y(q4) = 1 - k * (s(q4) - 3*pi/2);

    % Rotate
    a = cos(3 * pi / 4);
    b = sin(3 * pi / 4);
    tx = a * x - b * y;
    ty = b * x + a * y;

end


% Boundary circle
function [x, y] = boundary(s, a, m)
  x = cos(s) .* (1 + a * cos(m * s)) / sqrt(1 + (1 / 2) * a^2);
  y = sin(s) .* (1 + a * cos(m * s)) / sqrt(1 + (1 / 2) * a^2);
end