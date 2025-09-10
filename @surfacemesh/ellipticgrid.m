function dom = ellipticgrid(xb, yb, varargin)
%ELLIPTICGRID   Elliptic grid generation of a 2D domain from boundary data.
% xb, yb: boundary data (function handle or cell array).
% varargin: either (dom) or (n) or (n, nref), where dom is a surfacemesh
%   object (currently only surfacemesh.disk supported), n is the number of
%   Chebyshev points per element, and nref is the number of uniform
%   refinements to apply to the reference domain (default: 1).
  if nargin < 3
    error('Not enough input arguments.');
  end

  % Form reference domain
  if isa(varargin{1}, 'surfacemesh')
    
    % (WARNING: currently only stable when dom = surfacemesh.disk(varargin))
    dom = varargin{1};
    n = length(dom.x{1});

  elseif isa(varargin{1}, 'numeric')

    n = varargin{1};

    if nargin < 4
      nref = 1;
    else
      nref = varargin{2};
    end

    % Create reference domain
    dom = surfacemesh.disk(n, nref);

  end

  % Build Laplace-Beltrami operator
  pdo = struct('lap', 1);
  L = surfaceop(dom, pdo);
  L.build();
  xyz = L.patches{1}.xyz;

  num_elems = length(xyz) / (n - 2);

  % Arclength parameter
  s0 = chebpts(n - 2, [0 2 * pi / num_elems], 1);
  s = zeros(num_elems * (n - 2), 1);
  for elem = 1 : num_elems
    s((elem - 1) * (n - 2) + (1 : n - 2)) = s0 + (elem - 1) * 2 * pi / num_elems;
  end

  % Get boundary conditions
  if isa(xb, 'cell')
    
    n = length(xb{1});
    B = barymat(chebpts(n - 2, 1), chebpts(n, 2));
    for elem = 1 : length(xb)
      xb{elem} = B * xb{elem};
      yb{elem} = B * yb{elem};
    end
    xb = cell2mat(xb); xb = xb(:);
    yb = cell2mat(yb); yb = yb(:);

  elseif isa(xb, 'function_handle')
    xb = xb(s);
    yb = yb(s);
  end

  if length(xb) ~= size(xyz, 1)
    error('Incompatible resolution between boundary data and reference domain.');
  end

  % Initial guess
  x = L.solve(xb);
  y = L.solve(yb);

  % Picard conditions
  picard_tol = 1e-8;
  picard_iter = 100;
  
  % Start iteration counter
  iter = 0;
  err = 10 * picard_tol;

  % Picard iteration
  while err > picard_tol && iter < picard_iter

    % Gradient
    [xu, xv, ~] = grad(x);
    [yu, yv, ~] = grad(y);

    % Coefficients of inverse metric tensor
    a11 = xv .* xv + yv .* yv;
    a12 = -2 * (xu .* xv + yu .* yv);
    a22 = xu .* xu + yu .* yu;

    % Form operator
    pdo = struct('dxx', a11, 'dxy', a12, 'dyy', a22);
    L = surfaceop(dom, pdo);
    L.build();

    % Solve
    xp1 = L.solve(xb);
    yp1 = L.solve(yb);

    % Compute error
    err = max(norm(xp1 - x, 'inf'), norm(yp1 - y, 'inf'));

    % Prepare for next iteration
    x = xp1;
    y = yp1;

    % Update iteration
    iter = iter + 1;

  end

  if err > picard_tol
    error('Warning: Grid generation failed with tolerance %1.4e. Error was %1.4e\n', picard_tol, err)
  end
  
  % Output as cell array
  x = x.vals;
  y = y.vals;

  % Vertical coordinate (flat)
  z = cell(size(x));
  for elem = 1 : length(x)
    z{elem} = zeros(n);
  end

  [x, y] = surfacemesh.inherit_connectivity(x, y, dom);
  dom = surfacemesh(x, y, z);

end
