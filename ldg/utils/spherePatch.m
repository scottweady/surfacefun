
function [x, y, z] = spherePatch(xb, yb, zb)
% Constructs a mesh on the unit sphere with boundary xb, yb, zb

  % Number of grid points
  N = length(xb);

  % Vector size
  shp = [N N 3];
  
  % GMRES conditions
  gmres_tol = 1e-6;
  gmres_iter = 100;
  
  % Newton conditions
  newton_tol = 1e-6;
  newton_iter = 50;

  % Zero interior values
  b = cat(3, xb, yb, zb);
  b(2 : end - 1, 2 : end - 1, :) = 0;

  % Reshape
  b = b(:);

  % Compute Laplace precondition
  P = laplacePreconditioner(shp);

  % Initial guess for solution vector
  X = P(b);

  % PDE operator
  F = @(X) winslow(X, shp) - b;
  
  % Jacobian vector product
  JX = @(X, dX) jacobian(X, dX, shp);
  
  % Compute initial error
  err = max(abs(F(X)));
  
  % Start iteration counter
  iter = 0;
  
  % Newton-Krylov solver
  while err > newton_tol && iter < newton_iter
  
    % Solve for Newton update
    [dX, ~] = gmres(@(dX) JX(X, dX), F(X), [], gmres_tol, gmres_iter, P);

    % Update solution
    X = X - dX;

    % Compute error
    err = max(abs(F(X)));

    % Update iteration
    iter = iter + 1;
    
  end
  
  % Reformat
  X = reshape(X, shp);
  x = X(:, :, 1);
  y = X(:, :, 2);
  z = X(:, :, 3);

  % Ensure exactly on sphere
  xnorm = sqrt(x.^2 + y.^2 + z.^2);
  x = x ./ xnorm;
  y = y ./ xnorm;
  z = z ./ xnorm;

  if err > newton_tol
    x = []; y = []; z = [];
    fprintf('Grid generation failed.\n')
    return
  end

end

% Winslow PDE
function pde = winslow(X, shp)

  % Reshape
  X = reshape(X, shp);

  % Compute derivatives of solution vector
  [Xu, Xv, Xuu, Xuv, Xvv] = differentiate(X);
  
  % Compute metric tensor
  g11 = sum(Xu.^2, 3);
  g12 = sum(Xu .* Xv, 3);
  g22 = sum(Xv.^2, 3);
  detg = g11 .* g22 - g12.^2;

  % Compute coefficients
  a11 =  g22;
  a12 = -g12;
  a22 =  g11;

  % Evaluate PDE
  pde = a11 .* Xuu + 2 * a12 .* Xuv + a22 .* Xvv + 2 * detg .* X;

  % Apply Dirichlet boundary conditions
  pde([1 end], :, :) = X([1 end], :, :);
  pde(:, [1 end], :) = X(:, [1 end], :);

  % Reshape
  pde = pde(:);

end

% Jacobian of Winslow PDE
function dpde = jacobian(X, dX, shp)

  % Initialize
  dpde = zeros(shp);

  % Current solution
  X = reshape(X, shp);

  % Perturbation
  dX = reshape(dX, shp);

  % Compute derivatives
  [Xu, Xv, Xuu, Xuv, Xvv] = differentiate(X);
  [dXu, dXv, dXuu, dXuv, dXvv] = differentiate(dX);

  % Compute metric tensor
  g11 = sum(Xu.^2, 3);
  g12 = sum(Xu .* Xv, 3);
  g22 = sum(Xv.^2, 3);

  detg = g11 .* g22 - g12.^2;

  % Functional derivative of metric tensor
  dg11 = 2 * sum(Xu .* dXu, 3);
  dg12 = sum(Xu .* dXv + dXu .* Xv, 3);
  dg22 = 2 * sum(Xv .* dXv, 3);

  ddetg = dg11 .* g22 + g11 .* dg22 - 2 * g12 .* dg12;

  % Inverse metric tensor
  a11 =  g22;
  a12 = -g12;
  a22 =  g11;

  % Coefficients of perturbed PDE
  da11 = dg22;
  da12 = -dg12;
  da22 = dg11;

  % Derivative of PDE
  dpde(:, :, 1 : 3) =  a11 .* dXuu + 2 * a12 .* dXuv + a22 .* dXvv + ...
                       da11 .* Xuu + 2 * da12 .* Xuv + da22 .* Xvv + ...
                       2 * (detg .* dX + ddetg .* X);
  % Dirichlet boundary conditions
  dpde([1 end], :, 1 : 3) = dX([1 end], :, :);
  dpde(:, [1 end], 1 : 3) = dX(:, [1 end], :);

  % Reshape
  dpde = dpde(:);

end

function P = laplacePreconditioner(shp)

  % Initialize
  dof = prod(shp);
  Delta = zeros(dof);

  % Apply to basis vectors
  for m = 1 : dof
    e = zeros(dof, 1);
    e(m) = 1;
    Delta(:, m) = laplace(e, shp);
  end

  % Factor
  [L, U] = lu(Delta);

  % Assign function handle
  P = @(x) U \ (L \ x);

end


% Laplace operator
function pde = laplace(u, shp)

  % Reshape
  u = reshape(u, shp);
  [~, ~, uxx, ~, uyy] = differentiate(u);

  % Evaluate PDE
  pde = uxx + uyy + 2 * u;
  
  % Apply boundary conditions
  pde([1 end], :, :) = u([1 end], :, :);
  pde(:, [1 end], :) = u(:, [1 end], :);

  % Reshape
  pde = pde(:);

end
