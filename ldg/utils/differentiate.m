
function [ux, uy, uxx, uxy, uyy] = differentiate(u)
% Compute derivatives via Chebyshev differentiation

  % Number of grid points
  N = size(u, 1);

  % Chebyshev differentiation matrix
  D = diffmat(N);

  % Preallocate
  ux = zeros(size(u));
  uy = zeros(size(u));
  uxx = zeros(size(u));
  uxy = zeros(size(u));
  uyy = zeros(size(u));

  % First derivatives
  for i = 1 : 3
    ux(:, :, i) = D * u(:, :, i);
    uy(:, :, i) = u(:, :, i) * D';
  end

  % Second derivatives
  for i = 1 : 3
    uxx(:, :, i) = D * ux(:, :, i);
    uxy(:, :, i) = 0.5 * (D * uy(:, :, i) + ux(:, :, i) * D');
    uyy(:, :, i) = uy(:, :, i) * D';
  end

end