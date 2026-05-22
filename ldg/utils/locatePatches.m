
function Pk = locatePatches(v0, X, x, y, z)
% Locates existing patches that intersect the inclusion patch 
%   v0 - angle of inclusion
%   X - coordinate of inclusion
%   x, y, z - cell array of patches

  % Initialize
  Pk = [];

  % Distance threshold
  tol = 3 * (1 - cos(v0));

  for k = 1 : length(x) 

    % Squared distance to center
    xsq = (x{k} - X(1)).^2 + (y{k} - X(2)).^2 + (z{k} - X(3)).^2;

    % Find patches within tolerance
    if any(xsq(:) <= tol)
      Pk = [Pk k];
    end
    
  end

  % Get unique patches
  Pk = unique(Pk);

end