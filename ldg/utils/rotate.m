
function [xt, yt, zt] = rotate(x, y, z, P1, P2)
% Rotation matrix that takes point P1 to P2

  % Ensure the target point lies on the unit sphere
  P2 = P2 / norm(P2);

  % Compute rotation axis (cross product)
  v = cross(P1, P2);

  % If v is zero, P is either A or -A
  if norm(v) < 1e-6

    if all(abs(P2 - P1) < 1e-6)
      R = eye(3); % No rotation needed
    else
      % Rotate to other side of sphere
      R = eye(3);
      R(3, 3) = -1;
    end

  else

    v = v / norm(v); % Normalize rotation axis

    % Compute angle between A and P
    angle = acos(dot(P1, P2));

    % Cross product matrix
    K = [   0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0  ];

    % Rodrigues' rotation formula
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;

  end

  % Apply rotation
  
  if iscell(x)

    xt = x;
    yt = y;
    zt = z;

    for k = 1 : length(x)
      xt{k} = R(1, 1) * x{k} + R(1, 2) * y{k} + R(1, 3) * z{k};
      yt{k} = R(2, 1) * x{k} + R(2, 2) * y{k} + R(2, 3) * z{k};
      zt{k} = R(3, 1) * x{k} + R(3, 2) * y{k} + R(3, 3) * z{k};
    end

  else 

    xt = R(1, 1) * x + R(1, 2) * y + R(1, 3) * z;
    yt = R(2, 1) * x + R(2, 2) * y + R(2, 3) * z;
    zt = R(3, 1) * x + R(3, 2) * y + R(3, 3) * z;

  end

end
