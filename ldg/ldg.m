% Surface Landau-de Gennes via spin connection
%
% Q = q1*(e1 e1' - e2 e2') + q2*(e1 e2' + e2 e1')
%
% Spin connection:  omega_i = e1 . (d/dx_i) e2
%
% Euler-Lagrange (Coulomb gauge, div(omega) = 0):
%
%   -L*lap(q1) + 4L*(omega.grad(q2)) + 4L*|omega|^2*q1 + a*q1 + 2c*r^2*q1 = 0
%   -L*lap(q2) - 4L*(omega.grad(q1)) + 4L*|omega|^2*q2 + a*q2 + 2c*r^2*q2 = 0
%
% where r^2 = q1^2 + q2^2.
%
% Newton (solve for q^{k+1} directly):
%   DF(q^k) * q^{k+1} = DF(q^k)*q^k - F(q^k) = 4c*r^2 * q^k
%
%   DF = [-L*lap + C(q^k),  4L*(omega.grad)]
%        [-4L*(omega.grad), -L*lap + C(q^k)]
%
%   C_11 = 4L*|omega|^2 + a + 2c*(3q1^2 + q2^2)
%   C_12 = C_21 = 4c*q1*q2
%   C_22 = 4L*|omega|^2 + a + 2c*(q1^2 + 3q2^2)
%
% Frame: e1 = (-y, x, 0)/rho,  e2 = n x e1  (azimuthal/meridional on sphere)

%% Geometry: sphere with holes
addpath('utils')
close all

p   = 8;
N_i = 5;
v0  = 0.05 * pi;

[x, y, z] = cubedSphere(p, 4);
X = @(u, v) [cos(u) .* sin(v), sin(u) .* sin(v), cos(v)];

golden = (1 + sqrt(5)) / 2;
k   = (1:N_i)';
v_i = acos(1 - 2*k/(N_i+1));            % polar:    avoids poles
u_i = 2*pi * k / golden;                 % azimuthal: golden-angle spiral
X_i = X(u_i, v_i);

for n_i = 1:N_i
    P = locatePatches(v0, X_i(n_i,:), x, y, z);
    xedge = []; yedge = []; zedge = [];
    for p = P
        [xe, ye, ze] = locateEdges(v0, X_i(n_i,:), x{p}, y{p}, z{p});
        xedge = [xedge; xe]; yedge = [yedge; ye]; zedge = [zedge; ze];
    end
    [x_i, y_i, z_i] = inclusion(v0, X_i(n_i,:), xedge, yedge, zedge);
    x(P) = []; y(P) = []; z(P) = [];
    Np = length(x);
    for m_i = 1:length(x_i)
        x{Np+m_i} = x_i{m_i}; y{Np+m_i} = y_i{m_i}; z{Np+m_i} = z_i{m_i};
    end
end

dom = surfacemesh(x, y, z);

%% Parameters
L = 1.0;
a = -1.0;
c =  2.0;

maxit = 50;
tol   = 1e-10;

%% Quadrature-point coordinates and local frame
xyz_fun = surfacefunv(surfacefun(@(x,y,z) x, dom), ...
                      surfacefun(@(x,y,z) y, dom), ...
                      surfacefun(@(x,y,z) z, dom));
XYZ = sfunv_to_mat(xyz_fun);

xq = XYZ(:,1); yq = XYZ(:,2); zq = XYZ(:,3);
rho = max(sqrt(xq.^2 + yq.^2), eps);
E1 = [-yq, xq, zeros(size(xq))] ./ rho;   % azimuthal
E2 = cross(XYZ, E1, 2);                    % n x e1 (meridional)

%% Spin connection: omega_i = e1 . (d/dx_i) e2
% Build surfacefuns for e1, e2 components
e1x = sfun_from_vec(E1(:,1), dom);
e1y = sfun_from_vec(E1(:,2), dom);
e1z = sfun_from_vec(E1(:,3), dom);
e2x = sfun_from_vec(E2(:,1), dom);
e2y = sfun_from_vec(E2(:,2), dom);
e2z = sfun_from_vec(E2(:,3), dom);

% omega_x = e1 . (d e2 / dx), etc.
ox = e1x .* diff(e2x,1,1) + e1y .* diff(e2y,1,1) + e1z .* diff(e2z,1,1);
oy = e1x .* diff(e2x,1,2) + e1y .* diff(e2y,1,2) + e1z .* diff(e2z,1,2);
oz = e1x .* diff(e2x,1,3) + e1y .* diff(e2y,1,3) + e1z .* diff(e2z,1,3);

ox_vals = ox.vec();
oy_vals = oy.vec();
oz_vals = oz.vec();
om2_vals = ox_vals.^2 + oy_vals.^2 + oz_vals.^2;   % |omega|^2

%% Initial guess: director at each point oriented toward nearest hole
s0 = 0.5;
Np = size(XYZ, 1);
q_init = zeros(Np, 2);
for i = 1:Np
    p  = XYZ(i,:);
    [~, k] = max(X_i * p');
    Xi = X_i(k,:);
    nu = (p * Xi') * p - Xi;
    nu = nu / norm(nu);
    alpha = dot(E1(i,:), nu);
    beta  = dot(E2(i,:), nu);
    q_init(i,1) = s0 * (alpha^2 - beta^2) / 2;
    q_init(i,2) = s0 * alpha * beta;
end
q = sfunv_from_mat(q_init, dom);

%% Newton iteration
for iter = 1:maxit

    q_vals = sfunv_to_mat(q);
    q1 = q_vals(:,1);
    q2 = q_vals(:,2);
    r2 = q1.^2 + q2.^2;

    % Zero-order Jacobian block
    base = 4*L*om2_vals + a;
    c11 = base + 2*c*(3*q1.^2 + q2.^2);
    c12 = 4*c*q1.*q2;
    c22 = base + 2*c*(q1.^2 + 3*q2.^2);

    DG_cell = {sfun_from_vec(c11, dom), sfun_from_vec(c12, dom); ...
               sfun_from_vec(c12, dom), sfun_from_vec(c22, dom)};

    % First-derivative coupling: +4L*(omega.grad q2) in q1-eq,
    %                            -4L*(omega.grad q1) in q2-eq
    z0  = sfun_from_vec(zeros(Np, 1), dom);
    fox = sfun_from_vec( 4*L*ox_vals, dom);
    foy = sfun_from_vec( 4*L*oy_vals, dom);
    foz = sfun_from_vec( 4*L*oz_vals, dom);
    nox = sfun_from_vec(-4*L*ox_vals, dom);
    noy = sfun_from_vec(-4*L*oy_vals, dom);
    noz = sfun_from_vec(-4*L*oz_vals, dom);

    pdo     = [];
    pdo.lap = -L * eye(2);
    pdo.c   = DG_cell;
    pdo.dx  = {z0, fox;  nox, z0};
    pdo.dy  = {z0, foy;  noy, z0};
    pdo.dz  = {z0, foz;  noz, z0};

    % Newton RHS: DF(q^k)*q^k - F(q^k) = 4c*r^2 * q^k
    rhs_vals = 4*c*r2 .* q_vals;
    rhs_fun  = sfunv_from_mat(rhs_vals, dom);

    L_op = surfaceop(dom, pdo, rhs_fun);
    qnew = L_op.solve(@(x,y,z) bc_q(x, y, z, X_i, 1));

    dq  = qnew - q;
    res = max(cellfun(@(c) norm(c(:)), dq.components));
    fprintf('Newton iter %2d:  ||dq|| = %.3e\n', iter, res);

    q = qnew;

    if res < tol
        fprintf('Converged.\n');
        break
    end

end

if iter == maxit && res > tol
    warning('Newton did not converge after %d iterations.', maxit);
end

%% Scalar order parameter r = |q|
q_vals  = sfunv_to_mat(q);
r_fun   = sfun_from_vec(sqrt(sum(q_vals.^2, 2)), dom);

subplot(1,2,1)
plot(r_fun)
title('r = |q| (scalar order parameter)')
colorbar

%% Principal director field
subplot(1,2,2)
theta    = 0.5 * atan2(q_vals(:,2), q_vals(:,1));
dir_vals = cos(theta).*E1 + sin(theta).*E2;

m  = 6;
n  = order(dom) + 1;
B  = bary(linspace(-1,1,m).', eye(n));
np = length(dom);
XYZ_sub = zeros(m^2*np, 3);
dir_sub = zeros(m^2*np, 3);
XYZ_res = reshape(XYZ,     n, n, np, 3);
dir_res = reshape(dir_vals, n, n, np, 3);
for k = 1:np
    idx = (k-1)*m^2 + (1:m^2);
    for d = 1:3
        XYZ_sub(idx,d) = reshape(B * XYZ_res(:,:,k,d) * B.', [], 1);
        dir_sub(idx,d) = reshape(B * dir_res(:,:,k,d)  * B.', [], 1);
    end
end
dir_sub = dir_sub ./ max(sqrt(sum(dir_sub.^2, 2)), eps);

scl = 0.04;
x1 = XYZ_sub(:,1) - scl*dir_sub(:,1);  x2 = XYZ_sub(:,1) + scl*dir_sub(:,1);
y1 = XYZ_sub(:,2) - scl*dir_sub(:,2);  y2 = XYZ_sub(:,2) + scl*dir_sub(:,2);
z1 = XYZ_sub(:,3) - scl*dir_sub(:,3);  z2 = XYZ_sub(:,3) + scl*dir_sub(:,3);
nans = nan(size(x1));
Xp = [x1 x2 nans].';  Yp = [y1 y2 nans].';  Zp = [z1 z2 nans].';

plot3(Xp(:), Yp(:), Zp(:), 'k-', 'LineWidth', 0.5)
axis equal, grid on, view(3)
title('Principal director')

%% Helper functions

function Qmat = sfunv_to_mat(Q)
ncomp = numel(Q.components);
Qmat  = zeros(numel(Q.components{1}), ncomp);
for k = 1:ncomp
    Qmat(:,k) = Q.components{k}.vec();
end
end

function f = sfun_from_vec(v, dom)
n  = order(dom) + 1;
np = length(dom);
vals = cell(np, 1);
for k = 1:np
    vals{k} = reshape(v((k-1)*n^2 + (1:n^2)), n, n);
end
f = surfacefun(vals, dom);
end

function Q = sfunv_from_mat(Qmat, dom)
ncomp = size(Qmat, 2);
comps = cell(1, ncomp);
for k = 1:ncomp
    comps{k} = sfun_from_vec(Qmat(:,k), dom);
end
Q = surfacefunv(comps{:});
end

function qbc = bc_q(x, y, z, X_i, s)
% q1_bc = (s/2)*(alpha^2 - beta^2),  q2_bc = s*alpha*beta
% alpha = e1.nu,  beta = e2.nu,  nu = outward geodesic normal at hole boundary

x = x(:); y = y(:); z = z(:);
pts = [x y z];
N   = size(pts, 1);
qbc = zeros(N, 2);

rho = max(sqrt(x.^2 + y.^2), eps);
e1  = [-y, x, zeros(N,1)] ./ rho;
e2  = cross(pts, e1, 2);

for i = 1:N
    p = pts(i,:);
    [~, k] = max(X_i * p');
    Xi = X_i(k,:);
    nu = (p * Xi') * p - Xi;
    nu = nu / norm(nu);

    alpha = dot(e1(i,:), nu);
    beta  = dot(e2(i,:), nu);

    qbc(i,1) = (s/2) * (alpha^2 - beta^2);
    qbc(i,2) = s * alpha * beta;
end
end