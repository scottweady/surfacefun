% Fitzhugh-Nagumo on a surface (vector surfaceop solver)
%
% u_t = delta_u lap(u) + (1/alpha) u(1-u)(u-(v+b)/a)
% v_t = delta_v lap(v) + u - v

n = 16;
dom = surfacemesh.cyclide(n, 32, 16);

a = 0.75;
b = 0.02;
alpha = 0.02;
delta_u = 0.03;
delta_v = 0;
dt = 0.02;

Nu = @(u,v) 1/alpha * u.*(1-u).*(u-(v+b)/a);
Nv = @(u,v) u - v;

% IMEX: -dt*D*lap(U) + I*U = U^n + dt*N(U^n), D = diag([delta_u, delta_v])
pdo = [];
pdo.lap = diag([-dt*delta_u, -dt*delta_v]);
pdo.c   = eye(2);

L = surfaceop(dom, pdo);
build(L)

%% Initial conditions
u = surfacefun(@(x,y,z) 0.5*(1+tanh(2*x+y)), dom);
v = surfacefun(@(x,y,z) 0.5*(1-tanh(3*z)), dom);

%% Simulation
close all

doplot = @(u) chain(@()plot(u), @()view(60,40), @()material('dull'), ...
    @()lighting('gouraud'), @()camlight('headlight'), @()colorbar, @()caxis('auto'));

doplot(u), shg
t = 0;
for k = 1:10000
    L.rhs = surfacefunv(u + dt*Nu(u,v), v + dt*Nv(u,v));
    UV = solve(L);
    u = UV.components{1};
    v = UV.components{2};
    t = t + dt;
    if ( mod(k,10) == 0 )
        k
        clf, doplot(u), drawnow, shg
    end
end
