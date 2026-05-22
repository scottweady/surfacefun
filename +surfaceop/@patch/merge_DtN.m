function c = merge_DtN(a, b, rankdef)
%MERGE   Merge two patch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common edges will be
%   eliminated by enforcing continuity and continuity of the derivative
%   across the boundary.

% Parse inputs:
if ( nargin == 0 )
    c = [];
    return
elseif ( nargin == 1 )
    c = a;
    return
elseif ( nargin == 2 )
    rankdef = false;
end

% Compute the indices of intersecting points in a and b.
[i1, i2, s1, s2, flip1, flip2, scl1, scl2, D2N_scl, dom, edges] = intersect(a, b);
pi1 = i1; pi2 = i2;
ps1 = s1;
ncomp = size(a.BtB, 1)/size(a.xyz, 1);
parentIdx1 = [];
parentIdx2 = [];
if ( ncomp > 1 )
    [parentOrder, parentIdx1, parentIdx2] = boundaryOrder(ncomp, numel(i1), numel(i2));
    i1 = blockIdx(i1, ncomp, size(a.xyz, 1));
    s1 = blockIdx(s1, ncomp, size(a.xyz, 1));
    i2 = blockIdx(i2, ncomp, size(b.xyz, 1));
    s2 = blockIdx(s2, ncomp, size(b.xyz, 1));
    flip1 = kron(eye(ncomp), flip1);
    flip2 = kron(eye(ncomp), flip2);
    scl1 = repmat(scl1, ncomp, 1);
    scl2 = repmat(scl2, ncomp, 1);
end

% Extract D2N maps:
D2Na = a.BtB; D2Nb = b.BtB;

% Compute new solution operator:
% - The Dirichlet-to-Neumann maps on singular elements have their Jacobians
%   factored out. Therefore, we need to multiply the continuity conditions
%   by the multiplication matrices SCL1 and SCL2. The coordinate maps are
%   such that multiplying D2NA/B by SCL1/2 cancels out, so D2NA/B should
%   only be multiplied by SCL2/1, respectively.
A = -( scl2.*flip1*D2Na(s1,s1)*flip1.' + scl1.*flip2*D2Nb(s2,s2)*flip2.' );
z = [ scl2.*flip1*D2Na(s1,i1) scl1.*flip2*D2Nb(s2,i2) ];
z_part = scl2.*flip1*a.du_part(s1,:) + scl1.*flip2*b.du_part(s2,:);

% Check for a closed surface at the top level:
if ( rankdef && isempty(i1) && isempty(i2) )
    % Fix rank deficiency with Leslie's ones matrix trick:
    w = a.w(ps1);
    if ( ncomp > 1 )
        A = A + kron(eye(ncomp), w*w');
    else
        A = A + w*w'; % or is it sqrt(w)*sqrt(w)' ?
    end
end

% Store the decomposition for reuse in updateRHS():
dA = decomposition(A, 'CheckCondition', false);
S = dA \ z;
u_part = dA \ z_part;

if ( isIllConditioned(dA) )
    warning(['Schur complement linear system is nearly singular. ', ...
        'Did you forget to set rankdef = true?']);
end

% Compute new D2N maps:
M = [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ];
D2N = M*S;
b1 = 1:numel(i1);
b2 = numel(i1)+(1:numel(i2));
D2N(b1,b1) = D2N(b1,b1) + D2Na(i1,i1);
D2N(b2,b2) = D2N(b2,b2) + D2Nb(i2,i2);
du_part = [ a.du_part(i1,:) ; b.du_part(i2,:) ] + M * u_part;

if ( ncomp > 1 )
    S = S(:,parentOrder);
    D2N = D2N(parentOrder,parentOrder);
    du_part = du_part(parentOrder,:);
else
    parentIdx1 = 1:numel(i1);
    parentIdx2 = numel(i1)+(1:numel(i2));
end

% Construct the new patch:
xyz = [a.xyz(pi1,:) ; b.xyz(pi2,:)];
w = [a.w(pi1) ; b.w(pi2)];
id = [a.id ; b.id];
c = surfaceop.parent(dom, id, S, D2N, D2N_scl, u_part, du_part, A, dA, ...
    edges, xyz, w, a, b, {i1, s1, parentIdx1}, {i2, s2, parentIdx2}, flip1, flip2, scl1, scl2);

end

function [p, p1, p2] = boundaryOrder(ncomp, n1, n2)

p = zeros(ncomp*(n1+n2), 1);
p1 = zeros(ncomp*n1, 1);
p2 = zeros(ncomp*n2, 1);
pos = 0;
q1 = 0;
q2 = 0;
for k = 1:ncomp
    old1 = (k-1)*n1 + (1:n1);
    old2 = ncomp*n1 + (k-1)*n2 + (1:n2);
    rows1 = pos + (1:n1);
    rows2 = pos + n1 + (1:n2);

    p(rows1) = old1;
    p(rows2) = old2;
    p1(q1+(1:n1)) = rows1;
    p2(q2+(1:n2)) = rows2;

    pos = pos + n1 + n2;
    q1 = q1 + n1;
    q2 = q2 + n2;
end

end

function idx = blockIdx(idx, ncomp, nbase)

idx = idx(:);
idx = idx + (0:ncomp-1)*nbase;
idx = idx(:);

end
