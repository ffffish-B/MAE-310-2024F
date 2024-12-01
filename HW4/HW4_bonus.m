clear; clc; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 2;              % polynomial degree 拟合阶数（map的段数）
n_en = pp + 1;         % number of element or local nodes（map的节点数）
n_el = 10000;              % number of elements 单元数
n_np = n_el * pp + 1;  % number of nodal points 节点数
n_eq = n_np - 1;       % number of equations P Q
n_int = 10;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes 取等长单元h
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

IEN = zeros(n_el, n_en);

for ee = 1 : n_el      % 单元数
  for aa = 1 : n_en    % local节点数
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el   % 单元数
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix 建立Kab单元刚度阵
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee) % 局部节点X（ξ）

  % quadrature loop
  for qua = 1 : n_int       % map阶数
    dx_dxi = 0.0;           % dx/dξ
    x_l = 0.0;              % x(ξl)
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1); % Σxae Na,x (ξ)用高斯积分来积dx/dξ
    end 
    dxi_dx = 1.0 / dx_dxi;

% 求局部系数矩阵Kab

    for aa = 1 : n_en
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % Assembly of the matrix and vector based on the ID or LM data
  % check the ID(IEN(ee, aa)) and ID(IEN(ee, bb), if they are positive
  % put the element stiffness matrix into K
  % 把Kab往K里带

  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if(P > 0)
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;
tic
x = gmres(K, F, 10000, 1e-2, 10000);

x = [x; g];
toc

tic
d_temp = K \ F;

disp = [d_temp; g];
toc

delta = disp-x;
















