clear all; clc
run('tu.m')

%指定两个边界条件
r = @(x,y) sqrt((x+1)^2 + (y+1)^2);
th = @(x,y) atan2((y+1),(x+1));
Tx = 1e4;
R = 0.5;

xig_rr = @(x,y) Tx/2*(1-R^2/r(x,y)^2)+Tx/2*(1-4*R^2/r(x,y)^2+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_thth = @(x,y) Tx/2*(1+R^2/r(x,y)^2)-Tx/2*(1+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_rth = @(x,y) -Tx/2*(1+2*R^2/r(x,y)^2-3*R^4/r(x,y)^4)*sin(2*th(x,y));

xig_xx = @(x,y) xig_rr(x,y)*cos(-th(x,y))^2+xig_thth(x,y)*sin(-th(x,y))^2+2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_yy = @(x,y) xig_rr(x,y)*sin(-th(x,y))^2+xig_thth(x,y)*cos(-th(x,y))^2-2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_xy = @(x,y) -xig_rr(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_thth(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_rth(x,y)*(cos(-th(x,y))^2-sin(-th(x,y))^2);



E = 1E9;
v = 0.3;
% conductivity
lambda = v * E / (1 + v) / (1 - 2 * v);
miu = E / 2 / (1 + v);

D = zeros(3);

%plane strain
% D(1,1) = lambda + 2*miu;
% D(2,2) = D(1,1);
% D(1,2) = lambda;
% D(2,1) = D(1,2);
% D(3,3) = miu;

%plane stress
D(1,1) = E/(1-v^2);
D(2,2) = D(1,1);
D(1,2) = v*E/(1-v^2);
D(2,1) = D(1,2);
D(3,3) = E/2/(1+v);

%用户定义
DI1 = -1;%右边边界是迪利克雷为0，纽曼为-1
DI2 = -1;%上边边界是迪利克雷为0，纽曼为-1

pos = msh.POS;%点坐标
n_el = size(msh.QUADS,1);%元素数
n_en = 4;%局部节点数
n_np = msh.nbNod;%节点数
lines = msh.LINES;%线


x_coor = pos(:, 1);
y_coor = pos(:, 2);



%IEN array
IEN_tri = zeros(1,1);
IEN = msh.QUADS(:,1:4);
for ee = 1:size(IEN,1)
    IEN_tri(ee*2-1,1) = IEN(ee,1);
    IEN_tri(ee*2-1,2) = IEN(ee,2);
    IEN_tri(ee*2-1,3) = IEN(ee,3);
    IEN_tri(ee*2,1) = IEN(ee,1);
    IEN_tri(ee*2,2) = IEN(ee,3);
    IEN_tri(ee*2,3) = IEN(ee,4);
end

for i = 1 : n_el/2
    a = IEN(i, 1);
    b = IEN(i, 2);
    IEN(i, 1) = IEN(i, 4);
    IEN(i, 2) = IEN(i, 3);
    IEN(i, 4) = a;
    IEN(i, 3) = b;
end

ID = -1 * ones(n_np,2);
for nn = 1 : size(lines,1)
    if lines(nn,3) == 8
        ID(lines(nn,1), 1) = DI1;   %右边边的x方向是否是DI~条件
        ID(lines(nn,1), 2) = DI1;   %右边边的y方向是否是DI~条件
    elseif lines(nn,3) == 9
        ID(lines(nn,1), 1) = DI2;   %上边边的x方向是否是DI~条件
        ID(lines(nn,1), 2) = DI2;   %上边边的y方向是否是DI~条件
    elseif lines(nn,3) == 10        %左边对称线
        ID(lines(nn,1), 1) = 0;
    elseif lines(nn,3) == 11        %下边对称线
        ID(lines(nn,1), 2) = 0;
    end
end

ID(1,1) = DI1;  ID(1,2) = 0;        %右下
ID(2,1) = DI1;  ID(2,2) = DI2;      %右上
ID(3,1) = 0;    ID(3,2) = DI2;      %左上
ID(4,1) = DI2;  ID(4,2) = 0;        %左下下
ID(5,1) = 0;    ID(5,2) = DI1;      %左下上

counter = 1;
for nn = 1 : n_np
    if ID(nn, 1) ~= 0
        ID(nn, 1) = counter;
        counter = counter + 1;
    end
    if ID(nn, 2) ~= 0
        ID(nn, 2) = counter;
        counter = counter + 1;
    end
end


n_eq = counter - 1;

n_vector = lines; %向量
for nn = 1 : size(n_vector)
    if n_vector(nn, 3) == 13
        n_vector(nn, 1) = -(pos(lines(nn, 1), 1) + pos(lines(nn, 2), 1)) / 2 - 1;
        n_vector(nn, 2) = -(pos(lines(nn, 1), 2) + pos(lines(nn, 2), 2)) / 2 - 1;
        norm = sqrt(n_vector(nn, 1)^2+n_vector(nn, 2)^2);
        n_vector(nn, 1) = n_vector(nn, 1) / norm;
        n_vector(nn, 2) = n_vector(nn, 2) / norm;
    elseif n_vector(nn, 3) == 11
        n_vector(nn, 1) = 0;
        n_vector(nn, 2) = -1;
    elseif n_vector(nn, 3) == 10
        n_vector(nn, 1) = -1;
        n_vector(nn, 2) = 0;
    elseif n_vector(nn, 3) == 9
        n_vector(nn, 1) = 0;
        n_vector(nn, 2) = 1;
    elseif n_vector(nn, 3) == 8
        n_vector(nn, 1) = 1;
        n_vector(nn, 2) = 0;
    end
end

vector = n_vector(:, 1:2);



% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int2D     = n_int_xi * n_int_eta;
[xi2D, eta2D, weight2D] = Gauss2D(n_int_xi, n_int_eta);

n_int1D = 5;
[xi1D, weight1D] = Gauss(n_int1D, -1, 1);



% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
    x_ele = [x_coor( IEN(ee, 1:n_en) );x_coor( IEN(ee, 1) )];
    y_ele = [y_coor( IEN(ee, 1:n_en) );y_coor( IEN(ee, 1) )];

    k_ele = zeros(2*n_en, 2*n_en); % element stiffness matrix
    f_ele = zeros(2*n_en, 1);    % element load vector
    % quadrature loop
    for qua = 1 : n_int1D       % map阶数
        for aa = 1 : n_en
            pp = 2*(aa-1);
            if x_ele(aa) == 1 && x_ele(aa+1) == 1
                y_l    = y_ele(aa) * PolyShape(1, 1, xi1D(qua), 0) + y_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 0);
                dy_dxi = y_ele(aa) * PolyShape(1, 1, xi1D(qua), 1) + y_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 1); % Σxae Na,x (ξ)用高斯积分来积dx/dξ

                [h_x, h_y] = hh(1,y_l,1,0);
                f_ele(pp+1) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_x * dy_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_x * dy_dxi;
                f_ele(pp+2) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_y * dy_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_y * dy_dxi;
            end

            if x_ele(aa) == -1 && x_ele(aa+1) == -1
                y_l    = y_ele(aa) * PolyShape(1, 1, xi1D(qua), 0) + y_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 0);
                dy_dxi = y_ele(aa) * PolyShape(1, 1, xi1D(qua), 1) + y_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 1); % Σxae Na,x (ξ)用高斯积分来积dx/dξ

                [h_x, h_y] = hh(-1,y_l,-1,0);
                f_ele(pp+1) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_x * dy_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_x * dy_dxi;
                f_ele(pp+2) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_y * dy_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_y * dy_dxi;
            end

            if y_ele(aa) == 1 && y_ele(aa+1) == 1
                x_l    = x_ele(aa) * PolyShape(1, 1, xi1D(qua), 0) + x_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 0);
                dx_dxi = x_ele(aa) * PolyShape(1, 1, xi1D(qua), 1) + x_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 1); % Σxae Na,x (ξ)用高斯积分来积dx/dξ

                [h_x, h_y] = hh(x_l,1,0,1);
                f_ele(pp+1) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_x * dx_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_x * dx_dxi;
                f_ele(pp+2) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_y * dx_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_y * dx_dxi;
            end

            if y_ele(aa) == -1 && y_ele(aa+1) == -1
                x_l    = x_ele(aa) * PolyShape(1, 1, xi1D(qua), 0) + x_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 0);
                dx_dxi = x_ele(aa) * PolyShape(1, 1, xi1D(qua), 1) + x_ele(aa+1) * PolyShape(1, 2, xi1D(qua), 1); % Σxae Na,x (ξ)用高斯积分来积dx/dξ

                [h_x, h_y] = hh(x_l,-1,0,-1);
                f_ele(pp+1) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_x * dx_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_x * dx_dxi;
                f_ele(pp+2) = weight1D(qua) * PolyShape(1, 1, xi1D(qua), 0) * h_y * dx_dxi + weight1D(qua) * PolyShape(1, 2, xi1D(qua), 0) * h_y * dx_dxi;
            end
        end
    end
    for ll = 1 : n_int2D
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi2D(ll), eta2D(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi2D(ll), eta2D(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi2D(ll), eta2D(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1 : n_en
            Na = Quad(aa, xi2D(ll), eta2D(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi2D(ll), eta2D(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            B_a = zeros(3,2);
            B_a(1,1) = Na_x;
            B_a(2,2) = Na_y;
            B_a(3,1) = Na_y;
            B_a(3,2) = Na_x;


            for bb = 1 : n_en
                Nb = Quad(bb, xi2D(ll), eta2D(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi2D(ll), eta2D(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                B_b = zeros(3,2);
                B_b(1,1) = Nb_x;
                B_b(2,2) = Nb_y;
                B_b(3,1) = Nb_y;
                B_b(3,2) = Nb_x;

                for i = 1 : 2
                    for j = 1 : 2
                        p = 2*(aa-1)+i;
                        q = 2*(bb-1)+j;
                        k_ele(p, q) = k_ele(p, q) + weight2D(ll) * detJ...
                            * unit_vector(i)' * B_a' * D * B_b * unit_vector(j);
                    end
                end
            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for aa = 1 : n_en
        for i = 1 : 2
            PP = ID(IEN(ee, aa) ,i);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*(aa-1)+i);
                for bb = 1 : n_en
                    for j = 1 : 2
                        QQ = ID(IEN(ee, bb), j);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+i, 2*(bb-1)+j);

                            % else
                            %     % modify F with the boundary data
                            %     % here we do nothing because the boundary data g is zero or
                            %     % homogeneous
                        end
                    end
                end
            end
        end
    end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);

for ii = 1 : n_np
    for k = 1 : 2
        index = ID(ii, k);
        if index > 0
            disp(ii, k) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end
end





figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, disp(:,1));
axis equal;
colormap jet
shading interp
title('x - direction displacement (d_x)');
xlabel('x - coordinate');
ylabel('y - coordinate');

figure;
hold on;
trisurf(IEN_tri, x_coor, y_coor, disp(:,2));
axis equal;
colormap jet
shading interp
title('y - direction displacement (d_y)');
xlabel('x - coordinate');
ylabel('y - coordinate');




strain = zeros(n_el, 3);
for ee = 1:n_el
    x_ele = [x_coor( IEN(ee, 1:n_en) );x_coor( IEN(ee, 1) )];
    y_ele = [y_coor( IEN(ee, 1:n_en) );y_coor( IEN(ee, 1) )];
    for ll = 1:n_int2D
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1:n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi2D(ll), eta2D(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi2D(ll), eta2D(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi2D(ll), eta2D(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        epsilon = zeros(3, 1);
        for aa = 1:n_en
            Na = Quad(aa, xi2D(ll), eta2D(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi2D(ll), eta2D(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
            B_a = zeros(3, 2);
            B_a(1, 1) = Na_x;
            B_a(2, 2) = Na_y;
            B_a(3, 1) = Na_y;
            B_a(3, 2) = Na_x;
            for i = 1:2
                p = 2 * (aa - 1)+i;
                node_disp = disp(IEN(ee, aa), i);
                epsilon = epsilon + B_a(:, i) * node_disp;
            end
        end
        % 这里简单取高斯积分点的平均应变作为单元应变
        for qq = 1:3
            strain(ee, qq) = strain(ee, qq) + weight2D(ll) * epsilon(qq);
        end
    end
    strain(ee, :) = strain(ee, :) / n_int2D;
end

stress = zeros(n_el, 3);
for ee = 1:n_el
    stress(ee, :) = D * strain(ee, :)';
end

% 初始化节点应变数组，每个节点有3个应变分量（xx, yy, xy）
node_strain = zeros(n_np, 3);
% 计算每个单元的面积（假设单元为四边形，这里简单用叉积法估算）
element_area = zeros(n_el, 1);
for ee = 1:n_el
    x1 = x_coor(IEN(ee, 1));
    y1 = y_coor(IEN(ee, 1));
    x2 = x_coor(IEN(ee, 2));
    y2 = y_coor(IEN(ee, 2));
    x3 = x_coor(IEN(ee, 3));
    y3 = y_coor(IEN(ee, 3));
    x4 = x_coor(IEN(ee, 4));
    y4 = y_coor(IEN(ee, 4));
    area1 = 0.5 * abs((x1 - x3) * (y2 - y4)-(x2 - x4) * (y1 - y3));
    element_area(ee) = area1;
end
for ee = 1:n_el
    for aa = 1:n_en
        node_index = IEN(ee, aa);
        % 用单元面积加权
        node_strain(node_index, :) = node_strain(node_index, :) + element_area(ee) * strain(ee, :);
    end
end
% 对每个节点的应变进行平均
for nn = 1:n_np
    total_area = 0;
    for ee = 1:n_el
        if any(IEN(ee, :) == nn)
            total_area = total_area + element_area(ee);
        end
    end
    if total_area > 0
        node_strain(nn, :) = node_strain(nn, :) / total_area;
    end
end

% 初始化节点应力数组，每个节点有3个应力分量（xx, yy, xy）
node_stress = zeros(n_np, 3);
for nn = 1:n_np
    node_stress(nn, :) = D * node_strain(nn, :)';
end




figure;                     %x方向应力
hold on
trisurf(IEN_tri, x_coor, y_coor, node_stress(:, 1));
axis equal;
colormap jet
shading interp
title('x - direction stress (\sigma_{xx})');
xlabel('x - coordinate');
ylabel('y - coordinate');

figure;                     %y方向应力
hold on
trisurf(IEN_tri, x_coor, y_coor, node_stress(:, 2));
axis equal;
colormap jet
shading interp
title('y - direction stress (\sigma_{yy})');
xlabel('x - coordinate');
ylabel('y - coordinate');

figure;                     %xy扭矩
hold on
trisurf(IEN_tri, x_coor, y_coor, node_stress(:, 3));
axis equal;
colormap jet
shading interp
title('Shear stress (\sigma_{xy})');
xlabel('x - coordinate');
ylabel('y - coordinate');

%%%%%%%%%%%%%

figure;                     %x方向应变
hold on
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 1));
axis equal;
colormap jet
shading interp
title('x - direction strain (\sigma_{xx})');
xlabel('x - coordinate');
ylabel('y - coordinate');

figure;                     %y方向应变
hold on
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 2));
axis equal;
colormap jet
shading interp
title('y - direction strain (\sigma_{yy})');
xlabel('x - coordinate');
ylabel('y - coordinate');

figure;                     %xy扭转
hold on
trisurf(IEN_tri, x_coor, y_coor, node_strain(:, 3));
axis equal;
colormap jet
shading interp
title('Shear strain (\sigma_{xy})');
xlabel('x - coordinate');
ylabel('y - coordinate');
















function [h_x, h_y] = hh (x, y, n_x, n_y)

r = @(x,y) sqrt((x+1)^2 + (y+1)^2);
th = @(x,y) atan2((y+1),(x+1));
Tx = 1e4;
R = 0.5;

xig_rr = @(x,y) Tx/2*(1-R^2/r(x,y)^2)+Tx/2*(1-4*R^2/r(x,y)^2+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_thth = @(x,y) Tx/2*(1+R^2/r(x,y)^2)-Tx/2*(1+3*R^4/r(x,y)^4)*cos(2*th(x,y));
xig_rth = @(x,y) -Tx/2*(1+2*R^2/r(x,y)^2-3*R^4/r(x,y)^4)*sin(2*th(x,y));

xig_xx = @(x,y) xig_rr(x,y)*cos(-th(x,y))^2+xig_thth(x,y)*sin(-th(x,y))^2+2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_yy = @(x,y) xig_rr(x,y)*sin(-th(x,y))^2+xig_thth(x,y)*cos(-th(x,y))^2-2*xig_rth(x,y)*sin(-th(x,y))*cos(-th(x,y));
xig_xy = @(x,y) -xig_rr(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_thth(x,y)*sin(-th(x,y))*cos(-th(x,y))+xig_rth(x,y)*(cos(-th(x,y))^2-sin(-th(x,y))^2);

xig_ij = zeros(2,2);
h = zeros(2,1);

vector = [n_x, n_y]';
xig_ij(1,1) = xig_xx(x,y);
xig_ij(1,2) = xig_xy(x,y);
xig_ij(2,1) = xig_xy(x,y);
xig_ij(2,2) = xig_yy(x,y);
for i = 1 : 2
    for j = 1 : 2
        h(i) = h(i) + xig_ij(i,j) * vector(j);
    end
end
h_x = h(1);
h_y = h(2);
end
function e_i = unit_vector(ii)
if ii == 1
    e_i = [1,0]';
elseif ii == 2
    e_i = [0,1]';
end
end


%EOF