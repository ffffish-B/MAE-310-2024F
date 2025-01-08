clear all; clc
run('tu.m')

% quadrature rule
n_int_xi  = 10;
n_int_eta = 10;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);


%指定两个边界条件




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

pos = msh.POS;
n_en = 4;
n_np = msh.nbNod;
lines = msh.LINES;

x_coor = pos(:, 1);
y_coor = pos(:, 2);


IEN = msh.QUADS(:, 1:4);
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

n_vector = lines; %向量
for nn = 1 : size(n_vector)
    if n_vector(nn, 3) == 13
        n_vector(nn, 1) = -(pos(lines(nn, 1), 1) + pos(lines(nn, 2), 1)) / 2 - 1;
        n_vector(nn, 2) = -(pos(lines(nn, 1), 2) + pos(lines(nn, 2), 2)) / 2 - 1;
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






%EOF