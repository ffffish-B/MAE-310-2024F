clear all; clc;

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

% Force function in different direction for plane stress
% f_1 = @(x,y)(E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/((2*v - 1)*(v + 1))...
%     - (E*(v - 1/2)*((x - 1)*(y - 1) + x*y + 2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/((2*v - 1)*(v + 1)) - (2*E*y*(v - 1)*(y - 1))/((2*v - 1)*(v + 1));
% f_2 = @(x,y)(E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/((2*v - 1)*(v + 1))...
%     - (E*(v - 1/2)*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/((2*v - 1)*(v + 1)) - (2*E*x*(v - 1)*(x - 1))/((2*v - 1)*(v + 1));
%plane stress
D(1,1) = E/(1-v^2);
D(2,2) = D(1,1);
D(1,2) = v*E/(1-v^2);
D(2,1) = D(1,2);
D(3,3) = E/2/(1+v);

% Force function in different direction for plane stress
f_1 = @(x,y)(2*E*y*(y - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) ...
    + x*y + 2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);
f_2 = @(x,y)(2*E*x*(x - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) ...
    + x*y + x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);


% quadrature rule
n_int_xi  = 10;
n_int_eta = 10;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

errors_L2 = [];
errors_H1 = [];
mesh_sizes = [];
% mesh generation
for n_el_x = [10, 20, 40, 80] % 网格密度
    n_en   = 4;               % number of nodes in an element
    n_el_y = n_el_x;               % number of elements in y-dir
    n_el   = n_el_x * n_el_y; % total number of elements

    n_np_x = n_el_x + 1;      % number of nodal points in x-dir
    n_np_y = n_el_y + 1;      % number of nodal points in y-dir
    n_np   = n_np_x * n_np_y; % total number of nodal points

    x_coor = zeros(n_np, 1);
    y_coor = x_coor;

    hx = 1.0 / n_el_x;        % mesh size in x-dir
    hy = 1.0 / n_el_y;        % mesh size in y-dir
    h = max(hx ,hy);

    % generate the nodal coordinates
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx; % nodal index
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end

    % IEN array
    IEN = zeros(n_el, n_en);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index
            IEN(ee, 1) = (ey-1) * n_np_x + ex;
            IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
            IEN(ee, 3) =  ey    * n_np_x + ex + 1;
            IEN(ee, 4) =  ey    * n_np_x + ex;
        end
    end

    % ID array
    ID = zeros(n_np,2);
    counter = 1;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            ID(index, 1) = counter;
            ID(index, 2) = counter + 1;
            counter = counter + 2;
        end
    end

    n_eq = counter - 1;

    LM1 = zeros(size(IEN));
    LM2 = zeros(size(IEN));
    for AA = 1 : n_el
        for BB = 1 : n_en
            index = IEN(AA, BB);
            LM1(AA, BB)= ID(index, 1);
        end
    end


    for AA = 1 : n_el
        for BB = 1 : n_en
            index = IEN(AA, BB);
            LM2(AA, BB)= ID(index, 2);
        end
    end

    % allocate the stiffness matrix and load vector
    K = zeros(n_eq, n_eq);
    F = zeros(n_eq, 1);

    % loop over element to assembly the matrix and vector
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );

        k_ele = zeros(2*n_en, 2*n_en); % element stiffness matrix
        f_ele = zeros(2*n_en, 1);    % element load vector

        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

            for aa = 1 : n_en
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                B_a = zeros(3,2);
                B_a(1,1) = Na_x;
                B_a(2,2) = Na_y;
                B_a(3,1) = Na_y;
                B_a(3,2) = Na_x;

                pp = 2*(aa-1);

                f_ele(pp+1) = f_ele(pp+1) + weight(ll) * detJ * f_1(x_l, y_l) * Na;
                f_ele(pp+2) = f_ele(pp+2) + weight(ll) * detJ * f_2(x_l, y_l) * Na;

                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
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
                            k_ele(p, q) = k_ele(p, q) + weight(ll) * detJ...
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

    %error
    error_L2 = 0;
    error_H1 = 0;

    for ee = 1:n_el
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );
        u_h = [disp(IEN(ee, 1:n_en),1), disp(IEN(ee, 1:n_en),2)]';

        for ll = 1 : n_int
            x_l = 0.0; y_l = 0.0; dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

            u_exact = [exact(x_l, y_l), exact(x_l, y_l)]';
            ux_exact = [exact_x(x_l, y_l), exact_x(x_l, y_l)]';
            uy_exact = [exact_y(x_l, y_l), exact_y(x_l, y_l)]';

            uh = zeros(2,1); uh_x = zeros(2,1); uh_y = zeros(2,1);
            for aa = 1:n_en
                for i = 1:2
                    Na = Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] =Quad_grad(aa, xi(ll), eta(ll));
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                    uh(i) = uh(i) + u_h(i,aa) * Na;
                    uh_x(i) = uh_x(i) + u_h(i,aa) * Na_x;
                    uh_y(i) = uh_y(i) + u_h(i,aa) * Na_y;
                end
            end

            e = uh - u_exact;
            ex = uh_x - ux_exact;
            ey = uh_y - uy_exact;

            error_L2 = error_L2 + weight(ll) * detJ * e.^2;
            error_H1 = error_H1 + weight(ll) * detJ * (ex.^2 + ey.^2);
        end
    end
    error_L2 = sqrt(error_L2);
    error_H1 = sqrt(error_H1);

    errors_L2 = [errors_L2, error_L2];
    errors_H1 = [errors_H1, error_H1];
    mesh_sizes = [mesh_sizes, h];
    % save the solution vector and number of elements to disp with name
    % HEAT.mat
    save("HEAT", "disp", "n_el_x", "n_el_y");
end


% close all;
% [X, Y] = meshgrid(0 : hx : 1, 0 : hy : 1);
% hold on
% Z = reshape(disp(:,1), n_np_x, n_np_y)';
% surf(X, Y, Z);
%
% shading interp;
% axis equal;

figure;
loglog(mesh_sizes, errors_L2(1,:), '-o', 'DisplayName', 'L2 Norm');
hold on;
loglog(mesh_sizes, errors_H1(1,:), '-o', 'DisplayName', 'H1 Norm');

% Linear fit to find slopes
log_h = log(mesh_sizes);
log_L2 = log(errors_L2(1,:));
log_H1 = log(errors_H1(1,:));

coeff_L2 = polyfit(log_h, log_L2, 1);
coeff_H1 = polyfit(log_h, log_H1, 1);

slope_L2 = coeff_L2(1);
slope_H1 = coeff_H1(1);

% Add fitted lines
fit_L2 = exp(polyval(coeff_L2, log_h));
fit_H1 = exp(polyval(coeff_H1, log_h));

loglog(mesh_sizes, fit_L2, '--', 'DisplayName', sprintf('L2 Fit (Slope: %.2f)', slope_L2));
loglog(mesh_sizes, fit_H1, '--', 'DisplayName', sprintf('H1 Fit (Slope: %.2f)', slope_H1));

xlabel('Mesh Size (h)');
ylabel('Error');
legend;
grid on;
title('Error Convergence and Fitted Slopes');


function e_i = unit_vector(ii)
if ii == 1
    e_i = [1,0]';
elseif ii == 2
    e_i = [0,1]';
end
end


    % EOF