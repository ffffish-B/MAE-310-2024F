clear all; clc;clf;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

errors_L2 = [];
errors_H1 = [];
mesh_sizes = [];
% mesh generation
for n_el_x = [10, 20, 40, 80] % 网格密度
    n_en   = 4;               % number of nodes in an element
    n_el_y = n_el_x;               % number of elements in xy-dir
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
    IEN_tri = zeros(n_el, n_en); % 两倍的单元数，每个单元三个节点
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index
            IEN(ee, 1) = (ey-1) * n_np_x + ex;
            IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
            IEN(ee, 3) =  ey    * n_np_x + ex + 1;
            IEN(ee, 4) =  ey    * n_np_x + ex;
        end
    end


    % 上边所有的都在建立IEN数组
    % 建立ID数组，让边界的P=0

    % ID array
    ID = zeros(n_np,1);
    counter = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter = counter + 1;
            ID(index) = counter;
        end
    end

    n_eq = counter;

    LM = ID(IEN);
    % 打开K和F的空矩阵
    % allocate the stiffness matrix and load vector
    K = spalloc(n_eq, n_eq, 9 * n_eq);
    F = zeros(n_eq, 1);

    % loop over element to assembly the matrix and vector
    for ee = 1 : n_el %所有的节点循环，，建立ele
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );

        k_ele = zeros(n_en, n_en); % element stiffness matrix
        f_ele = zeros(n_en, 1);    % element load vector

        for ll = 1 : n_int %2D高斯点循环
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
            for aa = 1 : n_en %在map中用雅各比矩阵换元
                x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
                y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                dx_deta = dx_deta + x_ele(aa) * Na_eta;
                dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                dy_deta = dy_deta + y_ele(aa) * Na_eta;
            end

            detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

            for aa = 1 : n_en %用雅各比矩阵进行高斯积分-Na
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;

                for bb = 1 : n_en %''-Nb
                    Nb = Quad(bb, xi(ll), eta(ll));
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                    k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end % end of bb loop
            end % end of aa loop
        end % end of quadrature loop
        %结束单元内的组装
        %找P,Q进行整个大刚度阵的组装
        for aa = 1 : n_en
            PP = LM(ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(aa);

                for bb = 1 : n_en
                    QQ = LM(ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % modify F with the boundary data
                        % here we do nothing because the boundary data g is zero or
                        % homogeneous
                    end
                end
            end
        end
    end

    % solve the stiffness matrix
    dn = K \ F;

    % insert dn back into the vector for all nodes
    disp = zeros(n_np, 1);

    for ii = 1 : n_np
        index = ID(ii);
        if index > 0
            disp(ii) = dn(index);
        else
            %作业这里要改成非零的g
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end

    %error
    error_L2 = 0;
    error_H1 = 0;

    for ee = 1:n_el
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );
        u_h = disp(IEN(ee, 1:n_en));

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

            u_exact = exact(x_l, y_l);
            ux_exact = exact_x(x_l, y_l);
            uy_exact = exact_y(x_l, y_l);

            uh = 0.0; uh_x = 0.0; uh_y = 0.0;
            for aa = 1:n_en
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] =Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                uh = uh + u_h(aa) * Na;
                uh_x = uh_x + u_h(aa) * Na_x;
                uh_y = uh_y + u_h(aa) * Na_y;
            end

            e = uh - u_exact;
            ex = uh_x - ux_exact;
            ey = uh_y - uy_exact;

            error_L2 = error_L2 + weight(ll) * detJ * e^2;
            error_H1 = error_H1 + weight(ll) * detJ * (ex^2 + ey^2);
        end
    end
        error_L2 = sqrt(error_L2);
    error_H1 = sqrt(error_H1);

    errors_L2 = [errors_L2, error_L2];
    errors_H1 = [errors_H1, error_H1];
    mesh_sizes = [mesh_sizes, h];
    % save the solution vector and number of elements to disp with name
    % HEAT.mat
    % save("HEAT", "disp", "n_el_x", "n_el_y");
end


figure;
loglog(mesh_sizes, errors_L2, '-o', 'DisplayName', 'L2 Norm');
hold on;
loglog(mesh_sizes, errors_H1, '-o', 'DisplayName', 'H1 Norm');

% Linear fit to find slopes
log_h = log(mesh_sizes);
log_L2 = log(errors_L2);
log_H1 = log(errors_H1);

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







% EOF