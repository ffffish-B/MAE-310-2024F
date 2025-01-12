
clear; clc;close all

errors_L2_whole = zeros(4,3);
mesh_sizes = zeros(4,3);
run('tu_5.m')
run("get_error.m")
for ii = 1:3
    errors_L2_whole(1,ii) = errors_L2(ii);
    mesh_sizes(1,ii) = h;
end

run('tu_10.m')
run("get_error.m")
for ii = 1:3
    errors_L2_whole(2,ii) = errors_L2(ii);
    mesh_sizes(2,ii) = h;
end

run('tu_20.m')
run("get_error.m")
for ii = 1:3
    errors_L2_whole(3,ii) = errors_L2(ii);
    mesh_sizes(3,ii) = h;
end

run('tu_40.m')
run("get_error.m")
for ii = 1:3
    errors_L2_whole(4,ii) = errors_L2(ii);
    mesh_sizes(4,ii) = h;
end



figure;
loglog(mesh_sizes, errors_L2_whole(:,1), '-o', 'DisplayName', 'L2 Norm');

% Linear fit to find slopes
log_h = log(mesh_sizes(:,1));
log_L2 = log(errors_L2_whole(:,1));

coeff_L2 = polyfit(log_h, log_L2, 1);

slope_L2 = coeff_L2(1);
hold on
% Add fitted lines
fit_L2 = exp(polyval(coeff_L2, log_h));

loglog(mesh_sizes, fit_L2, '--', 'DisplayName', sprintf('L2 Fit (Slope: %.2f)', slope_L2));

xlabel('Mesh Size (h)');
ylabel('Error');
legend;
grid on;
title('Error Convergence and Fitted Slopes');
