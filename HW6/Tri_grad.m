function [N_xi, N_eta] = Tri_grad(node, xi, eta)
    if node == 1
        N_xi=-1;
        N_eta=-1;
    elseif node == 2
        N_xi = 1;
        N_eta = 0;
    elseif node == 3
        N_xi = 0;
        N_eta = 1;
    end
end