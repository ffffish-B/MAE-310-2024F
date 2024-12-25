function N = Tri(node, xi, eta)
    if node == 1
        N = 1 - xi - eta;
    elseif node == 2
        N = xi;
    elseif node == 3
        N = eta;
    end
end