function [K1 K1_dim K2 K2_dim] = splitter(dim,K)
    % Initialise K1 ceil(n/K)-dimensional PSO
    % Initialise K2 floor(n/K)-dimensional PSO
    n = dim;
    K1 = mod(n,K);
    K2 = K - K1;
    K1_dim = ceil(n/K);
    K2_dim = floor(n/K);
end
