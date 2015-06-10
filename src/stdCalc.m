function [x_std] = stdCalc(x,f,D)
    [nop dim] = size(x);
    if mod(dim-1,D) == 0
        dim = dim -1;
    end
    for i=1:nop
        x_fit(i) = ff({f},x(i,1:dim));
    end
        x_std = std(x_fit);
end
