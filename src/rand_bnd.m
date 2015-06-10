function [x] = rand_bnd(nop,dim,bnd)
    lbound = bnd(1);
    ubound = bnd(2)/4;
    x = [];
    for i =1:nop
        for j = 1:dim
            if(rand>0.5)
                x(i,j) = lbound*rand;
            else
                x(i,j) = ubound*rand;
            end
        end
    end
end
