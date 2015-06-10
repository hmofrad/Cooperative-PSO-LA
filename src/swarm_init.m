function [x_Q spd_Q pbest_Q gbest_Q lbest_Q] = swarm_init(nos,nop,dim,f,bnd)
    % number of swarms (nos)
    n = 'rn';
    nn = 5;
    x_Q = []; spd_Q = [];
    for i=1:nos
        x_Q = [x_Q; rand_bnd(nop,dim,bnd)];
        spd_Q = [spd_Q; rand([nop dim])];  
    end
    x_Q(:,end+1)=0;
    for i=1:nos*nop
        x_Q(i,end) = ff({f},x_Q(i,1:end-1));
    end
     gbest_Q = [];
    for i=1:nos
        [bst ind] = min(x_Q((i-1)*nop+1:i*nop,end));
        gbest_Q = [gbest_Q; x_Q(ind,:)];
    end
    pbest_Q = x_Q;
    lbest_Q = zeros(size(x_Q));
    for i=1:nos
        for j=1:nop
            lbest_Q((i-1)*nop+j,:) = neighborhood({n},x_Q(((i-1)*nop)+1:i*nop,:),j,nn);
        end
    end
end
