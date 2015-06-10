function [x_std sbest_fit sbest_hist] = cpsos6(f,bnd,dim,nop,endgen)
    vm = bnd(2);
    K = 6;
    c1 = 1.49;
    c2 = 1.49;
    w_max = 0.9;
    w_min = 0.4;
    x = rand_bnd(nop,dim,bnd); % initialize position
    spd = rand([nop dim]); % initial velocity
    pbest = x; %initialize Best Particle Position
    %initialize Swarm Best Position
    sbest = sbest_init_s(x,f,K);
    sbest_hist = []; % swarm best history
    for i=1:endgen
        w = w_max-(w_max-w_min)*i/endgen;
        for j=1:dim
            for k=1:nop  % inter swarm dimentional index
                if(ff({f},b(sbest,j,x(k,j))) < ff({f},b(sbest,j,pbest(k,j))))
                    pbest(k,j) = x(k,j);
                end
                if(ff({f},b(sbest,j,pbest(k,j))) < ff({f},sbest))
                    sbest(j) = pbest(k,j);
                end
            end
            if mod(j,5) == 0
                sl = j-4:j; %swarm length
                spd(:,sl) = w.* spd(:,sl)+c1.*rand(nop,length(sl)).*(pbest(:,sl)-x(:,sl))...
                                       +c2.*rand(nop,length(sl)).*(repmat(sbest(sl),nop,1) - x(:,sl));
                spdtmp = spd(:,sl);
                ind = find(abs(spdtmp>vm));
                spdtmp(ind) = vm.*sign(spdtmp(ind));
                spd(:,sl) = spdtmp;
                x(:,sl) = x(:,sl)+spd(:,sl);
            end
        end
        sbest_fit = ff({f},sbest);
        sbest_hist = [sbest_hist sbest_fit];
        fprintf('iteration=%u,sbest=%e\n',i,sbest_fit)
    end
    x_std = stdCalc(x,f);
end
