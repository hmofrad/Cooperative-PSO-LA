function [x_std sbest_fit sbest_hist] = cpsos(f,bnd,dim,nop,endgen)
    % dim as numbre of swarms
    vm = bnd(2);
    c1 = 1.49;      c2 = 1.49;
    w_max = 0.9;    w_min = 0.4;
    x = rand_bnd(nop,dim,bnd); % initialize position
    spd = rand([nop dim]); % initial velocity
    pbest = x; %initialize Best Particle Position
    for i=1:dim
        for j=1:nop
            B(j) = ff({f},x(j,i));
        end
        [mn ind] = min(B);
        sbest(i) = x(ind,i);
    end
    sbest_hist = []; % swarm best history
    for i=1:endgen
        w = w_max-(w_max-w_min)*i/endgen;
        for j=1:dim
            for k=1:nop
                if(ff({f},b(sbest,j,x(k,j))) < ff({f},b(sbest,j,pbest(k,j))))
                    pbest(k,j) = x(k,j);
                end
                if(ff({f},b(sbest,j,pbest(k,j))) < ff({f},sbest))
                    sbest(j) = pbest(k,j);
                end
            end            
            spd(:,j) = w.* spd(:,j)+c1.*rand(nop,1).*(pbest(:,j)-x(:,j))...
                                   +c2.*rand(nop,1).*(sbest(j) - x(:,j));
            spdtmp = spd(:,j);
            ind = find(abs(spdtmp>vm));
            spdtmp(ind) = vm.*sign(spdtmp(ind));
            spd(:,j) = spdtmp;
            x(:,j) = x(:,j)+spd(:,j);
        end
        sbest_fit = ff({f},sbest);
        sbest_hist = [sbest_hist sbest_fit];
%         fprintf('iteration=%u,sbest=%e\n',i,sbest_fit)
    end
    x_std = stdCalc(x,f);
end
