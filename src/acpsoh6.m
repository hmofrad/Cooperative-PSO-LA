function [x_std sbest_fit sbest_hist p_hist] = acpsoh6(f,bnd,dim,nop,endgen,alpha,beta)    
% Adaptive Cooperative Particle Swarm Optimization - HYBRID with Six Swarm
    vm = bnd(2);
    K = 6;
    c1 = 1.49;    c2 = 1.49;
    w_max = 0.9;  w_min = 0.4;
        w_max = 0.9;  w_min = 0.4;
    x = rand_bnd(nop,dim,bnd); % initialize position
    spd = rand([nop dim]); % initial velocity
    pbest = x; %initialize Best Particle Position
    sbest = sbest_init_s(x,f,K);     %initialize Swarm Best Position
    sbest_fit = ff({f},sbest);
    sbest_hist = []; % swarm best history
    % ========================
    % initilise Q Swarm parameters
    % ========================
    x_Q = rand_bnd(nop,dim,bnd);
    spd_Q = rand([nop dim]);
    x_Q(:,end+1)=0;
    for i=1:nop
        x_Q(i,end) = ff({f},x_Q(i,1:end-1));
    end
    pbest_Q = x_Q; %initialise Best Particle Position
    [bst ind] = min(x_Q(:,end));
    gbest_Q = x_Q(ind,:); % initialise global best position
        % =============================
        % Cooperative Learning Automata
        % =============================
    act1 = 1; % action 1 : do not start next swarm
    act2 = 2; % action 2 : start next swarm and perform cooperation
    action = [act1 act2];
    r = size(action,2);
    p = [.5 .5];
    p_hist = [];
    for i=1:endgen
        w = w_max-(w_max-w_min)*i/endgen;
        sbest_tmp = ff({f},sbest);
        rw = randsrc(1000,1,[action; p]); % roller wheel 
        act = rw(ceil(1000*rand));  % selected action
        switch (act)
            % ===== cooperative action number ONE ======
            case {1} % do not start next swarm
                for j=1:dim
                    for k=1:nop
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
            % ===== cooperative action numbet two =====
            case {2} % start next swarm
                for j=1:dim
                    for k=1:nop
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
                
                % ==== Interval TWO ====
                ind_Q = randi([1 ceil(nop/2)]);
                c = 0;
                while (pbest_Q(ind_Q,end) == gbest_Q(end))
                    ind_Q = randi([1 ceil(nop/2)]);
                    c = c+1;
                    if c == 7;
                        break
                    end
                end
                x_Q(ind_Q,1:end-1) =  sbest;
                x_Q(ind_Q,end) = sbest_fit; 
                
                % Q SWARM
                for j=1:nop
                    if (x_Q(j,end) < pbest_Q(j,end))
                        pbest_Q(j,:) = x_Q(j,:);
                    end
                    if (pbest_Q(j,end) < gbest_Q(end))
                        gbest_Q = pbest_Q(j,:);
                    end
                    % PSO update after each iteration
                    spd_Q(j,:) = w.* spd_Q(j,:)+c1.*rand(1,dim).*(pbest_Q(j,1:end-1)-x_Q(j,1:end-1))...
                                               +c2.*rand(1,dim).*(gbest_Q(1:end-1) - x_Q(j,1:end-1));
                    spdtmp_Q = spd_Q(j,:);
                    ind_Q = find(abs(spdtmp_Q>vm));
                    spdtmp_Q(ind_Q) = vm.*sign(spdtmp_Q(ind_Q));
                    sspd_Q(j,:) = spdtmp_Q;
                    x_Q(j,1:end-1) = x_Q(j,1:end-1)+spd_Q(j,:);
                    x_Q(j,end) = ff({f},x_Q(j,1:end-1));
                end
            % ==== Interval THREE ====
                for j=1:5:dim
                    % swarm length j:j+4
                    c = 0;
                    ind = randi([1 ceil(nop/2)]);
                    while (pbest(ind,j:j+4) == sbest(j:j+4))
                        ind = randi([1 ceil(nop/2)]);
                        c = c+1;
                        if c == 7
                            break
                        end
                    end
                    x(ind,j:j+4) =  gbest_Q(j:j+4); 
                end
        end
        sbest_fit = ff({f},sbest);
        sbest_hist = [sbest_hist sbest_fit];
        fprintf('iteration=%u,  sbest=%e\n',i,sbest_fit)        
        if sbest_fit < sbest_tmp
           p(act) = p(act) + alpha.*(1 - p(act)); % desired action
           p(action ~=act) = (1-alpha)*p(action ~=act);
        else
            p(act) = (1 - beta).*p(act); % non-desired action
            p(action ~=act) = (beta/(r-1))+(1-beta).*p(action ~=act);            
        end
        p_hist = [p_hist;p];
    end
    x_std = stdCalc(x,f);
end
