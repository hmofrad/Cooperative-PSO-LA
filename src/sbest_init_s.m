function [sbest] = sbest_init_s(x,f,K)
    % initilise swarm best by K1 and K2 swarms
    nop = size(x,1);
    dim = size(x,2);
    [K1 K1_d K2 K2_d] = splitter(dim,K); %split factor matrix [K1 K1dim;K2 K2dim]
    %initialize Swarm Best Position
    B =[];
    if K1 == 0;
        for i=K1+1:K2_d:dim  % i norm enforce number of swarms
            for j=i:i+K2_d-1 % j enforce as number of columns in each swarm  
                for k=1:nop  % k enforce as number of particles in each column of swarm
                B(k,j) = ff({f},x(k,j));
                end              
            end
            [x1 i1]  = min(B(:,i:j));
            [x2 i2]  = min(x1);
            sbest(i:j) = x(i1(i2),i2+i-1);
        end
    else
        for i=1:K1_d:K1*K1_d
            for j=i:i+K1_d-1
                for k=1:nop
                B(k,j) = ff({f},x(k,j));
                end              
            end
            [x1 i1]  = min(B(:,i:j));
            [x2 i2]  = min(x1);
            sbest(i:j) = x(i1(i2),i2+i-1);
        end
        for i=K1*K1_d+1:K2_d:dim
            for j=i:i+K2_d-1
                for k=1:nop
                B(k,j) = ff({f},x(k,j));
                end              
            end
            [x1 i1]  = min(B(:,i:j));
            [x2 i2]  = min(x1);
            sbest(i:j) = x(i1(i2),i2+i-1);
        end 
    end     
end
