function fit_val = ff(c,x)
% c as the fitness function number
    fit_func = c{1};
    fit_val = 0;
    
    switch (fit_func)
        case {'Rosenbrock',0}
%             xind = 1:length(x)/2;
%             fit_val = sum(100.*((x(2*xind) - (x(2*xind-1).^2))).^2 + (1 - x(2*xind-1)).^2);
            xind = 1:length(x)-1;
            fit_val = sum(100.*(x(xind).^2-x(xind+1)).^2+(x(xind)-1).^2);
        case {'Quadric',1}  
            for i=1:length(x)
                xtmp = (sum(x(1:i)))^2;
                fit_val  = fit_val+xtmp;
            end
        case {'Ackley',2}
        	fit_val = -20*exp(-0.2*sqrt(sum(x.^2)/length(x))) - exp(sum(cos(2*pi*x))/length(x)) + 20 + exp(1);
            
        case {'Rastrigin',3}
            fit_val = sum(x.^2 - 10*cos(2*pi*x) + 10);
        case {'Griewank',4}  
            xind = 1:length(x);
            fit_val = (1/4000)*sum(x.^2) - prod(cos(x./sqrt(xind))) + 1;
    end
end
