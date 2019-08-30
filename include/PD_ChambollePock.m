%**************************************************************************
% authors: Nelly Pustelnik                                                *
% institution: Laboratoire de Physique de l'ENS de Lyon                   *
% date: October 18 2018                                                   *
% License CeCILL-B                                                        *
%**************************************************************************

function [x,crit]=PD_ChambollePock(data, param, op, prox, objective)
    
    
    %% Fixing Proximal Parameters
    gamma = 0.99;
    %mu_g=2;
    tau = gamma/param.normL;
    sig = gamma/param.normL;
    if tau*sig*param.normL^2>1
        disp('ERROR');
    end
    theta=1;
    %% Initializing variables
    x = zeros(size(data));
    y = op.direct(x);
    x0 = x;
    bx = x;
    %% Criterion of convergence
    crit=zeros(1,param.iter);
    i=2;
    crit(1) = 0;
    crit(i) = objective.fidelity(x,data) + objective.regularization(op.direct(x),param.lambda);
    
    %% Algorithm
    while abs(crit(i)-crit(i-1))/crit(i) > 1e-10
        i = i+1;
        %for i=1:param.iter
        %Update of primal variable
        tmp = y + sig*op.direct(bx);
        y = tmp - sig*prox.regularization(tmp/sig, param.lambda/sig);

        %Update of dual variable
        x = prox.fidelity(x0 - tau * op.adjoint(y),data,tau);
       
        %Update of the descent steps
        if param.mu>=1
            theta = (1+2*param.mu*tau)^(-1/2);
            tau = theta*tau;
            sig=sig/theta;
        end
        
        %Update dual auxiliary variable
        bx = x + theta*(x - x0);
        x0 = x;
        crit(i) = objective.fidelity(x,data) + objective.regularization(op.direct(x),param.lambda);
        
    end
end