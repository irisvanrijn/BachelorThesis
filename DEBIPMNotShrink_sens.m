function MostSensPar  = DEBIPMNotShrink_sens(E_Y,MatrixSize,Lmin,Lp,rb,Rm,Lm,E_Ystdev,kappa,Lb,mu_juv,mu_ad);

% INSTRUCTIONS 
% to run this code, you have to give the values of the parameters that are
% in between brackets for each species. For example (based on the reef manta ray M. alfredi):

% MostSensPar  = DEBIPMNotShrink_sens(0.755102041,200,4,84.7,0.769,76.5,130,0.1,0.80,5,0.29,0.35)

% pars perturbed are in order: Lb, Lp, Lm, Rm, rb, mu_adult, mu_juv

% MostSensPar is an index value idenfitying which parameter is the most
% influential to the population growth rate. For example, if MostSensPar
% = 1, then Lb is the parameter that has the highest influence on lambda.
% If MostSensPar = 7, then juvenile mortality rate has the highest influence on lambda.

% INSTRUCTIONS END

% START CODE

fracperturb = 0.01; % perturbation of DEB pars

% pars perturbed are: Lb, Lp, Lm, Rm, rb, mu_adult, mu_juv
pars_orig = [Lb, Lp, Lm, Rm, rb, mu_ad, mu_juv]; % create vector with default parameter values

% calculate kernel and log population growth rate with default values
[S_d, R_d, G_d, D_d, y] = BigMatrixNotShrink(MatrixSize,Lmin,Lm,Lp,Rm,Lm,E_Y,E_Ystdev,rb,kappa,Lb,mu_juv,mu_ad);
kernel_b = G_d*S_d + D_d*R_d; 

[W,d] = eig(kernel_b); lambda1 = diag(d); imax = find(lambda1==max(lambda1)); 
V=conj(inv(W)); loglambda_orig = lambda1(imax); % population growth rate
            
% do sensitivity analysis whereby each model parameter is perturbed by fracperturb
% parameters are perturbed sequentially

       for i=1:length(pars_orig) % Lb, Lp, Lm, Rm, rb, mu_adult, mu_juv
           pars_perturb = pars_orig; % reset parameter values
           pars_perturb(i) = pars_orig(i)*(1+fracperturb); %perturb pars value 
           Lb = pars_perturb(1); Lp = pars_perturb(2); Lm = pars_perturb(3); Rm = pars_perturb(4); 
           rb = pars_perturb(5); mu_ad = pars_perturb(6); mu_juv = pars_perturb(7);

           % calculate kernel with perturbed parameter value
            [S_p, R_p, G_p, D_p, y] = BigMatrixNotShrink(MatrixSize,Lmin,Lm,Lp,Rm,Lm,E_Y,E_Ystdev,rb,kappa,Lb,mu_juv,mu_ad);
            kernel_p = G_p*S_p + D_p*R_p; 
           % calculate log population growth rate with perturbed parameter value 
            [W,d] = eig(kernel_p); lambda1 = diag(d); imax = find(lambda1==max(lambda1)); 
            V=conj(inv(W)); loglambda_pert = lambda1(imax); % population growth rate
            
            LogLambda_perturb(i,1) = loglambda_pert; % vector with perturbed log lambda values
         end
                
        %calculate highest sensitivity value and idenfity which model
        %parameter has highest sensitivity value
        
        LogLambda_Sens = (LogLambda_perturb - loglambda_orig)./fracperturb; 
        LogLambda_Sens = abs(LogLambda_Sens); % use absolute elasticity values
        [a b] = max(LogLambda_Sens);
        if b==1, MostSensPar = 1;     
        elseif b==2; MostSensPar=2; 
        elseif b==3, MostSensPar=3; 
        elseif b==4, MostSensPar=4; 
        elseif b==5, MostSensPar=5; 
        elseif b==6, MostSensPar=6; 
        else MostSensPar=7;
        end                
        
end

