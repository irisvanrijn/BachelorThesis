function [E_Y, SSD, PopQs] = DEBIPMNotShrink_YellowMudTurtle(E_Ymin,E_Ymax,step,MatrixSize)

% input: 
% E_Ymin: minimum feeding level; suggestion: 0.5
% E_Ymax: maximum feeding level; suggestion: 1
% step: step size between E_Ymin and E_Ymax: suggestion: 50
% MatrixSize: number of size bins of matrix approximation of DEBIPM. suggestion: 200

% to run code with above suggested parameter values: [E_Y, SSD, PopQs] = DEBIPMNotShrink_MantaRay(0.5,1,50,200);

% DEB parameters 

Lmin = 17.52; % cm used to define limits of matrix approximation
Lb = 21.9; % disc width at birth in cm 
Lp = 90; % disc width at puberty in cm
Lm = 117; % maximum disc width in cm 
kappa = 0.80; % fraction energy allocation to respiration (as opposed to reproduction)
mu_juv = 0.53; % juvenile mortality rate (per year)
mu_ad = 0.055; % adult mortality rate (per year)
rb = 0.11; % von Bertalanffy growth rate (per year)
Rm = 6; % maximum number of offspring produced per year at L_m

E_Ystdev=0.1; % sigma(Y)

E_Y=linspace(E_Ymin,E_Ymax,step); % Feeding level
       
        for n=1:length(E_Y); n/length(E_Y)                       
                        
            [S, R, G, D, meshpts] = BigMatrixNotShrink(MatrixSize,Lmin,Lm,Lp,Rm,Lm,E_Y(n),E_Ystdev,rb,kappa,Lb,mu_juv,mu_ad);
            kernelDEB = G*S + D*R; 
           % calculate lambda, ssd, and rv
            [W,d] = eig(kernelDEB); lambda1 = diag(d); imax = find(lambda1==max(lambda1)); 
            V=conj(inv(W)); lambda = lambda1(imax); % population growth rate
            w1=W(:,imax); v1 = real(V(imax,:))';
            ssd = w1/sum(w1); ssd = ssd'; % stable stage distribution

            % calculation of LRS
            [m,p]=size(kernelDEB); TS = G*S; Fmat = D*R; 
            R0mat = Fmat*inv((eye(m)-TS)); meanLRS = max(eig(R0mat)); 

            %  generation time
            GT = log(meanLRS)/log(lambda);

            PopQstemp = [E_Y(n),Lm,Lp,Rm,lambda,meanLRS,GT];
            PopQs(n,:) = PopQstemp;
            SSD(:,n) = ssd;
        end
        
        subplot(3,1,1); plot(E_Y,PopQs(:,5),'k-','LineWidth',2);  axis square;  box off
        ylabel('Population growth rate','FontSize',14); xlabel(' ');  
        subplot(3,1,2); plot(E_Y,PopQs(:,6),'k-','LineWidth',2);   axis square;  box off
        ylabel('R_0','FontSize',14);  xlabel(' ');  set(gca,'FontSize',12); 
        subplot(3,1,3); plot(E_Y,PopQs(:,7),'k-','LineWidth',2); axis square;  box off; hold on
        ylabel('Generation time','FontSize',14);  xlabel('Feeding level');  set(gca,'FontSize',12);