%%%%%%%%%%%%%%%% The 'big matrix' M of size n x n %%%%%%%%%%%%%%%%%%%%%%%%%
function [S, R, G, D, y] = BigMatrixNotShrink(n,lowval,hival,Lp,Rm,Lm,E_Y,Ystdev,rb,kappa,Lb,mu_juv,mu_ad)
	% upper and lower integration limits
	L = lowval; U = 1.1*hival;
	% boundary points b and mesh points y
	c = 0:1:n; b = L+c.*(U-L)./n;
    y = 0.5.*(b(1:n)+b(2:n+1));
	% create S, R, G and D matrices
	S = SurvivalFunction(y,E_Y,Lm,mu_juv,mu_ad,kappa,Lp); S=repmat(S,n,1); S=S.*eye(n);
	R = ReproductionFunction(y,Lp,Rm,Lm,kappa,E_Y); R=repmat(R,n,1); R=R.*eye(n);
    G = GrowthFunction(y,E_Y,Ystdev,Lm,rb); % plot(y,mutemp)
	D = DevelopmentFunction(y,Lb);
	% scale D and G so columns sum to 1
	[m,n]=size(G); G = G./repmat(sum(G),m,1); 
    [m,n]=size(D); D = D./repmat(sum(D),m,1); 

% Compute the kernel component functions from the fitted models

% survival
function [sx] = SurvivalFunction(x,E_Y,Lm,mu_juv,mu_ad,kappa,Lp)
    for i=1:length(x)
        if x(i) <= Lm*E_Y/kappa; % only survivors 
            if x(i)<=Lp, sx(i) = exp(-mu_juv);
            else sx(i) = exp(-mu_ad); end
        else sx(i)=0; end
    end
    

% growth
function gxy = GrowthFunction(x,E_Y,Ystdev,Lm,rb) 
    for i=1:length(x)
        if (x(i)>0 && x(i)<=Lm*E_Y); % when respiration energy is sufficient to cover maintenance
        mux = x(i).*exp(-rb) + (1 - exp(-rb))*Lm*E_Y;
        sigmax2 = (1 - exp(-rb))*Lm*Ystdev;
        sigmax = max(0.00001,sigmax2);
        fac1 = sqrt(2.*pi).*sigmax; fac2 = ((x'-mux).^2)./(2.*sigmax.^2);
        gxy(:,i) = (exp(-fac2)./fac1); 
        else gxy(:,i)=x; 
        end; 
     end

% reproduction
function [rx] = ReproductionFunction(x,Lp,Rm,Lm,kappa,E_Y)
    rx2 = x>=Lp; % matrix with zeros and ones for inds that are above the threshold size for reproduction
    rx1 = zeros(1,length(x));
    for i=1:length(x)
        if x(i) <= Lm*E_Y; % when respiration energy is sufficient to cover maintenance
        rx1(i) =  Rm/(Lm^2).*E_Y.*x(i).^2; 
        elseif x(i) <= Lm*E_Y/kappa; % only survivors reproduce
            rx(i)=Rm/(1-kappa)*(E_Y*x(i).^2-kappa*x(i).^3/Lm); end
        rx=rx1.*rx2; % only adults reproduce
    end

% development
function [dxy] = DevelopmentFunction(x,Lb)
    xtemp = x<=Lb; 
    dxy=repmat(xtemp',1,length(xtemp));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%