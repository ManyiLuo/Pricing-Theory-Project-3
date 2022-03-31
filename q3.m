%% q3
r0 = 0.02;
alpha = 3;
sigma= 0.01;
theta0 = 0.03;
beta = 1;
phi = 0.05;
eta = 0.005;

T = linspace(0,10,100);


%% Change in alpha
valpha = [2, 3, 4, 5, 10];
%subplot(3,3,1);
for i = 1: size(valpha,2)
    yield = bondyield(T, valpha(i), beta, eta, sigma, phi, r0, theta0);
    plotfig(T, yield);   
end
legend({'alpha=2','alpha=3','alpha=4','alpha=5','alpha=10'},'Location','best','FontSize',12);


%% Change in beta
vbeta = [0.1, 0.5, 1, 2, 5];
%subplot(3,3,2);
for i = 1: size(vbeta,2)
    yield = bondyield(T, alpha, vbeta(i), eta, sigma, phi, r0, theta0);
    plotfig(T, yield); 
end
legend({'beta=0.1','beta=0.5','beta=1','beta=2', 'beta=5'},'Location','best','FontSize',12);


%% Change in eta
%subplot(3,3,3);
veta = [0.001, 0.01, 0.05, 0.1, 0.15];
for i = 1: size(veta,2)
    yield = bondyield(T, alpha, beta, veta(i), sigma, phi, r0, theta0);
    plotfig(T, yield); 
end
legend({'eta=0.001','eta=0.01','eta=0.05','eta=0.1', 'eta=0.15'},'Location','best','FontSize',12);


%% Change in sigma
vsigma = [0.01, 0.05, 0.1, 0.2, 0.3];
%subplot(3,3,4);
for i = 1: size(vsigma,2)
    yield = bondyield(T, alpha, beta, eta, vsigma(i), phi, r0, theta0);
    plotfig(T, yield); 
end
legend({'sigma=0.01','sigma=0.05','sigma=0.1','sigma=0.2', 'sigma=0.3'},'Location','best','FontSize',12);


%% Change in phi
vphi = [0.01, 0.05, 0.1, 0.5, 1];
%subplot(3,3,5);
for i = 1: size(vsigma,2)
    yield = bondyield(T, alpha, beta, eta, sigma, vphi(i), r0, theta0);
    plotfig(T, yield); 
end
legend({'phi=0.01','phi=0.05','phi=0.1','phi=0.5', 'phi=1'},'Location','best','FontSize',12);


%% Change in r0
vr0 = [0.02, 0.05, 0.08, 0.1, 0.2];
%subplot(3,3,6);
for i = 1: size(vsigma,2)
    yield = bondyield(T, alpha, beta, eta, sigma, phi, vr0(i), theta0);
    plotfig(T, yield); 
end
legend({'r0=0.02','r0=0.05','r0=0.08','r0=0.1', 'r0=0.2'},'Location','best','FontSize',12);


%% change in theta0

vtheta0 = [0.03, 0.04, 0.05, 0.08, 0.1];
%subplot(3,3,7);
for i = 1: size(vsigma,2)
    yield = bondyield(T, alpha, beta, eta, sigma, phi, r0, vtheta0(i));
    plotfig(T, yield); 
end
legend({'theta0=0.03','theta0=0.04','theta0=0.05','theta0=0.08', 'theta0=0.1'},'Location','best','FontSize',12);




%%
function plotfig (T, yield)
    plot(T, yield,'LineWidth',1);
    title('Term Structure of Interest Rates','FontSize',12)
    xlabel('maturity time T','FontSize',12);
    ylabel('Bond Yield Y(T)','FontSize',12);
    hold on;
end

function yield_anal = bondyield (T, alpha, beta, eta, sigma, phi, r0, theta0)

    yield_anal = zeros(size(T,1), 1);

    for i = 1:size(T,2)
        B = (1 - exp(-alpha * T(i)))/alpha;
        C = exp(-alpha * T(i))/(alpha-beta) + 1/beta - alpha * exp(-beta*T(i))/(beta*(alpha-beta)); % - exp(-beta * T(i))/(alpha-beta) - exp(-beta * T(i))/beta + 1/beta;


        fun1 = @(s) ((1 - exp(-alpha * (T(i)-s)))/alpha).^2;

        omega1 = integral(fun1,0,T(i));

        fun2 = @(s) (exp(-alpha * (T(i)-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T(i)-s))/(beta*(alpha-beta))).^2;

        omega2 = integral(fun2,0,T(i));

        A = phi * (alpha * (1-exp(-beta*T(i)))/((alpha-beta)*beta) ...                  
                       - beta * (1-exp(-alpha*T(i)))/(alpha*(alpha-beta))...
                       - T(i))...
                       + 0.5 * sigma^2 * omega1 + 0.5 * eta^2 * omega2;

        price = exp(A - B * r0 - C * theta0);
        yield = -log(price)/T(i);

        yield_anal(i) = yield;
    end
    
   
   
end