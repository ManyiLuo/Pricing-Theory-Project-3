r0 = 0.02;
alpha = 3;
sigma= 0.01;
theta0 = 0.03;
beta = 1;
phi = 0.05;
eta = 0.005;

scale = [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05];
T1 = 3;
T2 = 5;
dt = 0.01;
Nsims = 1000;
option_value = zeros(size(scale,2),1);
option_value2 = zeros(size(scale,2),1);

[T1_A, T1_B, T1_C] = bondABC (T1, 0, alpha, beta, eta, sigma, phi);
[T2_A, T2_B, T2_C] = bondABC (T2, 0, alpha, beta, eta, sigma, phi);

P0_T1 = exp(T1_A - T1_B * r0 - T1_C * theta0);
P0_T2 = exp(T2_A - T2_B * r0 - T2_C * theta0);
K_ori = P0_T2/P0_T1;


%% Under the risk-neutral measure
for i = 1: size(scale,2)
    
    K = K_ori * scale(i);
   
    theta = zeros(Nsims, T2/dt);
    r = zeros(Nsims, T2/dt);
        
    theta(:,1) = theta0;
    r(:,1) = r0;

    %simulate theta and r using Euler Scheme
    for j = 1:T2/dt-1
          
       change_theta = beta * (phi - theta(:,j)) * dt + eta * sqrt(dt) * randn(Nsims,1);
       theta(:,j+1) = theta(:,j) + change_theta;

       change_r = alpha * (theta(:,j) - r(:,j)) * dt + sigma * sqrt(dt) * randn(Nsims,1);
       r(:,j+1) = r(:,j) + change_r;
    end

    PT1_T2 = exp(-dt*sum(r(:,(T1/dt:end-1)), 2));
    
    value = max(PT1_T2 - K, 0).*exp(-dt*sum(r(:,(1:T1/dt)), 2));

    option_value(i) = mean(value);
end
   
plot(scale*K_ori, option_value);
hold on;

%% Under forward neutral

for i = 1: size(scale,2)
    
    K = K_ori * scale(i);
   
    theta = zeros(Nsims, T2/dt);
    r = zeros(Nsims, T2/dt);
        
    theta(:,1) = theta0;
    r(:,1) = r0;
    %simulate theta and r using Euler Scheme
    for j = 1:T2/dt-1
       [T1t_A, T1t_B, T1t_C] = bondABC (T1, j*dt, alpha, beta, eta, sigma, phi);

       change_theta = (beta * (phi - theta(:,j)) - eta^2 * T1t_C) * dt + eta * sqrt(dt) * randn(Nsims,1);
       theta(:,j+1) = theta(:,j) + change_theta;

       change_r = (alpha * (theta(:,j) - r(:,j)) - sigma^2 * T1t_B) * dt + sigma * sqrt(dt) * randn(Nsims,1);
       r(:,j+1) = r(:,j) + change_r;
        
    end
    
    PT1_T22 = exp(-dt*sum(r(:,(T1/dt:end-1)), 2));
    
    value2 = max(PT1_T22 - K, 0).*exp(-dt*sum(r(:,(1:T1/dt)), 2));

    option_value2(i) = mean(value2);
end

plot(scale*K_ori, option_value2);
hold on;

%% Analytical
fun3 = @(s) sigma^2 * ((1 - exp(-alpha * (T1-s)))/alpha - (1 - exp(-alpha * (T2-s)))/alpha).^2;
fun4 = @(s) (eta^2 * (exp(-alpha * (T1-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T1-s))/(beta*(alpha-beta)) - (exp(-alpha * (T2-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T2-s))/(beta*(alpha-beta)))).^2);
fun5 = @(s) (sigma * ((1 - exp(-alpha * (T1-s)))/alpha - (1 - exp(-alpha * (T2-s)))/alpha));
fun6 = @(s) eta * (exp(-alpha * (T1-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T1-s))/(beta*(alpha-beta)) ...
                 -(exp(-alpha * (T2-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T2-s))/(beta*(alpha-beta))));
omega1 = integral(fun3, 0, T1);
omega2 = integral(fun4, 0, T1);
omegasq = omega1 + omega2;
omega = sqrt(omegasq);

dp = (-log(scale) + 0.5 * omegasq)./omega;
dm = (-log(scale) - 0.5 * omegasq)./omega;
anal_price = P0_T2 * (normcdf(dp) - scale .* normcdf(dm));

plot(scale*K_ori,anal_price);
xlabel('Strike');
ylabel('Bond Option Price');
legend('Risk Neutral','Forward Neutral','Analytical');
hold off;

%%
function [A,B,C] = bondABC (T, t, alpha, beta, eta, sigma, phi)
    B = (1 - exp(-alpha * (T-t)))/alpha;
    C = exp(-alpha * (T-t))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T-t))/(beta*(alpha-beta)); % - exp(-beta * T(i))/(alpha-beta) - exp(-beta * T(i))/beta + 1/beta;

    fun1 = @(s) ((1 - exp(-alpha * (T-s)))/alpha).^2;
    omega1 = integral(fun1,t,T);

    fun2 = @(s) (exp(-alpha * (T-s))/(alpha-beta) + 1/beta - alpha * exp(-beta*(T-s))/(beta*(alpha-beta))).^2;
    omega2 = integral(fun2,t,T);

    A = phi * (alpha * (1-exp(-beta*(T-t)))/((alpha-beta)*beta) ...                  
                       - beta * (1-exp(-alpha*(T-t)))/(alpha*(alpha-beta))...
                       - (T-t))...
                       + 0.5 * sigma^2 * omega1 + 0.5 * eta^2 * omega2;
        
end
