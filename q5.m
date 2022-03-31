r0 = 0.02;
alpha = 3;
sigma= 0.01;
theta0 = 0.03;
beta = 1;
phi = 0.05;
eta = 0.005;

scale = linspace(0.95,1.05,9);

dt = 0.01;
Nsims = 100;

tau = [3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6];
dtau = 0.25;
T0_A = size(tau,2);
T0_B = size(tau,2);
T0_C = size(tau,2);
P0_T = size(tau,2);

for i = 1:size(tau,2)
    [T0_A(i), T0_B(i), T0_C(i)] = bondABC (tau(i), 0, alpha, beta, eta, sigma, phi);
    P0_T(i) = exp(T0_A(i) - T0_B(i) * r0 - T0_C(i) * theta0);
end

%today's swap rate
F = (P0_T(1) - P0_T(end))/(sum(P0_T(2:end))*dtau);

%% Calculate V0
imp_vol = size(scale,2);
for k = 1:size(scale,2)
    theta = zeros(Nsims, tau(end)/dt);
    r = zeros(Nsims, tau(end)/dt);

    theta(:,1) = theta0;
    r(:,1) = r0;

    %simulate theta and r using Euler Scheme
    for j = 1:tau(end)/dt-1

       change_theta = beta * (phi - theta(:,j)) * dt + eta * sqrt(dt) * randn(Nsims,1);
       theta(:,j+1) = theta(:,j) + change_theta;

       change_r = alpha * (theta(:,j) - r(:,j)) * dt + sigma * sqrt(dt) * randn(Nsims,1);
       r(:,j+1) = r(:,j) + change_r;
    end

    T3_A = size(tau,2);
    T3_B = size(tau,2);
    T3_C = size(tau,2);
    P3_T = zeros(Nsims,size(tau,2));

    for i = 1:size(tau,2)
        [T3_A(i), T3_B(i), T3_C(i)] = bondABC (tau(i), 3, alpha, beta, eta, sigma, phi);
        P3_T(:,i) = exp(T3_A(i) - T3_B(i) * r(:,3/dt) - T3_C(i) * theta(:,3/dt));
    end

    S3 = (P3_T(:,1) - P3_T(:,end))./(sum(P3_T(:,(2:end)),2)*dtau);

    A0 = sum(P0_T(2:end))*dtau;
    
    A3 = sum(P3_T(:,(2:end)),2)*dtau;

    expectation = A3.*max(S3-scale(k)*F,0).* exp(-dt*sum(r(:,(1:3/dt)), 2));
    
    
    V0 = mean(expectation); %A0* mean(max(S3-scale(k)*F,0));

    syms i
    eqn = A0*F*normcdf((-log(scale(k))+0.5*i^2)/i) - scale(k)* normcdf((-log(scale(k))-0.5*i^2)/i) == V0;
    omega = vpasolve(eqn,i);
    imp_vol(k) = omega/sqrt(3);
end

plot(F*scale,imp_vol/100);
xlabel('Strike');
ylabel('Black Implied Volatility');


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

