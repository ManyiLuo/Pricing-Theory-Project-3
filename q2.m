%% Question 2 (by group 3)
r0 = 0.02;
alpha = 3;
sigma= 0.01;
theta0 = 0.03;
beta = 1;
phi = 0.05;
eta = 0.005;

T = linspace(0,10,100);
t = 0.01;
Nsims = 1000;

yield = zeros(size(T,1), 1);
error = zeros(size(T,1), 1);

%% Simulation

%simulate yield for each T
for i = 1:size(T,2)
    price = zeros(Nsims,1);
    steps = round(T(i)/t);
    %simulate different paths of r and get the average bond price
    for k = 1:Nsims
        theta = zeros(steps,1);
        r = zeros(steps,1);
        theta(1) = theta0;
        r(1) = r0;
        integ = r(1) * t;
       % integ = 0;

        %simulate theta and r using Euler Scheme
        for j = 2:steps

           change_theta = beta * (phi - theta(j-1)) * t + eta * sqrt(t) * normrnd(0, 1);
           theta(j) = theta(j-1) + change_theta;

           change_r = alpha * (theta(j-1) - r(j-1)) * t + sigma * sqrt(t) * normrnd(0, 1);
           r(j) = r(j-1) + change_r;
        
          % integ = integ + (r(j-1) + r(j)) * t/2;
           integ = integ + r(j) * t;
        end
    
        price(k) = exp(-integ);
        
    end
    
    
    yield(i) = -log(mean(price))/T(i);
    
    error(i) = std(-log(price*T(i)))/sqrt(Nsims);
    
end
    
plot(T,yield);
errorbar(T, yield,error);

hold on;

%% Q2 Analytical

yield_anal = zeros(size(T,1), 1);

for i = 1:size(T,2)
    B = (1 - exp(-alpha * T(i)))/alpha;
    C = exp(-alpha * T(i))/(alpha-beta) + 1/beta - alpha * exp(-beta*T(i))/(beta*(alpha-beta));
    
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


plot(T, yield_anal,'LineWidth',1);
title('Term Structure of Interest Rates','FontSize',12)
xlabel('maturity time T','FontSize',12);
ylabel('Bond Yield Y(T)','FontSize',12);
legend({'Euler-scheme Simulation','Analytical Formula'},'Location','southeast','FontSize',12);


