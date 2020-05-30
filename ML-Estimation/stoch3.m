% Crystal Wang Stoch3 project

%% Part 2
clear;
O = 20; % #of observations
N = 1000; % # of iterations
mu = 2; % Let the mean of the exponential distribution be 5
lamb = 1/mu;
expRV = exprnd(mu,[O,N]); % rv for exponential distribution
rayRV = raylrnd(lamb,O,N); % rv for rayleigh distribution
lambHat = 0.1:0.05:1; % Different values of lambda hat, given lambda is 0.5
lh = length(lambHat); % Te length of the chosen lambda hat vector
expM = zeros(O,N);    % The initialization of exponential distribution lambda maximum likelihood function for each iteration
raylM = zeros(O,N);   % The initialization of rayleigh distribution lambda maximum likelihood function for each iteration
for i = 1:1:N % For each iteration 
    expL = ones(O,lh); % The initialization of exponential distribution likelihood function
    rayL = ones(O,lh); % The initialization of rayleigh distribution likelihood function
    for j = 1:1:lh  % Loop for each lambda hat chosen
        for k = 1:1:O  % Loop for each observation
            for p = 1:1:k % Loop to multiply each pdf before the current observation
                expL(k,j) =  expL(p,j) * exppdf(expRV(p,i), 1/lambHat(j)); % Using the formula to calculate the likelihood function for each observation
                rayL(k,j) =  rayL(p,j) * raylpdf(rayRV(p,i), lambHat(j));% Using the formula to calculate the likelihood function for each observation
            end
        end
    end
    [~,expI] = max(expL, [], 2); %The index of exponential distribution maximum
    [~,rayI] = max(rayL, [], 2); %The index of Rayleigh distribution maximum
    expM(:,i) = lambHat(expI); % Assigning the lambda hat that maximizes the function 
    raylM(:,i) = lambHat(rayI); % Assigning the lambda hat that maximizes the function    
end

expMSE = mean((expM-lamb).^2,2); %Calculating MSE for exponential distribution
rayMSE = mean((raylM-lamb).^2,2); % Calculating MSE for rayleigh distribution
expBias = mean(expM - lamb,2); % Calculating Bias for exponential distribution
rayBias = mean(raylM - lamb,2); %Calculating Bias for rayleigh distribution 
expVar = mean((expM - mean(expM)).^2,2);%Calculating variance for exponential distribution
rayVar = mean((raylM - mean(raylM)).^2,2);% Calculating variance for rayleigh distribution
figure;
plot(1:1:O,expMSE);
hold on;
plot(1:1:O,rayMSE);
title('MSE with respect to number of Observations');
xlabel('Number of Observations');
ylabel('MSE');
legend('MSE for Exponential Distribution','MSE for Rayleigh Distribution');

figure;
plot(1:1:O, expBias);
hold on;
plot(1:1:O, expVar);
hold on;
plot(1:1:O, rayBias);
hold on;
plot(1:1:O, rayVar);
title('Bias and Variance for Exponential and Rayleigh Distributions');
xlabel('Number of Observations');
ylabel('Bias/Variance Value');
legend('Bias for Exponential Distribution','Variance for Exponential Distribution', 'Bias for Rayleigh Distribution','Variance for Rayleigh Distribution');

%% Part 3
clear;
X = importdata('data.mat');
L =length(X);
expLamb = L / sum(X); 
rayLamb = (1/(2*L) * sum(X.^2))^0.5;
%Log likelihood
expL = L*log(expLamb) - expLamb * sum(X)

rayL = prod(X) - 2*L*log(rayLamb) - 0.5*sum((X.*X)/rayLamb^2)

figure;
histogram(X);
xlabel('Histogram of the data');
ylabel('value');

%Thus we know that rayleigh is better in this case 

