%Crystal's Stoch Project 2

%% Scenario 1
clear;
N = 100000; % Number of iterations

Y = -1 + (1 - (-1)).* rand(1,N); % Generating N Ys that has a uniform distribution U~(-1,1)
W = -2 + (2 - (-2)).* rand(1,N); % Generating N Ws that has a uniform distribution U~(-2,2)

X = Y + W; %Definition of X

bEstimator = zeros(1,N); % The initialization of the bayes estimator
bSquareError = zeros(1,N); % The initialization of the square error.
fx = zeros(1,N); % The initialization of the marginal pdf
for i = 1:1:N % For each set of x and y
    if -3 <= X(i) < -1
        bEstimator(i) = 0.5 + 0.5 * X(i); %bayes Estimator equations
        fx(i) = (3+ X(i))/8;
    end
    if -1 <= X(i) < 1
        bEstimator(i) = 0;
        fx(i) = 0.25;
    end
    if 1 <= X(i) <= 3
        bEstimator(i) = -0.5 + 0.5 * X(i);
        fx(i) = (3- X(i))/8;        
    end
    % all the if statements contains the corresponding formula of
    % bayesEstimator according to its range of x
    bSquareError(i) = (Y(i) - bEstimator(i)).^ 2 ; % Computing the bayes square error here
    
end

meanbSquareError = sum(bSquareError)/N; % Summing up all the square error, taking the average of it 
result = zeros(1,N);% Initialization of each individual result

% For Linear MMSE
meanY = sum(Y)/N; % Expected Value of Y
meanX = sum(X)/N;  % Expected Value of X
coVarXY = sum(Y.* X)/N - (meanY * meanX); % The covariance of XY
varX = sum(X .* X)/N - (meanX).^2 ; % variance of X

lEstimator = zeros(1,N); % Initialization of linear Estimator
lSquareError = zeros(1,N); %Initialization of the linear square error
for i = 1:1:N % Calculating MMSE for each Y value for Bayes and calulating the LInear MMSE
    result(i) = meanbSquareError * fx(i);
    lEstimator(i) = meanY + coVarXY/varX * (X(i) - meanX);
    lSquareError(i) = (Y(i) - lEstimator(i)).^2;
end
bMSE = sum(result) / N;% The average of all the results is the expected value of the MMSE of
% overall x
lMSE = sum(lSquareError)/N;
TheoreticalMSE = 0.25; % Theoretical value of MSE
T = table(TheoreticalMSE, bMSE, lMSE) % View All the results in a table

%% Scenario 2
clear;
O = 40; % Number of observations
thstdY = [0.2,0.5,0.8]; % Theoretical std of Y
thstdR = [0.3,0.4,0.7]; % Teoretical std of R
N = 10000; % Number of Iterations
thvarY = zeros(1,3);% Initialization for theoretical variance of Y
thvarR = zeros(1,3);% Initialization for theoretical variance of R
thmeanY = zeros(1,3);% Initialzation for theoretical mean of Y
varY = zeros(1,3);% Initialization for variance of Y
varR = zeros(1,3);% Initialization for variance of R
meanY = zeros(1,3);% Initialzation for mean of Y
thYl = zeros(O,N); % Initialization of  theoreticalYl for different numbers of ovservations
Yl = zeros(O,N);% Initialization of Yl for different numbers of ovservations
% Initializaiton fo  theoretical MMSE for different numbers of observations
thMMSE = zeros(3,O); % Initialization of theoreticalMSE
MMSE = zeros(3,O);  % Initialization of MSE
thSquareError = zeros(1,N); %Initialzation for theoretical Square Error for each set of Y
SquareError = zeros(1,N); %Initialzation for Square Error for each set of Y
X = zeros(O,N); % Initializatin of the output X
CXY = zeros(O,1);   % Initialization of the CXY
for a = 1:1:3 % For loop for different variance of Y and R
Y = normrnd(1, thstdY(a),1,N); % Generating Y with normal distribution mean of 1, standard deviation of 0.5
R = normrnd(0, thstdR(a), O, N); % Generating R with normal distribution wiht mean of 0, standard deviation of 0.8
%The first row represents R1, the second row represents R2, the third row
%represents R3....
thmeanY(a) = 1; % The theoretical mean of Y
%meanY(a) =  sum(Y)/N;
thvarY(a) = thstdY(a) * thstdY(a); % Variance of theoretical Y
varY(a) = sum(Y.*Y,'All')/N - (meanY(a)).^2; % Variance of actual Y
thvarR(a) = thstdR(a) * thstdR(a); % Variance of theoretical R
varR(a) = sum(R.*R,'All')/(N*O); % mean of R is 0
    for i = 1:1:O %For loop for each iteration of observation
        X(i,:) = Y + R(i,:);% How to find each Xi
        thYl(i,:) = (1/((i*thvarY(a)) + thvarR(a))) * (thvarR(a)*thmeanY(a)); % Theoretical Yl using theoretical values
        covXY = cov(Y,X(i,:)); % Find the covariance of CXY
        CXY(i) = covXY(2); % CXY is the diagnal of covXY
       for k = 1:1:i %For loop for each varY*Xi term
            for j = 1:1:N % For loop for each Y  
                thYl(i,j) = thYl(i,j)+ (1/((i*thvarY(a))+thvarR(a))) * (thvarY(a)*X(k,j)); % Final theoretical yl value
                thSquareError(i,j) = (Y(j) - thYl(i,j)).^2; % Theoretical Square Error
            end 
        end
         thMMSE(a,i) = sum(thSquareError(i,:))/N;    % Theoretical MMSE   
    end
        meanX = mean(X,2); % The mean of X
        meanY(a) = mean(meanX); % taking the mean of X as the mean of y
        a0 = zeros(1,O) + meanY(a); % Initializaiont of ao
      
    for z = 1:1:O   % For loop for each observations
        CXX = cov(transpose(X(1:z,:))); % CXX 
        A = (CXX^(-1))*CXY(1:z); % A
            for f = 1:1:z % To get Yl for each iteration for each observation
            Yl(z,:) = Yl(z,:) + A(f)*X(f,:);
            a0(z) = a0(z) - meanX(z) * A(f);% Final values for a0
            end
            Yl(z,:)= Yl(z,:) + a0(z); % Final value for Yl
     end
    
    for c = 1:1:O % For loop for each observations
         for d = 1:1:N % For each iterations
            SquareError(c,d) = (Y(d) -  Yl(c,d)).^2;   % Square Error
         end
        MMSE(a,c) = sum(SquareError(c,:))/N; % MMSE
        
     end
Yl = zeros(O,N); % Make sure Yl is 0 after each variance
end

figure
for e = 1:1:3
    plot(1:1:O, MMSE(e,:));
    hold on;
    plot(1:1:O, thMMSE(e,:));
    hold on;
end
title("MMSE and Number of Observations using meanY = 1, mean R = 0");
xlabel("Number of Observatons");
ylabel("MMSE");
legend('varY = 0.04, varR = 0.09 ','varY = 0.25, varR = 0.16','varY = 0.64, varR = 0.49','empirical varY = 0.04,, varR = 0.09 ','empirical varY = 0.16, varR = 0.16','empirical varY = 0.64, varR = 0.49');







