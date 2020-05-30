% Crystal Wang Stoch Project4

%% Part1
clear;
 N = 1000; % Number of Iterations
 p0 = 0.8;
 p1 = 0.2;
 A = 1;
 target = randi(5,1,N); 
 target(find(target == 1)) = A;
 target(find(target  > 1)) = 0;
 SNR = [1,2,5,8];
 etaM = 0.01:0.01:10;
 sigZ = [30,40,50,40];
 Yz = zeros(1,N);
 figure
for k = 1:1:length(SNR)
 sigma = sqrt(A/SNR(k));
 X = normrnd(0, sigma, 1,N);
 Z = normrnd(0,sigZ(k),1,N);
 Y = X + target;
 index1 = find(target == A);
 index2 = find(target == 0);
 Yz(index1) = target(index1) + X(index1);
 Yz(index2) = target(index2) + Z(index2);
 
 p_ygivenH1 = normpdf(Y,A,sigma);
 p_ygivenH0 = normpdf(Y,0,sigma);
 
 pz_ygivenH1 = normpdf(Y,A,sigma);
 pz_ygivenH0 = normpdf(Y,A,sigZ(k));
 resultZ = zeros(1,N);
 resultZ(find((pz_ygivenH1./pz_ygivenH0) > (p0/p1))) = 1;
 resultZ(find((pz_ygivenH1./pz_ygivenH0) < (p0/p1))) = 0;
 
 result = zeros(1,N);
 result(find((p_ygivenH1./p_ygivenH0) > (p0/p1))) = 1;
 result(find((p_ygivenH1./p_ygivenH0) < (p0/p1))) = 0;
 
 count = 0;
 countz = 0;
 c0 = 0;
 c1 = 0;
 for i = 1:1:N
     a = result(i) - target(i);
     az = resultZ(i) - target(i);
     if abs(a) == 1 
         count = count + 1;
        
     end   
     if abs(az) == 1
         countz = countz + 1;
     end
     if target(i) == 0 
            c0 = c0 + 1;
     end
     if target(i) == 1
            c1 = c1 + 1;
     end
 end
pError = count/N;
pErrorZ = countz/N;
%b)
etas = (p0/p1) * etaM;
cD = zeros(1,length(etaM));
cF = zeros(1,length(etaM));
cDz = zeros(1,length(etaM));
cFz = zeros(1,length(etaM));
resultB = zeros(1,N);
resultBz = zeros(1,N);
for j = 1:1:length(etaM)
    resultB(find(p_ygivenH1./p_ygivenH0 > etas(j))) = 1;
    resultB(find(p_ygivenH1./p_ygivenH0 < etas(j))) = 0;
    resultBz(find(pz_ygivenH1./pz_ygivenH0 > etas(j))) = 1;
    resultBz(find(pz_ygivenH1./pz_ygivenH0 < etas(j))) = 0;
    for i = 1:1:N
        if resultB(i) - target(i) == 1
            cF(j) = cF(j) + 1;
        end
        if resultBz(i) - target(i) == 1
            cFz(j) = cFz(j) + 1;
        end
        if resultB(i) + target(i) == 2
            cD(j) = cD(j) + 1;
        end
        if resultBz(i) + target(i) == 2
            cDz(j) = cDz(j) + 1;
        end
        
    end
end
f = cF/c0;
d = cD/c1;
fz = cFz/c0;
dz = cDz/c1;
plot(f, d);
hold on;
plot(fz,dz);
end
plot(f(10),d(10), '*');
title("ROC with different SNR");
xlabel("pF");
ylabel("pD");
legend("SNR=1", "ratio = 1/30", "SNR=2", "ratio = 1/20", "SNR = 5",  "ratio = 1/10", "SNR = 8", "ration = 1/3","when missing target is 10 times worse than false alarm");

%d)
p1d = linspace(0,1,100);
resultD = zeros(1,N);
cM = zeros(1,N);
cFf = zeros(1,N);
expCost = zeros(1,length(p1d));
sigma = sqrt(A/8);

for a = 1:1:length(p1d)
    target = randi(100,1,N);
    target(find(target <= p1d(a)*100)) = A;
    target(find(target > p1d(a)*100)) = 0;
    X = normrnd(0, sigma, 1, N);
    Y = target + X;
    p_ygivenH1 = normpdf(Y,A,sigma);
    p_ygivenH0 = normpdf(Y,0,sigma);
    resultD(find((p_ygivenH1./p_ygivenH0) > 0.1*(1-p1d(a))/p1d(a))) = 1; 
    resultD(find((p_ygivenH1./p_ygivenH0) < 0.1*(1-p1d(a))/p1d(a))) = 0;
     for i = 1:1:N
        if resultD(i) - target(i) == -1
            cM(a) = cM(a) + 1;
        end
        if resultD(i) - target(i) == 1
            cFf(a) = cFf(a) + 1;
        end
     end
     expCost(a) = 10*cM(a)/N + 1*cFf(a)/N;
end

figure
plot(p1d,expCost);
xlabel('p1');
ylabel('Expected Cost');
title('When SNR = 8 the Expected Cost vs Priori Target Present Probabilities');


%% part 2
clear;
data = load('Iris.mat');

trainingF = zeros(75,4);
testF = zeros(75,4);
shuffledInd = randperm(size(data.features,1));%Shuffle the data so that the training and the testing set can be selected at random
shuffledF = data.features(shuffledInd,:);
shuffledL = data.labels(shuffledInd,:);
trainingF = shuffledF(1:75, :);
testF = shuffledF(76:end,:);
trainingL = shuffledL(1:75);
testL = shuffledL(76:end);

p1 = trainingF(trainingL == 1,:);
p2 = trainingF(trainingL == 2,:);
p3 = trainingF(trainingL == 3,:);

%Variances
cov1 = cov(p1);
cov2 = cov(p2);
cov3 = cov(p3);

%Means
mu1 = [mean(p1(:,1)),mean(p1(:,2)),mean(p1(:,3)),mean(p1(:,4))];
mu2 = [mean(p2(:,1)),mean(p2(:,2)),mean(p2(:,3)),mean(p2(:,4))];
mu3 = [mean(p3(:,1)),mean(p3(:,2)),mean(p3(:,3)),mean(p3(:,4))];
likelihood = [mvnpdf(testF,mu1,cov1),mvnpdf(testF,mu2,cov2),mvnpdf(testF,mu3,cov3)];

[~,result] = max(likelihood, [], 2);
error = sum(testL~=result);
pError = error/size(testL,1)
confutionM = confusionmat(testL,result)