%Crystal Wang
%% Question 1
clear;
%simulations:
N = 10000;  %Number of Iterations
get18 = zeros(1,N);     %Initializing for whether single dice can get a 18 array
get18in3 = zeros(1 , N); %Initializing for whether 1 in 3 dice can get a 18 array
stats = zeros(6,N);   %Initializing for the stats for a character          
for i = 1:1:N         %start of the for loop: N iterations
R = unidrnd(6,3,3);    % randomly generate a 3*3 matrix to represent all the values in 3d6
singleR = unidrnd(6,1,3);   %randomly generate 1d6 for part a)
for j = 1:1:6               %start of the for loops to calculate the stats for different character(fred and keene)
rAbility = unidrnd(6,3,3,N);    %using a 3D matrix to generate 6 abilities scores
sr = sum(rAbility,2);           %sum each colum to get the result of 3d6
abilities =(max(sr, [],1));     %each ability takes the highest score among the rolls
stats(j,:) = abilities;         %each column represent a 6ability charater
end                % The end of the loop for generating a character with 6 abilities scores
single = sum(singleR,2);    %sum the row to get the final score for a single die
S = sum(R,2);               %sum the row to get the final score for three dice
score = (max(S,[],1));      %for partb), need to take the maximum score out of the three dice value- for fun method
get18in3(i) = score;        % the resulting array for part b)
get18(i) = single;          % the resulting array for part a)
end                         % the end of the simulation 
%a)
p18 = length(find(get18 == 18))/N;  %probability of any one roll of 3dice, my value is around 0.005
%b)
p18in3 = length(find(get18in3 == 18))/N; %probability of generating 3cores and keeping the highest, using the fun method, my value is around 0.015
%c)
totalFred = sum(stats,1);
pfred = length(find(totalFred == 108))/N;
%if a character is fred, it has a sum of 108 scores, getting a zero
%everytime for my simulation since the probability of becoming a fred
%character is really really unlikely
%d)
count = 0;
for k = 1:1:N
if (stats(:,k)-9 == 0)
        count = count + 1;
end   
end
pkeene = count / N;
%if a character is keene, every single score is 9, the for loop is used to
%make sure each ability is a 9. Even though there are 25 ways of forming a
%9, getting six 9s is still very unlikely.


%% Question 2
N = 1000000;     %Number of iteration
HP = randi([1,4], 1, N);        %Generating hitting points for the trolls
keene = randi([1,2], 2, N); 
FB = sum(keene, 1);             %Generating randomly for keenes fireball spell
%a)
Av_HP = sum(HP)/N;          % The expected (average) value of hitting points is 2.50 around
Av_FB = sum(FB)/N;          % The expected (average) value of FireBall is around 3
kcount = 0;                 %There are only one possibile way to get more than 3 damage for fireball: getting a 4 (two 2s)
for a = 1:1:N
    if FB(a) == 4
        kcount = kcount + 1;
    end
end
pFBgreaterthan3 = kcount / N;% The probability is around 0.25

%b)
%pmfs are shown using a pmf diagram
nHP = [1, 2, 3, 4];
nFB = [2,3,4];
pmf_HP = zeros(1,4);
pmf_FB = zeros(1,3);
for i = 1:1:4
pmf_HP(i) = length(find(HP == i))/N;
end
for j = 1:1:3
pmf_FB(j) = length(find(FB == (j+1)))/N;
end
figure;
subplot(2,2,[1 2]);
stem(nHP, pmf_HP);
title("Probability mass functions of the hit points trolls have")
xlabel("The amount of hit points");
ylabel("Probability");
subplot(2,2,[3 4]);
stem(nFB, pmf_FB);
title("Probability mass functions of the amount of damage the FIREBALL does");
xlabel("The amount of damage");
ylabel("Probability");

%c)
HP6 = randi([1,4], 6, N);   %generating 6 trolls with N values
highest = zeros(1,N);
count6 = 0;
for b = 1:1:N
    highest(b) = max(HP6(:,b));
if (FB(1,b) - highest(b) >= 0)      %if FIREBALL can beat the max hitting points of the troll, it can beat any other trolls
   count6 = count6 + 1;
end
end
p_keeneslays = count6/ N;   % The probability of keenslays is around 0.344

%d)
trollDamage = HP6 - FB;     %The remaining hitting points troll has
countAlive = 0;
nDamage = 0;
oneSurvive = 0;
totalHP = 0;
for i = 1:1:N
    for j = 1:1:6       %This loop is to make sure that only one troll survive
        if (trollDamage(j,i) >=1)
            countAlive = countAlive + 1;    %if troll has higher points than fireball, it survives, and the count should increment by 1
            nDamage = nDamage + trollDamage(j,i);   % this is the sum of the total amount of damage the simulation has 
        end
    end
     if (countAlive == 1)
         oneSurvive = oneSurvive + 1;   % if there is only one troll survive, the total amount of HP is calculated
         totalHP = totalHP + nDamage;
     end
        countAlive = 0;
        nDamage = 0;
end
    expectedHP = totalHP/oneSurvive;    % expected value of hit points that the remainig troll has is equal to the total remaining amound / the count of 1survivor
    % it is around 1.15

%e)
shedjTotDamage = 0;
tuitionSword = sum(randi([1,6],2, N), 1);   %generating the sword
tenureHam = randi([1,4], 1,N);          %generating the hammer 
for k = 1:1:N
    if(randi([1,20],1,1) >= 11)     % if it is greater than 11, it can proceed to attack
        shedjTotDamage = shedjTotDamage + tuitionSword(k);  %adding the total amount of damage shedjam did to keene for the first round of attack
        if (randi([1,20],1,1) >= 11) % after being attacked by sword, if its creater than 11, one can use the hammer now 
           shedjTotDamage = shedjTotDamage + tenureHam(k); % adding up to the total amount of damage
        end      
    end
    
end
   expectedShedjDMG = shedjTotDamage / N;  %the expected value is around 4.12

