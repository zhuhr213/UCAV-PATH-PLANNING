clear
clc



global radar1
global radar2
global R   % 这个好像是每个rader的半径
% radar1 = [350 105 305 105 175 245 415 480 40 470];
% baili = size(radar1,2);
% radar2 = [200 0 -150 110 110 110 0 110 100 -50];
% R = [140 70 150 30 25 25 25 90 40 30];

% 破案了！rader1是横坐标，2是纵坐标
% 我是铁憨憨
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
baili = size(radar1,2);  % 第二个参数表示第几维度，matlab的下标好像都是从1开始的
                         % size(A)会返回A的每个维的size，组成一个向量
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];
 % R(1:5) = 10

NP=40; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
maxCycle=1000; %/*The number of cycles for foraging {a stopping criteria}*/
limit=0.1*maxCycle; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
D= 60; % 维度，分隔的块数
ub=ones(1,D).*50; %/*lower bounds of the parameters. */
lb=ones(1,D).*-50;%/*upper bound of the parameters.*/ //注释写反了
runtime=1;%/*Algorithm can be run many times in order to see its robustness*/
GlobalMins=zeros(1,runtime);
Range = repmat((ub-lb),[FoodNumber 1]);  % 以该矩阵为单位，生成[m n]的新矩阵
Lower = repmat(lb, [FoodNumber 1]);
Foods = rand(FoodNumber,D) .* Range + Lower;  % 每一行代表一个食物，列是表示维数
            % 生成一个这么大的随机数矩阵，代表随机的待定食物

for r=1:runtime  % 循环执行算法
for i = 1:FoodNumber
    ObjVal(i) = calcu(Foods(i,:));  % 每次它会自动扩充一列
end
Fitness = calculateFitness(ObjVal); % 利用objval计算可行度
trial=zeros(1,FoodNumber);
BestInd=find(ObjVal==min(ObjVal)); % 找到目前最好的解的位置
BestInd=BestInd(end);  % 如果存在多个最优解，选最后一个
GlobalMin = ObjVal(BestInd);  % 确定一个全局最优
GlobalParams=Foods(BestInd,:);  % 第BestInd个食物最厉害
iter=1;



while ((iter <= maxCycle)),  % 开始迭代，迭代次数为maxCycle
    for i=1:(FoodNumber)
        Param2Change=fix(rand*D)+1; % fix函数用来向零方向取证，rand生成（0，1）的随机数
        % 控制该变量在（1，30）内，根据论文，每次只修改一个seg（随机）
        neighbour=fix(rand*(FoodNumber))+1;     
            while(neighbour==i)  % 利用p622的公式（利用隔壁蜜蜂的位置）更新位置，保证k和i不等
                neighbour=fix(rand*(FoodNumber))+1;
            end;  
       sol=Foods(i,:); % 本次处理第i个食物，一个食物代表一种解法
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2; 
                                                                                                   % 控制这个区间在（-1，1）
        %  修正
        ind=find(sol<lb);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        ind=find(sol>ub);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        %
        ObjValSol = calcu(sol);
        FitnessSol=calculateFitness(ObjValSol); % 计算新位置的适应值，并于先前的比较
       if (FitnessSol>Fitness(i))   % 如果新的值更好，就替换
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1;    % 如果原来的更好，就另trial自增
       end;
     end;
     
prob=(0.9.*Fitness./max(Fitness))+0.1; % 修正每种解的概率
i=1;
t=0;

while(t<FoodNumber)  % 这一段没有的话，也能有不错的效果。这一段是模拟观察蜂，决定是否利用该解。
                     % 变量设置中，观察蜂跟食物的数量是一样的
    if(rand<prob(i))  % 使用该解，并再次寻找更优解
        t=t+1;
        Param2Change=fix(rand*D)+1;
        neighbour=fix(rand*(FoodNumber))+1;     
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        %
        ind=find(sol<lb);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        ind=find(sol>ub);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        %
        ObjValSol = calcu(sol);
        FitnessSol=calculateFitness(ObjValSol);
       if (FitnessSol>Fitness(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) % 如果已经遍历了一圈所有食物，但还有未修正的，就继续
        i=1;
    end;   
end; 

         ind=find(ObjVal==min(ObjVal));
         ind=ind(end);
         if (ObjVal(ind)<GlobalMin)
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:); 
         end;
         % 更新最小值
         
         

ind=find(trial==max(trial));  % 如果trial有超过limit，就修复
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    sol=(ub-lb).*rand(1,D)+lb;
    ObjValSol = calcu(sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;
fprintf('iteration = %d ObjVal=%g\n',iter,GlobalMin);  % 记录迭代次数
aaaaa(iter) = GlobalMin;
iter=iter+1;
end % End of ABC
storeer(r,:) = aaaaa;
end  % 这个是循环执行算法的最外层循环，默认只执行一次


% 
% figure (1)
% a = mean(storeer);
% plot(a)

figure (2)
hold on
plot(0,0,'k*');
plot(500,0,'ks');

for i = 1:baili  % 威胁的个数
hold on
plot(radar1(i),radar2(i),'ko');
cir_plot([radar1(i),radar2(i)],R(i)); % 参数方程画圆
end
legend('starting point','target point','threat center')

axis equal  % 横轴纵轴长度单位相同
axis([-100,620,-400,400])
for i = 1:(D-1)
    hold on
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalParams(i),GlobalParams(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,GlobalParams(1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[GlobalParams(D),0],'LineWidth',2);