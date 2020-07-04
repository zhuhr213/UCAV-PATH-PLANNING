% 蜘蛛猴算法（Spider Monkey Optimization）
% 它有点像是ABC与GWO的结合
% 当发生停滞的时候（ABC中的trial），种群就会分裂成小种群，小种群分散寻找
% 每个小种群都有一个LocalBest，全局有一个GlobalBest
% pr，perturbation rate，控制目前为止的变异率，一般取[0.1,0.8]

% 算法流程：
% 1. 初始化
% 2. 计算适应值
% 3. 选择全局最优和局部最优
% while 迭代次数或结果条件：
%     1. 更新每个解的每个维，利用每个小种群的局部最优和自身，计算适应值，贪心法更新
%     （如果概率大于pr，就依局部最优更新，否则保留原解）
%     2. 利用fitness计算每个解被选择的概率
%     3. 对每个解，依概率选择是否更新，随机更新某个维，依照全局最优和自身
%     4. 更新全局最优和每个小群体的局部最优，贪心法
%     5. 如果小群体的局部最优迟迟未更新，将该小群体所有成员全部更新位置
%     （看概率是否大于pr，若大于则直接生成新随机解，若小于则根据已有最优更新）
%     （变异率）
%     6. 如果全局最优迟迟没有更新，就将现有群体群划分为更多的小群体（2,3,4....）
%     （群体数有上限，当达到上限时就把全部小群体重新合并为一个），然后更新各局部最优
%     （不要忘记将更新计数器置零，如果迟迟没有更新）
% end

clear
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

maxCycle = 200;
D = 60;
N = 60;
MG = N/4;  % 最大小群体数
GlobalLeader = zeros(1,D);
LocalLeader = zeros(MG,D);  % 存放每个小群体的局部最优
GlobalLeaderLimit = N;  % 它应该在[N/2,2*N]之间
LocalLeaderLimit = D*N; % 文章里这么设置的
GlobalLeaderCount = 0;
LocalLeaderCount = zeros(1,MG);
GroupNumber = 1;  % 他也没说最开始整几个猴群，默认一个

ub=ones(1,D).*50;  % 上下界
lb=ones(1,D).*-50;

pr = 0.8;  % 一般取[0.1,0.8]，经实验这东西非常影响运行速度

SM = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;  % 生成每个猴的位置

for i = 1:N
    ObjVal(i) = calcu(SM(i,:));  % 计算objval值（代价函数），越小越好
end
Fitness = calculateFitness(ObjVal); % 计算Fitness，Fitness是拟合度，越高越拟合

ind = find(Fitness == max(Fitness));
ind = ind(end); 
GlobalLeader = SM(ind,:);  % 找到全局最优
GlobalLeaderObjVal = ObjVal(ind);

LocalLeaderObjVal = zeros(1,MG);
for i = 1:GroupNumber  % 暂时先不处理除不尽的情况
    ind = find(Fitness == max(Fitness(N/GroupNumber*(i-1)+1:N/GroupNumber*i)));
    ind = ind(end);
    LocalLeader(i,:) = SM(ind,:);  % 为每一个小群体找到局部最优
    LocalLeaderObjVal(i) = ObjVal(ind);
end

for iter = 1:maxCycle  % 开始迭代
    
    for k = 1:GroupNumber  % 遍历每个小群
        for i = N/GroupNumber*(k-1)+1:N/GroupNumber*k  % 遍历群里的每个个体
            for j = 1:D  % 遍历每一维
                if (rand>=pr)
                    range = -(N/GroupNumber*(k-1)+1 - N/GroupNumber*k);
                    neighbour = fix(N/GroupNumber*(k-1)+1 + rand*range); % 在该群体内生成随机数
                    while(neighbour==i) 
                        neighbour = fix(N/GroupNumber*(k-1)+1 + rand*range);
                    end
                    temp = SM(i,:);
                    temp(j) = SM(i,j) + rand*(LocalLeader(k,j)-SM(i,j)) + (rand-0.5)*2*(SM(neighbour,j)-SM(i,j));
                                          %
               if(temp(j)<lb(j)||temp(j)>ub(j))
                   temp(j) = rand * (ub(j)-lb(j)) + lb(j);
               end
                     %
                    tempObjVal = calcu(temp);
                    if (tempObjVal<ObjVal(i))
                        SM(i,j) = temp(j);
                        ObjVal(i) = tempObjVal;
                        Fitness(i) = calculateFitness(ObjVal(i));
                    end
                end
            end
        end
    end
    
    % 计算概率
    fitnessSum = sum(Fitness);
    for i = 1:N
        prob(i) = Fitness(i)/fitnessSum;
    end
    
    % 依概率prob和全局最优更新每一个
    
    for k = 1:GroupNumber  % 遍历每个小群体
        for i = N/GroupNumber*(k-1)+1:N/GroupNumber*k % 遍历每个个体
           if (prob(i)>rand)
               dim2change = fix(rand*D)+1;
               range = -(N/GroupNumber*(k-1)+1 - N/GroupNumber*k);
               neighbour = fix(N/GroupNumber*(k-1)+1 + rand*range); % 在该群体内生成随机数
               while(neighbour==i) 
                    neighbour = fix(N/GroupNumber*(k-1)+1 + rand*range);
               end
               temp = SM(i,:);
               temp(j) = SM(i,j) + rand*(GlobalLeader(j) - SM(i,j)) + (rand-0.5)*2*(SM(neighbour,j)-SM(i,j));
                       %
               if(temp(j)<lb(j)||temp(j)>ub(j))
                   temp(j) = rand * (ub(j)-lb(j)) + lb(j);
               end
                     %
               tempObjVal = calcu(temp);
               if (tempObjVal<ObjVal(i))
                   SM(i,j) = temp(j);
                   ObjVal(i) = tempObjVal;
                   Fitness(i) = calculateFitness(ObjVal(i));
               end
           end
        end
    end
    
    % 更新全局最优和各个小群体的最优
    ind = find(Fitness == max(Fitness));
    ind = ind(end); 
    if (ObjVal(ind) < GlobalLeaderObjVal)
        GlobalLeader = SM(ind,:);  % 新全局最优符合条件
        GlobalLeaderObjVal = ObjVal(ind);
        GlobalLeaderCount = 0;
    else
        GlobalLeaderCount = GlobalLeaderCount+1;
    end

    for i = 1:GroupNumber  % 更新局部最优
        ind = find(Fitness == max(Fitness(N/GroupNumber*(i-1)+1:N/GroupNumber*i)));
        ind = ind(end);
        if (ObjVal(ind) < LocalLeaderObjVal(i))
            LocalLeader(i,:) = SM(ind,:);  % 新全局最优符合条件
            LocalLeaderObjVal(i) = ObjVal(ind);
            LocalLeaderCount(i) = 0;
        else
            LocalLeaderCount(i) = GlobalLeaderCount(i)+1;
        end
    end
    
    % 如果有一个小群体长时间没更新，就重新生成他们，要么完全随机，要么根据全局最优，并不是回归正体
    for k = 1:GroupNumber
        if (LocalLeaderCount>LocalLeaderLimit)
            LocalLeaderCount(k) = 0;
        for i = N/GroupNumber*(k-1)+1:N/GroupNumber*k % 遍历每个小群个体
            for j = 1:D
                if(rand>pr)
                    SM(i,j) = lb(j) + rand*(ub(j)-lb(j));
                else
                    SM(i,j)=SM(i,j)+rand*(GlobalLeader(j)-SM(i,j))+rand*(SM(i,j)-LocalLeader(k,j));
                end
                                      %
               if(SM(i,j)<lb(j)||SM(i,j)>ub(j))
                   SM(i,j) = rand * (ub(j)-lb(j)) + lb(j);
               end
                     %
            end
            ObjVal(i) = calcu(SM(i,:));
            Fitness(i) = calculateFitness(ObjVal(i));
        end
        ind = find(Fitness == max(Fitness(N/GroupNumber*(k-1)+1:N/GroupNumber*k)));
        ind = ind(end);
        LocalLeader(k,:) = SM(ind,:);  % 为小群体找到局部最优
        LocalLeaderObjVal(k) = ObjVal(ind);
        
            ind = find(Fitness == max(Fitness));
    ind = ind(end); 
    if (ObjVal(ind) < GlobalLeaderObjVal)
        GlobalLeader = SM(ind,:);  % 新全局最优符合条件
        GlobalLeaderObjVal = ObjVal(ind);
        GlobalLeaderCount = 0;
    else
        GlobalLeaderCount = GlobalLeaderCount+1;
    end
        end
    end
    
    if (GlobalLeaderCount>GlobalLeaderLimit)
        GlobalLeaderCount = 0;
        if (GroupNumber < MG)
            GroupNumber = GroupNumber + 1;
        else
            GroupNumber = 1;
        end
        for i = 1:GroupNumber  % 更新局部
            ind = find(Fitness == max(Fitness(N/GroupNumber*(i-1)+1:N/GroupNumber*i)));
            ind = ind(end);
            LocalLeader(i,:) = SM(ind,:);  % 为每一个小群体找到局部最优
            LocalLeaderObjVal(i) = ObjVal(ind);
        end
    end
    fprintf('iteration = %d ObjVal=%g\n',iter,GlobalLeaderObjVal);  % 记录迭代次数
    aaaaa(iter) = GlobalLeaderObjVal;    
end

storeer = aaaaa;

% figure (1)
% a = mean(storeer);
% plot(a)

figure (2)
hold on
plot(0,0,'k*');
plot(500,0,'ks');

for i = 1:size(radar1,2)  % 威胁的个数
    hold on
    plot(radar1(i),radar2(i),'ko');
    cir_plot([radar1(i),radar2(i)],R(i)); % 参数方程画圆
end
legend('starting point','target point','threat center')

axis equal  % 横轴纵轴长度单位相同
axis([-100,620,-400,400])
for i = 1:(D-1)
    hold on
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalLeader(i),GlobalLeader(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,GlobalLeader(1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[GlobalLeader(D),0],'LineWidth',2);



% 结论：效果凑合，但速度极慢，大约在200次迭代收敛
    















        