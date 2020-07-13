% 布谷鸟算法（CS），模拟的是布谷鸟把蛋下别的鸟那里
% 他们的蛋出来的早，会把巢里的其他蛋弄死，以便自己吃更多
% 有一定几率被发现，越厉害的越不容易被发现

% 算法流程
% 步骤１ 定义目标函数ｆ（Ｘ），Ｘ ＝ （ｘ１，…ｘｄ）Ｔ ，函数初始化，并随机生成ｎ个鸟窝的初始位置
% 设置种群规模、问题维数、最大发现概率Ρ和最大迭代次数等参数；
% 步骤２ 选择适应度函数并计算每个鸟窝位置的目标函数值，得到当前的最优函数值；
% 步骤３ 记录上一代最优函数值，利 用 式 （１）对 其 他鸟窝的位置和状态进行更新；
% 步骤４ 现有位置函数值与上一代最优函数值进行比较，若较好，则改变当前最优值；
% 步骤５ 通过位置更新后，用随机数ｒ［０，１］与Ｐ 对比，若ｒ ＞Ｐ ，则对ｘ（ｔ＋１）ｉ 进行随机改变，反之则不变。最后保留最好的一组鸟窝位置ｙ（ｔ＋１）
% 步骤６ 若未达到最大迭代次数或最小误差要求，则返回步骤２，否则，继续下一步；
% 步骤７ 输出全局最优位置。

clear
runtime = 1;

for time = 1:runtime
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

maxCycle = 1000;

D = 30;
N = 60;

ub=ones(1,D).*50;  % 上下界 
lb=ones(1,D).*-50;

beta = 1.5;  % 又是莱维飞行

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;   % 为每个鸟窝生成随机位置

for i = 1:N
    ObjVal(i) = calcu(NP(i,:));
end

Fitness = calculateFitness(ObjVal);
BestInd = find(ObjVal == min(ObjVal));
BestInd = BestInd(end);
GlobalBest = NP(BestInd,:);
GlobalObj = ObjVal(BestInd);  % 找到最好的那个

% for r = 1:runtime

% https://ww2.mathworks.cn/matlabcentral/fileexchange/29809-cuckoo-search-cs-algorithm
    for iter = 1:maxCycle
        Ls = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
        for i = 1:N
            s = NP(i,:);
            % 模拟莱维飞行
            u = randn(size(s))*Ls;
            v = randn(size(s));
            step = u./abs(v).^(1/beta);
            stepsize = 0.5*step.*(s-GlobalBest);  % 保留最好的解
            % 这个系数影响代价函数曲线的顺滑程度和收敛速度
            % 太低的话收敛慢，太高的话收敛值不够低
            % 综合看，0.5最合适
            s = s+stepsize.*randn(size(s));
            ind = find(s>ub|s<lb);  % ||用于标量，|用于矩阵
            % 可能生成新随机解这个地方需要重新修改
            s(ind) = rand(size(ind)).*(ub(ind)-lb(ind))+lb(ind);  % 随机生成越界
            tempObj = calcu(s);
            if (tempObj < ObjVal(i))
                NP(i,:) = s;
                ObjVal(i) = tempObj;
                Fitness(i) = calculateFitness(ObjVal(i));
            end
        end
        
        % 以fitness值作为概率，在现实中，越合适的应该越不容易被发现
        % randperm：https://blog.csdn.net/xuxinrk/article/details/80961657
        
        tempfitness = Fitness';
        tempfitness = repmat(tempfitness,[1,D]);
        K = rand(size(NP))>tempfitness; % 这是个逻辑表达式，返回的是矩阵，对应位置是0或1
        % 非常神奇的是！如果，每个鸟都是弟弟鸟，不会发现不是自己的鸟蛋
        % 结果就很好了
        % 所以推测，是被发现后的调整函数不好
        % stepsize = 0.1*stepsize; % 被发现了，稍微移动一下
            u = randn(size(NP))*Ls;
            v = randn(size(NP));
            step = u./abs(v).^(1/beta);
            stepsize = 0.01*step;  % 由于是被发现，微调，所以系数选择更小
            % s = s+stepsize.*randn(size(s));
        NP = NP + stepsize.*K.*randn(size(NP));  % 这里乘不乘randn效果都很好
        % 越界检查
        for i = 1:N
            ind = find(NP(i,:)>ub|NP(i,:)<lb);
            NP(i,ind) = rand(size(ind)).*(ub(ind)-lb(ind))+lb(ind);  % 不随机生成越界
        end
                        
        tempBestInd = find(ObjVal == min(ObjVal));
        tempBestInd = tempBestInd(end);
        tempGlobalBest = NP(tempBestInd,:);
        tempGlobalObj = ObjVal(tempBestInd);  % 找到最好的那个
        
        if (tempGlobalObj < GlobalObj)
            GlobalBest = tempGlobalBest;
            GlobalObj = tempGlobalObj;
        end
            
        aaaaa(iter) = GlobalObj; 
        fprintf('iteration = %d ObjVal=%g Fitness=%g\n',iter,GlobalObj,calculateFitness(GlobalObj));
    end
        storeer(time,:) = aaaaa;
    LoopBest(time,:) = GlobalBest;
end

figure (1)
a = mean(storeer,1);
% save CSdata a;
plot([1:maxCycle],a);
ylim([0,10]);

figure (2)
hold on
plot(0,0,'k*');
plot(500,0,'ks');

for i = 1:size(radar1,2)  % 威胁的个数
hold on
plot(radar1(i),radar2(i),'ko');
cir_plot([radar1(i),radar2(i)],R(i)); % 参数方程画圆
end
% legend('starting point','target point','threat center')


axis equal  % 横轴纵轴长度单位相同
axis([-100,620,-400,400])

hold off

for i = 1:(D-1)
    hold on
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalBest(i),GlobalBest(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,GlobalBest(1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[GlobalBest(D),0],'LineWidth',2);

        
        
            
            








