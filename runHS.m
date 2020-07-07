% 和声搜索
% 将一个解看成D个乐器在演奏
% 初始化N个解。依概率选择是随机生成新解还是从N个解中随机选一个
% 如果是从N个解中选出来的，依概率选择是否进行微调
% 每一维，均如此
% 然后与最差的比较，决定是否替换

clear
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

runtime=1;
maxCycle = 2000;

D = 30;
N = 60;  

ub=ones(1,D).*50;  % 上下界
lb=ones(1,D).*-50;

L = max(ub)-min(lb);

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;   

for i = 1:N
    ObjVal(i) = calcu(NP(i,:));
end

HMCR = 0.875; % 记忆库取值概率，这个是影响效果的关键
PAR = 0.25;% 微调概率
BAR = L*0.01; % 微调贷款

Fitness = calculateFitness(ObjVal);

BestInd = find(ObjVal == min(ObjVal));
BestInd = BestInd(end);
GlobalBest = NP(BestInd,:);
GlobalObj = ObjVal(BestInd);

for r = 1:runtime
    for iter = 1:maxCycle
        % 产生一个新和声
        temp = zeros(1,D);
        for j = 1:D
            if (rand<HMCR)
                temp(j) = NP(fix(rand*N)+1,j);
            else
                temp(j) = lb(j) + rand*(ub(j)-lb(j)); % 随机一个
            end
            if (rand<PAR) % 微调概率
                temp(j) = temp(j) + 2*BAR*rand-BAR; % 微调
            end
        end
        tempObj = calcu(temp);
        
        BadInd = find(ObjVal == max(ObjVal));
        BadInd = BadInd(end);
        if (tempObj<ObjVal(BadInd)) % 新的和声更好
            NP(BadInd,:) = temp;
            ObjVal(BadInd) = tempObj;
        end
        
        BestInd = find(ObjVal == min(ObjVal));
        BestInd = BestInd(end);
        GlobalBest = NP(BestInd,:);
        GlobalObj = ObjVal(BestInd);
        
        aaaaa(iter) = GlobalObj; 
        fprintf('iteration = %d ObjVal=%g\n',iter,GlobalObj);
    end
end

% storeer(r,:) = aaaaa;

% 
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

hold off

for i = 1:(D-1)
    hold on
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalBest(i),GlobalBest(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,GlobalBest(1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[GlobalBest(D),0],'LineWidth',2);


% 运行极快，快的飞起
% 由于每次替换的是最差的，所以能推测是阶梯式下降
% 本以为表现糟糕，但降低变异概率，即，提升选取库中的概率后，效果强了很多很多
% 但，需要更多的迭代次数。但它实在太快了
% 但扩大战场范围（上下限后），表现一般，不稳定
                    







