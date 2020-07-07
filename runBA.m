% 蝙蝠算法，模拟蝙蝠通过声纳寻找猎物并区分障碍物
% 所有蝙蝠利用回声定位的方法感知距离,并且它们采用一种巧妙的方式来区别猎物和背景障碍物之间的不同。
% % 蝙蝠在位置xi以速度vi随机飞行,以固定的频率fmin、可变的波长λ和音量A0来搜索猎物。
% 蝙蝠根据自身与目标的邻近程度来自动调整发射的脉冲波长（或频率）和调整脉冲发射率r属于[0,1]
% 虽然音量的变化方式有多种但在蝙蝠算法中, 假定音量A是从一个最大值A0(整数)变化到固定最小值Amin

% https://www.cnblogs.com/caoer/p/12641369.html
% 它非常像差分进化！！！只不过没有交叉，且利用的是最优解来更新

% 源代码参考
% https://ww2.mathworks.cn/matlabcentral/fileexchange/74768-the-standard-bat-algorithm-ba

runtime=10;

for r = 1:runtime
clear
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];


maxCycle = 500;

D = 30;
N = 60;

Fmin = 1;
Fmax = 4;  % 频率上下限
f = zeros(N,1);  % 每个蝙蝠的频率
v = zeros(N,D);  % 速度，设每一维的都不一样

A = 5;  % 噪声
r0 = 10; % 脉冲
alpha = 0.97;
gamma = 1.5;  % 两个参数


ub=ones(1,D).*50;
lb=ones(1,D).*-50;

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;   

for i = 1:N
    ObjVal(i) = calcu(NP(i,:));
end
Fitness = calculateFitness(ObjVal);

BestInd = find(ObjVal == min(ObjVal));
BestInd = BestInd(end);
GlobalBest = NP(BestInd,:);
GlobalObj = ObjVal(BestInd);


    for iter = 1:maxCycle
        r = r0*(1-exp(-gamma*iter));
        A = alpha*A;  % 初始化脉冲率和噪声
        for i = 1:N
        	f(i) = Fmin + (Fmax-Fmin)*rand;
        	v(i,:) = v(i,:) + (NP(i,:) - GlobalBest).*f(i);  % 速度
        	temp = NP(i,:) + v(i,:);
            
            if (rand<r)
                temp = GlobalBest + 0.1*randn(1,D)*A; % randn是正态分布
            end
            
            for j = 1:D  % 越界检查
                if (temp(j)>ub(j)||temp(j)<lb(j))
                    temp(j) = rand*(ub(j)-lb(j))+lb(j);
                end
            end
            
            tempObj = calcu(temp);
            
            if (tempObj<ObjVal(i) && rand>A)  % 如果更优且噪声不大(与当前解相比）
                NP(i,:) = temp;
                ObjVal(i) = tempObj;
            end
            
            if (tempObj<=GlobalObj) % 更新全局
                GlobalBest = temp;
                GlobalObj = tempObj;
            end
        end
        aaaaa(iter) = GlobalObj; 
        fprintf('iteration = %d ObjVal=%g\n',iter,GlobalObj);
    end
    storeer(r,:) = aaaaa;
end



 
figure (1)
a = mean(storeer);
plot([1:maxCycle],a);

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


% 蝙蝠算法
% 挺快的，参数好多，好复杂，效果一般，调调参数
% 结果不稳定，总是陷入局部最优







