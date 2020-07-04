% 鲸群算法。模拟的是鲸群捕食鱼群的过程
% 为方便抓捕，鲸群会吐泡泡将鱼群包围，并驱赶至水面
% 吐出的泡泡是螺旋式上升的，旨在提高鱼群的密度
% 与狼群算法有类似的参数
% 但与狼群算法略有不同：通过概率（0.5）来选择是使用螺旋式包围，还是直接包围
% 并通过参数|A|决定，是在局部寻优，还是向全局搜索

% 算法流程：
% 初始化，计算适应值，选择全局最优
% for 迭代次数：
%     for 每个个体：
%         更新a,A,C,l,p参数
%         if (p<0.5)
%             if (|A|<1)
%                 使用直接包围更新
%             elseif (|A|>=1)
%                     选择一个随机个体
%                     利用该个体更新
%             end
%         else
%             利用螺旋包围更新
%         end
%     end
%     检查越界，更新适应值，更新全局最优
% end

clear
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

r=1;

maxCycle = 1000;
D = 30;
N = 60;

ub=ones(1,D).*50;  % 上下界
lb=ones(1,D).*-50;

% WG for whales group
WG = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;  

for i = 1:N
    ObjVal(i) = calcu(WG(i,:));  % 计算objval值（代价函数），越小越好
end
Fitness = calculateFitness(ObjVal); % 计算Fitness，Fitness是拟合度，越高越拟合

ind = find(Fitness == max(Fitness));
ind = ind(end); 
GlobalBest = WG(ind,:);  % 找到全局最优
GlobalBestObjVal = ObjVal(ind);

for iter = 1:maxCycle
    a=2-iter*((2)/maxCycle); 
    for i = 1:N
        b = 1;
        A = 2*a*rand - a;
        C = 2*rand;
        l = rand*2-1;
        p = rand;  % 对每个个体，采用一样的概率，来决定如何更新。且对同一个个体，使用同一参数
        for j = 1:D
            if (p<0.5)
                if (abs(A)<1)
                    DP = abs(C*GlobalBest(j)-WG(i,j));
                    WG(i,j) = GlobalBest(j) - A*DP;
                else
                    neighbour = fix(N*rand)+1;
                    while(neighbour == i)
                        neighbour = fix(N*rand)+1;
                    end
                    DP = abs(C*WG(neighbour,j)-WG(i,j));
                    WG(i,j) = WG(neighbour,j) - A*DP;
                end
            else
                DP = abs(GlobalBest(j)-WG(i,j));
                WG(i,j) = DP*exp(b*l)*cos(2*pi*l) + GlobalBest(j);
            end
        end
    end
    for i = 1:N
    	% 修正
    	repair = (ub-lb).*(rand(1,D))+lb;
    	ind = find(WG(i,:)>ub);
    	WG(i,ind) = repair(ind);
    	ind = find(WG(i,:)<lb);
    	WG(i,ind) = repair(ind); 
    	ObjVal(i) = calcu(WG(i,:));  % 计算objval值（代价函数），越小越好
    	Fitness = calculateFitness(ObjVal(i)); % 计算Fitness，Fitness是拟合度，越高越拟合
    end    
    ind = find(Fitness == max(Fitness));
    ind = ind(end); 
    if (ObjVal(ind) < GlobalBestObjVal)
        GlobalBest = WG(ind,:);  % 找到全局最优
        GlobalBestObjVal = ObjVal(ind);
    end
    aaaaa(iter) = GlobalBestObjVal;  
    fprintf('iteration = %d ObjVal=%g\n',iter,GlobalBestObjVal);  % 记录迭代次数
end
storeer(r,:) = aaaaa;  % r是算法的运行次数，用来画平均值图，比较优劣
% 应该还有一层循环，用来循环执行算法

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
for i = 1:(D-1)
    hold on
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalBest(i),GlobalBest(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,GlobalBest(1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[GlobalBest(D),0],'LineWidth',2);













