% 萤火虫算法
% 萤火虫不分性别，这样一个萤火虫将会吸引到所有其他的萤火虫;
% 吸引力与它们的亮度成正比，对于任何两个萤火虫，不那么明亮的萤火虫被吸引，因此移向更亮的一个(全局搜索)
% 然而，亮度又随着其距离的增加而减少;如果没有比一个给定的萤火虫更亮的萤火虫
% 它会随机移动（局部搜索）。亮度应与目标函数联系起来。
% 所以，重要的便是计算亮度和吸引度，他们不仅与ObjVal有关，也与笛卡尔距离有关

clear
clc

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

runtime=1;
maxCycle = 1000;

D = 30;
N = 60;  

ub=ones(1,D).*50;  % 上下界
lb=ones(1,D).*-50;

L = max(ub)-min(lb);
gama = 1/sqrt(L);  % 吸引度衰减参数，取0说明可被任意萤火虫发现，所有的均朝向最优（PSO）
                                 % 取无穷说明全都近视，退化为全部随机飞
alpha = 0.01*L;  % 控制随机运动步长，太大震荡，太小陷入局部，或者可以通过迭代次数来动态变化
% 它用来更新最亮的

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;   

for i = 1:N
    ObjVal(i) = calcu(NP(i,:));
end

Fitness = calculateFitness(ObjVal);

BestInd = find(ObjVal == min(ObjVal));
BestInd = BestInd(end);
GlobalBest = NP(BestInd,:);
GlobalObj = ObjVal(BestInd);

for r = 1:runtime
    for iter = 1:maxCycle
        for i = 1:N
            for j = 1:N
                if (i~=j)  % 取反为~
                    if (ObjVal(j)<ObjVal(i))  % j比i要更好，i就向j移动
                        dist = sqrt(sum(NP(i,:)-NP(j,:)).^2);
                        attr = exp(-(gama*dist^2));  % 吸引度
                        NP(i,:) = NP(i,:) + attr*(NP(j,:)-NP(i,:)) + rand*(rand-0.5)*alpha;
                        for k=1:D
                            if (NP(i,k)>ub(k)||NP(i,k)<lb(k))
                                NP(i,k) = rand;
                            end
                        end
                        ObjVal(i) = calcu(NP(i,:));
                    end
                end
            end
        end % 理论上，最亮的那个应该不会向别的萤火虫移动，所以放在下面更新最亮的
%         BestInd = find(ObjVal == min(ObjVal));
%         BestInd = BestInd(end);
        NP(BestInd,:) = NP(BestInd,:) + rand(1,D).*(rand(1,D)-0.5)*alpha;
        for k=1:D
        	if (NP(BestInd,k)>ub(k)||NP(BestInd,k)<lb(k))
            	NP(BestInd,k) = rand;
        	end
        end
        ObjVal(BestInd) = calcu(NP(BestInd,:));    
        
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
                        
                        
                        
% 极慢，比蜘蛛猴还慢   
% 没有贪心选择，可能不单调
                        
                        
                        
                        
                    

