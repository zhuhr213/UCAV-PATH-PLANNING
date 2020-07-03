% PSO，粒子群优化算法，是一种基于种群的随机优化算法，每个成员学习自身经验和伙伴经验来改变搜索模式
% 算法思想：种群中每个个体单独搜寻最优解，记为个体极值，并于其他个体分享，找到个体极值最优的那个
% 作为全局最优解。然后，所有个体根据自己的个体极值和全局极值来调整自己的位置。
% 粒子的两个属性：速度和位置，位置也包含了移动方向

% 非常像梯度下降！不过每次需要自己寻找方向

% 公式：
% vi = vi + c1*rand()*(pbest-xi)+c2*rand()*(gbest-xi) 更新速度（记忆、自我认知、整体认知）
% xi = xi + vi 更新位置。随机数取值（0，1）

% 由这两个式子，得到标准式：
% vi = Ω（小写）*vi + c1*rand()*(pbest-xi)+c2*rand()*(gbest-xi)
% 其中，Ω（小写）为惯性因子，非负。，w大，全局寻优能力强；w小，局部寻优能力强
% 可以这么理解：w小，vi的值就小，xi的变化就小->在局部活动

clear
clc

global radar1
global radar2
global R
% radar1 = [350 105 305 105 175 245 415 480 40 470];
% baili = size(radar1,2);
% radar2 = [200 0 -150 110 110 110 0 110 100 -50];
% R = [140 70 150 30 25 25 25 90 40 30];


radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
baili = size(radar1,2);
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];
%------给定初始化条件----------------------------------------------
c1=1.4962;             %学习因子1
c2=1.4962;             %学习因子2
w=0.7298;              %惯性权重
runtime =1;
MaxDT=1000;            %最大迭代次数
storeer = zeros(runtime,MaxDT);
Pbest = zeros(MaxDT,1);
D=60;                  %搜索空间维数（未知数个数）
N=20;                  %初始化群体个体数目
ub=ones(1,D).*100; %/*lower bounds of the parameters. */
lb=ones(1,D).*-50;%/*upper bound of the parameters.*/
for r=1:runtime
    
%------初始化种群的个体(可以在这里限定位置和速度的范围)------------
Range = repmat((ub-lb),[N 1]);
Lower = repmat(lb, [N 1]);
Foods = rand(N,D) .* Range + Lower;  % 随机生成初始位置和速度（Range只在这里使用了）
x = Foods;
Foods = rand(N,D) .* Range + Lower;
v = Foods;
%------先计算各个粒子的适应度，并初始化Pi和pg----------------------
for i=1:N
    p(i)=fitness(x(i,:)); % 求解适应值，利用新概率模型
    y(i,:)=x(i,:);  % D已经对平面划分了区域
end
pg=x(1,:);             %pg为全局最优
for i=2:N
    if fitness(x(i,:))<fitness(pg)
        pg=x(i,:);
    end
end  % 寻找初始全局最优
%------进入主要循环，按照公式依次迭代，直到满足精度要求------------
for t=1:MaxDT
    for i=1:N
        v(i,:)=w*v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg-x(i,:));
        x(i,:)=x(i,:)+v(i,:);
        sol=x(i,:);
        % 修正
        ind=find(sol<lb);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        ind=find(sol>ub);
        libai = rand(1,D).*(ub-lb)+lb;
        sol(ind)=libai(ind);
        %
        x(i,:)=sol;
        if fitness(x(i,:))<p(i) % 将新的适应值与原来的比较
            p(i)=fitness(x(i,:)); % 如果更优，就更新
            y(i,:)=x(i,:);
        end
        if p(i)<fitness(pg) % 更新全局最优
            pg=y(i,:);
        end
    end
    Pbest(t)=fitness(pg);
    fprintf('iteration = %d ObjVal=%g\n',t,Pbest(t));
end
%------最后给出计算结果
GlobalParams(r,:) = pg; % 这个是值
storeer(r,:) = Pbest;   % 这个是适应值
end




% figure (1)
% apso = mean(storeer);
% plot(apso)


figure (2)
hold on
plot(0,0,'k*')
plot(500,0,'ks')

for i = 1:baili
hold on
plot(radar1(i),radar2(i),'ko');
cir_plot([radar1(i),radar2(i)],R(i));
end
legend('starting point','target point','threat center')

axis equal
%axis([-100,620,-300,200])
for i = 1:(D-1)
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[GlobalParams(i),GlobalParams(i+1)],'r:','LineWidth',2);
    hold on
end
plot([0,500/(D+1)*(1)],[0,GlobalParams(1)],'r:','LineWidth',2);
plot([500/(D+1)*D,500],[GlobalParams(D),0],'r:','LineWidth',2);