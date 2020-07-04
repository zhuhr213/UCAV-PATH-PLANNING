% Moth Search algorithm，蛾子群算法，模拟的是蛾子的趋光性
% 将蛾子群分成了两组，第一组离光近，第二组离光远
% 每次都需要对蛾子群排序，这有点杀时间
% 离光远的，用跨大步算法，随机选择一个，这两个一个是越过best，一个是不越过best
% 离光近的，用随机游走法，绕着光源转圈圈

% 算法流程：
% 初始化，生成群解NP，设置参数Smax（跨大步每次跨的最大步长），β（一个指数算子，通常为1.5），加速因子phy
% 计算适应值fitness，对NP根据fitness排序
% for 迭代次数：
%     for 1 to N/2： 对厉害组进行随机游走，N是总数
%         随机游走
%     end
%     for N/2+1:N： 对不厉害组进行大步跨
%         if rand>0.5:
%             大大步跨
%         else
%             大步跨
%         end
%     end
%     更新fitness，基于fitness排序
% end

% 对于排序，可用sortrows()，将每一行看作整体，以某一列为基准
% 这样，需要多一列，存fitness或ObjVal

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

Smax = 15;  % 5太小，50太大，10~20是最好的选择
beta = 2;   % beta是2最合适
phy = (1+sqrt(5))/2;
% 对于gamma函数，直接y = gamma(x)即可，很厉害！

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;  

for i = 1:N
    NP(i,D+1) = calcu(NP(i,1:D));  % 计算objval值（代价函数），越小越好
end

% Fitness = calculateFitness(ObjVal); % 计算Fitness，Fitness是拟合度，越高越拟合
% 直接用ObjVal来排序吧

NP = sortrows(NP,D+1);  % 按照D+1维进行排序

for r = 1:runtime
    for iter = 1:maxCycle
        a = Smax/(iter^2);
        for i = 1:N/2
            s = rand+10;  % 算式里的一个参数，只说了它应该大于0，这里取rand+(5,10)都凑合
            Ls = (beta-1)*gamma(beta - 1)*sin(pi*(beta-1)/2)/pi*(s^beta);
            % 好多常数啊！可以提出去，减少计算量
            for j = 1:D
                NP(i,j) = NP(i,j) + a*Ls;
            end
        end
        for i = N/2+1:N
            if (rand>0.5)
                for j = 1:D
                    NP(i,j) = rand * (NP(i,j) + phy * (NP(1,j) - NP(i,j)));
                end
            else
                for j = 1:D
                    NP(i,j) = rand * (NP(i,j) + (1/phy) * (NP(1,j) - NP(i,j)));
                end
            end
        end
        for i = 1:N
    	% 修正
            repair = (ub-lb).*(rand(1,D))+lb;
            ind = find(NP(i,1:D)>ub);
            NP(i,ind) = repair(ind);
            ind = find(NP(i,1:D)<lb);
            NP(i,ind) = repair(ind); 
            NP(i,D+1) = calcu(NP(i,1:D));
        end  
        NP = sortrows(NP,D+1);  % 按照D+1维进行排序
        aaaaa(iter) = NP(1,D+1); 
        fprintf('iteration = %d ObjVal=%g\n',iter,NP(1,D+1));
    end
end

storeer(r,:) = aaaaa;

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
    plot([500/(D+1)*i,500/(D+1)*(i+1)],[NP(1,i),NP(1,i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(D+1)*(1)],[0,NP(1,1)],'LineWidth',2); % 起点和终点
plot([500/(D+1)*D,500],[NP(1,D),0],'LineWidth',2);
            
        








