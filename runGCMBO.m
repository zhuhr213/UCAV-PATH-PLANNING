% 帝王蝶优化算法（Monarch butterfly optimization）
% 现在的版本是引入了贪婪机制的改进算法
% 该算法模拟的是两个land的帝王蝶迁徙过程，秋天从land1到land2过冬，春天在回来
% land1中有ceil(p*N)个个体，ceil是向上取整，p是一个率，决定每个区域数量（这里就取0.5吧）
% land2中就有N-N1个
% 算法同时模拟了产卵过程。为了保持总数不变，出生一个小的就死去一个老的
% 所以在原算法中，都是不管新出生的好坏，统统保留。改进后的算法采用了贪心法
% 同时添加了自适应法，通过计算自适应值来繁衍

% 算法流程：
% 初始化种群P，总数为N，设置land1和land2的蝴蝶数量（N1和N2），计算适应值
% for 迭代次数：
%     对种群排序，参考蛾子算法，将整个群体分成两部分
%     for land1的所有个体：
%         for 对每一维：
%             更新每一维，更新方法参考Monarch butterfly optimization（原始算法）
%         end
%         贪心法决定是否保留（整体）
%     end
%     for land2的所有个体：
%         计算dx，alpha，其中dx是levy飞行函数，蛾子群算法出现过，就是下面这个东西
%         Ls = (beta-1)*gamma(beta - 1)*sin(pi*(beta-1)/2)/pi*(s^beta);
%         只是近似，真正的表达式是一个积分
%         for 对每一维：
%             更新每一维，同样参考上篇文章
%         end
%         对该个体，以自适应方式，生成一个新的值（模拟繁衍），并以贪心法选择是否保留
%     end
%     计算自适应值
% end

global radar1
global radar2
global R
radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

runtime=1;

maxCycle = 200;
D = 30;
N = 60;
p = 0.5;  %暂定两个land的个体数是一样的
N1 = ceil(p*N);
N2 = N-N1;
peri = 1;
beta = 2;
Smax = 15;
BAR = 0.75;  % butterfly adjusting rate，并不知道怎么设置

ub=ones(1,D).*50;  % 上下界
lb=ones(1,D).*-50;

Smax = 15;
beta = 1.5;

NP = rand(N,D).*(repmat((ub-lb),[N 1])) + lb;  

for i = 1:N
    NP(i,D+1) = calcu(NP(i,1:D));  % 计算objval值（代价函数），越小越好
end

NP = sortrows(NP,D+1);  % 按照D+1维进行排序

for r = 1:runtime

for iter = 1:maxCycle
    for i = 1:N1  % 遍历land1 （模拟迁徙）
        tempSol = NP(i,:);
        for j = 1:D  % 遍历每一维
            if (rand*peri <= p)  % 使用land1的个体更新
                neibour = fix(rand*(N1-1))+1;
                while (neibour == i)
                    neibour = fix(rand*(N1-1))+1;
                end
                tempSol(j) = NP(neibour,j);
            else
                neibour = fix(rand*(N-(N1+1)))+1;
                tempSol(j) = NP(neibour,j);
            end
        end
        tempSol(D+1) = calcu(tempSol(1:D));
        if (tempSol(D+1) < NP(i,D+1))
            NP(i,:) = tempSol;  % 贪心法  这里只是替换，理论上不会出现越界
        end
    end
    alpha = Smax/(iter^2);
    for i = N1+1:N
        s = abs(normrnd(0,1));
        dx = (beta-1)*gamma(beta - 1)*sin(pi*(beta-1)/2)/pi*(s^beta);
        
        tempSol1 = NP(i,:);
        for j = 1:D
            temprand = rand;
            if (temprand<=p)
                tempSol1(j) = NP(1,j);
            else
                neibour = fix(rand*(N-(N1+1)))+1;
                while (neibour == i)
                    neibour = fix(rand*(N-(N1+1)))+1;
                end
                tempSol1(j) = NP(neibour,j);
                if (temprand>BAR)
                    tempSol1(j) = tempSol1(j)+alpha*(dx-0.5);
                end
            end
        end
        
        Cr = 0.8 + 0.2*(NP(i,D+1)-NP(1,D+1))/(NP(N,D+1)-NP(1,D+1));  % 自适应值
        tempSol2(1:D) = tempSol1(1:D).*(1-Cr) + NP(i,1:D).*Cr;
        tempSol1(D+1) = calcu(tempSol1(1:D));
        tempSol2(D+1) = calcu(tempSol2(1:D));
        if (tempSol1(D+1)<tempSol2(D+1))
            NP(i,:) = tempSol1;
        else
            NP(i,:) = tempSol2;
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
            









