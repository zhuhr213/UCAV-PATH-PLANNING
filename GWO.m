% GWO，灰狼算法
% 算法中设置了alpha狼（头狼，最优候选解），beta狼（第二最优解），
% delta狼（第三最优解），omega狼（紧跟前三只狼的步伐）
% 算法模拟了狼群包围、狩猎、攻击、搜索的过程
% 迭代结束，头狼就是最优解

% 算法流程：
% 初始化狼群（随机解）
% 初始化a,A和C（相关参数）
% 计算每个个体的适应值
% Xa=最优解
% Xb=第二优解
% Xc=第三优解
% while(迭代次数）
%     for 每一个个体
%         更新每个个体的目前位置（式3.7）
%     end
%     更新a,A和C
%     计算每个个体的适应值
%     更新三个最优解
%     迭代次数++
% end
% 头狼就是最优解

% 似乎与ABC没有什么不同：多了两个待定解，每次更新所有位置而不只是一个维
% 更新的机制不一样
% 缺点：每头狼之间没有交流，很容易陷入局部最优
% 在dim = 30 ,swarm_number = 30的情况下效果相对最好

clear
clc

global radar1
global radar2
global R

radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];

% 画地图

dim = 30;  % 维数
max_iter = 1000;  % 迭代次数
swarm_no = 30;    % 狼的数量
lb_no = ones(1,dim).*-50;  % 上下界
ub_no = ones(1,dim).*50;

Alpha_pos=zeros(1,dim); % 初始化3狼。如果是求最大值问题，就改成-inf
Alpha_score= inf; 

Beta_pos=zeros(1,dim);
Beta_score= inf; 

Delta_pos=zeros(1,dim);
Delta_score= inf; 

Positions = repmat((ub_no-lb_no),[swarm_no,1]);
Positions = rand(swarm_no,dim).*Positions + lb_no;  % 初始化每只辣鸡狼的位置

for iter = 1:max_iter  % 开始迭代
    for i = 1:swarm_no  % 计算自适应值
        ObjVal(i) = calcu(Positions(i,:)); 
        % Fitness(i) = calculateFitness(ObjVal(i));
        
        % 更新三狼
        
        if ObjVal(i)<Alpha_score 
            Alpha_score=ObjVal(i); % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if ObjVal(i)>Alpha_score && ObjVal(i)<Beta_score 
            Beta_score=ObjVal(i); % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if ObjVal(i)>Alpha_score && ObjVal(i)>Beta_score && ObjVal(i)<Delta_score 
            Delta_score=ObjVal(i); % Update delta
            Delta_pos=Positions(i,:);
        end
        
    end
    
    a=2-iter*((2)/max_iter);  % 更新参数a
    
    for i=1:size(Positions,1)  % 更新每个辣鸡狼
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
        
        % 处理一下越界的
        Flag4ub=Positions(i,:)>ub_no;
        Flag4lb=Positions(i,:)<lb_no;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub_no.*Flag4ub+lb_no.*Flag4lb;
    end
    fprintf('iteration = %d ObjVal=%g\n',iter,Alpha_score);  % 记录迭代次数
end

% 迭代结束，此时头狼就是最优解

%画图

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
for i = 1:(dim-1)
    hold on
    plot([500/(dim+1)*i,500/(dim+1)*(i+1)],[Alpha_pos(i),Alpha_pos(i+1)],'LineWidth',2);
    % 画出路径
end
hold on
plot([0,500/(dim+1)*(1)],[0,Alpha_pos(1)],'LineWidth',2); % 起点和终点
plot([500/(dim+1)*dim,500],[Alpha_pos(dim),0],'LineWidth',2);












