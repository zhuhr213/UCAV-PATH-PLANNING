function name1 = calcu(path)
global radar1
global radar2
global R

% radar1 = [100 200 300 400 150 250 350 150 250 350 0 466 250 250 466 30];
% radar2 = [0 0 0 0 50 50 50 -50 -50 -50 40 40 -300 300 -40 -20];
% R = [40 40 40 40 40 40 40 40 40 40 20 20 260 277 20 30];
% 已在runABC中定义（好像是半径）
nimama = size(radar2,2);  % 这里是16
fenmu = R./log(20);
node_number = size(path,2); % 维数，path是传入参数，这里应该是30
counter = zeros(node_number,1);
countor = zeros(node_number,1);
d_x = 500./(node_number+1);  % 是个数，但是它代表啥？（就当是初始位置吧）
% 我好像懂了它是什么了，参考论文621页
% 参考论文第621页，将平面地图用D根直线切割，则划分成D+1块，每块的长度即为d_x
for i = 1:node_number
    for j = 1:nimama % 第一个威胁的第一维，第二个威胁的第一维.....
        d = sqrt((radar1(j) - d_x.*i)^2 + (radar2(j) - path(i))^2); % 距离？？？
        % 每一块，相当于每一步。需要计算的是每一步的代价函数
        % 一共由node_number步，需要对每一步，对每个威胁，计算各个威胁到该步的距离
        % 传入参数是一个food行向量，那么它的每一列，应该就是纵坐标
        counter(i) = counter(i) + exp(-1.*d./fenmu(j));   % 概率密度模型
    end
end
for i = 1:(node_number-1)  % 这里是完整代价函数的第二部分
    countor(i) = sqrt((d_x)^2 + (path(i)-path(i+1))^2);
end
bufen1 = sum(countor)./500;
bufen2 = sum(counter);
name1 = bufen1 + bufen2;





