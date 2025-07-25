clear; clc;

N = 100;    % 群体粒子个数
D = 4;     % 粒子维度
T = 200;   % 最大迭代次数
c1 = 2.05;    % 学习因子1
c2 = 2.05;    % 学习因子2
wmax = 0.8;   % 惯性权重最大值
wmin = 0.5;  % 惯性权重最小值
data = xlsread('Data/data_experiment_Liang.xlsx','ctrl');
par = setparameter();
%%%%%%%%%     初始化种群位置与种群速度    %%%%%%%%%
xmax =[1, 1e-10, 1, 1]; % 位置最大值
xmin = [0.1, 1e-12, 0.01, 0]; % 位置最小值

lambdaC_LHS = LHS_Call(xmax(1),(xmax(1)+xmin(1))/2,xmin(1), 0,N,'unif');
alpha1_LHS = LHS_Call( xmax(2),(xmax(2)+xmin(2))/2,xmin(2), 0,N,'unif');
lambdaT_LHS = LHS_Call(xmax(3),(xmax(3)+xmin(3))/2,xmin(3), 0,N,'unif');
rho_LHS = LHS_Call( xmax(4),(xmax(4)+xmin(4))/2,xmin(4), 0,N,'unif');

vmax = 0.01*xmax; % 速度最大值
vmin = -vmax;   % 速度最小值
x = [lambdaC_LHS,alpha1_LHS,lambdaT_LHS,rho_LHS];
v = rand(N, D) .* (vmax - vmin) + vmin;%初始种群速度

%%%%%%%%%     初始全局最优位置与个体最优位置    %%%%%%%%%
p = x;%初始化每个个体最优位置
pbest = ones(1,N);%初始化所有个体最优解
options = odeset('RelTol', 1e-4, 'AbsTol', [1e-6,1e-6,1e-6,1e-6]);% 定义ODE求解器的选项
t_span = [0,5,10,14,19,22,26,30];% 定义t间隔
y0 =[data(1,2),50,1.3*10^7,10^7];%C、S、E、T cells细胞的初值
soly = ones(length(t_span),N);%初始一个矩阵

for i=1:N
    par.lambdaC = x(i,1); par.alpha1 = x(i,2); par.lambdaT = x(i,3); par.rho = x(i,4);
    [~, jie] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);%解微分方程
    soly(:,i) = jie(:,1);%将每次肿瘤细胞数量计入soly矩阵中
    pbest(:,i) = ObjFun(data,soly(:,i));%将每一个个体的误差计入
end

[minimum,index] = min(pbest);%记录个体解误差最小的位置以及最小值
gbest = minimum;%所有个体误差最小值为初始全局最优值
g = p(index,:);%全局最优值代表的粒子为初始全局最优位置

%%
%%%%%%%%%%  按照公式依次迭代直到满足精度或者迭代次数 %%%%%%%%%
zhi=ones(1,N);%用来记录每一次迭代中个体的误差
for i = 1 : T
    %%%%%%%%%%  动态计算惯性权重值 %%%%%%%%%
    w = wmax - (wmax - wmin)*i/T;
    for j = 1 : N
        %%%%%%%%%%  更新位置和速度值 %%%%%%%%%
        v(j, :) = w*v(j, :) + c1*(p(j,:)- x(j, :)) + c2*(g - x(j, :));
        x(j, :) = x(j, :) + v(j, :);
        %%%%%%%%%%  边界条件处理 %%%%%%%%%
        for m = 1:D
            if (v(j, m) > vmax(m)) || (v(j, m) < vmin(m))
                v(j,:) = rand*(vmax - vmin) + vmin;
            end
            if (x(j, m) > xmax(m))||(x(j, m) < xmin(m))
                x(j,:) = rand*(xmax - xmin) + xmin;
            end
        end
        %%%%%%%%%%  更新个体最优位置和最优值 %%%%%%%%%
        par.lambdaC = x(j,1); par.alpha1 = x(j,2); par.lambdaT = x(j,3); par.rho = x(j,4);
        [~, jie] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);%解微分方程
        soly(:,j) = jie(:,1);%将每次肿瘤细胞数量计入soly矩阵中
        zhi(1,j) = ObjFun(data,soly(:,j));%将每一个个体的误差计入
        
        if  zhi(:,j)<pbest(:,j)
            p(j, :) = x(j, :);%更新个体最优位置
            pbest(:,j) = zhi(:,j);%更新个体最优值
        end
        
        [zho,ind] = min(pbest);%zho为每一次迭代更新后的最优解
        gbest = zho;%更新全局最优值
        g = p(ind,:);%更新全局最优位置
    end
end

%%
%%%%%%%%%%%%%%%%  画肿瘤细胞图   %%%%%%%%%%

dlmwrite('Data/data_BestFitPara_control.csv', g, 'delimiter', ',');
dlmwrite('Data/data_FitGood_R2_control.csv', gbest, 'delimiter', ',');

par.lambdaC = g(1);  par.alpha1 = g(2);  par.lambdaT = g(3);  par.rho = g(4);

t1_span = (0:0.5:40);
[t, R] = ode45(@(t,y) ODE_control(t,y,par), t1_span,y0, options);%解微分方程
C = R(:,1); E = R(:,2); T_1 = R(:,3); T_2 = R(:,4);
x2 = data(:,1); y2 = data(:,2);

plot(t1_span,C,'red-',x2,y2,'black o')
xlabel('Time(days)');
ylabel('Tumor cell number')



