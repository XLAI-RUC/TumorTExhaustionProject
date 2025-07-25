clear; clc;

N = 100;    % Ⱥ�����Ӹ���
D = 4;     % ����ά��
T = 200;   % ����������
c1 = 2.05;    % ѧϰ����1
c2 = 2.05;    % ѧϰ����2
wmax = 0.8;   % ����Ȩ�����ֵ
wmin = 0.5;  % ����Ȩ����Сֵ
data = xlsread('Data/data_experiment_Liang.xlsx','ctrl');
par = setparameter();
%%%%%%%%%     ��ʼ����Ⱥλ������Ⱥ�ٶ�    %%%%%%%%%
xmax =[1, 1e-10, 1, 1]; % λ�����ֵ
xmin = [0.1, 1e-12, 0.01, 0]; % λ����Сֵ

lambdaC_LHS = LHS_Call(xmax(1),(xmax(1)+xmin(1))/2,xmin(1), 0,N,'unif');
alpha1_LHS = LHS_Call( xmax(2),(xmax(2)+xmin(2))/2,xmin(2), 0,N,'unif');
lambdaT_LHS = LHS_Call(xmax(3),(xmax(3)+xmin(3))/2,xmin(3), 0,N,'unif');
rho_LHS = LHS_Call( xmax(4),(xmax(4)+xmin(4))/2,xmin(4), 0,N,'unif');

vmax = 0.01*xmax; % �ٶ����ֵ
vmin = -vmax;   % �ٶ���Сֵ
x = [lambdaC_LHS,alpha1_LHS,lambdaT_LHS,rho_LHS];
v = rand(N, D) .* (vmax - vmin) + vmin;%��ʼ��Ⱥ�ٶ�

%%%%%%%%%     ��ʼȫ������λ�����������λ��    %%%%%%%%%
p = x;%��ʼ��ÿ����������λ��
pbest = ones(1,N);%��ʼ�����и������Ž�
options = odeset('RelTol', 1e-4, 'AbsTol', [1e-6,1e-6,1e-6,1e-6]);% ����ODE�������ѡ��
t_span = [0,5,10,14,19,22,26,30];% ����t���
y0 =[data(1,2),50,1.3*10^7,10^7];%C��S��E��T cellsϸ���ĳ�ֵ
soly = ones(length(t_span),N);%��ʼһ������

for i=1:N
    par.lambdaC = x(i,1); par.alpha1 = x(i,2); par.lambdaT = x(i,3); par.rho = x(i,4);
    [~, jie] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);%��΢�ַ���
    soly(:,i) = jie(:,1);%��ÿ������ϸ����������soly������
    pbest(:,i) = ObjFun(data,soly(:,i));%��ÿһ�������������
end

[minimum,index] = min(pbest);%��¼����������С��λ���Լ���Сֵ
gbest = minimum;%���и��������СֵΪ��ʼȫ������ֵ
g = p(index,:);%ȫ������ֵ���������Ϊ��ʼȫ������λ��

%%
%%%%%%%%%%  ���չ�ʽ���ε���ֱ�����㾫�Ȼ��ߵ������� %%%%%%%%%
zhi=ones(1,N);%������¼ÿһ�ε����и�������
for i = 1 : T
    %%%%%%%%%%  ��̬�������Ȩ��ֵ %%%%%%%%%
    w = wmax - (wmax - wmin)*i/T;
    for j = 1 : N
        %%%%%%%%%%  ����λ�ú��ٶ�ֵ %%%%%%%%%
        v(j, :) = w*v(j, :) + c1*(p(j,:)- x(j, :)) + c2*(g - x(j, :));
        x(j, :) = x(j, :) + v(j, :);
        %%%%%%%%%%  �߽��������� %%%%%%%%%
        for m = 1:D
            if (v(j, m) > vmax(m)) || (v(j, m) < vmin(m))
                v(j,:) = rand*(vmax - vmin) + vmin;
            end
            if (x(j, m) > xmax(m))||(x(j, m) < xmin(m))
                x(j,:) = rand*(xmax - xmin) + xmin;
            end
        end
        %%%%%%%%%%  ���¸�������λ�ú�����ֵ %%%%%%%%%
        par.lambdaC = x(j,1); par.alpha1 = x(j,2); par.lambdaT = x(j,3); par.rho = x(j,4);
        [~, jie] = ode45(@(t,y) ODE_control(t,y,par), t_span,y0, options);%��΢�ַ���
        soly(:,j) = jie(:,1);%��ÿ������ϸ����������soly������
        zhi(1,j) = ObjFun(data,soly(:,j));%��ÿһ�������������
        
        if  zhi(:,j)<pbest(:,j)
            p(j, :) = x(j, :);%���¸�������λ��
            pbest(:,j) = zhi(:,j);%���¸�������ֵ
        end
        
        [zho,ind] = min(pbest);%zhoΪÿһ�ε������º�����Ž�
        gbest = zho;%����ȫ������ֵ
        g = p(ind,:);%����ȫ������λ��
    end
end

%%
%%%%%%%%%%%%%%%%  ������ϸ��ͼ   %%%%%%%%%%

dlmwrite('Data/data_BestFitPara_control.csv', g, 'delimiter', ',');
dlmwrite('Data/data_FitGood_R2_control.csv', gbest, 'delimiter', ',');

par.lambdaC = g(1);  par.alpha1 = g(2);  par.lambdaT = g(3);  par.rho = g(4);

t1_span = (0:0.5:40);
[t, R] = ode45(@(t,y) ODE_control(t,y,par), t1_span,y0, options);%��΢�ַ���
C = R(:,1); E = R(:,2); T_1 = R(:,3); T_2 = R(:,4);
x2 = data(:,1); y2 = data(:,2);

plot(t1_span,C,'red-',x2,y2,'black o')
xlabel('Time(days)');
ylabel('Tumor cell number')



