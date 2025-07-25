clear; clc;

N = 200;    % Ⱥ�����Ӹ���
D = 5;     % ����ά��
T = 200;   % ����������
c1 = 2.05;    % ѧϰ����1
c2 = 2.05;    % ѧϰ����2
wmax = 0.8;   % ����Ȩ�����ֵ
wmin = 0.5;  % ����Ȩ����Сֵ

data = xlsread('Data/data_exp_PD1_IL2v.xlsx');
par = setparameter();

%%%%%%%%%     ��ʼ����Ⱥλ������Ⱥ�ٶ�    %%%%%%%%%
xmax = [5e-5, 1.7e-8,3e-6,0.1,0.3]; % λ�����ֵ
xmin = [2e-5,1.5688e-8,1e-6,0.05,0.1]; % λ����Сֵ

Kb_LHS = LHS_Call(xmax(1),(xmax(1)+xmin(1))/2,xmin(1), 0,N,'unif');
beta2_LHS = LHS_Call( xmax(2),(xmax(2)+xmin(2))/2,xmin(2), 0,N,'unif');
zeta_LHS = LHS_Call( xmax(3),(xmax(3)+xmin(3))/2,xmin(3), 0,N,'unif');
rho_LHS = LHS_Call( xmax(4),(xmax(4)+xmin(4))/2,xmin(4), 0,N,'unif');
dB_LHS = LHS_Call( xmax(5),(xmax(5)+xmin(5))/2,xmin(5), 0,N,'unif');

vmax = 0.01*xmax; % �ٶ����ֵ
vmin = -vmax;   % �ٶ���Сֵ
x = [Kb_LHS,beta2_LHS,zeta_LHS,rho_LHS,dB_LHS];
v = rand(N, D) .* (vmax - vmin) + vmin;%��ʼ��Ⱥ�ٶ�

%%%%%%%%%     ��ʼȫ������λ�����������λ��    %%%%%%%%%
p = x;%��ʼ��ÿ����������λ��
pbest = ones(1,N);%��ʼ�����и������Ž�
options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% ����ODE�������ѡ��
t_span = [10,12,14,17,19,21,24,26,28];% ����t���
y0 = [7.2e7,5000,4*10^7,10^8,0];%��ֵ
soly = ones(length(t_span),N);%��ʼһ������

for i=1:N
    par.Kb = x(i,1); par.beta2 = x(i,2); par.zeta=x(i,3);  par.rho=x(i,4);par.dB=x(i,5);
    [~, jie] = ode45(@(t,y) ODE_fit_PD1IL2v(t,y,par), t_span,y0, options);%��΢�ַ���
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
        par.Kb = x(j,1); par.beta2 = x(j,2); par.zeta=x(j,3);  par.rho=x(j,4);par.dB=x(j,5);
        [~, jie] = ode45(@(t,y) ODE_fit_PD1IL2v(t,y,par), t_span,y0, options);%��΢�ַ���
        soly(:,j) = jie(:,1);%��ÿ������ϸ����������soly������
        zhi(1,j) = ObjFun(data,soly(:,j));%��ÿһ�������������
        
        if  zhi(:,j)<pbest(:,j)
            p(j, :) = x(j, :);%���¸�������λ��
            pbest(:,j)= zhi(:,j);%���¸�������ֵ
        end
        
        [zho,ind] = min(pbest);%zhoΪÿһ�ε������º�����Ž�
        gbest = zho;%����ȫ������ֵ
        g = p(ind,:);%����ȫ������λ��
    end
end
%%
%%%%%%%%%%%%%%%%  ������ϸ��ͼ   %%%%%%%%%%
% Save data
dlmwrite('Data/data_BestFitPara_PD1IL2V.csv', g, 'delimiter', ',');
dlmwrite('Data/data_FitGood_R2_PD1IL2V.csv', gbest, 'delimiter', ',');

par.Kb = g(1);
par.beta2 = g(2);
par.zeta = g(3);
par.rho = g(4);
par.dB = g(5);

options = odeset('RelTol', 1e-6, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]); 
t1_span=(0:0.5:40);
[t, R] = ode45(@(t,y) ODE_fit_PD1IL2v(t,y,par), t1_span,y0, options); 
C = R(:,1);  % Tumor cells
x2 = data(:,1);
y2 = data(:,2);
plot(t1_span,C,'red-',x2,y2,'black o-')
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')
legend('fitted curve','data point');
