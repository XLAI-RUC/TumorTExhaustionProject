clear; clc;

N = 400;    % Ⱥ�����Ӹ���
D = 4;     % ����ά��
T = 200;   % ����������
c1 = 2.05;    % ѧϰ����1
c2 = 2.05;    % ѧϰ����2
wmax = 0.8;   % ����Ȩ�����ֵ
wmin = 0.5;  % ����Ȩ����Сֵ

data = xlsread('Data/data_exp_IFNalphaAntiPDL1.xlsx');
par = setparameter();
par.gamma = 6.25e-6;

%%%%%%%%%     ��ʼ����Ⱥλ������Ⱥ�ٶ�    %%%%%%%%%
xmax =[1e-6,1,1e-6, 1.3e7]; % λ�����ֵ
xmin = [1e-8,0.1,3e-7,1.5e7]; % λ����Сֵ

Ka_LHS = LHS_Call(xmax(1),(xmax(1)+xmin(1))/2,xmin(1), 0,N,'unif');
dA_LHS = LHS_Call( xmax(2),(xmax(2)+xmin(2))/2,xmin(2), 0,N,'unif');
alphaA_LHS = LHS_Call( xmax(3),(xmax(3)+xmin(3))/2,xmin(3), 0,N,'unif');
Tint_LHS = LHS_Call( xmax(4),(xmax(4)+xmin(4))/2,xmin(4), 0,N,'unif');

vmax = 0.01*xmax; % �ٶ����ֵ
vmin = -vmax;   % �ٶ���Сֵ
x = [Ka_LHS,dA_LHS,alphaA_LHS,Tint_LHS];
v = rand(N, D) .* (vmax - vmin) + vmin;%��ʼ��Ⱥ�ٶ�

%%%%%%%%%     ��ʼȫ������λ�����������λ��    %%%%%%%%%
p = x;%��ʼ��ÿ����������λ��
pbest = ones(1,N);%��ʼ�����и������Ž�
options = odeset('RelTol', 1e-3, 'AbsTol', [1e-6,1e-6,1e-6,1e-6,1e-6]);% ����ODE�������ѡ��
t_span = [0,5,10,14,19,22,26,30];% ����t���
soly = ones(length(t_span),N);  %��ʼһ������

for i=1:N
    y0 =[data(1,2),5000,x(i,4),10^7,0]; %C��S��E��T cellsϸ���ĳ�ֵ
    par.Ka = x(i,1); par.dA = x(i,2); par.alphaA = x(i,3);
    [~, jie] = ode45(@(t,y) ODE_fit_INFalphaAntiPDL1(t,y,par), t_span,y0, options);%��΢�ַ���
    soly(:,i) = jie(:,1);  %��ÿ������ϸ����������soly������
    pbest(:,i) = ObjFun(data,soly(:,i)); %��ÿһ�������������
end

[minimum,index] = min(pbest);%��¼����������С��λ���Լ���Сֵ
gbest = minimum;%���и��������СֵΪ��ʼȫ������ֵ
g = p(index,:);%ȫ������ֵ���������Ϊ��ʼȫ������λ��

%%
%%%%%%%%%%  ���չ�ʽ���ε���ֱ�����㾫�Ȼ��ߵ������� %%%%%%%%%
zhi = ones(1,N);%������¼ÿһ�ε����и�������
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
        par.Ka = x(j,1);par.dA = x(j,2);par.alphaA=x(j,3);
        y0 =[data(1,2),5000,x(j,4),1e7,0]; %C��S��E��T cellsϸ���ĳ�ֵ
        [~, jie] = ode45(@(t,y) ODE_fit_INFalphaAntiPDL1(t,y,par), t_span,y0, options);%��΢�ַ���
        soly(:,j)= jie(:,1);%��ÿ������ϸ����������soly������
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
dlmwrite('Data/data_BestFitPara_INFalphaAntiPDL1.csv', g, 'delimiter', ',');
dlmwrite('Data/data_FitGood_R2_INFalphaAntiPDL1.csv', gbest, 'delimiter', ',');

par.Ka = g(1); par.dA = g(2);  par.alphaA = g(3);
 
y0 =[data(1,2),5000,g(4),1e7,0]; %C��S��E��T cellsϸ���ĳ�ֵ
options = odeset('RelTol', 1e-3, 'AbsTol', [1e-3,1e-3,1e-3,1e-3,1e-3]);% ����ODE�������ѡ��
t1_span = (0:0.5:40);
[t, R] = ode45(@(t,y) ODE_fit_INFalphaAntiPDL1(t,y,par), t1_span,y0, options);%��΢�ַ���
C = R(:,1);  E = R(:,2); T_1 = R(:,3); T_2 = R(:,4); x2 = data(:,1); y2 = data(:,2);

plot(t1_span,C,'red-',x2,y2,'black o')
xlabel('Time (days)')
ylabel('Tumor cell volume (mm^3)')
legend('fitted curve','data point');
