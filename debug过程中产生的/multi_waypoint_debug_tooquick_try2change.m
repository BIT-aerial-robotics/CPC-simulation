clc
clear all
close all

%%  �����Ż�����ʵ�����Ż�����
opti = casadi.Opti();
waypoints = [0 5 15];% ������ͷ�����ǲ�������β
M=size(waypoints,2);% M��waypoint�ĸ���
N =50;%����΢�ֶ�������õ���һ�㣬���ڵ���10��Ϊ����Ŀ�һ��
tN=opti.variable(1);
dt=tN/N;
u = opti.variable(N,1);
x = opti.variable(N+1,2);
mu = opti.variable(N+1,M);
v = opti.variable(N+1,M);%�ж��Ƿ񿿽�waypoint��Լ�����ƺ�Ҫ����(N+1)*M��
lamda = opti.variable(N+1,M);
%% �Ż�Ŀ��
opti.minimize(tN);
%% ��ʼ��

opti.set_initial(tN, 2);% tN��ʼ��Ϊ0/2���Խ��������ʼ����1��5�ⲻ����

x_ini = 0.5*ones(N+1,2);
x_ini(:,1)=[0:1:N];
opti.set_initial(x, x_ini);

% һֱ����0��ż������1
opti.set_initial(mu, 0);
opti.set_initial(u, 0);

% lamda ��ȫ��1���ȫ��0
% lamda_ini = 0.5*ones(N+1,1);
% lamda_ini(1:10)=1;
% lamda_ini(11:end)=0;
opti.set_initial(lamda, 1);
opti.set_initial(v, 0);
%% ����G����constraints
f = @(x,u) [x(2),u-x(2)];
for k=1:N-1 % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(x(k,:),         u(k,:));
   k2 = f(x(k,:)+dt/2*k1, u(k,:));
   k3 = f(x(k,:)+dt/2*k2, u(k,:));
   k4 = f(x(k,:)+dt*k3,   u(k,:));
   x_next = x(k,:) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(x(k+1,:)==x_next); % close the gaps
end

%��ʽ13��һ��ʽ�Ӷ�Ӧ��Լ��
for i =1:N  
    for j =1:M  
        opti.subject_to(mu(i,j) - lamda(i,j) + lamda(i,j) == 0);
    end
end

%��ʽ13���һ�еڶ���ʽ�Ӷ�Ӧ��Լ��
for i =1:N  
    for j =1:M  
        opti.subject_to(lamda(i,j) - lamda(i+1,j) >= 0);
    end
end


%���������
for j =1:N  
  opti.subject_to(-5 <= u(j) <= 5);
end
opti.subject_to( u(1) == 0);
%��ʽ13���һ�е�һ��ʽ�Ӷ�Ӧ��Լ��
% ���׼ȷ��Ӧ�øĳ�0/1
for j =1:N  
  opti.subject_to(0<=mu(j));
end
%��ʽ14��Ӧ��Լ��
% for j =1:N  
%   opti.subject_to(mu(j)*((x(j)-10)*(x(j)-10)-v(j))==0);
% end
for i =2:N+1  %i��1��ʼ�����Ǵ�2��ʼ�����ʱ����ȫ��һ��
    for j =1:M  
        opti.subject_to(mu(i,j)*((x(i,1)-waypoints(j))*(x(i,1)-waypoints(j))-v(i,j))==0);
    end
end

% for i =1:N+1  %i��1��ʼ�����Ǵ�2��ʼ�����ʱ����ȫ��һ��
%     for j =1:M  
%         opti.subject_to((mu(i,j)-1)*mu(i,j)==0);
%     end
% end


% for i =1:N+1  %i��1��ʼ�����Ǵ�2��ʼ�����ʱ����ȫ��һ��
%     for j =1:M  
%         opti.subject_to(0<=lamda(i,j)<=1);
%     end
% end

d = opti.parameter();

for j =1:N  
  opti.subject_to(0 <= v(j) <= d);
end

opti.set_value(d,0.02);
opti.subject_to([x(1,1);x(1,2)] == [0;0]);
opti.subject_to([x(N,1);x(N,2)] == [waypoints(end);0]);
opti.subject_to([lamda(1,:);lamda(N+1,:)] == [1;0]);
opti.subject_to(0 <= tN); 
% opti.solver('ipopt',struct('print_time',false),struct('print_level',0,'max_iter',10000)); 
opti.solver('ipopt',struct('print_time',false),struct('max_iter',10000)); 
sol = opti.solve(); 

%% ����
tN2=sol.value(tN);
fprintf('the value of tN is%6.2f\n',tN2)

m=tN2/(N-1);
t = 0: m:tN2;

subplot(2,3,1);
hold on
plot(t,sol.value(u),'-o');
legend('u');
title('���ٶ�');


m2=tN2/(N);
t = 0: m2:tN2;

subplot(2,3,2);
plot(t,sol.value(mu(:,1)),'-d');
hold on;
plot(t,sol.value(mu(:,2)),'-d');
hold on;
plot(t,sol.value(mu(:,3)),'-d');
legend('mu1','mu2','mu3');
title('һֱ����0��ż������1');

subplot(2,3,3);
mid=sol.value(x(:,1));
mid=mid';
plot(t,mid,'-o');
legend('x1');
title('λ��');

subplot(2,3,4);
plot(t,sol.value(lamda(:,1)),':*');
hold on;
plot(t,sol.value(lamda(:,2)),':*');
hold on;
plot(t,sol.value(lamda(:,3)),':*');
legend('lamda1','lamda2','lamda3');
title('��ȫ��1���ȫ��0');

subplot(2,3,5);
plot(t,sol.value(v),'-kX');
legend('v');
title('Ӧ��һֱͦС��');

subplot(2,3,6);
plot(t,sol.value(x(:,2)),'-rX');
legend('x2');
title('�ٶ�');

