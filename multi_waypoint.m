clc
clear all
close all
% prac
%%  定义优化问题实例和优化变量
opti = casadi.Opti();
waypoints = [0 20 100];
M=size(waypoints,2);
N =200;
tN=opti.variable(1);
dt=tN/N;
u = opti.variable(N,1);
x = opti.variable(N+1,2);
mu = opti.variable(N+1,M);%=N+1 判断靠近的公式好计算
v = opti.variable(N+1,M);
lamda = opti.variable(N+1,M);
%% 优化目标
opti.minimize(tN);
%% 初始化

opti.set_initial(tN, 0);% tN初始化为0/2可以解出来，初始化成1或5解不出来

x_ini = 0.5*ones(N+1,2);
x_ini(:,2) = 0.5*zeros(N+1,1);
% x_ini(:,1)=[0:1:N];

x_ini(:,1)= waypoints(1):(waypoints(end)-waypoints(1))/N:waypoints(end);

opti.set_initial(x, x_ini); 
% opti.set_initial(x, 0); %如果全初始化成0，会求解上千次
opti.set_initial(mu, 0);
opti.set_initial(u, 0);
% 
% lamda_ini = 0.5*ones(N+1,1);
% lamda_ini(1:10)=1;
% lamda_ini(11:end)=0;
opti.set_initial(lamda, 1);
opti.set_initial(v, 0);
%% 定义G，即constraints
f = @(x,u) [x(2),u];
for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   k1 = f(x(k,:),         u(k,:));
   k2 = f(x(k,:)+dt/2*k1, u(k,:));
   k3 = f(x(k,:)+dt/2*k2, u(k,:));
   k4 = f(x(k,:)+dt*k3,   u(k,:));
   x_next = x(k,:) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(x(k+1,:)==x_next); % close the gaps
end

%公式13第一个式子对应的约束
for j =1:N  
  opti.subject_to(mu(j,:) - lamda(j,:) + lamda(j+1,:) == 0);
end

%公式13最后一行第二个式子对应的约束
for i =1:N+1 
    for j =1:M-1 
        opti.subject_to(lamda(i,j) - lamda(i,j+1) <= 0);
    end
end


%输入的限制
for j =1:N  
  opti.subject_to(-5 <= u(j) <= 5);
end
opti.subject_to( u(1) == 0);
%公式13最后一行第一个式子对应的约束
% 这个准确讲应该改成0/1
% for j =1:N  
%   opti.subject_to(0<=mu(j));
% end
for i =1:N  
    for j =1:M  
        opti.subject_to(0<=mu(i,j));
    end
end

%公式14对应的约束
% for j =1:N  
%   opti.subject_to(mu(j)*((x(j)-10)*(x(j)-10)-v(j))==0);
% end
for i =1:N+1  %从2开始还是1开始？
    for j =1:M  
        opti.subject_to(mu(i,j)*((x(i,1)-waypoints(j))*(x(i,1)-waypoints(j))-v(i,j))==0);
    end
end
% for j =1:N  
%   opti.subject_to(lamda(j)>=0); 
% end

d = opti.parameter();
% for j =1:N  
%   opti.subject_to(0 <= v(j) <= d);
% end

opti.subject_to(0 <= v <= d);
opti.set_value(d,0.2);

  opti.subject_to([x(1,1);x(1,2)] == [waypoints(1);0]);
% opti.subject_to(x(1,1)==waypoints(1)); 
% opti.subject_to(x(1,2)==0); 
 opti.subject_to([x(N+1,1);x(N+1,2)] == [waypoints(end);0]);
% opti.subject_to(x(N,1)==waypoints(end)); 
% opti.subject_to(x(N,2)==0); 
opti.subject_to([lamda(1,:);lamda(N+1,:)] == [1*ones(1,M);zeros(1,M)]);
opti.subject_to(0 <= tN); 
opti.solver('ipopt',struct('print_time',false),struct('max_iter',20000)); 
sol = opti.solve(); 

%% 后处理
tN2=sol.value(tN);
fprintf('the value of tN is%6.2f\n',tN2)

m=tN2/(N-1);
t = 0: m:tN2;

subplot(2,3,1);
hold on
plot(t,sol.value(u),'-o');
legend('u');
title('加速度');


m2=tN2/(N);
t = 0: m2:tN2;

subplot(2,3,2);
plot(t,sol.value(mu(:,1)),'-d');
hold on;
plot(t,sol.value(mu(:,2)),'-d');
hold on;
plot(t,sol.value(mu(:,3)),'-d');
legend('mu1','mu2','mu3');
title('一直等于0，偶尔等于1');

subplot(2,3,3);
mid=sol.value(x(:,1));
mid=mid';
plot(t,mid,'-o');
legend('x1');
title(['位置是',num2str(waypoints)]);


subplot(2,3,4);
plot(t,sol.value(lamda(:,1)),':*');
hold on;
plot(t,sol.value(lamda(:,2)),':*');
hold on;
plot(t,sol.value(lamda(:,3)),':*');
legend('lamda1','lamda2','lamda3');
title('从全是1变成全是0');

subplot(2,3,5);
% plot(t,sol.value(v),'-kX');
plot(t,sol.value(v(:,1)),'-kX');
hold on;
plot(t,sol.value(v(:,2)),'-kX');
hold on;
plot(t,sol.value(v(:,3)),'-kX');
legend('v');
title('应该一直挺小的');

subplot(2,3,6);
plot(t,sol.value(x(:,2)),'-rX');
legend('x2');
title('速度');
