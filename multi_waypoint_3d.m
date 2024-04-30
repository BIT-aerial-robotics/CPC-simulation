clc
clear all;
close all;
% prac
%% optimization problem instances and optimization variables
opti = casadi.Opti();
waypoints = [5 40 60];
waypoints = [[0;0;0], [1;1;1], [0; 20; 100]];  %each column is a point 
start_point = 0 ;
end_point = 100;
M=size(waypoints,2);
N =50;
tN=opti.variable(1);
dt=tN/N;
u = opti.variable(N,3);
x = opti.variable(N+1,2*3);
mu = opti.variable(N+1,M);%=N+1 �жϿ���Ĺ�ʽ�ü���?
v = opti.variable(N+1,M);
lamda = opti.variable(N+1,M);
%% cost function
opti.minimize(tN);
%% initialization

opti.set_initial(tN, 0);% tN��ʼ��Ϊ0/2���Խ��������ʼ����?��5�ⲻ����

x_ini = 0.5*ones(N+1,2);
x_ini(:,2) = 0.5*zeros(N+1,1);
% x_ini(:,1)=[0:1:N];

x_ini(:,1)= waypoints(1):(waypoints(end)-waypoints(1))/N:waypoints(end);

opti.set_initial(x, x_ini); 
% opti.set_initial(x, 0); %���ȫ��ʼ����?���������ǧ��?opti.set_initial(mu, 0);
opti.set_initial(u, 0);
% 
lamda_ini =  ones(N+1,M);
% lamda_ini(1:10)=1;
lamda_ini(:,end)=0;
opti.set_initial(lamda, lamda_ini);
opti.set_initial(v, 0);
%% ����G����constraints
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

% equation 13
for j =1:N  
  opti.subject_to(mu(j,:) - lamda(j,:) + lamda(j+1,:) == 0);
end

% equation 13 
for i =1:N+1 
    for j =1:M-1 
        opti.subject_to(lamda(i,j) - lamda(i,j+1) <= 0);
    end
end


% input boundedness
for j =1:N  
  opti.subject_to(-5 <= u(j) <= 5);
end
opti.subject_to( u(1) == 0);
%��ʽ13���һ�е�һ��ʽ�Ӷ�Ӧ��Լ��?% ���׼ȷ��Ӧ�øĳ�?/1
% for j =1:N  
%   opti.subject_to(0<=mu(j));
% end
for i =1:N  
    for j =1:M  
        opti.subject_to(0<=mu(i,j));
    end
end

%equation 14
% for j =1:N  
%   opti.subject_to(mu(j)*((x(j)-10)*(x(j)-10)-v(j))==0);
% end
for i =1:N+1  %��2��ʼ����1��ʼ��
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
opti.set_value(d,0.4);


  opti.subject_to([x(1,1);x(1,2)] == [start_point;0]);
% opti.subject_to(x(1,1)==waypoints(1)); 
% opti.subject_to(x(1,2)==0); 
 opti.subject_to([x(N+1,1);x(N+1,2)] == [end_point;0]);
% opti.subject_to(x(N,1)==waypoints(end)); 
% opti.subject_to(x(N,2)==0); 
opti.subject_to([lamda(1,:);lamda(N+1,:)] == [1*ones(1,M);zeros(1,M)]);
opti.subject_to(0 <= tN); 

% opts["ipopt.tol"] = 1e-4;
% opti.solver('ipopt',struct('print_time',false),struct('max_iter',20000),struct('tol', 1e-4)); 

opti.solver('ipopt',struct('print_time',false),  struct('max_iter',80000)); 
% opti.solver('ipopt',struct('print_time',false) ); 
sol = opti.solve(); 

%% plot
tN2=sol.value(tN);
fprintf('the value of tN is%6.2f\n',tN2)

m=tN2/(N-1);
t = 0: m:tN2;

subplot(2,3,1);
hold on
plot(t,sol.value(u),'-o');
legend('u');
title('input');


m2=tN2/(N);
t = 0: m2:tN2;

subplot(2,3,2);
plot(t,sol.value(mu(:,1)),'-d');
hold on;
plot(t,sol.value(mu(:,2)),'-d');
hold on;
plot(t,sol.value(mu(:,3)),'-d');
legend('mu1','mu2','mu3');
title('Always equal to 0, occasionally equal to 1.');

subplot(2,3,3);
mid=sol.value(x(:,1));
mid=mid';
plot(t,mid,'-o');
legend('x1');
title(['Waypoint',num2str(waypoints)]);


subplot(2,3,4);
plot(t,sol.value(lamda(:,1)),':*');
hold on;
plot(t,sol.value(lamda(:,2)),':*');
hold on;
plot(t,sol.value(lamda(:,3)),':*');
legend('lamda1','lamda2','lamda3');
title('Transitioning from all 1s to all 0s.');

subplot(2,3,5);
% plot(t,sol.value(v),'-kX');
plot(t,sol.value(v(:,1)),'-kX');
hold on;
plot(t,sol.value(v(:,2)),'-kX');
hold on;
plot(t,sol.value(v(:,3)),'-kX');
legend('v');
title('small positive numbers');

subplot(2,3,6);
plot(t,sol.value(x(:,2)),'-rX');
legend('x2');
title('x_2, velocity');
