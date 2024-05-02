clc
clear all;
close all;
% prac
%% optimization problem instances and optimization variables
opti = casadi.Opti();

waypoints = [ [1;1;2],[2;3;2]];  %each column is a point 
start_point = [0;0;0];
end_point = [1; 2; 3];
M=size(waypoints,2);
N =100;
tN=opti.variable(1);
dt=tN/N;
u = opti.variable(N,3);
x = opti.variable(N+1,2*3);
mu = opti.variable(N,M); 
v = opti.variable(N,M);
lamda = opti.variable(N+1,M);
%% cost function
opti.minimize(tN);
%% initialization
opti.set_initial(tN, 0);  

x_ini = 0.5*ones(N+1,2*3);
x_ini(:,4:6) = 0.5*zeros(N+1,3); 

x_ini_col = zeros(3,N+1);
for aa= 1:3
 x_ini_col(aa,:) = start_point(aa):(end_point(aa)-start_point(aa))/N:end_point(aa);
end
x_ini(:,1:3) = x_ini_col';

opti.set_initial(x, x_ini); 
% opti.set_initial(x, 0); 
opti.set_initial(mu, 0);
opti.set_initial(u, 0);
% 
lamda_ini =  ones(N+1,M);
% lamda_ini(1:10)=1;
lamda_ini(:,end)=0;
opti.set_initial(lamda, lamda_ini);
opti.set_initial(v, 0);
%%constraints
f = @(x,u) [x(4:6),u];
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
% for i =1:N+1  
%     for j =1:3
%        opti.subject_to(-5 <= u(i,j) <= 5);
%     end
% end

opti.subject_to(-5 <= u  <= 5);
% opti.subject_to( u(1,:) == 0);

for i =1:N  
    for j =1:M  
        opti.subject_to(0<=mu(i,j));
    end
end
% opti.subject_to(0<=mu );

%equation 14
% for j =1:N  
%   opti.subject_to(mu(j)*((x(j)-10)*(x(j)-10)-v(j))==0);
% end
for i =1:N   
    for j =1:M  
        opti.subject_to(mu(i,j)*((x(i,1:3)-waypoints(:,j)')*(x(i,1:3)-waypoints(:,j)')'-v(i,j))==0);
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
opti.set_value(d,0.4^2);


  opti.subject_to([x(1,1:3);x(1,4:6)] == [start_point';[0,0,0]]);
% opti.subject_to(x(1,1)==waypoints(1)); 
% opti.subject_to(x(1,2)==0); 
 opti.subject_to([x(N+1,1:3);x(N+1,4:6)] == [end_point';[0,0,0]]);
% opti.subject_to(x(N+1,4:6)== [0,0,0]); 
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
t1 = 0: m:tN2;

m2=tN2/(N);
t2 = 0: m2:tN2;

subplot(2,3,1);
hold on
plot(t1,sol.value(u),'-o');
legend('u');
title('input');


subplot(2,3,2);
for aa = 1:M
    plot(t1,sol.value(mu(:,aa)),'-d');
    hold on;    
    legend(['mu',num2str(aa)]);
end
% legend('mu1','mu2','mu3');
title('Always  0, occasionally 1.');

subplot(2,3,3);
for aa=1:3
    mid=sol.value(x(:,aa));
    mid=mid';
    plot(t2,mid,'-o');
    hold on;
    legend(['x', num2str(aa)]);hold on;
    title(['Waypoint',num2str(waypoints(aa,:))]);
end


subplot(2,3,4);
for aa = 1:M
    plot(t2,sol.value(lamda(:,M)),':*');
    hold on;
    legend(['\lambda',num2str(aa)]);
end
% plot(t2,sol.value(lamda(:,2)),':*');
% hold on;
% plot(t2,sol.value(lamda(:,3)),':*');
% legend('lamda1','lamda2','lamda3');
title('From all 1s to all 0s.');

subplot(2,3,5);
% plot(t,sol.value(v),'-kX');
for aa = 1:M
    plot(t1,sol.value(v(:,M)),'-kX');
    hold on;
	legend(['v',num2str(aa)]);
% plot(t1,sol.value(v(:,2)),'-kX');
% hold on;
% plot(t1,sol.value(v(:,3)),'-kX');
% legend('v');
end
title('small positive numbers');

subplot(2,3,6);
plot(t2,sol.value(x(:,4:6)),'-rX');
legend('x(:,4,6)');
title('velocity');


for aa=1:3
    figure();
    mid=sol.value(x(:,aa));
    mid=mid';
    plot(t2,mid,'-o');
    hold on;
    legend(['x', num2str(aa)]);
    title(['Waypoint',num2str(waypoints(aa,:))]);    
    xlabel('time'); ylabel(['pos, ' num2str(aa)]);
end

figure();
pos = sol.value(x(:,1:3));
plot3(pos(:,1), pos(:,2), pos(:,3));
hold on; 
for bb=1:M
    plot3(waypoints(1,bb), waypoints(2,bb), waypoints(3,bb),'*');
end
xlabel('x'); ylabel('y');zlabel('z');
