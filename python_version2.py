
import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

def runge_kutta(y, t, u, dt, f):
    """ y is the initial value for y
        x is the initial value for x
        dx is the time step in x
        f is derivative of function y(t)
    """
    k1 = dt * f(y, u, t)
    k2 = dt * f(y + 0.5 * k1, u, t + 0.5 * dt)
    k3 = dt * f(y + 0.5 * k2, u, t + 0.5 * dt)
    k4 = dt * f(y + k3, u, t + dt)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.

opti = ca.Opti()

waypoints = [0, 5, 60]
start_point = -1
end_point = 80

waypoints = [0, 2, 10]
start_point = -1
end_point = 11

waypoints = [0, 20, 50]
start_point = -1
end_point = 59
N = 40


waypoints = [0, 10, 100]
start_point = -1
end_point = 130
N = 70
N = 200
N = 100

waypoints = [0, 2, 10]
start_point = -1
end_point = 11
N = 40

M = len(waypoints)

tN = opti.variable()
u = opti.variable(N)
x = opti.variable(N+1, 2)
mu = opti.variable(N, M)
v = opti.variable(N, M)
lamda = opti.variable(N+1, M)

opti.minimize(tN)
dt=tN/N

opti.subject_to(tN >= 0)
opti.subject_to(tN <= 20)

opti.subject_to(u >= -5)
opti.subject_to(u <= 5)
opti.subject_to(u[0] == 0)
opti.subject_to(x[0, :] == ca.horzcat(start_point, 0))
opti.subject_to(x[N, :] == ca.horzcat(end_point,0))
opti.subject_to(lamda[0, :] == 1)
opti.subject_to(lamda[N, :] == 0)



def func(y, u, t):
    return ca.vertcat(y[1], u)

#k = 0
#y = runge_kutta(ca.transpose(x[k, :]), 0, u[k], dt, func)

#for k in range(N):
#     # for j in ([0,1]):
#         # k1 = ca.horzcat(x[k, 1], u[k])
#         # k2 = [x[k, j] + dt / 2 * k1[0], u[k]]
#         # k3 = [x[k, j] + dt / 2 * k2[0], u[k]]
#         # k4 = [x[k, 1] + dt * k3[0], u[k]]
#         # #x_next = x[k, :] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
          # x_next = x[k, :] + k1
#         opti.subject_to(x[k+1, 1] == x[k, 1] + dt*u[k])
#         opti.subject_to(x[k+1, 0] == x[k, 0] + dt*x[k+1, 1] )

for k in range(N):
   x_next = runge_kutta(ca.transpose(x[k, :]), 0, u[k], dt, func)
   opti.subject_to(x[k+1, :] == ca.transpose(x_next))
 

for j in range(N):
    opti.subject_to(mu[j, :] - lamda[j, :] + lamda[j+1, :] == 0)

for i in range(N+1):
    for j in range(M-1):
        opti.subject_to(lamda[i, j] - lamda[i, j+1] <= 0)

for i in range(N):
    for j in range(M):
        opti.subject_to(mu[i, j] * ((x[i, 0] - waypoints[j])**2 - v[i, j]) == 0)

for i in range(N):
    for j in range(M):
        opti.subject_to(mu[i,j]>= 0)
        opti.subject_to(v[i,j] >= 0)
        opti.subject_to(v[i,j]<= 0.6)

opts = {"expand":True, "ipopt.max_iter":200000000}
#opts = {"ipopt.max_iter":200000000}
opts = {"ipopt.tol":1e-40, "expand":True, "ipopt.max_iter":200000000}
opti.solver('ipopt', opts)
#opti.solver('ipopt', {'print_time': False}, {'max_iter': 100000})

sol = opti.solve()

tN2 = sol.value(tN)
print(f'The value of tN is {tN2}')

t = ca.linspace(0, tN2, N)
plt.subplot(2, 3, 1)
plt.plot(t, sol.value(u), '-o')
plt.legend('u')
plt.title('Input')

 
plt.subplot(2, 3, 2)
for i in range(M):
    plt.plot(t, sol.value(mu[:, i]), '-d')
plt.legend(['mu1', 'mu2', 'mu3'])
plt.title('Always equal to 0, occasionally equal to 1.')

t = ca.linspace(0, tN2, N+1)
plt.subplot(2, 3, 3)
mid = sol.value(x[:, 0]).T
plt.plot(t, mid, '-o')
plt.legend('x1')
plt.title(f'Waypoint {waypoints}')

plt.subplot(2, 3, 4)
for i in range(M):
    plt.plot(t, sol.value(lamda[:, i]), ':*')
plt.legend(['lamda1', 'lamda2', 'lamda3'])
plt.title('Transitioning from all 1s to all 0s.')

t = ca.linspace(0, tN2, N)
plt.subplot(2, 3, 5)
for i in range(M):
    plt.plot(t, sol.value(v[:, i]), '-kX')
plt.legend('v')
plt.title('Small positive numbers')

t = ca.linspace(0, tN2, N+1)
plt.subplot(2, 3, 6)
plt.plot(t, sol.value(x[:, 1]), '-rX')
plt.legend('x2')
plt.title('x_2, velocity')

plt.show()


