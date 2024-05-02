
import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

def runge_kutta(y, t, u, dt, f):
    """ y is the initial value for y
        x is the initial value for x
        dx is the time step in x
        f is derivative of function y(t)
    """
    k1 =   f(y, u, t)
    k2 =  f(y + 0.5 *dt* k1, u, t + 0.5 * dt)
    k3 =  f(y + 0.5 *dt* k2, u, t + 0.5 * dt)
    k4 =  f(y + dt*k3, u, t + dt)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) *dt / 6.

opti = ca.Opti()

waypoints = np.array([[1, 1, 2], [2, 3, 2], [2, 3, 3]])  #each row is a point

print(waypoints[1,0:3])

start_point = np.array([0, 0, 0])
end_point = np.array([1, 2, 3])
N = 200  #seems that N should increments if v decrements   

M = waypoints.shape[0] #waypoint number

tN = opti.variable(1)
u = opti.variable(N, 3)
x = opti.variable(N + 1, 2 * 3)
mu = opti.variable(N, M)
v = opti.variable(N, M)
lamda = opti.variable(N+1, M)

 

#initialization
opti.set_initial(tN, 0)
x_ini = 0.5 * np.ones((N + 1, 2 * 3))
x_ini[:, 3:6] = 0.5 * np.zeros((N + 1, 3))

for aa in range(3):
    x_ini[:, aa] = np.linspace(start_point[aa], end_point[aa], N + 1)

opti.set_initial(x, x_ini)
opti.set_initial(mu, 0)
opti.set_initial(u, 0)
lamda_ini = np.ones((N + 1, M))
lamda_ini[:, -1] = 0
opti.set_initial(lamda, lamda_ini)
opti.set_initial(v, 0)


dt=tN/N

#add constraints:
opti.subject_to(tN >= 0)
#opti.subject_to(tN <= 100)


#opti.subject_to(u[0] == 0)
opti.subject_to(x[0, 0:3] == ca.horzcat(start_point[0], start_point[1], start_point[2]))
opti.subject_to(x[0, 3:6] == 0)
opti.subject_to(x[N, 0:3] == ca.horzcat(end_point[0], end_point[1], end_point[2]))
opti.subject_to(x[N, 3:6] == 0)
#opti.subject_to(x[N, 1] == 0)
opti.subject_to(lamda[0, :] == 1)
opti.subject_to(lamda[N, :] == 0)



def func(y, u, t):
    return ca.horzcat(y[3:6], u)

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
   x_next = runge_kutta(x[k, :], 0, u[k,:], dt, func)
   opti.subject_to(x[k+1, :] == x_next)
 

for j in range(N):
    opti.subject_to(mu[j, :] - lamda[j, :] + lamda[j+1, :] == 0)

for i in range(N+1):
    for j in range(M-1):
        opti.subject_to(lamda[i, j] - lamda[i, j+1] <= 0)

for i in range(N):
    for j in range(M):
       # opti.subject_to(mu[i, j] * ((x[i, 0] - waypoints[j])**2 - v[i, j]) == 0)
        opti.subject_to(mu[i, j] * (ca.norm_2(x[i, :3] - ca.transpose(waypoints[j, :]) )** 2- v[i, j]) == 0)

for i in range(N):
    for j in range(3):
        opti.subject_to(-5 <= u[i,j])
        opti.subject_to(u[i,j] <= 5)


for i in range(N):
    for j in range(M):
        opti.subject_to(mu[i,j]>= 0)
        #opti.subject_to(mu[i,j]<= 1)
        opti.subject_to(v[i,j] >= 0)
        opti.subject_to(v[i,j]<= 0.1**2)

opti.minimize(tN)
opts = {"expand":True, "ipopt.max_iter":200000000}
# opts = {"ipopt.max_iter":200000000}
# opts = {"ipopt.tol":1e-40, "expand":True, "ipopt.max_iter":200000000}
opti.solver('ipopt', opts)
#opti.solver('ipopt', {'print_time': False}, {'max_iter': 100000})

sol = opti.solve()

# tN2 = sol.value(tN)
# print(f'The value of tN is {tN2}')

# t = ca.linspace(0, tN2, N)
# plt.subplot(2, 3, 1)
# plt.plot(t, sol.value(u), '-o')
# plt.legend('u')
# plt.title('Input')

 
# plt.subplot(2, 3, 2)
# for i in range(M):
#     plt.plot(t, sol.value(mu[:, i]), '-d')
# plt.legend(['mu1', 'mu2', 'mu3'])
# plt.title('Always equal to 0, occasionally equal to 1.')

# t = ca.linspace(0, tN2, N+1)
# plt.subplot(2, 3, 3)
# mid = sol.value(x[:, 0]).T
# plt.plot(t, mid, '-o')
# plt.legend('x1')
# plt.title(f'Waypoint {waypoints}')

# plt.subplot(2, 3, 4)
# for i in range(M):
#     plt.plot(t, sol.value(lamda[:, i]), ':*')
# plt.legend(['lamda1', 'lamda2', 'lamda3'])
# plt.title('Transitioning from all 1s to all 0s.')

# t = ca.linspace(0, tN2, N)
# plt.subplot(2, 3, 5)
# for i in range(M):
#     plt.plot(t, sol.value(v[:, i]), '-kX')
# plt.legend('v')
# plt.title('Small positive numbers')

# t = ca.linspace(0, tN2, N+1)
# plt.subplot(2, 3, 6)
# plt.plot(t, sol.value(x[:, 1]), '-rX')
# plt.legend('x2')
# plt.title('x_2, velocity')

# plt.show()


# 绘制结果
tN2 = sol.value(tN)
print(f'The value of tN is {tN2}')

t1 = np.linspace(0, tN2, N)
t2 = np.linspace(0, tN2, N + 1)

plt.figure(figsize=(12, 8))

# 绘制 u 和 mu 的子图
plt.subplot(2, 3, 1)
for i in range(3):
    plt.plot(t1, sol.value(u[:, i]), '-o', label=f'u{i+1}')
plt.legend()
plt.title('Input')

plt.subplot(2, 3, 2)
for i in range(M):
    plt.plot(t1, sol.value(mu[:, i]), '-d', label=f'mu{i+1}')
plt.legend()
plt.title('Always 0, occasionally 1.')

#绘制 x 和 lamda 的子图
# for i in range(M):
#     plt.subplot(2, 3, i + 3)
#     plt.plot(t2, sol.value(x[:, i*3:(i+1)*3]), '-o')
#     plt.legend([f'x{i*3+1}', f'x{i*3+2}', f'x{i*3+3}'])
#     plt.title(f'Waypoint {waypoints[:, i]}')

     
plt.subplot(2, 3, 3)
for i in range(3):
    mid = sol.value(x[:, i]).T
    plt.plot(t2, mid, '-o')
    plt.legend('x1')
plt.title(f'Waypoint {waypoints}')

plt.subplot(2, 3, 4)
for i in range(M):
    plt.plot(t2, sol.value(lamda[:, i]), ':*', label=f'lamda{i+1}')
plt.legend()
plt.title('From all 1s to all 0s.')
 
plt.subplot(2, 3, 5)
for i in range(M):
    plt.plot(t1, sol.value(v[:, i]), '-kX')
plt.legend('v')
plt.title('Small positive numbers')

plt.subplot(2, 3, 6)
for i in range(3):
    mid = sol.value(x[:, i+3]).T
    plt.plot(t2, mid, '-o')
    plt.legend('velocity')
plt.title(f'velocity')

plt.tight_layout()
plt.show()

# 绘制 x 的三维轨迹
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
pos = sol.value(x[:, :3])
ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
for i in range(M):
    ax.scatter(waypoints[i,0], waypoints[i,1], waypoints[i,2], marker='*', label=f'Waypoint {i+1}')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.title('3D trajectory')
plt.show()