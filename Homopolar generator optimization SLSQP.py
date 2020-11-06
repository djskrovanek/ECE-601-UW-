from scipy.optimize import minimize
import numpy as np


def objective(x):
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]

    return x[0] + x[1] + x[2] + x[4]  # try to minimize all of the geometric variables


def constraint1(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    sigma_cu = 5.96 * 10 ** 7  # conductivity of copper
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]

    B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
    K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant

    R_brushes = 0.01 / x[9]  # resistance of brushes
    R_disk = np.log(x[4] / x[5]) / (2 * np.pi * x[3] * sigma_cu)  # resistance of disk
    return 6.1 - K * omega - 500 * (R_disk + R_brushes)


def constraint2(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    sigma_cu = 5.96 * 10 ** 7  # conductivity of copper
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]

    B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
    K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant

    R_brushes = 0.01 / x[9]  # resistance of brushes
    R_disk = np.log(x[4] / x[5]) / (2 * np.pi * x[3] * sigma_cu)  # resistance of disk
    return K * omega - 500 * (R_disk + R_brushes) - 5.9


def constraint3(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]
    return 1.5 - N_f * I_f * u_0 / (g - w)


def constraint4(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]
    return N_f * I_f * u_0 / (g - w) - 0.0001


def constraint5(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    sigma_cu = 5.96 * 10 ** 7  # conductivity of copper
    I_0 = 500
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]

    B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
    K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant

    r_w = np.sqrt(I_f / (4 * 10 ** 6 * np.pi))  # radius of winding wire
    N_layer = 2 * N_f * r_w / (c - 2 * b - g)  # number of winding layers needed
    l_wire = sum([2 * np.pi * N_f * (r_o + ((i - 1) * 2 * r_w))
                  for i in range(1, round(N_layer))])  # length of field winding

    R_coil = l_wire / (np.pi * r_w ** 2 * sigma_cu)  # resistance of coil winding wire
    R_brushes = 0.01 / N_br  # resistance of brushes
    R_disk = np.log(r_o / r_i) / (2 * np.pi * w * sigma_cu)  # resistance of disk

    P_loss = I_0 ** 2 * (R_disk + R_brushes) + I_f ** 2 * (R_coil)  # Ohmic losses
    P_mech = K * omega * I_0  # mechanical power input
    P_in = P_mech + P_loss
    P_out = 3000
    eta = P_out / P_in
    return eta - 0.85  # eta > 85%


def constraint6(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    I_0 = 500
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]
    B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
    K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant
    P_mech = K * omega * I_0  # mechanical power input
    return 3001 - P_mech  # mechanical power input must be 3000 W


def constraint7(x):
    u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
    I_0 = 500
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]
    B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
    K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant
    P_mech = K * omega * I_0  # mechanical power input
    return P_mech - 2999  # mechanical power input must be 3000 W


def constraint8(x):
    return x[6] - x[5]  # r_o > r_i


def constraint9(x):
    return x[2] - 2 * x[1] - x[3]  # c > 2b + g


def constraint10(x):
    return x[3] - x[7] - 10 ** -3  # g > w + 10^-3


def constraint11(x):
    return x[4] / 2 - x[6] - x[1]  # h/2 > r_o + b


def constraint12(x):
    return x[4] - x[2]  # h > c


def constraint13(x):
    return x[4] - x[0]  # h > a


def constraint14(x):
    return x[0] - x[1]  # a > b


def constraint15(x):
    return x[2] - x[5]  # c > r_i


def constraint16(x):
    return x[2] - x[3]  # c > g


def constraint17(x):
    return x[0] - 2 * x[6] - 0.05  # a > 2r_o + 0.05


def constraint18(x):
    return x[1] - x[3]  # b > g


def constraint19(x):
    a = x[0]
    b = x[1]
    c = x[2]
    g = x[3]
    h = x[4]
    r_i = x[5]
    r_o = x[6]
    w = x[7]
    I_f = x[8]
    N_br = x[9]
    N_f = x[10]
    omega = x[11]
    return ((8 * 10 ** 6) * 0.4 * (h / 2 - r_o - b) * 0.5 * (c - 2 * b - g)) - (N_f * I_f)


# initial guesses
n = 12
x0 = np.zeros(n)
x0[0] = 0.5  # a
x0[1] = 0.06  # b
x0[2] = 0.5  # c
x0[3] = 0.01  # g
x0[4] = 0.8  # h
x0[5] = 0.008  # r_i
x0[6] = 0.2  # r_o
x0[7] = 0.008  # w
x0[8] = 3.8  # I_f
x0[9] = 90  # N_brushes
x0[10] = 300  # N_f
x0[11] = 250  # omega

# show initial objective
print('Initial SSE objective: ' + str(objective(x0)))

# bounds for variables
b0 = (0.001, 1)  # a
b1 = (0.0001, 0.1)  # b
b2 = (0.001, 1)  # c
b3 = (0.00001, 0.1)  # g
b4 = (0.001, 1)  # h
b5 = (0.00001, 0.01)  # r_i
b6 = (0.001, 1)  # r_o
b7 = (0.00001, 0.01)  # w
b8 = (0.01, 10)  # I_f
b9 = (1, 150)  # N_brushes
b10 = (1, 1000)  # N_f
b11 = (1, 314)  # omega
bnds = (b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11)

con1 = {'type': 'ineq', 'fun': constraint1}
con2 = {'type': 'ineq', 'fun': constraint2}
con3 = {'type': 'ineq', 'fun': constraint3}
con4 = {'type': 'ineq', 'fun': constraint4}
con5 = {'type': 'ineq', 'fun': constraint5}
con6 = {'type': 'ineq', 'fun': constraint6}
con7 = {'type': 'ineq', 'fun': constraint7}
con8 = {'type': 'ineq', 'fun': constraint8}
con9 = {'type': 'ineq', 'fun': constraint9}
con10 = {'type': 'ineq', 'fun': constraint10}
con11 = {'type': 'ineq', 'fun': constraint11}
con12 = {'type': 'ineq', 'fun': constraint12}
con13 = {'type': 'ineq', 'fun': constraint13}
con14 = {'type': 'ineq', 'fun': constraint14}
con15 = {'type': 'ineq', 'fun': constraint15}
con16 = {'type': 'ineq', 'fun': constraint16}
con17 = {'type': 'ineq', 'fun': constraint17}
con18 = {'type': 'ineq', 'fun': constraint18}
con19 = {'type': 'ineq', 'fun': constraint19}
cons = ([con1, con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12, con13, con14, con15, con16, con17, con18, con19])
print('Iteration in progress')
solution = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons, options={'maxiter':14})
x = solution.x
print('Final SSE Objective: ' + str(objective(x)))
print(str(solution))
print('Solution')
print('a = ' + str(x[0]))
print('b = ' + str(x[1]))
print('c = ' + str(x[2]))
print('g = ' + str(x[3]))
print('h = ' + str(x[4]))
print('r_i = ' + str(x[5]))
print('r_o = ' + str(x[6]))
print('w = ' + str(x[7]))
print('I_f = ' + str(x[8]))
print('N_brushes = ' + str(x[9]))
print('N_f = ' + str(x[10]))
print('omega = ' + str(x[11]))

#############################################
u_0 = 4 * np.pi * 10 ** -7  # permeability of free space
sigma_cu = 5.96 * 10 ** 7  # conductivity of copper
I_0 = 500
a = x[0]
b = x[1]
c = x[2]
g = x[3]
h = x[4]
r_i = x[5]
r_o = x[6]
w = x[7]
I_f = x[8]
N_br = x[9]
N_f = x[10]
omega = x[11]

B = u_0 * N_f * I_f / (g - w)  # magnetic field of excitation coil
K = B / 2 * (r_o ** 2 - r_i ** 2)  # machine constant
r_w = np.sqrt(I_f / (4 * 10 ** 6 * np.pi))  # radius of winding wire
N_layer = 2 * N_f * r_w / (c - 2 * b - g)  # number of winding layers needed
l_wire = sum([2 * np.pi * N_f * (r_o + ((i - 1) * 2 * r_w))
              for i in range(1, round(N_layer))])  # length of field winding

R_coil = l_wire / (np.pi * r_w ** 2 * sigma_cu)  # resistance of coil winding wire
R_brushes = 0.01 / N_br  # resistance of brushes
R_disk = np.log(r_o / r_i) / (2 * np.pi * w * sigma_cu)  # resistance of disk

P_loss = I_0 ** 2 * (R_disk + R_brushes) + I_f ** 2 * (R_coil)  # Ohmic losses
P_mech = K * omega * I_0  # mechanical power input
P_in = P_mech + P_loss
P_out = 3000
eta = P_out / P_in
r_d = r_o + 0.05
l_p = 0.5 * (c - 2 * b - g)
w_f = h / 2 - r_o - b
M_Fe = (7.87 * 100 ** 3) * ((2 * a * b * h) + (4 * (c / 2 - b) * a * b) + (2 * np.pi * r_o ** 2 * l_p) - (2 * np.pi * r_i ** 2 * (l_p + b)))
M_Cu = (8.92 * 100 ** 3) * ((np.pi * r_d ** 2 * w) - (np.pi * r_i ** 2 * w))

print('B = ' + str(B))
print('KVL loop: ' + str(K * omega - 500 * (R_disk + R_brushes)))
print('eta = ' + str(eta))
print('Mass Fe = ' + str(M_Fe / 1000) + ' kg')
print('Mass Cu = ' + str(M_Cu / 1000) + ' kg')
print('d_wire ' + str(r_w * 1000 * 2) + ' mm')
print('K = ' + str(K))
print('Torque = ' + str(500 * B * 0.5 * (r_o ** 2 - r_i ** 2)) + ' N-m')
print('R_brushes = ' + str(R_brushes))
print('R_disk = ' + str(R_disk))
print('N_layers = ' + str(round(N_layer)))
print('P_in = ' + str(K * omega * 500))
print('l_p = ' + str(l_p))
print('w_f = ' + str(w_f))
