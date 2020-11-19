clear
clc

%constants
u_0 = 4*pi*10^-7;
sigma = 5.96 * 10^7;
omega = 188.5;
P_out = 372.85;
alpha = pi/6;
J = 4 * 10^6;
Kw = 0.4;
rho_Fe = 7860;  %kg/m³
rho_Cu = 8960;  %kg/m³

%parameters
I_f = 10;
I_a = 10;
N_a = 60;
N_f = 260;
l = 0.05;
r_ao = 0.154;
r_ai = 0.07;
r_yo = 0.25;
r_yi = 0.2;
tip = 0.01;
V_f = 12;
V_a = 36.2;
w_f = 0.02;


%armature
r_wa = sqrt(I_a / (pi * J));
l_a = 2 * N_a * (l + r_ao - r_ai);
R_a = l_a / (sigma * pi * r_wa^2);

%geometry
g = 0.7 * 10^-3 + r_wa;
r_p = r_yi - g - r_ao - tip;
theta = atan(w_f / (2 * (r_ao + g)));

%excitation field
r_wf = sqrt(I_f / (pi * J));
turnlayer = r_p / r_wf;
N_layer = N_f / turnlayer;
for i = 0 : N_layer
    l_f = 2*l + 2*(w_f + 4*r_wf*i);
end
R_f = l_f / (sigma * pi * r_wf^2);
B = N_f * I_f * u_0 / (2 * g);
A_w = r_yi * (pi/2 - alpha/2 - theta) - ...
    (r_ao + g + tip) * (pi/2 - alpha/2 - theta);

%power conversion
K = N_a * (1 - alpha/pi) * N_f * u_0 * l * r_ao / (4 * g);
eta = P_out / (P_out + I_f^2 * R_f + I_a^2 * R_a);
torque = K * I_f * I_a;

%mass
Fe_a = rho_Fe * (pi*r_ao^2*l - pi*r_ai^2*l);
Fe_f = rho_Fe * (pi*r_yo^2*l - pi*r_yi^2*l + 2*w_f*l*r_p +...
    (pi-pi/3)*(r_ao+g+tip)^2*l - (pi-pi/3)*(r_ao+g)^2*l);
Cu_f = rho_Cu * pi*r_wf^2*l_f;
Cu_a = rho_Cu * pi*r_wa^2*l_a;
density = torque / (Fe_a + Fe_f + Cu_f + Cu_a);

if (J * A_w * Kw) < (N_f * I_f)
    disp('winding area violated')
    disp(num2str(J * A_w * Kw))
    disp(num2str(N_f * I_f))
end
if B > 1.2
    disp('magnetic field too big')
end
if torque > 1.988 
    disp('too much torque')
    elseif torque < 1.968
        disp('not enough torque')
        disp(num2str(1.978 - torque))
end
if eta < 0.9
    disp('inefficient')
end
if (K*omega*I_f - I_a*R_a - V_a > 0.1) || (K*omega*I_f - I_a*R_a - V_a < -0.1)
    disp('KVL violated')
    disp(num2str(K*omega*I_f))
    disp(num2str(I_a*R_a + V_a))
end
if w_f > 2*(r_ao + g) * (pi/2 - alpha/2 - theta)
    disp('saturation at yoke tips')
    disp(num2str(2*(r_ao + g) * (pi/2 - alpha/2 - theta)))
end
if w_f > pi * r_ao
    disp('saturation at armature')
    disp(num2str(pi * r_ao))
end
if w_f > 0.5 * (r_yo - r_yi)
    disp('saturation in yoke body')
    disp(num2str(0.5 * (r_yo - r_yi)))
end
disp(num2str((I_a^2 * R_a)/(I_f^2 * R_f + I_a^2 * R_a)))