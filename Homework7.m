clear
clc

%constants
omega = 1500 * 2 * pi / 60;                                                 %angular speed [rad/s]
P_out = 150;                                                                %output power [W]
q = 3;                                                                      %slots/pole []
u_0 = 4 * pi * 10 ^ -7;                                                     %vacuum permeability [H/m]
sigma = 5.96 * 10^7;                                                        %conductivity of copper [S]
K_w = 0.4;                                                                  %form factor of winding []
J = 4 * 10 ^ 6;                                                             %current density [A/m²]
rho_Fe = 7860;                                                              %mass density of iron [kg/m³]
rho_Cu = 8960;                                                              %mass density of copper [kg/m³]

tau = P_out / omega;                                                        %torque [Nm]
R = 84.4976 / 1000;                                                         %resistance of 24-gauge wire [ohm/m]

%set parameters
r_yo = 50 * 10 ^ -3;                                                        %outer radius of yoke [m]
r_yi = 35 * 10 ^ -3;                                                        %inner radius of yoke [m]
r_a = 22 * 10 ^ -3;                                                         %radius of armature
g = 0.7 * 10 ^ -3;                                                          %air gap length [m]
l = 85 * 10 ^ -3;                                                           %depth of armature [m]
w_f = 24 * 10 ^ -3;                                                         %winding face [m]

d = 10.2 * 10 ^ -3;                                                            %slot depth [m]
b_0 = 2.5 * 10 ^ -3;                                                        %slot width [m]
t_s = (2 * pi * r_a - 8 * b_0) / 8 + b_0;                                   %slot + tooth arc length [m]
k_c = t_s / (t_s - b_0 + 4 * g / pi * log(1 + pi / 4 * b_0 / g));           %Carter coefficient []


%unknown parameters
I_a = 1.2;                                                                  %armature current [A]
I_f = 1.2;                                                                  %field current [A]
NaIa = 326.4;
B = tau / (2 * q * NaIa * l * (r_a - d / 2));
NfIf = B * 2 * k_c * g / u_0;
N_a = NaIa / I_a;
N_f = NfIf / I_f;
l_a = N_a * (2 * l + 2 * (2*r_a - 2*d));                                    %length of armature winding [m]
l_f = N_f * (2 * w_f + 2 * l);                                              %length of field winding [m]
R_a = R * l_a;                                                              %armature winding resistance [ohms]
R_f = R * l_f;                                                              %field winding resistance [ohms]

K = q * N_a * l * (r_a - d / 2) * N_f * u_0 / (k_c * g);                    %machine constant [Nm/A]
B = N_f * I_f * u_0 / (2 * k_c * g);                                        %magnetic field [T]
eta = P_out / (P_out + I_a ^ 2 * R_a + I_f ^ 2 * R_f);                      %efficiency [%]
torque = I_a * I_f * K;
V_a = I_a*R_a + K*omega*I_f;
V_f = I_f*R_f;

Fe_a = rho_Fe * (pi*r_a^2*l - 8*d*b_0*l);
Fe_f = rho_Fe * (pi*r_yo^2*l - pi*r_yi^2*l + 2*460*l / 10^6);
Cu_f = rho_Cu * pi*(0.51054/2 * 10^-3)^2*l_f;
Cu_a = rho_Cu * pi*(0.51054/2 * 10^-3)^2*l_a;
density = torque / (Fe_a + Fe_f + Cu_f + Cu_a);