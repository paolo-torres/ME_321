clc;

% Gear Design

% Initial Givens
P = 20.0; % hp
L = 12000; % hrs
Y = 22; % in
CW = 1.0; % in; clearances + wall thicknesses

% Train Value
TV = 20.25;
TV1 = 4.5;
TV2 = 4.5;

% Number of Teeth
N2 = 16;
N3 = TV1*N2;
N4 = 16;
N5 = TV2*N4;

% Rotational Speeds (rpm)
w2 = 1750;
w3 = (N2/N3)*w2;
w4 = w3;
w5 = w2/(TV1*TV2);
w = [w2 w3 w4 w5];

% Torques (lb*ft)
T2 = (P/w2)*550*(1/(2*pi))*60;
T3 = T2*(w2/w3);
T4 = T3;
T5 = T2*(w2/w5);
T = [T2 T3 T4 T5];

% Diameters (in; D2 = D4 and D3 = D5)
% D3 + D2 + (D2/2 + D3 - D3/2 - D2) + 2a + C < Y
% (3/2)*(N3/Pd)+(1/2)*(N2/Pd)+2*(1.00/Pd)+C < 22
% Pd > 5.62 teeth/in; round to 6 teeth/in
syms Pd;
eqn = (3/2)*(N3/Pd)+(1/2)*(N2/Pd)+2*(1.00/Pd)+CW == 22;
Pd = round(solve(eqn, Pd));
D2 = N2/Pd;
D3 = N3/Pd;
D4 = D2;
D5 = D3;
D = [D2 D3 D4 D5];

% Geometry Factors
J2 = 0.27;
J3 = 0.41;
J4 = J2;
J5 = J3;
J = [J2 J3 J4 J5]; % Fig. 9-10

% Quality Numbers
Av2 = 8;
Av3 = 8;
Av4 = 10;
Av5 = 10;
Av = [Av2 Av3 Av4 Av5];

% Gear Materials (App. 3)
% Gear 2: SAE 3140 OQT 1300; HB = 233
% Gear 3: SAE 1020 Hot-rolled; HB = 111
% Gear 4: SAE 4140 OQT 400; HB = 578
% Gear 5: SAE 1050 OQT 1300; HB = 192
HB_mat = [233 111 578 192];

for index = 1:4
    % Bending Stress (sigma_T)
    Wt = double(T(index)/((D(index)/2)/12)); % lbs
    F = 4; % in
    sigma_T = double((Wt*Pd)/(F*J(index))); % psi

    % Bending Stress Number (s_t)
    Ko = 2.00; % Table 9-1
    Ks = 1.00; % Table 9-2
    Cpf = double((F/(10*D(index)))-0.0375+0.0125*F); % Fig. 9-12
    Cma = double(0.127+0.0158*F-(0.0001093)*F^2); % Fig. 9-13
    Km = double(1.0+Cpf+Cma);
    Kb = 1.0; % Fig. 9-14
    vt = (double(pi*D(index)*w(index))/12);
    B = 0.25*(Av(index)-5.0)^0.667;
    C = 50+56*(1.0-B);
    Kv = double((C/(C+sqrt(vt)))^-B); % Table 9-6
    s_t = double(sigma_T*Ko*Ks*Km*Kb*Kv); % psi

    % Contact Stress Number (s_c)
    Cp = 2300; % Table 9-7
    I = 0.103; % Fig. 9-17
    s_c = double(Cp*sqrt((Wt*Ko*Ks*Km*Kv)/(F*D(index)*I))); % psi

    % Bending and Pitting (s_at, s_ac)
    Kr = 1.00; % Table 9-11
    SF = 1; % Assumed
    q = 1;
    Nc = 60*L*w(index)*q;
    Yn = 1.3558*Nc^-0.0178; % Fig. 9-21
    Zn = 1.4488*Nc^-0.023; % Fig. 9-22
    s_at = s_t*(SF*Kr)/Yn; % psi
    s_ac = s_c*(SF*Kr)/Zn; % psi
    HB_at = ((s_at/1000)-12.80)/0.0773;
    HB_ac = ((s_ac/1000)-29.10)/0.322;
    HB = max(HB_at, HB_ac);

    % Safety Factor
    s_at_mat = 0.0773*HB_mat(index)+12.80; % ksi
    s_ac_mat = 0.322*HB_mat(index)+29.10; % ksi 
    SF_t = ((s_at_mat*1000)/s_t)*(Yn/Kr);
    SF_c = ((s_ac_mat*1000)/s_c)*(Zn/Kr);

    % Power Transmitting Capacity
    P_cap_t = double((s_at_mat*1000*Yn*F*J(index)*w(index)*D(index))/(126000*Pd*SF*Kr*Ko*Ks*Km*Kb*Kv));
    P_cap_c = double(((w(index)*F*I)/(126000*Ko*Ks*Km*Kv))*((s_ac_mat*1000*D(index)*Zn)/(SF*Kr*Cp))^2);
end

% Shaft Design

% Torque
T_shaft = T3*12;

% Forces
syms RLx WT3 WT4 RRx % Top View
syms RLy WR3 WR4 RRy % Front View

% Top View Forces (lbs)
WT3 = double(T3*12/(D3/2));
WT4 = double((T4*12)/(D4/2));
RRx = double((-WT3*3.125+WT4*10.625)/13.75);
RLx = double(-WT3+WT4-RRx);

% Front View Forces (lbs)
WR3 = WT3*tan(20*pi/180);
WR4 = WT4*tan(20*pi/180);
RRy = (WR3*3.125+WR4*10.625)/13.75;
RLy = WR3+WR4-RRy;

% Moments (lb in)
M_front = [0 1103 2276 0];
M_top = [0 422.0 5487 0];
M_front_max = max(M_front);
M_top_max = max(M_top);
M_max = sqrt((M_front_max^2)+(M_top_max^2));

% Design Equation Parameters
N = 2.5; % Asssumed
Kt = 3.0; % Asssumed
sy = 251; % Appendix 3
sn = 96; % Table 5-2
Cr = 0.81; % Table 5-3
Cs = 0.81493128507815697523532994381114; % Table 5-4
sn_prime = (sn*1000)*Cr*Cs; % psi
D = ((32*N/pi)*sqrt(((Kt*M_max/sn_prime)^2)+((3/4)*(T_shaft/(sy*1000))^2)))^(1/3); % in
D_actual = D*1.06; % in
D_pref = round(D_actual); % Table A2-1
