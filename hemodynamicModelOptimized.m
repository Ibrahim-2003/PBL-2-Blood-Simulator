function hemodynamicModelOptimized
V_lv0 = 500;
V_rv0 = 400;

y0 = zeros(40,1);
y0(1) = 1;
y0(21) = V_lv0;
y0(22) = V_rv0;

t0 = 0;
tf = 1;

tspan = [t0 tf];

[ts,ys] = ode45(@bloodFlowModel,tspan,y0);
fprintf('y(T) = %g, y''(T) = %g\n',ys(end,1),ys(end,2))

disp('       t           y1          y2')
disp([ts,ys])                     % table of t and y values of solution

for i = [1,5,6,10,11,12]
plot(ts,ys(:,i),'b')              % plot solution y(t)
title('Solution y(t) of IVP')
xlabel('t'); grid on
figure;
end
end

function dy = bloodFlowModel(t,y)
%Define constants

%Displacement Parameters 
K_st_a = 2.5;
K_st_v = 20.0;
K_f = 0.0004;
K_st = 9000.0;
M = 0.0004;
A = 0.00047;
%Valve parameters 
K_p_ao = 5500;
K_f_ao = 50;
K_p_mi = 6000;
K_f_mi = 40;
K_p_po = 5000;
K_f_po = 60;
K_p_ti = 5800;
K_f_ti = 45;


%Parameters for the heart
CQ_ao = 350;
CQ_mi = 400;
P_lv_0 = 1.0;
V_lv_0 = 5.0;
P_la_0 = 1.0;
V_la_0 = 4.0;

CQ_po = 350;
CQ_ti = 400;
P_rv_0 = 1.0;
V_rv_0 = 10;
P_ra_0 = 1.0;
V_ra_0 = 4.0;

%Parameters for blood vessels 
C_sas = 0.08;
R_sas = 0.003;
L_sas = 0.000062;
C_sat = 1.6;
R_sat = 0.05;
L_sat = 0.0017;
R_sar = 0.5;
R_scp = 0.52;
R_svn = 0.075;
C_svn = 20.5;
V_lv0 = 500;

C_pas = 0.18;
R_pas = 0.002;
L_pas = 0.000052;
C_pat = 3.8;
R_pat = 0.01;
L_pat = 0.0017;
R_par = 0.05;
R_pcp = 0.25;
R_pvn = 0.006;
C_pvn = 20.5;
V_rv0 = 400;

dy = zeros(40,1);
%Calculate elastance values
e = elastance(t);

%Pressures in the heart
y(11) = P_lv_0 + e(1)*(y(21) - V_lv_0);
y(16) = P_rv_0 + e(2)*(y(23) - V_rv_0);
y(20) = P_la_0 + e(3)*(y(22) - V_la_0);
y(15) = P_ra_0 + e(4)*(y(24) - V_ra_0);

%Blood volumes derived from diaphram displacement 
dy(33) = y(34); 
dy(34) = (K_st_a*e(3))/M - (K_st_v*e(1))/M + ...
                    (y(11)-y(20))*A/M - (K_f/M)*y(34) - ...
                    (K_st/M)*y(33);

dy(35) = y(36);
dy(36) = (K_st_a*e(3))/M - (K_st_v*e(1))/M + ...
                    (y(16)-y(15))*A/M - (K_f/M)*y(36) - ...
                    (K_st/M)*y(35);


y(21) = y(21) + A*y(33);
y(22) = y(22) - A*y(33);
y(23) = y(23) + A*y(35);
y(24) = y(24) - A*y(35);

%Flow rates
y(4) = (y(14) - y(15))/R_svn;
y(9) = (y(19) - y(20))/R_pvn;

%Blood volume in each chamber of the heart
dy(21) = y(10) - y(1);
dy(22) = y(9) - y(10);
dy(23) = y(5) - y(6);
dy(24) = y(4) - y(5);

% Blood flow through the valves
dy(25) = y(26);
dy(26) = K_p_ao * (y(11)-y(12)) * cos(y(25)) - K_f_ao*y(26);
dy(27) = y(28);
dy(28) = K_p_mi * (y(20)-y(11)) * cos(y(27)) - K_f_mi*y(28);
dy(29) = y(30);
dy(30) = K_p_po * (y(16)-y(17)) * cos(y(29)) - K_f_po*y(30);
dy(31) = y(32);
dy(32) = K_p_ti * (y(15)-y(16)) * cos(y(31)) - K_f_ti*y(32);

y(37) = ((1-cos(y(25))).^2/(1-cos(pi/2)).^2);
y(38) = ((1-cos(y(27))).^2/(1-cos(pi/2)).^2);
y(39) = ((1-cos(y(29))).^2/(1-cos(pi/2)).^2);
y(40) = ((1-cos(y(31))).^2/(1-cos(pi/2)).^2);

y(1) = CQ_ao*y(37)*sqrt(abs(y(11) - y(12)));
y(10) = CQ_mi*y(38)*sqrt(abs(y(20) - y(11)));
y(6) = CQ_po*y(39)*sqrt(abs(y(16) - y(17)));
y(5) = CQ_ti*y(40)*sqrt(abs(y(15) - y(16)));

%Differential pressures and flow rates within systemic circulation 
dy(12) = (y(1) - y(2))/C_sas;
dy(13) = (y(2) - y(3))/C_sat;
dy(14) = (y(3) - y(4))/C_svn;

dy(2) = (y(12) - y(13) - R_sas *y(3))/L_sas;
dy(3) = (y(12) - y(14) - (R_sat + R_sar + R_scp)*y(3))/L_sat;

%Differential pressures within pulmonary circultion
dy(17) = (y(6) - y(7))/C_pas;
dy(18) = (y(7) - y(8))/C_pat;
dy(19) = (y(8) - y(9))/C_pvn;

dy(7) = (y(17) - y(18) - R_pas * y(7))/L_pas;
dy(8) = (y(17) - y(19) -(R_pat + R_par + R_pcp) * y(7))/L_pat;

end 

function theta = valveAngle(t,P_lv, P_sas,P_la,P_rv,P_pas,P_ra)


f_ao = @(t,y) [y(2);K_p_ao * (P_lv-P_sas) * cos(y(1)) - K_f_ao*y(2)];
f_mi = @(t,y) [y(2);K_p_mi * (P_la-P_lv) * cos(y(1)) - K_f_mi*y(2)];
f_po = @(t,y) [y(2);K_p_po * (P_rv-P_pas) * cos(y(1)) - K_f_po*y(2)];
f_ti = @(t,y) [y(2);K_p_ti * (P_ra-P_rv) * cos(y(1)) - K_f_ti*y(2)];

t0 = 0;

if (t == 0) 
    theta = [0;0;0;0];
    return
end 

tf = t;

y0 = [0;0];

[~,ys_ao] = ode45(f_ao,[t0 tf],y0);
[~,ys_mi] = ode45(f_mi,[t0 tf],y0);
[~,ys_po] = ode45(f_po,[t0 tf],y0);
[~,ys_ti] = ode45(f_ti,[t0 tf],y0);

theta = [ys_ao(end,1); ys_mi(end,1); ys_po(end,1);ys_ti(end,1)];
end

function [l_L,l_R] = displacement(t,P_lv,P_la,P_rv,P_ra)
%Displacement from the left and right side of the heart 

if (t == 0)
l_L = 0;
l_R = 0;
return

end
K_st_a = 2.5;
K_st_v = 20.0;
K_f = 0.0004;
K_st = 9000.0;
M = 0.0004;
A = 0.00047;

e = elastance(t);

f1 = @(t,y) [y(2); (K_st_a*e(3))/M - (K_st_v*e(1))/M + ...
                    (P_lv-P_la)*A/M - (K_f/M)*y(2) - ...
                    (K_st/M)*y(1)];

f2 = @(t,y) [y(2); (K_st_a*e(3))/M - (K_st_v*e(1))/M + ...
                    (P_rv-P_ra)*A/M - (K_f/M)*y(2) - ...
                    (K_st/M)*y(1)];
y0 = [0;0];
t0 = 0;
tf = t;

[~,ys1] = ode45(f1,[t0,tf],y0);
[~,ys2] = ode45(f2,[t0,tf],y0);

l_L = ys1(end,1);
l_R = ys2(end,1);

end

function e = elastance(t)


%Elastance(1,1) = left ventricle
%Elastance(2,1) = right ventricle
%Elastance(3,1) = left atrium
%Elastance(4,1) = right atrium

T_s1 = 0.3;
T_s2 = 0.45;
T_pwb = 0.92;
T_pww = 0.09;
T = 1.0;
E_lv_s = 2.5;
E_lv_d = 0.1;
E_la_max = 0.25;
E_la_min = 0.15;

E_rv_s = 1.15;
E_rv_d = 0.1;
E_ra_max = 0.25;
E_ra_min = 0.15;


t = mod(t,T);
if t < T_s1
    ebar_v = 1 - cos(t/T_s1 * pi);
elseif t < T_s2
    ebar_v = 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi);
else
    ebar_v = 0;
end

e(1,1) = E_lv_d + (E_lv_s - E_lv_d)/2 * ebar_v;
e(2,1) = E_rv_d + (E_rv_s - E_rv_d)/2 * ebar_v;

if t < T_pwb
    ebar_a = 0;
elseif t < T_pwb + T_pww
    ebar_a = 1 - cos((t-T_pwb)/T_pww * 2 * pi);
else 
    ebar_a = 0;
end

e(3,1) = E_la_min + (E_la_max - E_la_min)/2 * ebar_a;
e(4,1) = E_ra_min + (E_ra_max - E_ra_min)/2 * ebar_a;


end

