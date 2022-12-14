function BloodModel
%Additional Parameters 
dt = 0.0001;
T = 1;
T_s1 = 0.3;
T_s2 = 0.45;
T_pwb = 0.91;
T_pww = 0.09;

%Parameters for the heart
CQ_ao = 350;
CQ_mi = 400;
E_lv_s = 2.5;
E_lv_d = 0.1;
P_lv_0 = 1.0;
V_lv_0 = 5.0;
E_la_max = 0.25;
E_la_min = 0.15;
P_la_0 = 1.0;
V_la_0 = 4.0;

CQ_po = 350;
CQ_ti = 400;
E_rv_s = 1.15;
E_rv_d = 0.1;
P_rv_0 = 1.0;
V_rv_0 = 10;
E_ra_max = 0.25;
E_ra_min = 0.15;
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
V_lv0 = 800;

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
V_rv0 = 500;

%Parameters for variable valve opening model
K_p_ao = 5500;
K_f_ao = 50;
K_p_mi = 6000;
K_f_mi = 40;
K_p_po = 5000;
K_f_po = 60;
K_p_ti = 5800;
K_f_ti = 45;

t0 = 0;  y0 = [0;0];
[ts,ys] = ode45(@test,[t0,T],y0);
end

function y = test(t,y)
Plv_0 = 1.0;
Pla_0 = 1.0;
k_stla = 2.5;
k_stlv = 20.0;
Msav = 0.0004;
Asav = 0.00047;
k_fsav = 0.0004;
k_esav = 9000.0;
x = theta(t,1,0);
y = [y(2); (k_stla)/Msav - (k_stlv)/Msav + ...
                    (Plv_0-Pla_0)*Asav/Msav - (k_fsav/Msav)*y(2) - ...
                    (k_esav/Msav)*y(1)*x];
end


function dy = bloodFlowModel(t,y)
%Calculate elastance values
e = elastance(t);
%Pressures in the heart
P_lv = P_lv_0 + e(1)*(V_lv - V_lv_0);
P_rv = P_rv_0 + e(2)*(V_rv - V_rv_0);
P_la = P_la_0 + e(3)*(V_la - V_la_0);
P_ra = P_ra_0 + e(4)*(V_ra - V_ra_0);

%Flow rates
Q_svn = (P_svn - P_ra)/R_svn;
Q_pvn = (P_pvn - P_la)/R_pvn;

% Blood volume in each chamber of the heart
dV_lv = Q_mi - Q_ao;
dV_la = Q_pvn - Q_mi;
dV_rv = Q_ti - Q_po;
dV_ra = Q_svn - Q_ti;

%Differential pressures and flow rates within systemic circulation 
dP_sas = (Q_ao - Q_sas)/C_sas;
dP_sat = (Q_sas - Q_sat)/C_sat;
dP_svn = (Q_sat - Q_svn)/C_svn;

dQ_sas = (P_sas - P_sat - R_sas *Q_sat)/L_sas;
dQ_sat = (P_sas - P_svn - (R_sat + R_sar + R_scp)*Q_sat)/L_sat;

%Differential pressures within pulmonary circultion
dP_pas = (Q_po - Q_pas)/C_pas;
dP_pat = (Q_pas - Q_pat)/C_pat;
dP_pvn = (Q_pat - Q_pvn)/C_pvn;

dQ_pas = (P_pas - P_pat - R_pas * Q_pas)/L_pas;
dQ_pat = (P_pas - P_pvn -(R_pat + R_par + R_pcp) * Q_pas)/L_pat;


end 

function y = theta(t,P_lv, P_sas,P_la,P_rv,P_pas,P_ra)
K_p_ao = 5500;
K_f_ao = 50;
K_p_mi = 6000;
K_f_mi = 40;
K_p_po = 5000;
K_f_po = 60;
K_p_ti = 5800;
K_f_ti = 45;

f_ao = @(t,y) [y(2);K_p_ao * (P_lv-P_sas) * cos(y(1)) - K_f_ao*y(2)];
f_mi = @(t,y) [y(2);K_p_mi * (P_la-P_lv) * cos(y(1)) - K_f_mi*y(2)];
f_po = @(t,y) [y(2);K_p_po * (P_rv-P_pas) * cos(y(1)) - K_f_po*y(2)];
f_ti = @(t,y) [y(2);K_p_ti * (P_ra-P_rv) * cos(y(1)) - K_f_ti*y(2)];

t0 = 0;

if (t == 0) 
    y = 0;
    return
end 

tf = t;

y0 = [0;0];

[~,ys_ao] = ode45(f_ao,[t0 tf],y0);
[~,ys_mi] = ode45(f_mi,[t0 tf],y0);
[~,ys_po] = ode45(f_po,[t0 tf],y0);
[~,ys_ti] = ode45(f_ti,[t0 tf],y0);

y = [ys_ao(end,1); ys_mi(end,1), ys_po(end,1),ys_ti(end,1)];
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



% theta_max = 30;

% syms ebar_v(t) ebar_a(t) e_lv(t) e_ra(t) e_rv(t) e_la(t)
% syms V_lv(t) V_la(t) V_rv(t) V_ra(t)
% syms P_lv(t) P_rv(t) P_la(t) P_ra(t) P_sas(t) P_sat(t) P_svn(t)
% syms P_pas(t) P_pat(t) P_pvn(t)
% syms Q_mi(t) Q_sas(t) Q_ao(t) Q_svn(t) Q_sat(t) Q_pvn(t) Q_ti(t) Q_po(t) Q_pas(t) Q_pat(t)
% syms AR_ao(t) AR_po(t) AR_mi(t) AR_ti(t) theta_ao(t) theta_po(t) theta_mi(t) theta_ti(t)
% 
% %Activaiton funciton for time depedent elastance
% ebarV = piecewise(mod(t,T) < T_s1, ebar_v == 1 - cos(t/T_s1 * pi), (mod(t,T) >= T_s1) & ...
%     (mod(t,T) < T_s2), ebar_v == 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi), ebar_v == 0);
% 
% ebarA = piecewise(mod(t,T) < T_pwb, ebar_a == 0, (mod(t,T) >= T_pwb) & ...
%     (mod(t,T) < T_pwb + T_pww), ebar_v == 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi), ebar_v == 0);
% 
% %Elastance for each chamber of the heart
% elastanceLV = e_lv == E_lv_d + (E_lv_s - E_lv_d)/2 * ebar_v(t);
% elastanceLA = e_la == E_la_min + (E_la_max - E_la_min)/2 * ebar_a(t);
% elastanceRV = e_rv == E_rv_d + (E_rv_s - E_rv_d)/2 * ebar_v(t);
% elastanceRA = e_ra == E_ra_min + (E_ra_max - E_ra_min)/2 * ebar_a(t);
% 
% %Blood volume in each chamber of the heart
% volumeLV = diff(V_lv) == Q_mi - Q_ao;
% volumeLA = diff(V_la) == Q_pvn - Q_mi;
% volumeRV = diff(V_rv) == Q_ti - Q_po;
% volumeRA = diff(V_ra) == Q_svn - Q_ti;
% 
% %Pressures in each chamber of the heart
% pressureLV = P_lv == P_lv_0 + e_lv*(V_lv - V_lv_0);
% pressureRV = P_rv == P_rv_0 + e_rv*(V_rv - V_rv_0);
% pressureLA = P_la == P_la_0 + e_la*(V_la - V_la_0);
% pressureRA = P_ra == P_ra_0 + e_ra*(V_ra - V_ra_0);
% 
% %Pressures within systemic circulation 
% pressureSAS = diff(P_sas) == (Q_ao - Q_sas)/C_sas;
% pressureSAT = diff(P_sat) == (Q_sas - Q_sat)/C_sat;
% pressureSVN = diff(P_svn) == (Q_sat - Q_svn)/C_svn;
% 
% %Pressures within pulmonary circultion
% pressurePAS = diff(P_pas) == (Q_po - Q_pas)/C_pas;
% pressurePAT = diff(P_pat) == (Q_pas - Q_pat)/C_pat;
% pressurePVN = diff(P_pvn) == (Q_pat - Q_pvn)/C_pvn;
% 
% %Flow rates 
% flowSAS = Q_sas == (P_sas - P_svn - (R_sas + R_sat + R_sar + R_scp)*Q_sat)/L_sas;
% flowSAT = diff(Q_sat) == (P_sas - P_svn - (R_sas + R_sat + R_sar + R_scp)*Q_sat)/L_sat;
% flowSVN = Q_svn == (P_svn - P_ra)/R_svn;
% 
% flowPAS = diff(Q_pas) == (P_pas - P_pvn - (R_pas + R_pat + R_par + R_pcp + R_pvn) * Q_pas)/L_pas;
% flowPAT = diff(Q_pat) == (P_pas - P_pvn -(R_pas + R_pat + R_par + R_pcp + R_pvn) * Q_pas)/L_pat;
% flowPVN = Q_pvn == (P_pvn - P_la)/R_pvn;
% 
% 
% %Flow rates through heart valves 
% flowAO = piecewise(P_lv >= P_sas, Q_ao == CQ_ao * AR_ao * sqrt(P_lv - P_sas), ...
%     P_lv < P_sas, Q_ao == -CQ_ao * AR_ao * sqrt(P_sas - P_lv));
% valveAO = AR_ao == (1-cos(theta_ao))^2/(1-cos(theta_max))^2;
% [angleAO] = odeToVectorField(diff(theta_ao,2) == K_p_ao * (P_lv - P_sas)*cos(theta_ao)-K_f_ao*diff(theta_ao));
% 
% flowMI = piecewise(P_la >= P_lv, CQ_mi == CQ_mi * AR_mi * sqrt(P_la - P_lv), ...
%     P_la < P_lv, Q_mi == -CQ_mi * AR_mi * sqrt(P_lv - P_la));
% valveMI = AR_mi == (1-cos(theta_mi))^2/(1-cos(theta_max))^2;
% [angleMI] = odeToVectorField(diff(theta_mi,2) == K_p_mi * (P_la - P_lv)*cos(theta_mi)-K_f_mi*diff(theta_mi));
% 
% flowPO = piecewise(P_rv >= P_pas, CQ_po == CQ_po * AR_po * sqrt(P_rv - P_pas), ...
%     P_rv < P_pas, Q_po == -CQ_po * AR_po * sqrt(P_pas - P_rv));
% valvePO = AR_po == (1-cos(theta_po))^2/(1-cos(theta_max))^2;
% [anglePO] = odeToVectorField(diff(theta_po,2) == K_p_po * (P_rv - P_pas)*cos(theta_po)-K_f_po*diff(theta_po));
% 
% flowTI = piecewise(P_ra >= P_rv, CQ_ti == CQ_ti * AR_ti * sqrt(P_ra - P_rv), ...
%     P_ra < P_rv, Q_ti == -CQ_ti * AR_ti * sqrt(P_rv - P_ra));
% valveTI = AR_ti == (1-cos(theta_ti))^2/(1-cos(theta_max))^2;
% [angleTI] = odeToVectorField(diff(theta_ti,2) == K_p_ti * (P_ra - P_rv)*cos(theta_ti)-K_f_ti*diff(theta_ti));
% 
% allEq = [ebarV;ebarA;elastanceLV;elastanceLA;elastanceRV;elastanceRA;volumeLV;volumeLA
%         volumeRV;volumeRA;pressureLV;pressureRV;pressureLA;pressureRA;pressureSAS;
%         pressureSAT;pressureSVN;pressurePAS;pressurePAT;pressurePVN;flowSAS;flowSAT;
%         flowSVN;flowPAS;flowPAT;flowPVN;flowAO;valveAO;angleAO;flowMI;valveMI;angleMI;
%         flowPO;valvePO;anglePO;flowTI;valveTI;angleTI];
% 
% odes1 = [volumeLV;volumeLA;volumeRV;volumeRA;pressureSAS;pressureSAT;pressureSVN;pressurePAS;pressurePAT;pressurePVN;flowSAS;flowSAT;
%         flowPAS;flowPAT];
% 
% odes2 = [pressureSAS;pressurePAS;pressurePAT;flowSAS;flowSAT;flowPAT];
% 
% odes3 = [pressureSAT;pressurePAS;pressurePAT;flowSAT;flowPAS;flowPAT];
% 
% %M = matlabFunction(odes,'vars',{'t','Y'})
% cond1 = V_lv(0) == V_lv0;
% cond2 = V_rv(0) == V_rv0;
% conds = [cond1; cond2];
% e = dsolve(odes3,conds);
% 
% sol1 = e.Q_sat
