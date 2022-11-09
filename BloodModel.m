function BloodModel
%Additional Parameters 

T = 1;
T_s1 = 0.3;
T_s2 = 0.45;
T_pwb = 0.92;
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

%May need to change 
theta_max = 30;

syms ebar_v(t) ebar_a(t) e_lv(t) e_ra(t) e_rv(t) e_la(t)
syms V_lv(t) V_la(t) V_rv(t) V_ra(t)
syms P_lv(t) P_rv(t) P_la(t) P_ra(t) P_sas(t) P_sat(t) P_svn(t)
syms P_pas(t) P_pat(t) P_pvn(t)
syms Q_mi(t) Q_sas(t) Q_ao(t) Q_svn(t) Q_sat(t) Q_pvn(t) Q_ti(t) Q_po(t) Q_pas(t) Q_pat(t)
syms AR_ao(t) AR_po(t) AR_mi(t) AR_ti(t) theta_ao(t) theta_po(t) theta_mi(t) theta_ti(t)

%Activaiton funciton for time depedent elastance
ebarV = piecewise(mod(t,T) < T_s1, ebar_v == 1 - cos(t/T_s1 * pi), (mod(t,T) >= T_s1) & ...
    (mod(t,T) < T_s2), ebar_v == 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi), ebar_v == 0);

ebarA = piecewise(mod(t,T) < T_pwb, ebar_a == 0, (mod(t,T) >= T_pwb) & ...
    (mod(t,T) < T_pwb + T_pww), ebar_v == 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi), ebar_v == 0);

%Elastance for each chamber of the heart
elastanceLV = e_lv == E_lv_d + (E_lv_s - E_lv_d)/2 * ebar_v(t);
elastanceLA = e_la == E_la_min + (E_la_max - E_la_min)/2 * ebar_a(t);
elastanceRV = e_rv == E_rv_d + (E_rv_s - E_rv_d)/2 * ebar_v(t);
elastanceRA = e_ra == E_ra_min + (E_ra_max - E_ra_min)/2 * ebar_a(t);

%Blood volume in each chamber of the heart
volumeLV = diff(V_lv) == Q_mi - Q_ao;
volumeLA = diff(V_la) == Q_pvn - Q_mi;
volumeRV = diff(V_rv) == Q_ti - Q_po;
volumeRA = diff(V_ra) == Q_svn - Q_ti;

%Pressures in each chamber of the heart
pressureLV = P_lv == P_lv_0 + e_lv*(V_lv - V_lv_0);
pressureRV = P_rv == P_rv_0 + e_rv*(V_rv - V_rv_0);
pressureLA = P_la == P_la_0 + e_la*(V_la - V_la_0);
pressureRA = P_ra == P_ra_0 + e_ra*(V_ra - V_ra_0);

%Pressures within systemic circulation 
pressureSAS = diff(P_sas) == (Q_ao - Q_sas)/C_sas;
pressureSAT = diff(P_sat) == (Q_sas - Q_sat)/C_sat;
pressureSVN = diff(P_svn) == (Q_sat - Q_svn)/C_svn;

%Pressures within pulmonary circultion
pressurePAS = diff(P_pas) == (Q_po - Q_pas)/C_pas;
pressurePAT = diff(P_pat) == (Q_pas - Q_pat)/C_pat;
pressurePVN = diff(P_pvn) == (Q_pat - Q_pvn)/C_pvn;

%Flow rates 
flowSAS = Q_sas == (P_sas - P_svn - (R_sas + R_sat + R_sar + R_scp)*Q_sat)/L_sas;
flowSAT = diff(Q_sat) == (P_sas - P_svn - (R_sas + R_sat + R_sar + R_scp)*Q_sat)/L_sat;
flowSVN = Q_svn == (P_svn - P_ra)/R_svn;

flowPAS = diff(Q_pas) == (P_pas - P_pvn - (R_pas + R_pat + R_par + R_pcp + R_pvn) * Q_pas)/L_pas;
flowPAT = diff(Q_pat) == (P_pas - P_pvn -(R_pas + R_pat + R_par + R_pcp + R_pvn) * Q_pas)/L_pat;
flowPVN = Q_pvn == (P_pvn - P_la)/R_pvn;


%Flow rates through heart valves 
flowAO = piecewise(P_lv >= P_sas, Q_ao == CQ_ao * AR_ao * sqrt(P_lv - P_sas), ...
    P_lv < P_sas, Q_ao == -CQ_ao * AR_ao * sqrt(P_sas - P_lv));
valveAO = AR_ao == (1-cos(theta_ao))^2/(1-cos(theta_max))^2;
angleAO = diff(theta_ao,2) == K_p_ao * (P_lv - P_sas)*cos(theta_ao)-K_f_ao*diff(theta_ao);

flowMI = piecewise(P_la >= P_lv, CQ_mi == CQ_mi * AR_mi * sqrt(P_la - P_lv), ...
    P_la < P_lv, Q_mi == -CQ_mi * AR_mi * sqrt(P_lv - P_la));
valveMI = AR_mi == (1-cos(theta_mi))^2/(1-cos(theta_max))^2;
angleMI = diff(theta_mi,2) == K_p_mi * (P_la - P_lv)*cos(theta_mi)-K_f_mi*diff(theta_mi);

flowPO = piecewise(P_rv >= P_pas, CQ_po == CQ_po * AR_po * sqrt(P_rv - P_pas), ...
    P_rv < P_pas, Q_po == -CQ_po * AR_po * sqrt(P_pas - P_rv));
valvePO = AR_po == (1-cos(theta_po))^2/(1-cos(theta_max))^2;
anglePO = diff(theta_po,2) == K_p_po * (P_rv - P_pas)*cos(theta_po)-K_f_po*diff(theta_po);

flowTI = piecewise(P_ra >= P_rv, CQ_ti == CQ_ti * AR_ti * sqrt(P_ra - P_rv), ...
    P_ra < P_rv, Q_ti == -CQ_ti * AR_ti * sqrt(P_rv - P_ra));
valveTI = AR_ti == (1-cos(theta_ti))^2/(1-cos(theta_max))^2;
angleTI = diff(theta_ti,2) == K_p_ti * (P_ra - P_rv)*cos(theta_ti)-K_f_ti*diff(theta_ti);

odes = [ebarV;ebarA;elastanceLV;elastanceLA;elastanceRV;elastanceRA;volumeLV;volumeLA
        volumeRV;volumeRA;pressureLV;pressureRV;pressureLA;pressureRA;pressureSAS;
        pressureSAT;pressureSVN;pressurePAS;pressurePAT;pressurePVN;flowSAS;flowSAT;
        flowSVN;flowPAS;flowPAT;flowPVN;flowAO;valveAO;angleAO;flowMI;valveMI;angleMI;
        flowPO;valvePO;anglePO;flowTI;valveTI;angleTI];
S = dsolve(odes)

end
%Ignore the Functions
function Q_po = CardiacOutput(t)
  

end

function elastance = ActivationFunction(t,ventricle)
T_s1 = 0.3;
T_s2 = 0.45;
T_pwb = 0.92;
T_pww = 0.09;
E_lv_s = 2.5;
E_lv_d = 0.1;
E_la_max = 0.25;
E_la_min = 0.15;


% make t a 
while t > 1
    t = t - 1;
end

if ventricle 
    if t < T_s1
        ebar = 1 - cos(t/T_s1 * pi);
    elseif t < T_s2
        ebar = 1 + cos((t-T_s1)/(T_s2 - T_s1) * pi);
    else 
        ebar = 0;
    end
    
    elastance = E_lv_d + (E_lv_s - E_lv_d)/2 * ebar;

else %if ventricle is false (if atrium)
    if t < T_pwb
        ebar = 0;
    elseif t < T_pwb + T_pww
        ebar = 1 - cos((t-T_pwb)/T_pww * 2 * pi);
    end

    elastance = E_la_min + (E_la_max - E_la_min)/2 * ebar;
end

end

function AorticSinus
end