function BloodFlow
% This function models the blood flow and blood pressures of vessels and
% organs

% Defining unknown variables
syms Qao(t) Qsas(t) Qsat(t) Qo(t) Qsvn(t) Qmi(t) Qti(t) Qoi(t) Qof(t)...
    Qpat(t) Qpvn(t) Qpas(t) Qpo(t) Qlv(t) Qrv(t) Qra(t) Qla(t) Psas(t)...
    Psat(t) Poi(t) Ppat(t) Ppvn(t) Ppas(t) Pof(t) Pra(t) Pla(t) Plv(t)...
    Prv(t) ARmi(t) ARti(t) ARao(t) ARpo(t) Vra(t) Vrv(t) Vla(t) Vlv(t)...
    lsav(t) lpav(t) theta_ao(t) theta_po(t) theta_ti(t) theta_mi(t) Psvn(t)

% Defining known variables
syms Elv_s Elv_d Cpas Cpat Cpvn Erv_s Erv_d Csas Csat Csvn Rpas Lpas Rpat...
    Rpar Rpcp Rpvn Lsas Rsas Rsat Lsat Rsar Rscp Rsvn CQao CQmi Plv_0 Vlv_0...
    Ela_max Ela_min Pla_0 Vla_0 CQpo CQti Prv_0 Vrv_0 Era_max...
    Era_min Pra_0 Vra_0 Kst_la Kst_lv Kf_sav Ke_sav Msav Asav Kst_ra Kst_rv...
    Kf_pav Ke_pav Mpav Apav Kp_ao Kf_ao Kp_mi Kf_mi Kp_sv Kf_sv Kp_po Kp_ti...
    Kf_ti Kp_pv Kf_pv theta_ao_max theta_po_max theta_ti_max theta_mi_max...
    Rsm Rbrain Rkid Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat Kf_po 

% Defining System of Differential Algebraic Equations

% Flow and Pressure Equations
eqn1 = diff(Psas(t)) == (Qao(t) - Qsas(t))/Csas;
eqn2 = diff(Psat(t)) == (Qsas(t) - Qsat(t))/Csat;
eqn3 = diff(Psvn(t)) == (Qo(t) - Qsvn(t))/Csvn;
eqn4 = diff(Ppvn(t)) == (Qpat(t) - Qpvn(t))/Cpvn;
eqn5 = diff(Ppat(t)) == (Qpas(t) - Qpat(t))/Cpat;
eqn6 = diff(Ppas(t)) == (Qpo(t) - Qpas(t))/Cpas;
eqn7 = diff(Qsas(t)) == (Psas(t) - Psat(t) - Rsas*Qsas(t))/Lsas;
eqn8 = diff(Qsat(t)) == (Psat(t) - Poi(t) - (Rsat + Rsar + Rscp)...
                        *Qsat(t))/Lsat;
eqn9 = diff(Qpat(t)) == (Ppat(t) - Ppvn(t) - (Rpat + Rpcp + Rpar)...
                        *Qpat(t))/Lpat;
eqn10 = diff(Qpas(t)) == (Ppas(t) - Ppat(t) - Rpas*Qpas(t))/Lpas;

eqn11 = Qao(t) == (CQao*ARao(t)*sqrt(Plv(t) - Psas(t))).*(Plv(t) >= Psas(t))...
        + (CQao*ARao(t)*sqrt(Psas(t) - Plv(t))).*(Plv(t) < Psas(t));
eqn12 = ARao(t) == ((1-cos(theta_ao(t)))^2) / ((1-cos(theta_ao_max))^2);

eqn13 = Qpo(t) == (CQpo*ARpo(t)*sqrt(Prv(t) - Ppas(t))).*(Prv(t) >= Ppas(t))...
        + (CQpo*ARpo(t)*sqrt(Ppas(t) - Prv(t))).*(Prv(t) < Ppas(t));
eqn14 = ARpo(t) == ((1-cos(theta_po(t)))^2) / ((1-cos(theta_po_max))^2);

eqn15 = Qo(t) == Qsat(t);
eqn16 = Qsvn(t) == (Psvn(t) - Pra(t))/Rsvn;
eqn17 = Qmi(t) + Qlv(t) == Qao(t);

eqn18 = Qti(t) == (CQti*ARti(t)*sqrt(Pra(t) - Prv(t))).*(Pra(t) >= Prv(t))...
        + (CQti*ARti(t)*sqrt(Prv(t) - Pra(t))).*(Pra(t) < Prv(t));
eqn19 = ARti(t) == ((1-cos(theta_ti(t)))^2) / ((1-cos(theta_ti_max))^2);

eqn20 = Qmi(t) == (CQmi*ARmi(t)*sqrt(Pla(t) - Plv(t))).*(Pla(t) >= Plv(t))...
        + (CQmi*ARmi(t)*sqrt(Plv(t) - Pla(t))).*(Pla(t) < Plv(t));
eqn21 = ARmi(t) == ((1-cos(theta_mi(t)))^2) / ((1-cos(theta_mi_max))^2);
    
eqn22 = Poi(t) - Pof(t) == Qsat(t)/ ((1/Rsm) + (1/Rbrain) + (1/Rkid) + ...
                                    (1/Rliv) + (1/Rsp) + (1/Rint));
                                
eqn23 = Qof(t) == (Poi(t) - Pof(t))/Rsm;
eqn24 = Qof(t) + Qoi(t) == Qo(t);
eqn25 = Qsat(t) == Qoi(t) + Qof(t);
eqn26 = Qpvn(t) + Qla(t) == Qmi(t);
eqn27 = Qti(t) + Qrv(t) == Qpo(t);
eqn28 = Qsvn(t) + Qra(t) == Qti(t);


% Elastance Equations
% eqn29 = era(t) == Era_min.*((mod(t,T) >= 0) & (mod(t,T) < Tpwb)) + ...
%     (Era_min + (Era_max - Era_min)/2 * (1 - cos(((mod(t,T) - Tpwb)/Tpwb)*2*pi))).*...
%     ((mod(t,T) >= Tpwb) & (mod(t,T) < (Tpwb+Tpww))) + Era_min.*((mod(t,T) >= (Tpwb+Tpww))...
%     & (mod(t,T) <= T));
% 
% eqn30 = ela(t) == Ela_min.*((mod(t,T) >= 0) & (mod(t,T) < Tpwb)) + ...
%     (Ela_min + (Ela_max - Ela_min)/2 * (1 - cos(((mod(t,T) - Tpwb)/Tpwb)*2*pi))).*...
%     ((mod(t,T) >= Tpwb) & (mod(t,T) < (Tpwb+Tpww))) + Ela_min.*((mod(t,T) >= (Tpwb+Tpww))...
%     & (mod(t,T) <= T));
% 
% eqn31 = erv(t) == (Erv_d + (Erv_s - Erv_d)/2 * (1 - cos((mod(t,T)/Ts1)*pi)))...
%         .*((mod(t,T) >= 0) & (mod(t,T) < Ts1)) + ...
%         (Erv_d + (Erv_s - Erv_d)/2 * (1 - cos(((mod(t,T) - Ts1)/(Ts2 - Ts1))*pi)))...
%         .*((mod(t,T) >= Ts1) & (mod(t,T) < Ts2)) + ...
%         Erv_d.*((mod(t,T) >= Ts2) & (mod(t,T) <= T));
%     
% eqn32 = elv(t) == (Elv_d + (Elv_s - Elv_d)/2 * (1 - cos((mod(t,T)/Ts1)*pi)))...
%         .*((mod(t,T) >= 0) & (mod(t,T) < Ts1)) + ...
%         (Elv_d + (Elv_s - Elv_d)/2 * (1 - cos(((mod(t,T) - Ts1)/(Ts2 - Ts1))*pi)))...
%         .*((mod(t,T) >= Ts1) & (mod(t,T) < Ts2)) + ...
%         Elv_d.*((mod(t,T)) >= Ts2 & (mod(t,T) <= T));
    
% Valve Pressure Equations
eqn33 = Pra(t) == Pra_0 + era(t)*(Vra(t) - Vra_0);
eqn34 = Prv(t) == Prv_0 + erv(t)*(Vrv(t) - Vrv_0);
eqn35 = Pla(t) == Pla_0 + ela(t)*(Vla(t) - Vla_0);
eqn36 = Plv(t) == Plv_0 + elv(t)*(Vlv(t) - Vlv_0);

% Valve Volume Equations
eqn37 = Vlv(t) == Vlv_0 + Asav*lsav(t);
eqn38 = Vla(t) == Vla_0 - Asav*lsav(t);
eqn39 = Vrv(t) == Vrv_0 + Apav*lpav(t);
eqn40 = Vra(t) == Vra_0 - Apav*lpav(t);

% Displacement Equations
eqn41 = Msav * diff(lsav(t), 2) == Kst_la*ela(t) - Kst_lv*elv(t) + (Plv(t)...
        - Pla(t))*Asav - Kf_sav*diff(lsav(t)) - Ke_sav*lsav(t);
eqn42 = Mpav * diff(lpav(t), 2) == Kst_ra*era(t) - Kst_rv*erv(t) + (Prv(t)...
        - Pra(t))*Apav - Kf_pav*diff(lpav(t)) - Ke_pav*lpav(t);
    
% Leaflet Angle Equations
eqn43 = diff(theta_ao(t), 2) == Kp_ao*(Plv(t) - Psas(t))*cos(theta_ao(t))...
        - Kf_ao*diff(theta_ao(t));
eqn44 = diff(theta_po(t), 2) == Kp_po*(Prv(t) - Ppas(t))*cos(theta_po(t))...
        - Kf_po*diff(theta_po(t));
eqn45 = diff(theta_ti(t), 2) == Kp_ti*(Pra(t) - Prv(t))*cos(theta_ti(t))...
        - Kf_ti*diff(theta_ti(t));
eqn46 = diff(theta_mi(t), 2) == Kp_mi*(Pla(t) - Plv(t))*cos(theta_mi(t))...
        - Kf_mi*diff(theta_mi(t));
    
% (I found I was missing this equation)
eqn47 = Qpvn(t) == (Ppvn(t) - Pla(t))/Rpvn;

vars = [Qao(t) Qsas(t) Qsat(t) Qo(t) Qsvn(t) Qmi(t) Qti(t)...
    Qoi(t) Qof(t) Qpat(t) Qpvn(t) Qpas(t) Qpo(t) Qlv(t) Qrv(t)...
    Qra(t) Qla(t) Psas(t) Psat(t) Poi(t) Ppat(t) Ppvn(t) Ppas(t)...
    Pof(t) Pra(t) Pla(t) Plv(t) Prv(t) ARmi(t) ARti(t) ARao(t)...
    ARpo(t) Vra(t) Vrv(t) Vla(t) Vlv(t) lsav(t) lpav(t) theta_ao(t)...
    theta_po(t) theta_ti(t) theta_mi(t) Psvn(t)];

origVars = length(vars);

eqns = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10 eqn11 eqn12...
        eqn13 eqn14 eqn15 eqn16 eqn17 eqn18 eqn19 eqn20 eqn21 eqn22 eqn23...
        eqn24 eqn25 eqn26 eqn27 eqn28 eqn33 eqn34...
        eqn35 eqn36 eqn37 eqn38 eqn39 eqn40 eqn41 eqn42 eqn43 eqn44 eqn45...
        eqn46 eqn47];

[eqns, vars, newVars] = reduceDifferentialOrder(eqns,vars);

isLowIndexDAE(eqns, vars)

F = daeFunction(eqns, vars, [Elv_s Elv_d Cpas Cpat Cpvn Erv_s Erv_d Csas...
                            Csat Csvn Rpas Lpas Rpat Rpar Rpcp Rpvn Lsas...
                            Rsas Rsat Lsat Rsar Rscp Rsvn CQao CQmi Plv_0...
                            Vlv_0 Ela_max Ela_min Pla_0 Vla_0 CQpo CQti...
                            Prv_0 Vrv_0 Era_max Era_min Pra_0...
                            Vra_0 Kst_la Kst_lv Kf_sav Ke_sav Msav Asav...
                            Kst_ra Kst_rv Kf_pav Ke_pav Mpav Apav Kp_ao...
                            Kf_ao Kp_mi Kf_mi Kp_sv Kf_sv Kp_po Kp_ti...
                            Kf_ti Kp_pv Kf_pv theta_ao_max theta_po_max...
                            theta_ti_max theta_mi_max Rsm Rbrain Rkid...
                            Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat Kf_po]);

for i=1:length(eqns)
    disp(eqns(i));
end

% [ODEs,constraints] = reduceDAEToODE(eqns,vars);
% 
% [massM,f] = massMatrixForm(ODEs,vars)
% 
% pODEs = symvar(ODEs);
% pvars = symvar(vars);
% extraParams = setdiff(pODEs, pvars)

% massM = odeFunction(massM, vars, m, r, g);
% f = odeFunction(f, vars, m, r, g);
% 
% m = 1;
% r = 1;
% g = 9.81;
% ODEsNumeric = subs(ODEs);
% constraintsNumeric = subs(constraints);
% 
% M = @(t,Y) massM(t,Y,m,r,g);
% F = @(t,Y) f(t,Y,m,r,g);
% 
% y0est = [r*sin(pi/6); -r*cos(pi/6); 0; 0; 0];
% yp0est = zeros(5,1);
% 
% opt = odeset('Mass', M, 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));
% 
% [y0, yp0] = decic(ODEsNumeric, vars, constraintsNumeric, 0,...
%                 y0est, [1,0,0,0,1], yp0est, opt)
%             
% opt = odeset(opt, 'InitialSlope', yp0);
% 
% [tSol,ySol] = ode15s(F, [0, 0.5], y0, opt);
% plot(tSol,ySol(:,1:origVars),'-o')
% 
% for k = 1:origVars
%   S{k} = char(vars(k));
% end
% 
% legend(S, 'Location', 'Best')
% grid on

end

function elastance = ela(t)

T = 1;
T_pwb = 0.91;
T_pww = 0.09;
E_la_max = 0.25;
E_la_min = 0.15;

t = mod(t,T);

if t < T_pwb
    elastance = E_la_min;
elseif t < (T_pwb + T_pww) && t >= T_pwb
    elastance = E_la_min + (E_la_max - E_la_min)/2 * (1-cos(((t-T_pwb)/T_pww)*2*pi));
else
    elastance = E_la_min;
end

end

function elastance = elv(t)

T = 1;
T_s1 = 0.3;
T_s2 = 0.45;
E_lv_s = 2.5;
E_lv_d = 0.1;

t = mod(t,T);

if t < T_s1
    elastance = E_lv_d + (E_lv_s - E_lv_d)/2*(1-cos(t/T_s1 * pi));
elseif t >= T_s1 && t < T_s2
    elastance = E_lv_d + (E_lv_s - E_lv_d)/2*(1+cos((t-T_s1)/(T_s2-T_s1) * pi));
else
    elastance = E_lv_d;
end

end

function elastance = era(t)

T = 1;
T_pwb = 0.91;
T_pww = 0.09;
E_ra_max = 0.25;
E_ra_min = 0.15;

t = mod(t,T);

if t < T_pwb
    elastance = E_ra_min;
elseif t < (T_pwb + T_pww) && t >= T_pwb
    elastance = E_ra_min + (E_ra_max - E_ra_min)/2 * (1-cos(((t-T_pwb)/T_pww)*2*pi));
else
    elastance = E_ra_min;
end

end

function elastance = erv(t)

T = 1;
T_s1 = 0.3;
T_s2 = 0.45;
E_rv_s = 1.15;
E_rv_d = 0.1;

t = mod(t,T);

if t < T_s1
    elastance = E_rv_d + (E_rv_s - E_rv_d)/2*(1-cos(t/T_s1 * pi));
elseif t >= T_s1 && t < T_s2
    elastance = E_rv_d + (E_rv_s - E_rv_d)/2*(1+cos((t-T_s1)/(T_s2-T_s1) * pi));
else
    elastance = E_rv_d;
end

end