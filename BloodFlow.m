function BloodFlow
% This function models the blood flow and blood pressures of vessels and
% organs

% Defining unknown variables
syms Qao(t) Qsas(t) Qsat(t) Qsvn(t) Qmi(t) Qti(t)...
    Qpat(t) Qpvn(t) Qpas(t) Qlv(t) Qrv(t) Qra(t) Qla(t) Psas(t)...
    Psat(t) Ppat(t) Ppvn(t) Ppas(t) Pra(t) Pla(t) Plv(t) Qpo(t)...
    Prv(t) ARmi(t) ARti(t) ARao(t) ARpo(t) Vra(t) Vrv(t) Vla(t) Vlv(t)...
    lsav(t) lpav(t) theta_ao(t) theta_po(t) theta_ti(t) theta_mi(t) Psvn(t)

% Defining known variables
syms Elv_s Elv_d Cpas Cpat Cpvn Erv_s Erv_d Csas Csat Csvn Rpas Lpas Rpat...
    Rpar Rpcp Rpvn Lsas Rsas Rsat Lsat Rsar Rscp Rsvn CQao CQmi Plv_0 Vlv_0...
    Ela_max Ela_min Pla_0 Vla_0 CQpo CQti Prv_0 Vrv_0 Era_max...
    Era_min Pra_0 Vra_0 Kst_la Kst_lv Kf_sav Ke_sav Msav Asav Kst_ra Kst_rv...
    Kf_pav Ke_pav Mpav Apav Kp_ao Kf_ao Kp_mi Kf_mi Kp_sv Kf_sv Kp_po Kp_ti...
    Kf_ti Kp_pv Kf_pv theta_ao_max theta_po_max theta_ti_max theta_mi_max...
    Rsm Rbrain Rkid Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat Kf_po ela elv era erv

% Defining System of Differential Algebraic Equations

% Flow and Pressure Equations
eqn1 = diff(Psas(t)) == (Qao(t) - Qsas(t))/Csas;
eqn2 = diff(Psat(t)) == (Qsas(t) - Qsat(t))/Csat;
eqn3 = diff(Psvn(t)) == (Qsat(t) - Qsvn(t))/Csvn;
eqn4 = diff(Ppvn(t)) == (Qpat(t) - Qpvn(t))/Cpvn;
eqn5 = diff(Ppat(t)) == (Qpas(t) - Qpat(t))/Cpat;
eqn6 = diff(Ppas(t)) == (Qpo(t) - Qpas(t))/Cpas;
eqn7 = diff(Qsas(t)) == (Psas(t) - Psat(t) - Rsas*Qsas(t))/Lsas;
eqn8 = diff(Qsat(t)) == (Psat(t) - Psvn(t) - (Rsat + Rsar + Rscp)...
                        *Qsat(t))/Lsat;
eqn9 = diff(Qpat(t)) == (Ppat(t) - Ppvn(t) - (Rpat + Rpcp + Rpar)...
                        *Qpat(t))/Lpat;
eqn10 = diff(Qpas(t)) == (Ppas(t) - Ppat(t) - Rpas*Qpas(t))/Lpas;

eqn11 = Qao(t) == (CQao*ARao(t)*sqrt(abs(Plv(t) - Psas(t))));
eqn12 = ARao(t) == ((1-cos(theta_ao(t)))^2) / ((1-cos(theta_ao_max))^2);

eqn13 = Qpo(t) == (CQpo*ARpo(t)*sqrt(abs(Prv(t) - Ppas(t))));
eqn14 = ARpo(t) == ((1-cos(theta_po(t)))^2) / ((1-cos(theta_po_max))^2);

% eqn15 = Qo(t) == Qsat(t);
eqn16 = Qsvn(t) == (Psvn(t) - Pra(t))/Rsvn;
eqn17 = Qmi(t) + Qlv(t) == Qao(t);

eqn18 = Qti(t) == (CQti*ARti(t)*sqrt(abs(Pra(t) - Prv(t))));
eqn19 = ARti(t) == ((1-cos(theta_ti(t)))^2) / ((1-cos(theta_ti_max))^2);

eqn20 = Qmi(t) == (CQmi*ARmi(t)*sqrt(abs(Pla(t) - Plv(t))));
eqn21 = ARmi(t) == ((1-cos(theta_mi(t)))^2) / ((1-cos(theta_mi_max))^2);
    
% eqn22 = Poi(t) - Pof(t) == Qsat(t)/ ((1/Rsm) + (1/Rbrain) + (1/Rkid) + ...
%                                     (1/Rliv) + (1/Rsp) + (1/Rint));
                                
% eqn23 = Qof(t) == (Poi(t) - Pof(t))/Rsm;
% eqn24 = Qof(t) + Qoi(t) == Qo(t);
% eqn25 = Qsat(t) == Qoi(t) + Qof(t);
eqn26 = Qpvn(t) + Qla(t) == Qmi(t);
eqn27 = Qti(t) + Qrv(t) == Qpo(t);
eqn28 = Qsvn(t) + Qra(t) == Qti(t);
    
% Valve Pressure Equations
eqn33 = Pra(t) == Pra_0 + era*(Vra(t) - Vra_0);
eqn34 = Prv(t) == Prv_0 + erv*(Vrv(t) - Vrv_0);
eqn35 = Pla(t) == Pla_0 + ela*(Vla(t) - Vla_0);
eqn36 = Plv(t) == Plv_0 + elv*(Vlv(t) - Vlv_0);

% Valve Volume Equations
eqn37 = Vlv(t) == Vlv_0 + Asav*lsav(t);
eqn38 = Vla(t) == Vla_0 - Asav*lsav(t);
eqn39 = Vrv(t) == Vrv_0 + Apav*lpav(t);
eqn40 = Vra(t) == Vra_0 - Apav*lpav(t);

% Displacement Equations
eqn41 = Msav * diff(lsav(t), 2) == Kst_la*ela - Kst_lv*elv + (Plv(t)...
        - Pla(t))*Asav - Kf_sav*diff(lsav(t)) - Ke_sav*lsav(t);
eqn42 = Mpav * diff(lpav(t), 2) == Kst_ra*era - Kst_rv*erv + (Prv(t)...
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

vars = [Qao(t) Qsas(t) Qsat(t) Qsvn(t) Qmi(t) Qti(t)...
    Qpat(t) Qpvn(t) Qpas(t) Qlv(t) Qrv(t) Qpo(t)...
    Qra(t) Qla(t) Psas(t) Psat(t) Ppat(t) Ppvn(t) Ppas(t)...
    Pra(t) Pla(t) Plv(t) Prv(t) ARmi(t) ARti(t) ARao(t)...
    ARpo(t) Vra(t) Vrv(t) Vla(t) Vlv(t) lsav(t) lpav(t) theta_ao(t)...
    theta_po(t) theta_ti(t) theta_mi(t) Psvn(t)];

origVars = length(vars)

eqns = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9 eqn10 eqn11 eqn12...
        eqn13 eqn14 eqn16 eqn17 eqn18 eqn19 eqn20 eqn21 eqn26 eqn27 eqn28...
        eqn33 eqn34 eqn35 eqn36 eqn37 eqn38 eqn39 eqn40 eqn41 eqn42 eqn43...
        eqn44 eqn45 eqn46 eqn47];

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
                            Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat Kf_po...
                            ela era elv erv]);

for i=1:length(eqns)
    disp(eqns(i));
end

[eqs, vars, newVars, index] = reduceDAEIndex(eqns, vars)

disp(length(eqs));
disp(length(vars));
[eqs, vars, S] = reduceRedundancies(eqs, vars)

disp(length(eqs));
disp(length(vars));
S.solvedEquations
S.constantVariables
S.replacedVariables
S.otherEquations

isLowIndexDAE(eqs, vars)

[ODEs,constraints] = reduceDAEToODE(eqs,vars)
[massM,f] = massMatrixForm(ODEs,vars)

F = daeFunction(eqs, vars, [Elv_s Elv_d Cpas Cpat Cpvn Erv_s Erv_d Csas...
                            Csat Csvn Rpas Lpas Rpat Rpar Rpcp Rpvn Lsas...
                            Rsas Rsat Lsat Rsar Rscp Rsvn CQao CQmi Plv_0...
                            Vlv_0 Ela_max Ela_min Pla_0 Vla_0 CQpo CQti...
                            Prv_0 Vrv_0 Era_max Era_min Pra_0...
                            Vra_0 Kst_la Kst_lv Kf_sav Ke_sav Msav Asav...
                            Kst_ra Kst_rv Kf_pav Ke_pav Mpav Apav Kp_ao...
                            Kf_ao Kp_mi Kf_mi Kp_sv Kf_sv Kp_po Kf_po Kp_ti...
                            Kf_ti Kp_pv Kf_pv theta_ao_max theta_po_max...
                            theta_ti_max theta_mi_max Rsm Rbrain Rkid...
                            Rliv Rsp Rint T Tpwb Tpww Ts1 Ts2 Lpat ...
                            ela era elv erv])
                        
disp(vars);
                        
                      
f = @(t, y, yp)  F(t, y, yp, [2.5 0.1 0.18 3.8 20.5 1.15 0.1 0.08...
                           1.6 20.5 0.002 0.000052 0.01 0.05 0.25 0.006 0.000062...
                             0.003 0.05 0.0017 0.5 0.52 0.075 350 400 1.0...
                             5.0 0.25 0.15 1.0 4.0 350 400 ...
                             1.0 10 0.25 0.15 1.0...
                             4.0 2.5 20.0 0.0004 9000.0 0.0004 0.00047...
                            2.5 20.0 0.0004 9000.0 0.0004 0.00047 5500 ...
                            50 5500 50 5500 50 5500 50 5500 ...
                           50 5500 50 0.42*pi 0.42*pi...
                           0.42*pi 0.42*pi 0.1 0.1 0.1... %Note 0.1 is placeholder for organ resistances
                            0.1 0.1 0.1 1 0.91 0.09 0.30 0.45 0.0017...
                           ela1(t) era1(t) elv1(t) erv1(t)]);

opt = odeset('RelTol', 10.0^(-4), 'AbsTol' , 10.0^(-4));
t0 = 0;
tfinal = 1;
                        

% Note: Need to figure out decic function for accurate results
y0_est = [zeros(32,1)];
y0_est(1,1) = 15;
y0_est(5,1) = 15;
y0_est(6,1) = 45;
y0_est(12,1) = 0;
y0_est(15,1) = 76;
y0_est(19,1) = 13.9;
y0_est(20,1) = 7;
y0_est(21,1) = 9;
y0_est(23,1) = 0;
y0_est(24,1) = 0;
y0_est(25:26,1) = 0.42*pi;

fixedVars = [zeros(32,1)];
fixedVars(1,1) = 1;

[y0, yp0] = decic(f, t0, [ones(32,1)]*0.3, [], [zeros(32,1)], [], opt)
% [y0,yp0] = decic(eqs,vars,constraintEqs,t0,y0_est,fixedVars,yp0_est,opt);


[t,y] = ode15i(f, [t0, tfinal], y0, yp0, opt);
Qao = y(:,1);
Qsas = y(:,2);
Qsat = y(:,3);
[row, col] = size(y);
for i=1:col
    subplot(8,4,i)
    plot(t,y(:,i));
    tit = sprintf("Solution to %s",vars(i));
    title(tit);
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

function elastance = ela1(t)

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

function elastance = elv1(t)

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

function elastance = era1(t)

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

function elastance = erv1(t)

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