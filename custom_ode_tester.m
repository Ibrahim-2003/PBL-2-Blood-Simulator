function custom_ode_tester

Plv_0 = 1.0;
Pla_0 = 1.0;
k_stla = 2.5;
k_stlv = 20.0;
Msav = 0.0004;
Asav = 0.00047;
k_fsav = 0.0004;
k_esav = 9000.0;

T = 0.9999; % final t-value


f = @(t,y) [y(2); (k_stla*ela(t))/Msav - (k_stlv*elv(t))/Msav + ...
                    (Plv_0-Pla_0)*Asav/Msav - (k_fsav/Msav)*y(2) - ...
                    (k_esav/Msav)*y(1)]; % define function f(t,y)
t0 = 0;  y0 = [0.0032;0];             % initial condition with vector y0


[ts,ys] = ode23t(f,[t0,T],y0);     % solve IVP

% fprintf('y(T) = %g, y''(T) = %g\n',ys(end,1),ys(end,2))
% 
% disp('       t           y1          y2')
% disp([ts,ys])                     % table of t and y values of solution

plot(ts,ys(:,1),'b')              % plot solution y(t)
title('Solution y(t) of IVP')
xlabel('t'); grid on

plot_left_elastances();
plot_right_elastances();

end

function plot_left_elastances

ts = [0:0.0001:2];

for i=1:length(ts)
    elas(i) = ela(ts(i));
end

figure
plot(ts,elas);
hold on
for i=1:length(ts)
    elav(i) = elv(ts(i));
end
plot(ts,elav);
ylim([0 3]);
title('Elastances of Left Side of Heart');
xlabel('Time (s)');
legend('Left Atrium', 'Left Ventricle');
grid off

end

function plot_right_elastances

ts = [0:0.0001:2];

for i=1:length(ts)
    elas(i) = era(ts(i));
end

figure
plot(ts,elas);
hold on
for i=1:length(ts)
    elav(i) = erv(ts(i));
end
plot(ts,elav);
ylim([0 3]);
title('Elastances of Right Side of Heart');
xlabel('Time (s)');
legend('Right Atrium', 'Right Ventricle');
grid off

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

function dydt = RHS(t,y)
% Defines system of ODEs for problem

T = 1;
T_s1 = 0.3;
T_s2 = 0.45;
E_lv_s = 2.5;
E_lv_d = 0.1;
T_pwb = 0.92;
T_pww = 0.09;
E_la_max = 0.25;
E_la_min = 0.15;
E_rv_s = 1.15;
E_rv_d = 0.1;
E_ra_max = 0.25;
E_ra_min = 0.15;

theta_ao_max = 0.42*pi;
theta_po_max = 0.42*pi;

% Implements the piecewise function, ela
tmod = mod(t,T);
ela1 = E_la_min;
ela2 = E_la_min + (E_la_max - E_la_min)/2 * (1-cos((tmod-T_pwb)/T_pwb)*2*pi);
ela3 = E_la_min;
ela = ela1.*(0<tmod & tmod<T_pwb) +ela2.*(T_pwb<=tmod & tmod<(T_pwb + T_pww))+ela3.*...
        ((T_pwb + T_pww)<=tmod & tmod<T);
    
% Implements the piecewise function, elv
elv1 = E_lv_d + (E_lv_s - E_lv_d)/2*(1-cos(tmod/T_s1 * pi));
elv2 = E_lv_d + (E_lv_s - E_lv_d)/2*(1+cos((tmod-T_s1)/(T_s2-T_s1) * pi));
elv3 = E_lv_d;
elv = elv1.*(0<tmod & tmod<T_s1) +elv2.*(T_s1<=tmod & tmod<T_s2)+elv3.*...
        (T_s2<=tmod & tmod<T);
   
% Implements the piecewise function, era
era1 = E_ra_min;
era2 = E_ra_min + (E_ra_max - E_ra_min)/2 * (1-cos((tmod-T_pwb)/T_pwb)*2*pi);
era3 = E_ra_min;
era = era1.*(0<tmod & tmod<T_pwb) +era2.*(T_pwb<=tmod & tmod<(T_pwb + T_pww))+era3.*...
        ((T_pwb + T_pww)<=tmod & tmod<T);
    
% Implements the piecewise function, erv
erv1 = E_rv_d + (E_rv_s - E_rv_d)/2*(1-cos(tmod/T_s1 * pi));
erv2 = E_rv_d + (E_rv_s - E_rv_d)/2*(1+cos((tmod-T_s1)/(T_s2-T_s1) * pi));
erv3 = E_rv_d;
erv = erv1.*(0<tmod & tmod<T_s1) +erv2.*(T_s1<=tmod & tmod<T_s2)+erv3.*...
        (T_s2<=tmod & tmod<T);
    


    

    

end