function custom_ode_tester

Plv_0 = 1.0;
Pla_0 = 1.0;
k_stla = 2.5;
k_stlv = 20.0;
Msav = 0.0004;
Asav = 0.00047;
k_fsav = 0.0004;
k_esav = 9000.0;

T = 0.9999;                            % final t-value


f = @(t,y) [y(2); (k_stla*ela(t))/Msav - (k_stlv*elv(t))/Msav + ...
                    (Plv_0-Pla_0)*Asav/Msav - (k_fsav/Msav)*y(2) - ...
                    (k_esav/Msav)*y(1)]; % define function f(t,y)
t0 = 0;  y0 = [0;0];             % initial condition with vector y0


[ts,ys] = ode45(f,[t0,T],y0);     % solve IVP

fprintf('y(T) = %g, y''(T) = %g\n',ys(end,1),ys(end,2))

disp('       t           y1          y2')
disp([ts,ys])                     % table of t and y values of solution

plot(ts,ys(:,1),'b')              % plot solution y(t)
title('Solution y(t) of IVP')
xlabel('t'); grid on

end

function elastance = ela(t)

T_pwb = 0.92;
T_pww = 0.09;
E_la_max = 0.25;
E_la_min = 0.15;

if t < T_pwb
    elastance = E_la_min;
elseif t < (T_pwb + T_pww) && t >= T_pwb
    elastance = E_la_min + (E_la_max - E_la_min)/2 * (1-cos((t-T_pwb)/T_pwb)*2*pi);
else
    elastance = E_la_min;
end

end

function elastance = elv(t)

T_s1 = 0.3;
T_s2 = 0.45;
E_lv_s = 2.5;
E_lv_d = 0.1;

if t < T_s1
    elastance = E_lv_d + (E_lv_s - E_lv_d)/2*(1-cos(t/T_s1 * pi));
elseif t >= T_s1 && t < T_s2
    elastance = E_lv_d + (E_lv_s - E_lv_d)/2*(1+cos((t-T_s1)/(T_s2-T_s1) * pi));
else
    elastance = E_lv_d;
end

end