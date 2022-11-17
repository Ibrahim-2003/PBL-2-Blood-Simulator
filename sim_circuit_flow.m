t = [0:0.01:10];
elax = [zeros(length(t),1)];
elvx = [zeros(length(t),1)];
erax = [zeros(length(t),1)];
ervx = [zeros(length(t),1)];
for n=1:length(t)
    elax(n) = ela1(t(n));
    erax(n) = era1(t(n));
    elvx(n) = elv1(t(n));
    ervx(n) = erv1(t(n));
end

disp(size(elax'))
disp(size(t))

ela = [t elax'];
disp(size(ela))

% ela=timeseries(elax,t); % Create time series data
era=timeseries(erax,t);
elv=timeseries(elvx,t);
erv=timeseries(ervx,t);

dataSet=Simulink.SimulationData.Dataset(ela); % Create Simulink data set
save('ela.mat','dataSet');% Save it to MAT-file
dataSet=Simulink.SimulationData.Dataset(era); % Create Simulink data set
save('era.mat','dataSet');% Save it to MAT-file
dataSet=Simulink.SimulationData.Dataset(elv); % Create Simulink data set
save('elv.mat','dataSet');% Save it to MAT-file
dataSet=Simulink.SimulationData.Dataset(erv); % Create Simulink data set
save('erv.mat','dataSet');% Save it to MAT-file
% 
% simOut=sim('flow_circuit');
% voltage = simOut.voltage1.Data;
% time = simOut.voltage1.Time;

figure
subplot(2,2,1);
plot(era);
title('Elastance of Right Atrium');
ylabel('Elastance (mmHg/mL)');
xlabel('Time (s)');
subplot(2,2,2);
plot(ela(:,1:1001),ela(:,1002:2002));
title('Elastance of Left Atrium');
ylabel('Elastance (mmHg/mL)');
xlabel('Time (s)');
subplot(2,2,3);
plot(erv);
title('Elastance of Right Ventricle');
ylabel('Elastance (mmHg/mL)');
xlabel('Time (s)');
subplot(2,2,4);
plot(elv);
title('Elastance of Left Ventricle');
ylabel('Elastance (mmHg/mL)');
xlabel('Time (s)');

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