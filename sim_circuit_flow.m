 simulinkSimulator(["120" "10" "0" "0.5" "0.5" "0.25" "1.5"], ["25" "4" "0" "0.5" "0.5" "0.25" "1.5"],...
     '0.08', '0.000062','0.003', '0.0017', '1.6', '0.05', '0.5', '0.52', '1', '1', '1', '1', '1',...
     '0.075', '0.18', '0.002', '0.000052', '3.8', '0.01', '0.0017', '0.05', '0.25', '20.5', '0.006');

function simulinkSimulator(heart_l, heart_r, C_sas, L_sas, R_sas, L_sat,...
                            C_sat, R_sat, R_sar, R_scp, R_brain, R_liver,...
                            R_spleen, R_kidney, R_bone_marrow, R_svn, C_pas,...
                            R_pas, L_pas, C_pat, R_pat, L_pat, R_par, R_pcp,...
                            C_pvn, R_pvn)
% ----------------------------------------------------------
% Function Description:
% This function runs the simulation and plots the results. Takes in initial
% conditions as input defined below.
% ----------------------------------------------------------

% ----------------------------------------------------------
% Variables:
% heart_l = vector describing the pulsatile voltage source model of left
%           heart in the format ['V1', 'V2', 'TD', 'TR', 'TF', 'pW', 'PER']
% heart_r = vector describing the pulsatile voltage source model of right
%           heart in the format ['V1', 'V2', 'TD', 'TR', 'TF', 'pW', 'PER']
% C_sas = scalar capacitance value of the systemic aortic sinus modeling
%         the compliance of the vessel
% L_sas = scalar inductance value of the systemic aortic sinus modeling the
%         inertia of the vessel
% R_sas = scalar resistance value of the systemic aortic sinus modeling the
%         resistance of blood flow in the vessel
% L_sat = scalar inductance value of the systemic artery modeling the
%         inertia of the vessel
% C_sat = scalar capacitance value of the systemic artery modeling
%         the compliance of the vessel
% R_sat = scalar resistance value of the systemic artery modeling the
%         resistance of blood flow in the vessel
% R_sar = scalar resistance value of the systemic arterioles modeling the
%         resistance of blood flow in the vessel
% ...CONTINUE VARIABLE DEFINITIONS...
% ----------------------------------------------------------

% ----------------------------------------------------------
% PARAMETERS FOR MODEL
% ----------------------------------------------------------

% capacitor_param = c (in Farad)
% inductor_param = l (in Henry)
% resistor_param = R (in Ohms)

% heart_params = ['V1', 'V2', 'TD', 'TR', 'TF', 'pW', 'PER'];
% V1 = initial voltage 1 (in Volts)
% V2 = initial voltage 2 (in Volts)
% TD = pulse delay time (in seconds)
% TR = pulse rise time (in seconds)
% TF = pulse fall time (in seconds)
% pW = pulse width (in seconds)
% PER = pulse period (in seconds)

% Sets parameters of Simulink model
set_param('CirculationCircuitv2021_v3/R_scp', 'R', R_scp);
set_param('CirculationCircuitv2021_v3/R_sas', 'R', R_sas);
set_param('CirculationCircuitv2021_v3/R_sat', 'R', R_sat);
set_param('CirculationCircuitv2021_v3/R_sar', 'R', R_sar);
set_param('CirculationCircuitv2021_v3/R_brain', 'R', R_brain);
set_param('CirculationCircuitv2021_v3/R_kidney', 'R', R_kidney);
set_param('CirculationCircuitv2021_v3/R_spleen', 'R', R_spleen);
set_param('CirculationCircuitv2021_v3/R_liver', 'R', R_liver);
% set_param('CirculationCircuitv2021_v3/R_bone_marrow', 'R', R_bone_marrow);
set_param('CirculationCircuitv2021_v3/R_svn', 'R', R_svn);
set_param('CirculationCircuitv2021_v3/R_pas', 'R', R_pas);
set_param('CirculationCircuitv2021_v3/R_pat', 'R', R_pat);
set_param('CirculationCircuitv2021_v3/R_par', 'R', R_par);
set_param('CirculationCircuitv2021_v3/R_pcp', 'R', R_pcp);
set_param('CirculationCircuitv2021_v3/R_pvn', 'R', R_pvn);

set_param('CirculationCircuitv2021_v3/C_sas', 'c', C_sas);
set_param('CirculationCircuitv2021_v3/C_sat', 'c', C_sat);
set_param('CirculationCircuitv2021_v3/C_sas', 'c', C_sas);
set_param('CirculationCircuitv2021_v3/C_pas', 'c', C_pas);
set_param('CirculationCircuitv2021_v3/C_pat', 'c', C_pat);
set_param('CirculationCircuitv2021_v3/C_pvn', 'c', C_pvn);

set_param('CirculationCircuitv2021_v3/L_sas', 'l', L_sas);
set_param('CirculationCircuitv2021_v3/L_pas', 'l', L_pas);
set_param('CirculationCircuitv2021_v3/L_sat', 'l', L_sat);
set_param('CirculationCircuitv2021_v3/L_pat', 'l', L_pat);

% v1l = get_param('CirculationCircuitv2021/Heart_L', 'V1');
% v2l = get_param('CirculationCircuitv2021/Heart_L', 'V2');
% tdl = get_param('CirculationCircuitv2021/Heart_L', 'TD');
% tfl = get_param('CirculationCircuitv2021/Heart_L', 'TF');
% pwl = get_param('CirculationCircuitv2021/Heart_L', 'pW');
% perl = get_param('CirculationCircuitv2021/Heart_L', 'PER');
% 
% v1r = get_param('CirculationCircuitv2021/Heart_R', 'V1');
% v2r = get_param('CirculationCircuitv2021/Heart_R', 'V2');
% tdr = get_param('CirculationCircuitv2021/Heart_R', 'TD');
% tfr = get_param('CirculationCircuitv2021/Heart_R', 'TF');
% pwr = get_param('CirculationCircuitv2021/Heart_R', 'pW');
% perr = get_param('CirculationCircuitv2021/Heart_R', 'PER');

% set_param('CirculationCircuitv2021_v3/Heart_L', 'V1', v1l);
% set_param('CirculationCircuitv2021_v3/Heart_L', 'V2', v2l);
% set_param('CirculationCircuitv2021_v3/Heart_L', 'TD', tdl);
% set_param('CirculationCircuitv2021_v3/Heart_L', 'TF', tfl);
% set_param('CirculationCircuitv2021_v3/Heart_L', 'pW', pwl);
% set_param('CirculationCircuitv2021_v3/Heart_L', 'PER', perl);
% 
%  
% set_param('CirculationCircuitv2021_v3/Heart_R', 'V1', v1r);
% set_param('CirculationCircuitv2021_v3/Heart_R', 'V2', v2r);
% set_param('CirculationCircuitv2021_v3/Heart_R', 'TD', tdr);
% set_param('CirculationCircuitv2021_v3/Heart_R', 'TF', tfr);
% set_param('CirculationCircuitv2021_v3/Heart_R', 'pW', pwr);
% set_param('CirculationCircuitv2021_v3/Heart_R', 'PER', perr);

simOut=sim('CirculationCircuitv2021_v3'); % Loads and runs model

% Reads in results of Simulink model
Q_right_heart = simOut.Q_right_heart.signals.values;
time = simOut.Q_right_heart.time;
Q_left_heart = simOut.Q_left_heart.signals.values;
V_left_heart = simOut.V_left_heart.signals.values;
V_right_heart = simOut.V_right_heart.signals.values;
P_right_heart = simOut.P_right_heart.signals.values;
P_left_heart = simOut.P_left_heart.signals.values;

Q_sas = simOut.Q_sas.signals.values;
Q_sat = simOut.Q_sat.signals.values;
Q_sar = simOut.Q_sar.signals.values;
Q_scp = simOut.Q_scp.signals.values;

V_sas = simOut.V_sas.signals.values;
V_sat = simOut.V_sat.signals.values;
V_sar = simOut.V_sar.signals.values;
V_scp = simOut.V_scp.signals.values;

% ----------------------------------------------------------
% Heart Dynamics
% ----------------------------------------------------------
figure
subplot(2,3,1)
plot(time,Q_right_heart);
title('Q_{Right Heart}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(2,3,4)
plot(time,Q_left_heart);
title('Q_{Left Heart}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(2,3,2)
plot(time,V_right_heart);
title('V_{Right Heart}');
xlabel('Time (s)');
ylabel('Flow (mL)');

subplot(2,3,5)
plot(time,V_left_heart);
title('V_{Left Heart}');
xlabel('Time (s)');
ylabel('Flow (mL)');

subplot(2,3,6)
plot(time,P_left_heart);
title('P_{Left Heart}');
xlabel('Time (s)');
ylabel('Blood Pressure (mmHg)');

subplot(2,3,3)
plot(time,P_right_heart);
title('P_{Right Heart}');
xlabel('Time (s)');
ylabel('Blood Pressure (mmHg)');

sgtitle('Heart Blood Flow Dynamics');


% ----------------------------------------------------------
% Systemic Loop Dynamics
% ----------------------------------------------------------
figure
subplot(4,2,1)
plot(time, Q_sas);
title('Q_{sas}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(4,2,3)
plot(time, Q_sat);
title('Q_{sat}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(4,2,5)
plot(time, Q_sar);
title('Q_{sar}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(4,2,7)
plot(time, Q_scp);
title('Q_{scp}');
xlabel('Time (s)');
ylabel('Flow Rate (mL/s)');

subplot(4,2,2)
plot(time, V_sas);
title('V_{sas}');
xlabel('Time (s)');
ylabel('Flow (mL)');

subplot(4,2,4)
plot(time, V_sat);
title('V_{sat}');
xlabel('Time (s)');
ylabel('Flow (mL)');

subplot(4,2,6)
plot(time, V_sar);
title('V_{sar}');
xlabel('Time (s)');
ylabel('Flow (mL)');

subplot(4,2,8)
plot(time, V_scp);
title('V_{scp}');
xlabel('Time (s)');
ylabel('Flow (mL)');


sgtitle('Systemic Circuit Blood Flow Rates');

end