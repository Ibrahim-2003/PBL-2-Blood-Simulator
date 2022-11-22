% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel.m
% Models O2 concentration through the human body in different conditions
% Last modified: 11/22/22

function oxygenModel
    %the driver
    %setting up initial values, need to compare this with the flow model
    time = transpose(1:5);

    %init = table(time,flow,vol,pressure,oconca,oconc,oconcv,co2conc) %initial values vector
    %flow = init.flow
    %ho = heart(init);
    %bo = brain(ho);
    
    tspan = 0:0.0001:15;
    y0 = 975;
    [t,y] = ode45(@heartoxygencons,tspan,y0);
    plot(t,y)

    figure;
    load("Q_left_heart.mat");
    plot(t,Q_left_heart)
    %x = linspace(0,20,250);
    %y = deval(sol,x);

    %plot(x,y)
end

function ho = heart(values)
    %computes change in blood through heart
    %inputs: values - matrix of blood values coming in
    %outputs: ho - blood values matrix altered by heart
    
    ho = values;
    t = ho.time;

%     %Simple model
%     %oxygen consumption
%     ocons = ho.oconc - 10*t;
%     ho.oconc = ocons;
%     %co2 generation, assuming all o2 consumption goes towards making co2
%     co2gen = ho.co2conc + ocons;
%     ho.co2conc = co2gen;
%     ho
    
    %My attempt at the model from the paper
    %we might wanna specify these in the driver or just assume that they
    %are one value
    weight = 150;
    height = 180;
    age = 21;
    OP0 = 88.362 + 13.397*weight + 4.799*height - 5.677*age; %BMR

    %trying the ODE
    %syms O2(t) OP(t) Qdot(t) sigmaO2 t
    sigmaO2 = (ho.oconca(1,:) - (OP0/ho.flow(1,:)))/ho.oconc(1,:);
    OP(t) = -33.1*t;
    Qdot(t) = 150*t;
    ode = diff(O2,t) == (OP(t) + Qdot(t).*(ho.oconca-sigmaO2.*ho.oconc))./331;
    ode = matlabFunction(ode);
    O2sol(t) = dsolve(ode)

end

function dy = heartoxygencons(t,y)
    flow = 1;     
    vol = 5000; %mL
    pressure = 15; %mmHg
    oconca = 0.195*vol;
    oconcv = 0.145*vol;
    co2conc = 0.492*vol;

%     oconca = values.oconca;
%     flow = values.flow;
%     oconc = values.oconc;

    weight = 150;
    height = 180;
    age = 21;
    OP0 = 88.362 + 13.397*weight + 4.799*height - 5.677*age; %BMR
    
%     v1l = '120';
%     v2l = '10';
%     tdl = '0';
%     tfl = '0.5';
%     trl = '0.5';
%     pwl = '0.25';
%     perl = '1.5';
%     v1r = '25';
%     v2r = '4';
%     tdr = '0';
%     tfr = '0.5';
%     trr = '0.5';
%     pwr = '0.25';
%     perr = '1.5';
%     C_sas = '0.08';
%     L_sas = '0.000062';
%     R_sas = '0.003';
%     L_sat = '0.0017';
%     C_sat = '1.6';
%     R_sat = '0.05';
%     R_sar = '0.5';
%     R_scp = '0.52';
%     R_brain = '1';
%     R_liver = '1';
%     R_spleen = '1';
%     R_kidney = '1';
%     R_svn = '0.075';
%     C_pas = '0.18';
%     R_pas = '0.002';
%     L_pas = '0.000052';
%     C_pat = '3.8';
%     R_pat = '0.01';
%     L_pat = '0.0017';
%     R_par = '0.05';
%     R_pcp = '0.25';
%     C_pvn = '20.5';
%     R_pvn = '0.006';
% 
%     %sim_circuit_flow;
% 
%     [time, Q_right_heart, Q_left_heart, V_left_heart,...
%             V_right_heart, P_right_heart, P_left_heart, Q_sas, Q_sat,...
%             Q_sar, Q_scp, V_sas, V_sat, V_sar, V_scp, P_sas, P_sat,...
%             P_sar, P_scp, Q_svn, V_svn, P_svn, Q_organs, V_organs,...
%             P_organs, Q_brain, Q_kidney, Q_spleen, Q_liver, V_brain,...
%             V_kidney, V_spleen, V_liver, Q_pas, Q_pat, Q_par, Q_pcp,...
%             V_pas, V_pat, V_par, V_pcp, P_pas, P_pat, P_par, P_pcp,...
%             Q_pvn, V_pvn, P_pvn] = simulinkSimulator(...
%                             v1l, v2l, tdl, tfl, trl, pwl, perl, v1r, v2r, tdr,...
%                             tfr, trr, pwr, perr, C_sas, L_sas, R_sas, L_sat,...
%                             C_sat, R_sat, R_sar, R_scp, R_brain, R_liver,...
%                             R_spleen, R_kidney, R_svn, C_pas,...
%                             R_pas, L_pas, C_pat, R_pat, L_pat, R_par, R_pcp,...
%                             C_pvn, R_pvn);

    %trying the ODE
    %syms O2(t) OP(t) Qdot(t) sigmaO2 t
    
    OP = -0.3;
    load("Q_left_heart.mat");
    %disp(size(Q_left_heart))
    Qdot = Q_left_heart(round(t/0.0001,0)+1);
    sigmaO2 = (oconca - (OP0/Qdot))/975;
    dy = (OP + Qdot.*(oconca-sigmaO2.*y))./331;
end

function bo = brain(values)
    %computs change in blood through brain
    %inputs: values - blood values coming in
    %outputs: bo - blood values matrix altered by brain

    bo = values;
    
end