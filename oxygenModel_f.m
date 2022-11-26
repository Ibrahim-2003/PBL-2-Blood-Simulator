% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel_f.m
% Models O2 concentration through the human body in different conditions
% Last modified: 11/21/22

function oxygenModel_f
    %the driver
    %setting up initial values, need to compare this with the flow model
    tspan = 0:0.01:2497;
    %y0 = 0.195; %initial mL of O2
    y0 = 0.0068; %mmol/mL

    %running the ode
    [t,y] = ode45(@heartoxygencons,tspan,y0);
    plot(t,y)
    figure;

    %check Q_left_heart
    load("Q_left_heart.mat");
    t = 0:0.01:2500;
    plot(t,Q_left_heart)
end

function dy = heartoxygencons(t,y)
    %an ode modeling oxygen volume in the left heart
    vol = 214; %mL
    %oconca = 0.195; %mL O2/mL blood
    oconca = 0.017; %mmol/mL
    oconcv = 0.145;
    oconcavg = 0.17;
  
    %values to input in ode
    %OP = -0.43; %consumption rate mL/sec
    OP = -0.019; %consumption rate mmol/sec
    load("Q_left_heart.mat"); %this is from simulink

    Qdot = Q_left_heart(round(t/0.01,0)+3)*0.05; %flow rate
    %sigmaO2 = (oconca - (OP/250))/oconcavg;
    %sigmaO2 = oconcv/oconcavg;
    sigmaO2 = 0.3; %lol we figured this out

    %the actual ode
    dy = (OP + Qdot.*(oconca-sigmaO2.*y))./vol;
    %disp(Qdot.*(oconca-sigmaO2.*y));
end

function dy = simplified(t,y)
    %an ode modeling oxygen volume in the left heart
    vol = 214; %mL
    oconca = 0.195; %mL O2/mL blood
    oconcv = 0.145;
    oconcavg = 0.17;

    %values to input in ode
    %OP = -0.43; %consumption rate mL/min
    load("Q_left_heart.mat"); %this is from simulink

    Qdot = Q_left_heart(round(t/0.01,0)+3)*0.05; %flow rate
    OP = mean(Q_left_heart(round(20/0.001,0)+3:end));
    %sigmaO2 = (oconca - (OP/250))/oconcavg;
    %sigmaO2 = oconcv/oconcavg;
    sigmaO2 = 0.3; %lol we figured this out

    %the actual ode
    dy = (OP + Qdot.*(oconca-oconcv))./vol;
end