% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel_f.m
% Models O2 concentration through the human body in different conditions
% Last modified: 11/21/22

function oxygenModel_f
    %the driver
    %setting up initial values, need to compare this with the flow model
    tspan = 0:0.0001:15;
    y0 = 25.935;

    %running the ode
    [t,y] = ode45(@heartoxygencons,tspan,y0);
    plot(t,y)
    figure;

    %check Q_left_heart
    load("Q_left_heart.mat");
    plot(t,Q_left_heart)
end

function dy = heartoxygencons(t,y)
    %an ode modeling oxygen concentration in the left heart
    vol = 133; %mL
    oconca = 0.195*vol;
    oconcv = 0.145*vol;

    weight = 150; %you should probably check these too
    height = 180;
    age = 21;
    OP0 = 88.362 + 13.397*weight + 4.799*height - 5.677*age; %BMR
    
    %values to input in ode
    OP = -0.3; %consumption rate
    load("Q_left_heart.mat"); %this is from simulink

    Qdot = Q_left_heart(round(t/0.0001,0)+1); %flow rate
    sigmaO2 = (oconcv - (OP0/Qdot))/975; %partition coefficient?

    %the actual ode
    dy = (OP + Qdot.*(oconcv-sigmaO2.*y))./331;
end