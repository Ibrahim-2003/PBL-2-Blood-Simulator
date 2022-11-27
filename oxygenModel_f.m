% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel_f.m
% Models O2 concentration through the human body in different conditions
% Last modified: 11/21/22

function oxygenModel_f
    %the driver
    %setting up initial values, need to compare this with the flow model
    tspan = 0:0.0001:57;
    %y0 = 0.195; %initial mL of O2
    y0 = 0.0068; %mmol/mL

    %running the ode
    [t,y] = ode45(@heartoxygencons,tspan,y0);
    plot(t,y)

    %check Q_left_heart
    load("Q_left_heart.mat");
    t = 0:0.0001:60;
    figure;
    plot(t,Q_left_heart)

    %initializing second ode
    %O0 = 0.0068; % mmol/mL
    O0 = 0.195; %mL O2/mL blood
    %C0 = 0.0020; % mmol/mL
    C0 = 0.492;
    %G0 = 0.0061; % mmol/mL
    G0 = 0.705; %calculated from mg/dL using 1.56 g/mL density

    x0 = [O0;C0;G0];

    %running second ode
    [t,x] = ode15s(@combined,tspan,x0);

    %separate variables
    O = x(:,1);
    C = x(:,2);
    G = x(:,3);
    
    figure;
    subplot(1,3,1);
    plot(t,O);

    subplot(1,3,2);
    plot(t,C);

    subplot(1,3,3);
    plot(t,G);
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

    Qdot = Q_left_heart(round(t/0.0001,0)+3)*0.05; %flow rate
    %sigmaO2 = (oconca - (OP/250))/oconcavg;
    %sigmaO2 = oconcv/oconcavg;
    sigmaO2 = 0.75; %lol we figured this out

    %the actual ode
    dy = (OP + Qdot.*(oconca-sigmaO2.*y))./vol;
    %disp(Qdot.*(oconca-sigmaO2.*y));
end

function xdot = combined(t,x)
    %an ode modeling component concentration in the left heart
    vol = 214; %mL

    O = x(1,1);
    C = x(2,1);
    G = x(3,1);

    o_a = 0.195;
    c_a = 0.492; % mmol/mL
    g_a = 0.705; % mmol/mL

    %values to input in ode
    %OC = -0.019;
    OC = 0.3567; %mL/sec
    %CC = 0.0152;
    CC = 0.2853;
    %GC = -0.00316;
    GC = 0.0594;

    load("Q_left_heart.mat"); %this is from simulink
    Qdot = Q_left_heart(round(t/0.0001,0)+3)*0.05; %flow rate

    sigmaO2 = (o_a - (OC/250))/0.145;
    %disp(sigmaO2);
    %sigmaCO2 = 0.75;
    sigmaCO2 = (c_a - (CC/250))/0.532;
    %disp(sigmaCO2);
    %sigmaG = 1.5;
    sigmaG = (g_a - (GC/250))/0.5;
    %disp(sigmaG);

    %the actual odes
    Odot = (-OC + Qdot.*(o_a-sigmaO2.*O))./vol;
    Cdot = (CC + Qdot.*(c_a-sigmaCO2.*C))./vol;
    Gdot = (-GC + Qdot.*(g_a-sigmaG.*G))./vol;

    xdot = [Odot;Cdot;Gdot];
end