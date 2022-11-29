% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel_h.m
% Models O2 concentration through the left heart
% Last modified: 11/26/22

function oxygenModel_h
    %the driver
    %setting up initial values based on flow model
    tspan = 0:0.0001:57;
    
    %initial concentrations (arterial concentration)
    %O0 = 0.0068; % mmol/mL
    O0 = 0.05; %mL O2/mL blood
    %C0 = 0.0020; % mmol/mL
    C0 = 0.492;
    %G0 = 0.0061; % mmol/mL
    G0 = 0.705; %calculated from mg/dL using 1.56 g/mL density

    x0 = [O0;C0;G0];

    %running ode
    [t,x] = ode15s(@combined,tspan,x0);

    %separate variables
    O = x(:,1);
    C = x(:,2);
    G = x(:,3);
    
    %plot figure
    figure(1);
    subplot(1,3,1);
    plot(t,O);
    xlabel('Time (s)');
    ylabel('O2 Conc (mL O2/mL blood)');

    subplot(1,3,2);
    plot(t,C);
    xlabel('Time (s)');
    ylabel('CO2 Conc (mL CO2/mL blood)');

    subplot(1,3,3);
    plot(t,G);
    xlabel('Time (s)');
    ylabel('Glucose Conc (mL glucose/mL blood)');

    %check Q_left_heart
    load("Q_left_heart.mat");
    t = 0:0.0001:60;
    figure;
    plot(t,Q_left_heart)
    disp(mean(Q_left_heart(find(t>3))))
end

function xdot = combined(t,x)
    %an ode modeling component concentration in the left heart
    vol = 214; %mL of left heart

    O = x(1,1); %initial concentrations
    C = x(2,1);
    G = x(3,1);

    o_a = 0.195; % mL O2/mL blood
    c_a = 0.492; % mL CO2/mL blood
    g_a = 0.705; % mL glucose/mL blood

    %consumption to input in ode
    %OC = -0.019;
    OC = 0.3567; %mL/sec
    %CC = 0.0152;
    CC = 0.2853;
    %GC = -0.00316;
    GC = 0.0594;

    load("Q_left_heart.mat"); %this is from simulink
    Qdot = Q_left_heart(round(t/0.0001,0)+3)*0.05; %flow rate without the weirdness at the beginnin
    
    sigmaO2 = (o_a - (OC/250))/0.145; %partition coefficient
    sigmaO2 = 0.095/0.05;
    %disp(sigmaO2);
    %sigmaCO2 = 0.75;
    sigmaCO2 = (c_a - (CC/250))/0.532;
    %disp(sigmaCO2);
    %sigmaG = 1.5;
    sigmaG = (g_a - (GC/250))/0.5;
    %disp(sigmaG);

    %the actual odes
    %Odot = (-OC + (3.567)*(cos(t)+1)*(o_a-sigmaO2.*O))./vol;
    Odot = (-OC + Qdot*(o_a-sigmaO2.*O))./vol;
    Cdot = (CC + Qdot.*(c_a-sigmaCO2.*C))./vol;
    Gdot = (-GC + Qdot.*(g_a-sigmaG.*G))./vol;

    xdot = [Odot;Cdot;Gdot]; %return vector
end

function xdot = heart(t,x)
    %an ode modeling component concentration in the left heart
    vol = 214; %mL of left heart

    O = x(1,1); %initial concentrations
    C = x(2,1);
    G = x(3,1);

    o_a = 0.195; % mL O2/mL blood
    c_a = 0.492; % mL CO2/mL blood
    g_a = 0.705; % mL glucose/mL blood

    %consumption to input in ode
    %OC = -0.019;
    OC = 0.3567; %mL/sec
    %CC = 0.0152;
    CC = 0.2853;
    %GC = -0.00316;
    GC = 0.0594;

    load("Q_left_heart.mat"); %this is from simulink
    Qdot = Q_left_heart(round(t/0.0001,0)+3)*0.05; %flow rate without the weirdness at the beginnin
    
    sigmaO2 = (o_a - (OC/250))/0.145; %partition coefficient
    sigmaO2 = 0.095/0.05;
    %disp(sigmaO2);
    %sigmaCO2 = 0.75;
    sigmaCO2 = (c_a - (CC/250))/0.532;
    %disp(sigmaCO2);
    %sigmaG = 1.5;
    sigmaG = (g_a - (GC/250))/0.5;
    %disp(sigmaG);

    %the actual odes
    %Odot = (-OC + (3.567)*(cos(t)+1)*(o_a-sigmaO2.*O))./vol;
    Odot = (-OC + Qdot*(o_a-sigmaO2.*O))./vol;
    Cdot = (CC + Qdot.*(c_a-sigmaCO2.*C))./vol;
    Gdot = (-GC + Qdot.*(g_a-sigmaG.*G))./vol;

    xdot = [Odot;Cdot;Gdot]; %return vector
end

