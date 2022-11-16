% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% oxygenModel.m
% Models O2 concentration through the human body in different conditions
% Last modified: 10/31/22

function oxygenModel
    %the driver
    %setting up initial values, need to compare this with the flow model
    time = transpose(1:5);
    flow = ones(5,1);     
    vol = zeros(5,1); 
    vol(1,1) = 5000; %mL
    pressure = zeros(5,1);
    pressure(1,1) = 15; %mmHg
    oconca = 0.195*vol;
    oconc = 0.195*vol;
    oconcv = 0.145*vol;
    co2conc = 0.492*vol;

    init = table(time,flow,vol,pressure,oconca,oconc,oconcv,co2conc) %initial values vector
    ho = heart(init);
    %bo = brain(ho);
    
end

function ho = heart(values)
    %computes change in blood through heart
    %inputs: values - matrix of blood values coming in
    %outputs: ho - blood values matrix altered by heart
    
    ho = values;
    t = ho.time;

    %Simple model
    %oxygen consumption
    ocons = ho.oconc - 10*t;
    ho.oconc = ocons;
    %co2 generation, assuming all o2 consumption goes towards making co2
    co2gen = ho.co2conc + ocons;
    ho.co2conc = co2gen;
    ho
    
    %My attempt at the model from the paper
    %we might wanna specify these in the driver or just assume that they
    %are one value
    weight = 150;
    height = 180;
    age = 21;
    OP0 = 88.362 + 13.397*weight + 4.799*height - 5.677*age; %BMR

    sigmaO2 = (ho.oconca(1,:) - (OP0/ho.flow(1,:)))/ho.oconc(1,:);
    OP = -33.1*t;
    Qdot = 150*t;

    %trying the ODE
    syms O2(t) 
    ode = diff(O2,t) == (OP + Qdot.*(ho.oconca-sigmaO2.*ho.oconc))./331;
    O2sol(t) = dsolve(ode)

end

function bo = brain(values)
    %computs change in blood through brain
    %inputs: values - blood values coming in
    %outputs: bo - blood values matrix altered by brain

    bo = values;
    
end