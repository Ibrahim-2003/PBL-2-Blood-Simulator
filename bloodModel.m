% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% bloodModel.m
% Models blood flow through the human body in different conditions
% Last modified: 10/30/22

function bloodModel
    %the driver
    init = [0 0 0 0 0 0 0 0 0]; % Initial Conditions
    ho = heart(init);
    brain(ho);
    
    % Compliances and Blood Volumes
    Ca = 0.3; % mL/mmHg
    Cs = 3; % mL/mmHg
    Cv = 61.11; % mL/mmHg
    V0 = 5300; % mL
    Vu = 4700; % mL
    Vu_vent = 16.77; % mL
end

function healthy_simple_circulation(values)
    %computes the blood flow dynamics in circulatory system using "simple"
    %4-compartment 0D model from Cardiovascular Mathematics book
    %inputs: values - matrix of blood values coming in
    %outputs: None, output figures of results
    
    % System of ODEs
    % Ca * (dPa/dt) = Qa - (Pa-Ps)/Ra
    % Cs * (dPs/dt) = (Pa-Ps)/Ra - (Ps-Pv)/Rs
    % Cv * (dPv/dt) = Qa - Ca*(dPa/dt) - Cs*(dPs/dt)
    % Qsm = (Ps-Pv)/Rsm
    % Qsp = (Ps-Pv)/Rsp
    % Qo1 = (Ps-Pv)/Ro1
    % Qo2 = (Ps-Pv)/Ro2
    % Qo3 = (Ps-Pv)/Ro3
    % Qo4 = (Ps-Pv)/Ro4
    % Qo5 = (Ps-Pv)/Ro5
    % Qa = Vstr/T
    % Vstr = Ved(Pv) - Vu_vent - Pa/E
    
    % Variables
    % Vstr = Stroke Volume of Heart
    % T = Heart Period
    % Ved = End-of-Diastole Volume
    % Pv = Venous Pressure of Vu_vent
    % Vu_vent = Unstressed Ventricular Volume
    % E = Elastance
    
    
    

end

function ho = heart(values)
    %computes change in blood through heart
    %inputs: values - matrix of blood values coming in
    %outputs: ho - blood values matrix altered by heart
    
    ho = values;

end

function bo = brain(values)
    %computs change in blood through brain
    %inputs: values - blood values coming in
    %outputs: bo - blood values matrix altered by brain

    bo = values;
    
end

function luo = lung(values)
    %computs change in blood through lungs
    %inputs: values - blood values coming in
    %outputs: luo - blood values matrix altered by lungs
    
    luo = values;

end

function lio = liver(values)
    %computs change in blood through liver
    %inputs: values - blood values coming in
    %outputs: lio - blood values matrix altered by liver
    
    lio = values;

end

function ko = kidney(values)
    %computs change in blood through kidney
    %inputs: values - blood values coming in
    %outputs: ko - blood values matrix altered by kidney
    
    ko = values;

end

function mo = marrow(values)
    %computs change in blood through bone marrow
    %inputs: values - blood values coming in
    %outputs: mo - blood values matrix altered by bone marrow
    
    mo = values;

end

function so = spleen(values)
    %computs change in blood through spleen
    %inputs: values - blood values coming in
    %outputs: so - blood values matrix altered by spleen
    
    so = values;

end
