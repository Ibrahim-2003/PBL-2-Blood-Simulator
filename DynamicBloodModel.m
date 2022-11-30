function DynamicBloodModel
t = 0:(0.01):100;
disp(cardiacOutputFraction(3850))
bloodVolume = shock(0.3,3300,t);
%plot(t,bloodVolume)

%Initial arterial concentrations 
O2_a = 0.195; % mL O2/mL blood
CO2_a = 0.492; % mL CO2/mL blood
G_a = 1099.8; % mg glucose/mL blood
C_a = [O2_a CO2_a G_a];

%Heart Values 
RQ_h = 0.8;
perCO_h = 0.05;
OC_h = 0.3567;

%Liver Values
RQ_l = 0.8;
perCO_l = 0.25;
OC_l = 6.08/60;
GG_l0 = 30000; %mg/s
%GG_l0 = 2.04293;

%Kidney 
RQ_k = 0.88;
perCO_k = 0.2;
OC_k = 0.35;
GG_k = 0.045399312;

%Brain 
RQ_b = 0.95;
perCO_b = 0.15;
OC_b = 49/60;

%Lung value 
RQ_lung = 0.84;
OC_lung = 0.09166667 ; %negative to represent generation 
perCO_lung = 0.02;




for i = 1:length(t)
GG_l = GG_l0 * ((bloodVolume(i)/5500)*-3.333333+4.33333); %Scale glucose consumption 
[C_h(i,:),Q_h(i)] = organConc(C_a(i,:),bloodVolume(i),perCO_h,OC_h,RQ_h,0,@heartDynamic);
[C_l(i,:),Q_l(i)] = organConc(C_a(i,:),bloodVolume(i),perCO_l,OC_l,RQ_l,GG_l,@liverDyanmic);
[C_b(i,:),Q_b(i)] = organConc(C_a(i,:),bloodVolume(i),perCO_b,OC_b,RQ_b,0,@brainDynamic);
[C_k(i,:),Q_k(i)] = organConc(C_a(i,:),bloodVolume(i),perCO_k,OC_k,RQ_k,GG_k,@kidneyDynamic);
[C_lung(i,:),Q_lung(i)] = organConc(C_a(i,:),bloodVolume(i),perCO_lung,OC_lung,RQ_lung,0,@lungDynamic);

C_body = [0.035042 0.7109339 0];

%Calculate venous concentrations
C_ven(i,:) = (C_h(i,:).*Q_h(i) + C_l(i,:).*Q_l(i) + C_b(i,:).*Q_b(i) + C_k(i,:).*Q_k(i) + C_lung(i,:).*Q_lung(i))/(44.2660*0.67);
disp((C_h(i,:).*Q_h(i) + C_l(i,:).*Q_l(i) + C_b(i,:).*Q_b(i) + C_k(i,:).*Q_k(i) + C_lung(i,:).*Q_lung(i))/(44.2660*0.67))
C_ven(i,:) = (C_ven(i,:) + C_body)/2;
%disp(C_h(i,:).*Q_h(i) + C_l(i,:).*Q_l(i) + C_b(i,:).*Q_b(i) + C_k(i,:).*Q_k(i) + C_lung(i,:).*Q_lung(i))
%Pulmonary circulation
[C_pul(i,:),Q_pul(i)] = pulmonaryCirculation(C_ven(i,:),bloodVolume(i));

%Oxygenated blood from pulmonary veins enter systemic circulation
C_a(i+1,:) = C_pul(i,:);

end
figure;

plot(t,C_a(1:(end-1),1),'LineWidth',2);
hold on;
plot(t,C_a(1:(end-1),2),'LineWidth',2);
hold off;
legend("Oxygen Concentration","Carbon Dioxide Concentration")
xlabel("Time (h)")
ylabel("Concentration (mL/mL of blood)")

figure;
plot(t,C_a(1:(end-1),3),'LineWidth',2);
xlabel("Time (h)")
ylabel("Concentration (mg/mL of blood)")
figure
%y = C_a(1:(end-1),1).*Q_b';

plot(t,C_a(1:(end-1),1).*Q_l','LineWidth',2)
hold on
plot(t,C_a(1:(end-1),1).*Q_b','LineWidth',2)
plot(t,C_a(1:(end-1),1).*Q_k','LineWidth',2)
plot(t,C_a(1:(end-1),1).*Q_lung','LineWidth',2)
plot(t,C_a(1:(end-1),1).*Q_h','LineWidth',2)

xlabel("Time (h)")
ylabel("Oyxgen Delivery (mL/s)")
legend("Liver","Brain","Kidneys","Lungs","Myocardium")
 hold off



%plot(t,C_v(1,:))

end



function [C_v,Q] = organConc(C_a,volume,perc_CO,OC,RQ,GG,dynamicFun)
%Calculates the concentration of compoenents in the vein leading draining
%the organ
%Inputs:
    %time - the time period at which the model is evaluated 
    %C_a - a matrix of the arterial concentrations 
        %C_art(:,1) - oxygen
        %C_art(:,2) - carbon dioxide
        %C_art(:,3) - glucose 
    %volume - volume of blood inside the body
%Outputs:
    %C_v - a matrix of the venous concentrations 
        %C_art(1) - oxygen
        %C_art(2) - carbon dioxide
        %C_art(3) - glucose 
    %Q - blood flow into and out of the organ

%Calculate blood flow to the myocardium
[Q_frac, VO2_frac] = dynamicFun(volume);
Q = 44.2660 * perc_CO * Q_frac;

%Consumption rates (mL/sec)
VO2 =  OC*VO2_frac; 

VCO2 = RQ * VO2; %0.8 is respiratory quotient 
VG = -((VCO2*0.0425)/6*180.156) + GG; %stoichiometric ratio for glycolysis 

%Accounting equation
C_v(1) = C_a(1) - VO2/Q;
C_v(2) = C_a(2) - VCO2/Q;
C_v(3) = C_a(3) + VG/Q;

end





function [bloodVolume] = shock(bloodLoss,lossRate,t) 
%This function models blood volume changes with regards to an episode of
%hemorrhagic shock
%Input:
    %bloodLoss - the percentage of bood lost (%)
    %lossRate - the rate at which blood is lost (mL/h)
    %time - an array for the time interval (h)
%Output:
    %bloodVolume - an array for the volume of blood at different times    


bloodVolume = zeros(size(t));
%Initial blood volume (mL)
volume_0 = 5500;

%Steady-state
t_0 = 5; %time at which perturbation begins
for i = find(t<=t_0)
    bloodVolume(i) = volume_0;
end

%Active blood loss 
volume_f = 5500 * (1 - bloodLoss);
tspan_bleed = abs((volume_f - volume_0)/lossRate);

for i = find((t <= tspan_bleed+ t_0) & t > t_0)
    bloodVolume(i) = volume_0 - lossRate * (t(i)-t_0);
end

%Recovery
for i = find(t > tspan_bleed + t_0)
    bloodVolume(i) = bloodVolumeRecovery(t(t == tspan_bleed+t_0),t(i),volume_f);
end




end

function volume = bloodVolumeRecovery(t0, tf,V0)
%Calculates an array of volumes of blood at different time periods after an
%inital blood loss
%Output:
    %volume - an array containing the volume of blood in mL in entire body
%Input: 
    %t0 - initial time when blood loss occurs (h)
    %tf - time during which instantaenous blood volume should be calculated (h)
    %volume_i - the amount of blood lost at t0 (mL)

%Set up variables
t = (tf - t0)+0.1; %time (h) 
normalBV = 5500; %normal blood volume (mL)
 
k = 0.05; %blood replenishment rate

volume = normalBV./((normalBV-V0)/V0*exp(-k*t)+1);

end

function HR = heartRate(volume,HR0)

%Calculates the heart rate based on the current blood volume
%Outputs:
    %HR - heart rate in bpm
%Input:
    %volume - the volume of blood in mL
    %HR0 - the baseline heart rate


BL = 5500 - volume; %Calculate blood loss
LBNP = BL / 17.2619; %Convert blood loss to LBNP equivalent
HR = 0.18*LBNP + HR0;

end

function CO = cardiacOutputFraction(volume)
%Calculates fractional change in cardiac output adjusted for blood loss 
%Output:
    %CO - cardiac output fractional change (multiply with baseline cardiac
    %output for adjusted cardiac output
%Inputs: 
    
    %Volume - blood volume


BL = 5500 - volume; %Calculate blood loss
%LBNP = BL / 17.2619; %Convert blood loss to LBNP equivalent
%CO = -0.0075*LBNP + 1;

CO = -0.0001377 * BL + 1.01219;
end

function [C_v,Q] = pulmonaryCirculation(C_a,volume)

Q =  44.2660 *cardiacOutputFraction(volume);
%Pulmonary circulation values 
OG_0 = 4.73333; %oxygen entering the blood from the aveoli (mL/s)
CC_0 = 3.78333;
%fracBloodLoss = 1 - volume/5500;
fracBlood = volume /5500;
%24 breadth/min with 30% blood loss, 18 breadth/min normal
%respiratoryRate = (20*fracBloodLoss + 18)/60;

%OG = OG_0/(18/60)*respiratoryRate;
%CC = CC_0/(18/60)*respiratoryRate;

OG = (1.83333*fracBlood - 0.83333)*OG_0;
CC = (1.83333*fracBlood - 0.83333)*CC_0;
C_v(1) = C_a(1) + OG/Q;
C_v(2) = C_a(2) - CC/Q;
C_v(3) = C_a(3);

end

function [Q,VO2] = heartDynamic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

%Adjusts blood flow rate to the heart
CO = cardiacOutputFraction(volume);
Q = CO; %Blood flow change to the heart is 1:1 with to cardiac output 
VO2 = 0.382 * CO + 0.618; %oxygen consumption decreases during shock 
end

function [Q,VO2] = brainDynamic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

CO = cardiacOutputFraction(volume);


Q = 1; %Blood flow change does not change
%VO2 = (CO * -1.5 +2.5); %oxygen consumption increases during shock
BL = volume/5500;
VO2 = 0.83333*BL + 0.1666667;
end

function [Q,VO2] = lungDynamic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

CO = cardiacOutputFraction(volume);
fractionalVolume = volume / 5500;
Q = fractionalVolume * 2.1205 - 1.1259;
VO2 = 1; %assume lung fractional conversion stays the same


end

function [Q,VO2] = liverDyanmic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

CO = cardiacOutputFraction(volume);

Q = CO; 
VO2 = CO * 1.7 - 0.7;

end

function [Q,VO2] = kidneyDynamic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

CO = cardiacOutputFraction(volume);

Q = CO * 1.82 - 0.82; 
VO2 = CO * 1.26 - 0.26; 
end


function [C_v,Q] = myocardium(C_a,volume)
%Calculates the concentration of compoenents in the vein leading draining
%the organ
%Inputs:
    %time - the time period at which the model is evaluated 
    %C_a - a matrix of the arterial concentrations 
        %C_art(:,1) - oxygen
        %C_art(:,2) - carbon dioxide
        %C_art(:,3) - glucose 
    %volume - volume of blood inside the body
%Outputs:
    %C_v - a matrix of the venous concentrations 
        %C_art(1) - oxygen
        %C_art(2) - carbon dioxide
        %C_art(3) - glucose 
    %Q - blood flow into and out of the organ

%Calculate blood flow to the myocardium
[Q_frac, VO2_frac] = heartDynamic(volume);
Q = 44.2660 * 0.05 * Q_frac;

%Consumption rates (mL/sec)
VO2 =  0.3567*VO2_frac; 
VCO2 = 0.8 * VO2; %0.8 is respiratory quotient 
VG = VO2/6; %stoichiometric ratio for glycolysis 

%Accounting equation
C_v(1) = C_a(1) - VO2/Q;
C_v(2) = C_a(2) - VCO2/Q;
C_v(3) = C_a(3) - VG/Q;

end
