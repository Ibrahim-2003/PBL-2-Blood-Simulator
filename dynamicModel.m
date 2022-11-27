function dynamicModel
    t = 0:0.01:100;
    volume = bloodVolumeRecovery(0,t,1650);
    plot(t,volume)
    HR = heartRate(volume,60);
    figure;
    plot(t,HR)
   
   
end

function volume = bloodVolumeRecovery(t0, tf,bloodLoss)
%Calculates an array of volumes of blood at different time periods after an
%inital blood loss
%Output:
    %volume - an array containing the volume of blood in mL in entire body
%Input: 
    %t0 - initial time when blood loss occurs (h)
    %tf - time during which instantaenous blood volume should be calculated (h)
    %bloodLoss - the amount of blood lost at t0 (mL)

%Set up variables
t = (tf - t0); %time (h) 
normalBV = 5500; %normal blood volume (mL)
V0 = 5500 - bloodLoss; 
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

function CO = cardiacOutput(CO0,volume)
%Calculates cardiac output adjusted for blood loss 
%Output:
    %CO - adjusted cardiac output
%Inputs: 
    %CO0 - baseline cardiac output
    %Volume - blood volume
BL = 5500 - volume; %Calculate blood loss
LBNP = BL / 17.2619; %Convert blood loss to LBNP equivalent
CO = 0.18*LBNP + HR0;
end

