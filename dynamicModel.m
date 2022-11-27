function dynamicModel
    t = 0:0.01:100;
    volume = bloodVolumeRecovery(0,t,1650);
    plot(t,volume)
    HR = heartRate(volume,60);
    figure;
    plot(t,HR)

    
  
    %disp(mean(Q_left_heart(2500:end)))

   
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

function CO = cardiacOutputFraction(volume)
%Calculates fractional change in cardiac output adjusted for blood loss 
%Output:
    %CO - cardiac output fractional change (multiply with baseline cardiac
    %output for adjusted cardiac output
%Inputs: 
    
    %Volume - blood volume


BL = 5500 - volume; %Calculate blood loss
LBNP = BL / 17.2619; %Convert blood loss to LBNP equivalent
CO = -0.0075*LBNP + 1;
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
VO2 = (CO * -1.5 +2.5); %oxygen consumption increases during shock
end

function Q = lungDynamic(volume)
%Calcualates fractional changes for organ blood flow and oxygen consumption
%based on the current blood volume
%Outputs:
    %Q - the fractional adjustment coefficient for heart rate
    %VO2 - the fractional adjustment coefficeint for oxygen consumption
%Inputs: 
    %volume - current blood volume in the body (not blood lost!)

CO = cardiacOutputFraction(volume);
fractionalVolume = volume / 5500;
Q = fractionalVolume * 2.1205 - 1.1259; %Blood flow change does not change

%Did not include lung fractional oxygen consumption because assumes same
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




