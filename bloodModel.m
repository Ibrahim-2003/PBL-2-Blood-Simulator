% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% bloodModel.m
% Models blood flow through the human body in different conditions
% Last modified: 10/28/22

function bloodModel
    %the driver
    init = %initial conditions;
    heart(init);
    brain(ho);
end

function heart(values)
    %computes change in blood through heart
    %inputs: values - matrix of blood values coming in
    %outputs: ho - blood values matrix altered by heart

end

function brain(values)
    %computs change in blood through brain
    %inputs: values - blood values coming in
    %outputs: bo - blood values matrix altered by brain

end

function lung(values)
    %computs change in blood through lungs
    %inputs: values - blood values coming in
    %outputs: luo - blood values matrix altered by lungs

end

function liver(values)
    %computs change in blood through liver
    %inputs: values - blood values coming in
    %outputs: lio - blood values matrix altered by liver

end

function kidney(values)
    %computs change in blood through kidney
    %inputs: values - blood values coming in
    %outputs: ko - blood values matrix altered by kidney

end

function marrow(values)
    %computs change in blood through bone marrow
    %inputs: values - blood values coming in
    %outputs: mo - blood values matrix altered by bone marrow

end

function spleen(values)
    %computs change in blood through spleen
    %inputs: values - blood values coming in
    %outputs: so - blood values matrix altered by spleen

end
