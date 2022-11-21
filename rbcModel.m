% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% rbcModel.m
% Models rbc concentration through the human body in different conditions
% Last modified: 11/20/22

function rbcModel
    %The driver

    %setting up initial values
    time = (1:10000);
    n_rbc = zeros(10000,1);
    n_rbc(1,:) = 5.5E12;
    %flow = ones(86400);
    
    %flow_rbc(1,:) = 5.4E9*flow; %initial flow rate

%     n_rbc = marrow(time,n_rbc);
%     plot(time,n_rbc)
    n_rbc = spleen(time,n_rbc);
    figure;
    plot(time,n_rbc)

end

% function n_rbc = marrow(time,n_rbc)
%     for i = 2:length(time)
%         n_rbc(i,:) = n_rbc(i-1,:) + (5.555E-4).*n_rbc(i-1,:);
%     end
% end

function n_rbc = spleen(time,n_rbc)
    f = 5.555E-4;
    c = 0.347E9;

    for i = 2:length(time)
        n_rbc(i,:) = (1-f)^i .* (n_rbc(1,:) - c/f) + c/f;
    end
end
