% PBL Team 1, BIOE 252, Fall 2022, PBL 2: Modeling Blood Flow
% rbcModel.m
% Models rbc concentration through the human body in different conditions
% Last modified: 11/20/22

function plateletModel
    %The driver

    %setting up initial values
    time = (1:200);
    n_platelet = zeros(200,1);
    n_platelet(1,:) = 1.68E12; 
    %flow = ones(86400);
    
    %flow_rbc(1,:) = 5.4E9*flow; %initial flow rate

%     n_rbc = marrow(time,n_rbc);
%     plot(time,n_rbc)
    n_platelet = spleen_p(time,n_platelet);
    figure;
    plot(time,n_platelet)
    xlabel('Time (min)')
    ylabel('Number of Platelets')
    title('Platelet Recovery')
    saveas(gcf, 'platelet.png')

end

% function n_rbc = marrow(time,n_rbc)
%     for i = 2:length(time)
%         n_rbc(i,:) = n_rbc(i-1,:) + (5.555E-4).*n_rbc(i-1,:);
%     end
% end

function n_platelet = spleen_p(time,n_platelet)
    f = 0.020833;
    c = 35E9;

    for i = 2:length(time)
        n_platelet(i,:) = (1-f)^i .* (n_platelet(1,:) - c/f) + c/f;
    end
end
