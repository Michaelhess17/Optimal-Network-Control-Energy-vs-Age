clear all
close all
clc
% Selected Nodes as Drivers (SD)
%% load network data
load('../mat/nki_mats.mat'); % SC = interregional structural connectivity matrix; ci = system labels
ciprime = zeros(size(ci));
ciprime(ci >= 1 & ci <= 3) = 1;
ciprime(ci >= 4 & ci <= 6) = 2;
ciprime(ci >= 7 & ci <= 8) = 3;
ciprime(ci == 9) = 4;
ciprime(ci >= 10 & ci <= 11) = 5;
ciprime(ci >= 12 & ci <= 13) = 6;
ciprime(ci == 14) = 7;
ciprime(ci >= 15 & ci <= 16) = 8;
ciprime(ci == 17) = 9;
ci = ciprime;
%%
n = size(SC, 1);  % number of nodes;
rho = 100;      % balance parameter -- determines relative contribution of "distance from target state" to "total input energy"
energySD = nan(n, max(ci), max(ci), length(SC));
totalenergySD = zeros(max(ci), max(ci),length(SC));
keep = ci == 1 | ci == 3 | ci == 5;  % we'll use brain regions  in systems 1, 3 and 5 as drivers
T = 1;          % time horizon -- the duration of time in which the control action takes place, i.e. time to go from x_0 to x_T
%% scale matrix
for i = 1:length(SC)
    disp(i)
    A = SC(:,:,i);
    Aprime = (A/eigs(A,1)) - eye(n);    % need to transform connectivity matrix so all eigenvalues < 0
    for j = 1:max(ci)
        for k = 1:max(ci)
            x0 = ci == j;   % starting state -- here we set it so that nodes in system labeled "1" are all active
            xT = ci == k;  % target state -- here we set it so that nodes in system labeled "16" will be active
            %% choose input nodes
            B = eye(n);     % start with identity matrix
            B = B(:,keep);  % index just those regions
            %% run optimal control
            [x,u] = optimalControlContinuous(Aprime,B,rho,x0,xT,T); % x = trajectories; % u = inputs
            %% calculate total energy       
            energySD(keep,j,k,i) = sum(u.^2,2);   % for each node
            totalenergySD(j, k, i) = sum(u(:).^2);  % for brain as a whole
        end
    end
end
%%
save('../output/energySD.mat','energySD', 'totalenergySD')
%% All Nodes as Drivers (AD) 
n = size(SC, 1);  % number of nodes;
rho = 100;      % balance parameter -- determines relative contribution of "distance from target state" to "total input energy"
energyAD = nan(n, max(ci), max(ci), length(SC));
totalenergyAD = zeros(max(ci), max(ci),length(SC));
T = 1;          % time horizon -- the duration of time in which the control action takes place, i.e. time to go from x_0 to x_T
%% scale matrix
for i = 1:length(SC)
    disp(i)
    A = SC(:,:,i);
    Aprime = (A/eigs(A,1)) - eye(n);    % need to transform connectivity matrix so all eigenvalues < 0
    for j = 1:max(ci)
        for k = 1:max(ci)
            x0 = ci == j;   % starting state -- here we set it so that nodes in system labeled "1" are all active
            xT = ci == k;  % target state -- here we set it so that nodes in system labeled "16" will be active
            %% choose input nodes
            B = eye(n);     % start with identity matrix
%             B = B(:,keep);  % index just those regions
            %% run optimal control
            [x,u] = optimalControlContinuous(Aprime,B,rho,x0,xT,T); % x = trajectories; % u = inputs
            %% calculate total energy
%             energy{j, k, i} = sum(u.^2,2);       % for each node
            energyAD(:,j,k,i) = sum(u.^2,2);
            totalenergyAD(j, k, i) = sum(u(:).^2);  % for brain as a whole
        end
    end
end
save('../output/energyAD.mat','energyAD', 'totalenergyAD')