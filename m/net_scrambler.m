clc
clear all
close all
%% Load data
load('../mat/age.mat')
load('../mat/nki_mats.mat', 'coor', 'SC', 'ci')
keep = age <= 30 & age >= 25;
subjects = SC(:,:, keep);
for i = 1:size(subjects, 3)
   for j = 1:100
       scrambled_nets(:,:,i,j) = randmio_und(subjects(:,:,i), 16);
   end
end
save('../output/scrambled_nets.mat','scrambled_nets')

%% Run scrambled nets through control pipeline
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

n = size(SC, 1);  % number of nodes;
rho = 100;      % balance parameter -- determines relative contribution of "distance from target state" to "total input energy"
energy = nan(n, max(ci), max(ci), size(scrambled_nets, 3), size(scrambled_nets, 4));
totalenergy = zeros(max(ci), max(ci),size(scrambled_nets, 3), size(scrambled_nets, 4));
T = 1;          % time horizon -- the duration of time in which the control action takes place, i.e. time to go from x_0 to x_T


for i = 1:size(scrambled_nets, 3)
    for p = 1:size(scrambled_nets, 4)
    disp(i)
    disp(p)
    
    A = scrambled_nets(:,:,i,p);
    Aprime = (A/eigs(A,1)) - eye(n);    % need to transform connectivity matrix so all eigenvalues < 0
        for j = 1:max(ci)
            for k = 1:max(ci)
                x0 = ci == j;   % starting state -- here we set it so that nodes in system labeled "1" are all active
                xT = ci == k;  % target state -- here we set it so that nodes in system labeled "16" will be active
                %% choose input nodes
                B = eye(n);     % start with identity matrix
                keep = ci ~= j & ci ~= k;  % we'll use brain regions not in systems 1 and 16 as drivers
                B = B(:,keep);  % index just those regions
                %% run optimal control
                [x,u] = optimalControlContinuous(Aprime,B,rho,x0,xT,T); % x = trajectories; % u = inputs
                %% calculate energy
                energy(keep,j,k,i,p) = sum(u.^2,2);
                totalenergy(j,k,i,p) = sum(u(:).^2);  % for brain as a whole
            end
        end
    end
end