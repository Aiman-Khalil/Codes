% Clear the workspace
clc;
clear;
% Set the random number generator seed for reproducibility
seed = 42; % Replace 42 with the desired seed value
rng(seed);

% Parameters
De = 15;   % (μm^2/s; diffusion coefficient for exosomal tRNAs in Brain Extracellular space)
Dm = 0.3774;   % (μm^2/s; diffusion coefficient for exosomal tRNAs in PES membrane;)
Ds= 0.3668;    % (μm^2/s; effective diffusion coefficient for exosomal tRNAs in hydrogel scaffold)
epsilon = 0.9;    % ( Membrane porosity;)
totuosity = 2;  % ( Membrane Tortuosity;)
D_m = Ds * epsilon / totuosity;  % (μm^2/sec; effective diffusion coefficient for exosomal tRNAs in membrane)
Pm=0.001; %(μm/s; Permeability coefficient;)

Ps=0.00012;%(μm/s; Permeability coefficient_scaffold)


Re = 95;      % (μm, brain extracellular radius)
Rm = 85;      % (μm, membrane radius)
Rs = 150;      % (μm, scaffold radius)
lambda = 1.67; % (Tortuosity of the Brain Extracellular Space) ref:4]
N = 700;      % (molecules)
Vmax= 2.84e-6;  % (molecules/cell.s.;  uptake rate;)
N0= 4285.71;      % (3000 cells/0.7 micorL or 4285 cells/mm^3)
k=3.8e14; % Half_Saturation Constant

% Number of seeds
num_seeds = 10;

% Spatial domain
% Spatial domain
num_points_Re = 161;   % Number of grid points
num_points_Rm = 111;   % Number of grid points
num_points_Rs = 101;   % Number of grid points




dr_Re = Re / (num_points_Re-1);
dr_Rm = (Rm - Re) / num_points_Rm;
dr_Rs = (Rs - Rm) / num_points_Rs;


r_Re = 0:dr_Re:Re;           % Vector of x-axis (r_Re) values to be used for plotting
r_Rm = Re + dr_Rm:dr_Rm:Rm;  % Vector of x-axis (r_Rm) values to be used for plotting
r_Rs = Rm + dr_Rs:dr_Rs:Rs;  % Vector of x-axis (r_Rs) values to be used for plotting


%temporal grid

dt = 0.01;
t=0:dt:500;
numt=length(t);

% Arrays to store results for each seed
C_results = zeros(num_points_Re, numt, num_seeds);
Cm_results = zeros(num_points_Rm, numt, num_seeds);
Cs_results = zeros(num_points_Rs, numt, num_seeds);

for seed = 1:num_seeds
   
    % Initialize concentrations for each seed
   C = zeros(num_points_Re, numt);   % Initialize everything to zero
   Cm = zeros(num_points_Rm, numt);  % Initialize everything to zero
   Cs = zeros(num_points_Rs, numt);  % Initialize everything to zero





mu = 0  ;
sigma = 0.05;


  
    lamb = N;  % Mean number of molecules (lambda) is set to the initial release rate
    random_release_count = poissonRandom(lamb);
     C(:, 1) = random_release_count * exp(-(r_Re).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);

     
% Iterating the diffusion equation in the brain extracellular space (Re)
for j = 1:numt
    % Initialize source term (only for time = 0)
   source_term = zeros(num_points_Re, 1);
   if j == 1
    % Calculate the number of molecules to be released in this time step using a Poisson distribution
    lamb = N;  % Mean number of molecules (lambda) is set to the initial release rate
    random_release_count = poissonRandom(lamb);

    % Now, you can distribute these random_release_count molecules spatially
    source_term = random_release_count * exp(-(r_Re).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);
     
   end
    for i = 2:num_points_Re-1
        C(i, j+1) = C(i, j) + ((De/lambda^2) * dt / dr_Re^2) * (C(i+1, j) - 2 * C(i, j) + C(i-1, j))+source_term(i)*dt;
    end
    C(end, j) = C(end-1, j); % Set the boundary value
   
end
% Set the boundary value and initialize concnetration for Rm
Cm(1, 1) = C(end, numt);
Cm(1, j+1) = Pm*C(1, j+1);

% Iterating the diffusion equation in the membrane (Rm)
for j = 1:numt
    
    for i = 2:num_points_Rm-1
        d2Cm_dr2 = (Cm(i + 1, j) - 2 * Cm(i, j) + Cm(i - 1, j)) / dr_Rm^2;
        dCm_dr = (Cm(i + 1, j) - Cm(i - 1, j)) / (2 * dr_Rm);
        Cm(i, j+1) = Cm(i, j) + D_m * (d2Cm_dr2 + (1 / r_Rm(i)) * dCm_dr) * dt;
    end
    Cm(end, j) = Cm(end - 1, j); % Set the boundary value
end


Cs(1, 1) = Cm(end, numt);
Cs(1, j+1) = Pm*Cm(1, j+1);

for j = 1:numt
    for i = 2:num_points_Rs-1
        d2Cs_dr2 = (Cs(i + 1, j) - 2 * Cs(i, j) + Cs(i - 1,  j)) / dr_Rs^2;
        dCs_dr = (Cs(i + 1,  j) - Cs(i - 1,  j)) / (2 * dr_Rs);

        % Reaction term
        
        reaction_term = (Vmax* N0 * (Cs(i, j) / (k + Cs(i, j))));
        
      

         % Update Cs_Rs using the concentration from the previous time step
        Cs(i, j+1) = Cs(i, j) + Ds * (d2Cs_dr2 + (1/r_Rs(i)) * dCs_dr)*dt - reaction_term*dt;
    end
    Cs(end, j) = Cs(end - 1, j); % Set the boundary value
    
   
 end




 % Store results for this seed
    C_results(:, :, seed) = C(:, 1:numt);
    Cm_results(:, :, seed) = Cm(:, 1:numt);
    Cs_results(:, :, seed) = Cs(:, 1:numt);






% Save concentration data to CSV files with headers
    header = {'Time', 'Concentration'};

    % Brain ECS
     %data_BrainECS = [t', C(end, 1:numt)'];
     %csvwrite(['Concentration_BrainECS_without_source_Dif_Poisson_Seed_', num2str(seed), '.csv'], data_BrainECS);
     %addHeadersToFile(['Concentration_BrainECS_without_source_Dif_Poisson_Seed_', num2str(seed), '.csv'], header);

% Membrane
     %data_Membrane = [t', Cm(end, 1:numt)'];
     %csvwrite(['Concentration_Membrane_without_source_Dif_Poisson_Seed_', num2str(seed), '.csv'], data_Membrane);
     %addHeadersToFile(['Concentration_Membrane_without_source_Dif_Poisson_Seed_', num2str(seed), '.csv'], header);

    % Scaffold
     data_Scaffold = [t', Cs(end, 1:numt)'];
     csvwrite(['Concentration_Scaffold_without_source_Dif_Poisson_thick_Seed_', num2str(seed), '.csv'], data_Scaffold);
      addHeadersToFile(['Concentration_Scaffold_without_source_Dif_Poisson_thick_Seed_', num2str(seed), '.csv'], header);


end


% Function to add headers to a CSV file
function addHeadersToFile(filename, header)
    fid = fopen(filename, 'r+');
    if fid == -1
        error('Could not open the file for writing.');
    end
    line = fgetl(fid);
    data = textscan(line, '%s', 'Delimiter', ',');
    if isequal(data{1}, header)
        fclose(fid);
        return;
    end
    fseek(fid, 0, 'bof');
    fprintf(fid, '%s,%s\n', header{1}, header{2});
    fclose(fid);
end

function x = poissonRandom(lamb)
    x = 0;
    p = exp(-lamb);
    s = p;
    u = rand;
    while u > s
        x = x + 1;
        p = p * lamb / x;
        s = s + p;
    end
end   
