% Clear the workspace
clc;
clear;
% Set the random number generator seed for reproducibility
seed = 42; % Replace 42 with the desired seed value
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




for seed = 1:num_seeds
% Concentration matrix
C = zeros(num_points_Re, numt);   % Initialize everything to zero
Cm = zeros(num_points_Rm, numt);  % Initialize everything to zero
Cs = zeros(num_points_Rs, numt);  % Initialize everything to zero


mu = 0  ;
sigma = 0.05;


    % Now, you can distribute these random_release_count molecules spatially
    
        lamb = N;  % Mean number of molecules (lambda) is set to the initial release rate
    random_release_count = poissonRandom(lamb);

     C(:, 1) = random_release_count * exp(-(r_Re).^2 / (2 * sigma^2)) / sqrt(2 * pi * sigma^2);

     


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
            C(i, j+1) = C(i, j) + ((De/lambda^2) * dt / dr_Re^2) * (C(i+1, j) - 2 * C(i, j) + C(i-1, j)) + source_term(i) * dt;
        end
        C(end, j) = C(end - 1, j); % Set the boundary value
  end



% Set the boundary value and initialize concnetration for Rm
Cm(1, 1) = C(end, numt);
Cm(1, j+1) =Pm*C(1, j+1);

% Iterating the diffusion equation in the membrane (Rm)
for j = 1:numt
    for i = 2:num_points_Rm-1
        d2Cm_dr2 = (Cm(i + 1, j) - 2 * Cm(i, j) + Cm(i - 1, j)) / dr_Rm^2;
        dCm_dr = (Cm(i + 1, j) - Cm(i - 1, j)) / (2 * dr_Rm);
        Cm(i, j+1) = Cm(i, j) + D_m * (d2Cm_dr2 + (1 / r_Rm(i)) * dCm_dr) * dt;
    end
    Cm(end, j) = Cm(end - 1, j); % Set the boundary value
end


P_values = [0.0001, 0.001, 0.01]; % Permeability values at the interface between membrane and scaffold (Rm)
% Initialize arrays to store the results for each permeability value
results = cell(length(P_values), 1);

for p = 1:length(P_values)
    P = P_values(p); % Set the current permeability value
% Set the boundary value and initialize concnetration for Rm
Cs(1, 1) = Cm(end, numt);
Cs(1, j+1) = P*Cm(1, j+1);

for j = 1:numt
    for i = 2:num_points_Rs-1
        d2Cs_dr2 = (Cs(i + 1, j) - 2 * Cs(i, j) + Cs(i - 1,  j)) / dr_Rs^2;
        dCs_dr = (Cs(i + 1,  j) - Cs(i - 1,  j)) / (2 * dr_Rs);

        % Reaction term
        
        reaction_term = Vmax* N0 * ((Cs(i, j) / (k + Cs(i, j))));
    
      

         % Update Cs_Rs using the concentration from the previous time step
        Cs(i, j+1) = Cs(i, j) + Ds *dt* (d2Cs_dr2 + (1/r_Rs(i)) * dCs_dr) - reaction_term *dt;
    end
    Cs(end, j) = Cs(end - 1, j); % Set the boundary value
   
 end
results{p} = Cs(end, 1:numt);
end


% Define the base directory where you want to create subfolders
baseDirectory = '/Users/walton_user/Documents/MATLAB/Scaffold_Permeability_Diffusion'; % Update this path

% Create subfolders for each permeability value
P_values = [0.0001, 0.001, 0.01]; % Permeability values
for p = 1:length(P_values)
    folderName = sprintf('permeability_%d', p);
    folderPath = fullfile(baseDirectory, folderName);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
end
    
    for p = 1:length(P_values)
        % Create a matrix containing your data
        data = [t', results{p}'];

        % Create a cell array for headers
        header = {'Time', ['Concentration_permeability_', num2str(p)]};

        % Define the folder for this permeability value
        folderName = sprintf('permeability_%d', p);
        folderPath = fullfile(baseDirectory, folderName);

        % Define the filename within the folder
        filename = sprintf('concentration_permeability_seed%d.csv', seed);

        % Create the full file path
        filePath = fullfile(folderPath, filename);

        % Write the table to a CSV file within the appropriate folder
        table_data = array2table(data, 'VariableNames', header);
        writetable(table_data, filePath);
    end
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

