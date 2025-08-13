% Given values
thickness_values = 85; % Thickness in micrometers
Ds = 0.3668; % Diffusion coefficient (μm^2/s)
porosity_values = 0.9; % Porosity value
d = 85; % Distance in micrometers
d0=85;
% Permeability values
P = (0.01:0.001:0.5); % Permeability (μm/s)

% Initialize a matrix to store loss values for different permeability values
loss_values = zeros(size(P));

% Calculate loss
for i = 1:length(P)
    tortuosity = (porosity_values / P(i)) * (Ds / thickness_values);
    % Calculate the loss using a sigmoidal function
    n = 1 / (1 + exp(-tortuosity * d));
    loss_dB = n * 10 * log10(d); % Calculate loss in dB
    
    loss_values(i) = loss_dB;
    
end


% Store time and concentration loss percentage data in a table
data = table(P', loss_values', 'VariableNames', {'Permeability', 'Concentration_Loss_dB'});

% Write data to CSV file
writetable(data, 'Membrane_Permeability_concentration_loss_14OCT.csv');


