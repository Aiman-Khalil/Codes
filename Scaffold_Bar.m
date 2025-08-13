% Given Vmax values
Vmax_values = [0.84e-6 , 1.84e-6 , 2.84e-6];

% k_values
N0= [1000, 2000, 3000, 4000]
k_values = 3.8e14; % Half saturation constants
d = 150; % Constant value

% Initialize a matrix to store loss values for different Vmax values
loss_values = zeros(length(N0), length(Vmax_values));

% Calculate loss for each combination of k and Vmax
for j = 1:length(Vmax_values)
    Vmax = Vmax_values(j);
    for i = 1:length(N0)
        Absorption = N0(i) * Vmax / k_values;
        n=exp(Absorption);
        loss_dB = n * 10 * log10(d);
        loss_values(i, j) = loss_dB;
    end
end

% Convert loss to percentage for each Vmax value
max_loss = max(loss_values, [], 'all');
loss_percentage = (loss_values / max_loss) * 100;




% Define the file names for each porosity value
file_names = {'uptake_Rate_0.84e-6.csv', 'uptake_Rate_1.84e-6.csv', 'uptake_Rate_2.84e-6.csv'};

% Write data to CSV files
for j = 1:length(Vmax_values)
    % Create file name
    file_name = file_names{j};
    
    % Create header
    headers = {'N0', ['Concentration_Loss_for_uptake_rate', num2str(Vmax_values(j))]};

    % Combine thickness values and concentration losses for current porosity
    data = [N0', loss_values(:, j)];

    % Write data to CSV
    % Write headers first
    fid = fopen(file_name, 'w');
    fprintf(fid, '%s,%s\n', headers{:});
    fclose(fid);
    % Append data
    dlmwrite(file_name, data, '-append', 'precision', '%.6f', 'delimiter', ',');
end
