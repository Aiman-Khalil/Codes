Ds=0.3668
P=0.001; %(μm/s)
d = 85; %(μm) or d = 150, 
% Given values
porosity_values = [0.5, 0.7, 0.9];  % Provide three porosity values
thickness_values = [50, 100, 150, 200]; %(μm)

% Initialize a matrix to store loss values for different porosity values
loss_values = zeros(length(thickness_values), length(porosity_values));

% Calculate loss for each combination of thickness and porosity
for j = 1:length(porosity_values)
    for i = 1:length(thickness_values)
        tortuosity  = (porosity_values(j) / P) * (Ds ./ thickness_values(i));
        n = 1 / (1 + exp(tortuosity * d));

        % Calculate loss in dB
        loss_dB = n * 10 * log10(d);

        % Store the loss in the matrix
        loss_values(i, j) = loss_dB;
    end
end

max_linear_loss = max(abs(loss_values), [], 'all');


% Convert losses to percentages
percentage_loss_values = (abs(loss_values) / max_linear_loss);


% Create file names for each porosity value
file_names = { ...
    'loss_porosity_0.5.csv', 
    'loss_porosity_0.7.csv', 
    'loss_porosity_0.9.csv' 
};

% Write data to CSV files
for j = 1:length(porosity_values)
    file_name = file_names{j};
    
    % Headers
    headers = {'Thickness_um',
               ['Concentration_Loss_percent_for_porosity_', num2str(porosity_values(j))]};
    
    % Combine thickness and loss
    data = [thickness_values', percentage_loss_values(:, j)];
    
    % Write headers
    fid = fopen(file_name, 'w');
    fprintf(fid, '%s,%s\n', headers{:});
    fclose(fid);
    
    % Append numeric data
    dlmwrite(file_name, data, '-append', 'precision', '%.6f', 'delimiter', ',');
end

