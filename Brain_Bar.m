
% Thickness values
d_values= [20, 40, 60, 80];  

d = 95; %(μm)
% Given values

lambda_values = [1.40, 1.67, 1.80];

% 
loss_values = zeros(length(d_values), length(lambda_values));


for j = 1:length(lambda_values)
    for i = 1:length(d_values)
       n = exp(lambda_values(j));

        % Calculate loss in dB
        loss_dB = n * 10 * log10(d_values(i));

        % Store the loss in the matrix
        loss_values(i, j) = loss_dB;
    end
end

% Convert loss to percentage for each porosity value
max_loss = max(loss_values, [], 'all');
loss_percentage = (loss_values / max_loss) * 100;





% Define the file names for each porosity value
file_names = {'tortuosity_1.40_oct.csv', 'tortuosity_1.67_oct.csv', 'tortuosity_1.80_oct.csv'};

% Write data to CSV files
for j = 1:length(lambda_values)
    % Create file name
    file_name = file_names{j};
    
    % Create header
    headers = {'Thickness', ['Concentration_Loss_for_Tortuosity ', num2str(lambda_values(j))]};

    % Combine thickness values and concentration losses for current porosity
    data = [d_values', loss_values(:, j)];

    % Write data to CSV
    % Write headers first
    fid = fopen(file_name, 'w');
    fprintf(fid, '%s,%s\n', headers{:});
    fclose(fid);
    % Append data
    dlmwrite(file_name, data, '-append', 'precision', '%.6f', 'delimiter', ',');
end

