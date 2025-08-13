% Define the fixed lambda value
lambda_val = 0:0.01:2; % Lambda values
d = 95; % Distance in micrometers
d0=95;
% Initialize loss values for each lambda value
loss_values = zeros(size(lambda_val)); 

% Calculate  the loss for each lambda value
for j = 1:length(lambda_val)
    lamb = lambda_val(j);
    n = exp(lamb); % Calculate n based on lambda
    loss= n * 10 * log10(d);
    loss_values(j) = loss; 
end


% Store time and concentration loss percentage data in a table
data = table(lambda_val', loss_values', 'VariableNames', {'Tortuosity', 'Concentration_Loss_dB'});

% Write data to CSV file
writetable(data, 'Tortuosity_concentration_loss_updated_14OCT.csv');

