N0 = 4285;
k_values=3.8e14;

Vmax=(0.8:0.01:2)*10^-6;
d = 150;

loss_values = zeros(size(Vmax));

for i = 1:length(Vmax)
    Absorption =N0*(Vmax(i)) / k_values;  
    n = (exp(Absorption));
    loss_dB = n * 10 * log10(d);
    loss_values(i) = loss_dB;
end

max_loss = max(loss_values);
loss_percentage = (loss_values / max_loss) * 100;


data = table(Vmax', loss_values','VariableNames', {'Uptake_Rate', 'Concentration_Loss_dB'});

% Write table to CSV
writetable(data, 'Scaffold_Uptake_Rate_loss_percent_15OCT.csv');


