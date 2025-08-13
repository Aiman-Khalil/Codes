clc;
clear;

% Given parameters
transmitted_molecules = 10*log10(700);
distance_range = 1:330;

% Calculate LA losses
tortuosity_LA = 1.67;
distances_LA = 1:95;

loss_values_Brain = zeros(size(distances_LA));

for i = 1:length(distances_LA)
    d = distances_LA(i);
    n = exp(tortuosity_LA);
    loss = n * 10 * log10(d);
    loss_values_Brain(i) = loss;
end

% Calculate LM losses
Ds = 0.3668;
P = 0.001;
porosity_values = 0.9;
thickness=85;
d_values_LM = 96:180;
loss_values_membrane = zeros(size(d_values_LM));

for i = 1:length(d_values_LM)
    tortuosity_LM = (porosity_values/P) * (Ds/thickness);
    n = 1 / (1 + exp(tortuosity_LM * d_values_LM(i)));
    loss_dB_LM = n * 10 * log10(d_values_LM(i));
    loss_values_membrane(i) = loss_dB_LM;
end
disp('Loss Values - Membrane:');
disp(loss_values_membrane);

% Calculate LS losses
Vmax = 2.84e-6;
d_values_LS = 181:330;

k_value = 3.8e14;
N0=4285;

loss_values_scaffold = zeros(size(d_values_LS));

for i = 1:length(d_values_LS)
    Absorption = N0*(Vmax / k_value);
    n=exp(Absorption);
    loss_dB_LS = n * 10 * log10(d_values_LS(i));
    loss_values_scaffold(i) = loss_dB_LS;
end
disp('Loss Values - Scaffold:');
disp(loss_values_scaffold);



common_distances = 1:330; % Choose a common set of distances

% Interpolate values for Brain, Membrane, and Scaffold losses
interp_loss_values_Brain = interp1(distances_LA, loss_values_Brain, common_distances, 'linear', 'extrap');
interp_loss_values_membrane = interp1(d_values_LM, loss_values_membrane, common_distances, 'linear', 'extrap');
interp_loss_values_scaffold = interp1(d_values_LS, loss_values_scaffold, common_distances, 'linear', 'extrap');

% Calculate total losses using interpolated values
total_losses = interp_loss_values_Brain + interp_loss_values_membrane + interp_loss_values_scaffold;

% Calculate received molecules for the adjusted distance range
received_molecules = transmitted_molecules - interp_loss_values_Brain-interp_loss_values_membrane-interp_loss_values_scaffold;

% Create CSV files for header Distance and Received Molecules Percentage

% Write headers to CSV file
header = {'Distance', 'received_molecules'};
writecell(header, 'received_molecules_12JUNE.csv');

% Write distance and received molecules percentage data to CSV file
data = [common_distances' received_molecules'];
dlmwrite('received_molecules_12JUNE.csv', data, '-append');
