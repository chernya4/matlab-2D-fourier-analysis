close all;
clear;

% Define the path to the Excel file
filePath = './Ag111_360nm_Data.xlsx';

% Read the Excel file into a matrix
data = readmatrix(filePath);

[rows, cols] = size(data);

disp(rows);
disp(cols);

% Generate Energy Axis (electronVolts)
startValue = 6.26;
increment = 0.0121569;
numPoints = 256;
endValue = startValue + (numPoints - 1) * increment;

energyAxisEv = linspace(startValue, endValue, numPoints);

% Generate Time Axis (femptoseconds)
startValue = -73.1844;
increment = 0.034375;
numPoints = 4352;
endValue = startValue + (numPoints - 1) * increment;

timeAxisFs = linspace(startValue, endValue, numPoints);
timeAxisS = timeAxisFs * 1e-15; % Convert femptoseconds to seconds
L = length(timeAxisS);
disp("This is length of time axis" + L);

%% POLARIZATION AXIS FOR FFT %%
% Interferrometric equivalent of frequency axis. Related to frequency by formula: E = h * f. 
% where h is the planks constant. 

% 15 is the Points Per Period Sampling Rate of our Detection Device
% 1240 / 360 is the wavelength that we collected this data under
% L is simply the dimension of the time axis

PA = 15 .* 1240 ./ 360 ./ L .* (0:L-1)

%% POLARIZATION AXIS FOR SHIFTFFT %%

shiftFFT_PA = 15 .* 1240 ./ 360 ./ L .* (-L/2:L/2-1)

% Create a 2D matrix for the fast fourier transform data
fftData = zeros(size(data));


%%%%%%
% Create a new figure for the combined plots
figure;

% Plot the raw data with appropriate axes in the first subplot
ax = subplot(2, 2, 1);
imagesc(timeAxisFs, energyAxisEv, data);
c = colorbar;
c.Label.String = 'Intensity';
xlabel('Time (fs)');
ylabel('Energy (eV)');
title('Raw Data');

% Plot the first row of the data matrix in the second subplot
subplot(2, 2, 2);
plot(timeAxisFs, data(1,:));
xlabel('Time (fs)');
ylabel('Intensity');
title('First Row of Data Matrix');

% Plot the first row of the fast fourier transformed data matrix in the third subplot
subplot(2, 2, 3);
plot(PA, abs(fft(data(1,:))));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('FFT of First Row');

% Plot the first row of the fast fourier transformed and shifted data matrix in the fourth subplot
subplot(2, 2, 4);
plot(shiftFFT_PA, abs(fftshift(fft(data(1,:)))));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('FFT Shift of First Row');

