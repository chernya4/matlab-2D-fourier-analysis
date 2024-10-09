close all;
clear;

% Define the path to the Excel file
filePath = './Cu111_806nm_660kc';

% Read the Excel file into a matrix
data = readmatrix(filePath);

[rows, cols] = size(data);

disp(rows);
disp(cols);

% Generate Energy Axis (electronVolts)
startValue = 4.8125;
increment = 0.0157035;
numPoints = 200;
endValue = startValue + (numPoints - 1) * increment;

energyAxisEv = linspace(startValue, endValue, numPoints);

% Generate Time Axis (femptoseconds)
startValue = -224.58;
increment = 0.179234;
numPoints = 1801;
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

PA = 15 * 1240 / 806 / L * (0:L-1);

%% POLARIZATION AXIS FOR SHIFTFFT %%

shiftFFT_PA = 15 .* 1240 ./ 806 ./ L .* (-L/2:L/2-1);


% ###############################################################
% #                                                             #
% #                                                             #
% #                      FIRST ROW DATA                         #
% #                                                             #
% ###############################################################

% Create a new figure for the combined plots
figure;

% Set the title of the entire figure
sgtitle('Filtered Harmonics and FFT of First Row Data');

%%%%%%%%%%%%%%%%%%%%%%% Unfiltered Plots %%%%%%%%%%%%%%%%%%%%%%%

% Plot the first row of the data matrix in the second subplot
subplot(4, 2, 2);
plot(timeAxisFs, data(1,:));
xlabel('Time (fs)');
ylabel('Intensity');
title('First Row of Data Matrix');

% Plot the first row of the fast fourier transformed and shifted data matrix in the fourth subplot
subplot(4, 2, 1);
plot(shiftFFT_PA, abs(fftshift(fft(data(1,:)))));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('FFT Shift of First Row');

%%%%%%%%%%%%%%%%%%%%%%% 0w Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the first row of data to select the zeroth harmonic
filterRow = zeros(1, cols);
zerothHarmonicRange = [-0.1, 0.1]; % Define the range for the zeroth harmonic in eV

% Find the indices corresponding to the zeroth harmonic range
zerothHarmonicIndices = find(shiftFFT_PA >= zerothHarmonicRange(1) & shiftFFT_PA <= zerothHarmonicRange(2));

% Set the filter values to one in the zeroth harmonic range
filterRow(zerothHarmonicIndices) = 1;

% Apply the filter to the FFT data of the first row
filteredFFTDataRow = fftshift(fft(data(1,:))) .* filterRow;

% Perform the inverse FFT to get the filtered time-domain data for the first row
filteredDataRow = ifft(ifftshift(filteredFFTDataRow));

% Plot the filtered first row data in the time domain

subplot(4, 2, 3);
plot(shiftFFT_PA, abs(filteredFFTDataRow));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('Filtered FFT Shift of First Row (Zeroth Harmonic)');

subplot(4, 2, 4);
plot(timeAxisFs, abs(filteredDataRow));
xlabel('Time (fs)');
ylabel('Intensity');
title('Filtered First Row Data (Zeroth Harmonic)');

%%%%%%%%%%%%%%%%%%%%%%% 1w Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the first row of data to select the first harmonic
filterRow = zeros(1, cols);
firstHarmonicRange1 = [1.3, 1.6]; % Define the first range for the first harmonic in eV
firstHarmonicRange2 = [-1.7, -1.3]; % Define the second range for the first harmonic in eV

% Find the indices corresponding to the first harmonic ranges
firstHarmonicIndices1 = find(shiftFFT_PA >= firstHarmonicRange1(1) & shiftFFT_PA <= firstHarmonicRange1(2));
firstHarmonicIndices2 = find(shiftFFT_PA >= firstHarmonicRange2(1) & shiftFFT_PA <= firstHarmonicRange2(2));

% Set the filter values to one in the first harmonic ranges
filterRow(firstHarmonicIndices1) = 1;
filterRow(firstHarmonicIndices2) = 1;

% Apply the filter to the FFT data of the first row
filteredFFTDataRow = fftshift(fft(data(1,:))) .* filterRow;

% Perform the inverse FFT to get the filtered time-domain data for the first row
filteredDataRow = ifft(ifftshift(filteredFFTDataRow));

subplot(4, 2, 5);
plot(shiftFFT_PA, abs(filteredFFTDataRow));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('Filtered FFT Shift of First Row (First Harmonic)');

subplot(4, 2, 6);
plot(timeAxisFs, abs(filteredDataRow));
xlabel('Time (fs)');
ylabel('Intensity');
title('Filtered First Row Data (First Harmonic)');

%%%%%%%%%%%%%%%%%%%%%%% 2w Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the first row of data to select the second harmonic
filterRow = zeros(1, cols);
secondHarmonicRange1 = [2.8, 3.2]; % Define the first range for the second harmonic in eV
secondHarmonicRange2 = [-3.2, -2.8]; % Define the second range for the second harmonic in eV

% Find the indices corresponding to the second harmonic ranges
secondHarmonicIndices1 = find(shiftFFT_PA >= secondHarmonicRange1(1) & shiftFFT_PA <= secondHarmonicRange1(2));
secondHarmonicIndices2 = find(shiftFFT_PA >= secondHarmonicRange2(1) & shiftFFT_PA <= secondHarmonicRange2(2));

% Set the filter values to one in the second harmonic ranges
filterRow(secondHarmonicIndices1) = 1;
filterRow(secondHarmonicIndices2) = 1;

% Apply the filter to the FFT data of the first row
filteredFFTDataRow = fftshift(fft(data(1,:))) .* filterRow;

% Perform the inverse FFT to get the filtered time-domain data for the first row
filteredDataRow = ifft(ifftshift(filteredFFTDataRow));

% Plot the filtered first row data in the time domain

subplot(4, 2, 7);
plot(shiftFFT_PA, abs(filteredFFTDataRow));
xlabel('Polarization (eV)');
ylabel('Magnitude');
title('Filtered FFT Shift of First Row (Second Harmonic)');

subplot(4, 2, 8);
plot(timeAxisFs, abs(filteredDataRow));
xlabel('Time (fs)');
ylabel('Intensity');
title('Filtered First Row Data (Second Harmonic)');

% ###############################################################
% #                                                             #
% #                                                             #
% #                          2D DATA                            #
% #                                                             #
% ###############################################################



% Perform FFT on each row of the data matrix
fftData = fft(data, [], 2);

% Perform FFT shift on each row of the data matrix
shiftFFTData = fftshift(fft(data, [], 2), 2);

% Create a new figure for the combined plots
figure;

% Set the title of the entire figure
sgtitle('Filtered Harmonics and FFT of Data');


%%%%%%%%%%%%%%%%%%%%%%% 2D Shifted FFT Data %%%%%%%%%%%%%%%%%%%%%%%

% Plot the 2D shifted FFT data with appropriate axes in the second subplot
subplot(4, 2, 1);
imagesc(shiftFFT_PA, energyAxisEv, abs(shiftFFTData));
c = colorbar;
xlabel('Polarization (eV)');
ylabel('Energy (eV)');
title('2D Shifted FFT of Data');

%%%%%%%%%%%%%%%%%%%%%%% Raw Data %%%%%%%%%%%%%%%%%%%%%%%

% Plot the raw data with appropriate axes in the first subplot
subplot(4, 2, 2);
imagesc(timeAxisFs, energyAxisEv, data);
c = colorbar;
xlabel('Time (fs)');
ylabel('Energy (eV)');
title('Raw Data');

%%%%%%%%%%%%%%%%%%%%%%% Zeroth Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the entire data matrix to select the zeroth harmonic
filterMatrix = zeros(size(shiftFFTData));
zerothHarmonicRange = [-0.1, 0.1]; % Define the range for the zeroth harmonic in eV

% Find the indices corresponding to the zeroth harmonic range
zerothHarmonicIndices = find(shiftFFT_PA >= zerothHarmonicRange(1) & shiftFFT_PA <= zerothHarmonicRange(2));

% Set the filter values to one in the zeroth harmonic range for all rows
filterMatrix(:, zerothHarmonicIndices) = 1;

% Apply the filter to the shifted FFT data matrix
filteredShiftFFTData = shiftFFTData .* filterMatrix;

% Plot the filtered shift FFT data with appropriate axes
subplot(4, 2, 3);
imagesc(shiftFFT_PA, energyAxisEv, abs(filteredShiftFFTData));
c = colorbar;
xlabel('Polarization (eV)');
ylabel('Energy (eV)');
title('Filtered Shift FFT Data (Zeroth Harmonic)');

% Perform the inverse FFT to get the filtered time-domain data matrix for the zeroth harmonic
filteredDataMatrix = ifft(ifftshift(filteredShiftFFTData, 2), [], 2);

% Plot the filtered data matrix with appropriate axes
subplot(4, 2, 4);
imagesc(timeAxisFs, energyAxisEv, abs(filteredDataMatrix));
c = colorbar;
xlabel('Time (fs)');
ylabel('Energy (eV)');
title('Filtered Data Matrix (Zeroth Harmonic)');

%%%%%%%%%%%%%%%%%%%%%%% First Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the entire data matrix to select the first harmonic
filterMatrix = zeros(size(shiftFFTData));
firstHarmonicRange1 = [1.3, 1.6]; % Define the first range for the first harmonic in eV
firstHarmonicRange2 = [-1.7, -1.3]; % Define the second range for the first harmonic in eV

% Find the indices corresponding to the first harmonic ranges
firstHarmonicIndices1 = find(shiftFFT_PA >= firstHarmonicRange1(1) & shiftFFT_PA <= firstHarmonicRange1(2));
firstHarmonicIndices2 = find(shiftFFT_PA >= firstHarmonicRange2(1) & shiftFFT_PA <= firstHarmonicRange2(2));

% Set the filter values to one in the first harmonic ranges for all rows
filterMatrix(:, firstHarmonicIndices1) = 1;
filterMatrix(:, firstHarmonicIndices2) = 1;

% Apply the filter to the shifted FFT data matrix
filteredShiftFFTData = shiftFFTData .* filterMatrix;

% Plot the filtered shift FFT data with appropriate axes
subplot(4, 2, 5);
imagesc(shiftFFT_PA, energyAxisEv, abs(filteredShiftFFTData));
c = colorbar;
xlabel('Polarization (eV)');
ylabel('Energy (eV)');
title('Filtered Shift FFT Data (First Harmonic)');

% Perform the inverse FFT to get the filtered time-domain data matrix for the first harmonic
filteredDataMatrix = ifft(ifftshift(filteredShiftFFTData, 2), [], 2);

% Plot the filtered data matrix with appropriate axes
subplot(4, 2, 6);
imagesc(timeAxisFs, energyAxisEv, abs(filteredDataMatrix));
c = colorbar;
xlabel('Time (fs)');
ylabel('Energy (eV)');
title('Filtered Data Matrix (First Harmonic)');

%%%%%%%%%%%%%%%%%%%%%%% Second Harmonic %%%%%%%%%%%%%%%%%%%%%%%

% Construct a square filter for the entire data matrix to select the second harmonic
filterMatrix = zeros(size(shiftFFTData));
secondHarmonicRange1 = [2.8, 3.2]; % Define the first range for the second harmonic in eV
secondHarmonicRange2 = [-3.2, -2.8]; % Define the second range for the second harmonic in eV

% Find the indices corresponding to the second harmonic ranges
secondHarmonicIndices1 = find(shiftFFT_PA >= secondHarmonicRange1(1) & shiftFFT_PA <= secondHarmonicRange1(2));
secondHarmonicIndices2 = find(shiftFFT_PA >= secondHarmonicRange2(1) & shiftFFT_PA <= secondHarmonicRange2(2));

% Set the filter values to one in the second harmonic ranges for all rows
filterMatrix(:, secondHarmonicIndices1) = 1;
filterMatrix(:, secondHarmonicIndices2) = 1;

% Apply the filter to the shifted FFT data matrix
filteredShiftFFTData = shiftFFTData .* filterMatrix;

% Plot the filtered shift FFT data with appropriate axes
subplot(4, 2, 7);
imagesc(shiftFFT_PA, energyAxisEv, abs(filteredShiftFFTData));
c = colorbar;
xlabel('Polarization (eV)');
ylabel('Energy (eV)');
title('Filtered Shift FFT Data (Second Harmonic)');

% Perform the inverse FFT to get the filtered time-domain data matrix for the second harmonic
filteredDataMatrix = ifft(ifftshift(filteredShiftFFTData, 2), [], 2);

% Plot the filtered data matrix with appropriate axes
subplot(4, 2, 8);
imagesc(timeAxisFs, energyAxisEv, abs(filteredDataMatrix));
c = colorbar;
xlabel('Time (fs)');
ylabel('Energy (eV)');
title('Filtered Data Matrix (Second Harmonic)');






% % Create a 2D matrix for the fast fourier transform data
% nfft = nextpow2(L);
% fftData = zeros(size(data)); %%% NEED TO FIX THIS 


