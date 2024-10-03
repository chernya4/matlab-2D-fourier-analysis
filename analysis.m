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

% Create a 2D matrix for the fast fourier transform data
fftData = zeros(size(data));

% Plot the raw data with appropriate axes.
figure;
imagesc(timeAxisFs, energyAxisEv, data);
colorbar;

% Plot the first row of the data matrix
figure;
plot(timeAxisFs, data(1,:));

% Plot the first row of the fast fourier transformed data matrix
figure;
plot(abs(fft(data(1,:))));

% Plot the first row of the fast fourier transformed and shifted data matrix
figure;
plot(abs(fftshift(fft(data(1,:)))));









%%%%%%%%%%%%%   3D WORK %%%%%%%%%%%%%%%
% for loop runs through all the rows of the code. 
% argument of the for loop stores the fast fourier transformed row into a row of the fftData matrix.

% for i = 1:size(data, 1)
%     fftData(i,:) = fft(data(i,:));
% end



% figure;
% imagesc(data);
% colorbar;

% figure;
% imagesc(abs(fftData));
% colorbar;

% % % Just a sample that prints the upper left hand corner 10x10 cutout of the data matrix.
% % disp(data(1:10, 1:10))

% Just a Test for a Commit