% Define the path to the Excel file
filePath = './Ag111_360nm_Data.xlsx';

% Read the Excel file into a matrix
data = readmatrix(filePath);

% Create a 2D matrix for the fast fourier transform data
fftData = zeros(size(data));

figure;
plot(data(1,:));

% for loop runs through all the rows of the code. 
% argument of the for loop stores the fast fourier transformed row into a row of the fftData matrix.

for i = 1:size(data, 1)
    fftData(i,:) = fft(data(i,:));
end

figure;
plot(abs(fft(data(1,:))));

% figure;
% imagesc(data);
% colorbar;

% figure;
% imagesc(abs(fftData));
% colorbar;

% % % Just a sample that prints the upper left hand corner 10x10 cutout of the data matrix.
% % disp(data(1:10, 1:10))

% Just a Test for a Commit
