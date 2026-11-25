close all
clear
clc
format longEng

%%
RNo = 133;
loc = '30 mm (diode)';
deviceID = 'Dev2';
channelID = 'ai0';
fs = 10000;           % Sampling frequency in Hz
numSamples = 50000;   % Number of samples to acquire
duration = numSamples / fs;

% Create DAQ session
dq = daq("ni");
ch = addinput(dq, deviceID, channelID, "Voltage");
ch.Range = [-1 1];

% Configure session
dq.Rate = fs;
% dq.DurationInSeconds = duration;

% Display acquisition start
disp('Starting acquisition...');

% Acquire data
[data, timestamps] = read(dq, seconds(duration), "OutputFormat", "Matrix");

% Display acquisition complete
disp('Acquisition complete.');

% Save data with metadata
outputFile = 'daq_data.mat';
channelInfo = struct( ...
    'Device', deviceID, ...
    'Channel', channelID, ...
    'Rate', fs, ...
    'NumSamples', numSamples ...
);

% Compute RMS after subtracting mean
meanSubtracted = data - mean(data);
rmsValue = sqrt(mean(meanSubtracted.^2));

disp('Writing start.');
filenameD = ['t' num2str(RNo) '.mat'];
save(filenameD, 'data', 'timestamps', 'channelInfo', 'loc', 'fs', 'duration', 'numSamples', 'RNo', 'rmsValue');
disp('Writing end.');

% Display RMS
fprintf('Mean-subtracted RMS voltage: %.4f V\n', rmsValue);

% Plot the acquired data
timeVector = (0:numSamples-1) / fs;  % Time in seconds
figure;
plot(timeVector, data);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Acquired Voltage Signal');
grid on;
filename = ['t' num2str(RNo) '.jpg'];
saveas(gcf, filename);