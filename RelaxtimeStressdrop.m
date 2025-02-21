clc;
clear;

%% Define the input file paths
filePaths = {
    'X:\Research\ADMET\Processing data\09127-2-2025-01-10-12-54-45.csv', ...
    'X:\Research\ADMET\Processing data\09905-2-2025-01-09-16-00-11.csv', ...
    'X:\Research\ADMET\Processing data\08641-re-1-2025-01-07-13-31-06.csv', ...
    'X:\Research\ADMET\Processing data\04271-3-2025-02-10-16-34-14.csv', ...
    'X:\Research\ADMET\Processing data\04062-1-2025-02-19-15-24-54.csv', ...
    'X:\Research\ADMET\Processing data\03899-3-2025-02-10-15-22-02.csv', ...
    'X:\Research\ADMET\Processing data\08363-3-2024-12-20-16-32-03.csv', ...
    'X:\Research\ADMET\Processing data\08570-1-2025-02-19-14-28-38.csv', ...
    'X:\Research\ADMET\Processing data\08219-1-2025-01-09-15-24-04.csv', ...
    'X:\Research\ADMET\Processing data\04974-2-2024-12-23-10-15-31.csv', ...
    'X:\Research\ADMET\Processing data\08183-2025-02-05-15-25-57.csv', ...
    };
% Preallocate cell arrays for parameters
numFiles = numel(filePaths);
Time = cell(1, numFiles);
Load = cell(1, numFiles);
Position = cell(1, numFiles);
C = cell(1, numFiles);
Thickness = cell(1, numFiles);
L = cell(1, numFiles);
Stretch = cell(1, numFiles);
CStress = cell(1, numFiles);
NormalCS = cell(1, numFiles); % Cell array to store NormalCS

% Process each dataset
for i = 1:numFiles
    % Read the dataset
    dataset = xlsread(filePaths{i});
    
    % Extract parameters
    Time{i} = dataset(4:end, 1);
    Load{i} = dataset(4:end, 2);
    Position{i} = dataset(4:end, 3);
    C{i} = dataset(1, 11); % Circumference
    Thickness{i} = dataset(1, 12);
    L{i} = 6; % Long (constant value in this case)

    % Calculate the slope for Load adjustment
    Slope = (Load{i}(end) - Load{i}(1)) / Time{i}(end);
    
    % Adjust the Load values
    Load{i} = Load{i} - Slope * Time{i};
    Load{i} = Load{i} - Load{i}(1);

    
    % Calculate derived quantities
    Stretch{i} = (Position{i} + C{i}) / C{i};
    CStress{i} = Load{i} .* Stretch{i} / L{i} / Thickness{i};
end

%% Parameters
threshold = 5;       % Number of consecutive cells
tolerance = 1e-4;    % Range around one
minDistance = 15;    % Minimum number of indices between the first and second close-to-1 values

% Preallocate cell arrays for results
constantIndices = cell(1, numFiles);
zeroIndices = cell(1, numFiles);
secondZeroIndices = cell(1, numFiles);
LastCycleStretch = cell(1, numFiles);
LastCycleStress = cell(1, numFiles);

% Process each dataset
for i = 1:numFiles
    % Find the constant index and associated indices
    [constantIndex, zeroIndex, secondZeroIndex] = findConstantIndexWithDisplay(Stretch{i}, threshold, tolerance, minDistance);
    
    % Store results
    constantIndices{i} = constantIndex;
    zeroIndices{i} = zeroIndex;
    secondZeroIndices{i} = secondZeroIndex;

    % Check if valid indices were found
    if zeroIndex > 0 && secondZeroIndex > 0
        % Extract LastCycleStretch and LastCycleStress
        % Loading only!!!
        LastCycleStretch{i} = Stretch{i}(secondZeroIndex:(zeroIndex+secondZeroIndex)/2);
        LastCycleStress{i} = CStress{i}(secondZeroIndex:(zeroIndex+secondZeroIndex)/2);
        
        % Display the range of LastCycleStretch
        disp(['Dataset ', num2str(i), ': LastCycleStretch starts at index ', ...
              num2str(secondZeroIndex), ' and ends at index ', num2str(zeroIndex)]);
        
    else
        disp(['Dataset ', num2str(i), ': Valid zeroIndex or secondZeroIndex was not found. Unable to create LastCycleStretch.']);
    end
end

%% Slope fixer
for i = 1:numFiles
    % Calculate the slope for Load adjustment
    Slope = (Load{i}(end) - Load{i}(zeroIndices{i})) / (Time{i}(end)-Time{i}(zeroIndices{i}));
    
    % Adjust the Load values
    Load{i} = Load{i} - Slope * Time{i};
    Load{i} = Load{i} - Load{i}(1);

    
    % Calculate derived quantities
    Stretch{i} = (Position{i} + C{i}) / C{i};
    CStress{i} = Load{i} .* Stretch{i} / L{i} / Thickness{i};
end

%% Calculate RelaxTime for each dataset based on zeroIndex
RelaxTime = cell(1, numFiles); % Preallocate RelaxTime cell array

for i = 1:numFiles
    % Check if zeroIndex is valid for the current dataset
    if zeroIndices{i} > 0
        % Define RelaxTime as the difference between Time and the value at zeroIndex
        RelaxTime{i} = Time{i} - Time{i}(zeroIndices{i});
        
    else
        % If zeroIndex is invalid, display a warning and set RelaxTime to empty
        warning('No valid zeroIndex found for dataset %d. RelaxTime not calculated.', i);
        RelaxTime{i} = [];
    end
end

%% Fix spike for i=5
% Define the spike range
o = 9; %The 08219 one that has a noise spike!
spikeRange = RelaxTime{o} >= 177.987 & RelaxTime{o} <= 186.366;

% Exclude spike data
RelaxTime_clean = RelaxTime{o}(~spikeRange); % x-values without the spike
CS_clean = CStress{o}(~spikeRange);  % y-values without the spike

% Interpolate to fill the missing values
RelaxTime_spike = RelaxTime{o}(spikeRange); % x-values in the spike range
CS_interpolated = interp1(RelaxTime_clean, CS_clean, RelaxTime_spike, 'linear');

% Replace the spike in NormalCS with the interpolated values
CStress{o}(spikeRange) = CS_interpolated;

%Plot the cleaned data
% figure;
% plot(RelaxTime{5}, CStress{5}, 'b', 'LineWidth', 1.5);
% xlabel('Relaxation Time');
% ylabel('Cauchy Stress');
% title('Filtered Data with Spike Replaced by Interpolation');
% grid on;

%% Plot CS{i} against log RelaxTime{i} with RelaxTime >= 0
figure;
hold on; % Allow multiple lines on the same plot
% Use the provided exact colors
colors = {
    '#cd001a',  % Red
    '#1961ae',  % Blue
    '#f2cd00',  % Yellow
    '#61007d',  % Purple
    '#ef6a00',  % Orange
    '#79c300',  % Green
    '#5A5AFF',  % Bright Blue
    '#00BFFF',  % Dark blue
    '#D4A017',  % Gold
    '#008080',  % Teal
    '#404040',  % Dark Gray
    '#FF69B4',  % Hot Pink
    '#A52A2A',  % Brown
    '#00CED1',  % Dark Turquoise
    '#708090'   % Slate Gray
};

% Convert hex colors to RGB values for MATLAB
colors_rgb = cellfun(@(c) sscanf(c(2:end), '%2x%2x%2x', [1 3])/255, colors, 'UniformOutput', false);

% Darken the colors by scaling each RGB component by 0.7
colors_rgb_darker = cellfun(@(c) c * 0.7, colors_rgb, 'UniformOutput', false);
StressDropPerc = cell(1, numFiles);
maxStress = cell(1, numFiles);
FinishStress = cell(1, numFiles);
StartIndices = cell(1,numFiles);

for i = 1:numFiles
    % Check if RelaxTime{i} is valid (not empty)
    if ~isempty(RelaxTime{i})
        % Filter out values where RelaxTime < 0
        validIndices = RelaxTime{i} >= 0;
        
        % Plot only the data where RelaxTime >= 0
        plot(RelaxTime{i}(validIndices), CStress{i}(validIndices)*1e3, 'LineWidth', 1.5, 'Color', colors_rgb_darker{i}, ...
            'DisplayName', sprintf('NormalCS%d', i));
        % Calculate the stress drop percentage
        tolerance = 1; % Define a small tolerance value
        targetTime = 1180;

        
        StartIndices{i} = RelaxTime{i} < 5; % Find indices where RealxTime is less than 5
        maxStress{i} = max(CStress{i}(StartIndices{i})); % Find the max CStress value within those indices
        % This one I don't know what is the problem
        maxStress{o} = 3.60852*1e-3;
        FinishIndex = find(abs(RelaxTime{i} - targetTime) <= tolerance, 1); % Find the first index within the tolerance
        FinishStress{i} = CStress{i}(FinishIndex);
        StressDropPerc{i} = 100*(maxStress{i} - FinishStress{i}) / maxStress{i};
        
    else
        % Skip the plot for this dataset if RelaxTime is empty
        warning('Skipping plot for dataset %d due to invalid RelaxTime.', i);
    end
end

hold off;

% Add labels, legend, and title
legendEntries = {'27 years old', '39 years old', '47 years old', '59 years old', '62 years old', '64 years old', '70 years old', '72 years old', '74 years old', '87 years old'}; % Adjust according to numFiles
set(gca, 'XScale', 'log');
ylim([0, 19.5]);
xlabel('log Relaxation Time');
ylabel('Cauchy Stress');
title('Cauchy Stress vs. Relaxation Time');
legend(legendEntries, 'Location', 'best'); % Use the manually defined legend entries
grid on;

%% Plot stress drop percentage versus age

age = [26, 39, 47, 59, 62, 64, 70, 72, 74, 87, 90];
% Convert StressDropPerc to a numeric array if it's a cell array
StressDropPercNumeric = cell2mat(StressDropPerc); % Converts {a, b, c} to [a, b, c]
% Fit a linear regression line to the data
p = polyfit(age, StressDropPercNumeric, 1); % 1 for a linear fit (first degree polynomial)

% Generate fitted values using the polynomial coefficients
fittedValues = polyval(p, age);

% Plot StressDropPerc against age
figure;
plot(age, StressDropPercNumeric, '.', 'LineStyle', 'none', 'MarkerSize', 25, 'DisplayName', 'Experimental Data');
% Hold the plot to overlay the fitted line
hold on;

% Plot the fitted line as dotted
plot(age, fittedValues, ':r', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit'); % Dotted red line for fitted data

% Add labels and title
xlabel('Age (years)', 'FontSize', 14, 'FontWeight','bold');
ylabel('Stress Drop Percentage', 'FontSize', 14, 'FontWeight','bold');
 
% Set axis limits and hide specific ticks
ax = gca;
ax.FontSize = 14; % Adjust axis number font size
box off;
xlim([20, 100]);
ylim([5, 62]);   % Start y-axis from -1

% Add legend
legend('show', 'Location', 'best', 'Fontsize', 12);

% Set the plot aspect ratio
pbaspect([1 1 1]); % Set the plot aspect ratio

grid off;





%% Helper function to process a dataset
function outputMatrix = processDataset(dataset, L)
    % Extract relevant data
    Time = dataset(4:end, 1);
    Load = dataset(4:end, 2);
    Position = dataset(4:end, 3);
    C = dataset(1, 11); % Circumference
    Thickness = dataset(1, 12);
    
    % Calculate derived quantities
    Stretch = (Position + C) / C;
    CStress = Load .* Stretch / L / Thickness;
    
    % Combine into an output matrix
    outputMatrix = [Time, CStress, Load, Position, Stretch];
end


