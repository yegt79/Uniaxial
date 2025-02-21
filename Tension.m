clc;
clear;

% Define the input file paths
filePaths = {
    'X:\Research\ADMET\Processing data\08641-re-1-2025-01-07-13-31-06.csv', ...
    'X:\Research\ADMET\Processing data\04974-2-2024-12-23-10-15-31.csv', ...
    'X:\Research\ADMET\Processing data\08219-1-2025-01-09-15-24-04.csv', ...
    'X:\Research\ADMET\Processing data\09905-2-2025-01-09-16-00-11.csv', ...
    'X:\Research\ADMET\Processing data\09127-2-2025-01-10-12-54-45.csv', ...
    'X:\Research\ADMET\Processing data\08363-3-2024-12-20-16-32-03.csv'};

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

% Parameters
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
        LastCycleStretch{i} = Stretch{i}(secondZeroIndex:zeroIndex);
        LastCycleStress{i} = CStress{i}(secondZeroIndex:zeroIndex);
        
        % Display the range of LastCycleStretch
        disp(['Dataset ', num2str(i), ': LastCycleStretch starts at index ', ...
              num2str(secondZeroIndex), ' and ends at index ', num2str(zeroIndex)]);
        
    else
        disp(['Dataset ', num2str(i), ': Valid zeroIndex or secondZeroIndex was not found. Unable to create LastCycleStretch.']);
    end
end

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
% Calculate RelaxTime for each dataset based on zeroIndex
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

% Define NormalCS{i} after everything else has been defined
NormalCS = cell(1, numFiles); % Preallocate NormalCS cell array

for i = 1:numFiles
    if zeroIndices{i} > 0
        % Find max of CStress after zeroIndex
        maxCStress = max(CStress{i}(zeroIndices{i}:end));
        NormalCS{i} = CStress{i} / maxCStress; % Normalize after zeroIndex
    else
        NormalCS{i} = []; % Empty if zeroIndex is invalid
    end
end

% Plot all NormalCS{i} against RelaxTime{i} with RelaxTime >= 0
figure;
hold on; % Allow multiple lines on the same plot
colors = lines(numFiles); % Generate distinct colors for each line

for i = 1:numFiles
    % Check if RelaxTime{i} is valid (not empty)
    if ~isempty(RelaxTime{i})
        % Filter out values where RelaxTime < 0
        validIndices = RelaxTime{i} >= 0;
        
        % Plot only the data where RelaxTime >= 0
        plot(RelaxTime{i}(validIndices), NormalCS{i}(validIndices), 'LineWidth', 1.5, 'Color', colors(i, :), ...
            'DisplayName', sprintf('NormalCS%d', i));
    else
        % Skip the plot for this dataset if RelaxTime is empty
        warning('Skipping plot for dataset %d due to invalid RelaxTime.', i);
    end
end

hold off;

% Add labels, legend, and title
xlabel('Relaxation Time (s)');
ylabel('Normal Cauchy Stress');
title('Normal Cauchy Stress vs. Relaxation Time');
legend('show'); % Show legend based on DisplayName properties
grid on;

% Set y-axis limits from 0 to 1.2
ylim([0, 1.2]);

% Plot all CStress{i} against Stretch{i} 
figure;
hold on; % Allow multiple lines on the same plot
colors = lines(numFiles); % Generate distinct colors for each line

for i = 1:numFiles
        plot(LastCycleStretch{i}, (LastCycleStress{i}-LastCycleStress{i}(1))*Thickness{i}*1e3, 'LineWidth', 1.5, 'Color', colors(i, :), ...
            'DisplayName', sprintf('NormalCS%d', i));
end

hold off;

% Add labels, legend, and title
xlabel('Stretch');
ylabel('Tension (N/m)');
title('Tension vs. Stretch');
legend('show'); % Show legend based on DisplayName properties
grid on;





% Helper function to process a dataset
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


