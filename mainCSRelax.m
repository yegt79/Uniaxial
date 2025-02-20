clc;
clear;

%% Define the input file paths
filePaths = {
    'X:\Research\ADMET\Processing data\09127-2-2025-01-10-12-54-45.csv', ...
    'X:\Research\ADMET\Processing data\09905-2-2025-01-09-16-00-11.csv', ...
    'X:\Research\ADMET\Processing data\08641-re-1-2025-01-07-13-31-06.csv', ...
    'X:\Research\ADMET\Processing data\04271-3-2025-02-10-16-34-14.csv', ...
    'X:\Research\ADMET\Processing data\03899-3-2025-02-10-15-22-02.csv', ...
    'X:\Research\ADMET\Processing data\08363-3-2024-12-20-16-32-03.csv', ...
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

%% Define NormalCS{i} after everything else has been defined
NormalCS = cell(1, numFiles); % Preallocate NormalCS cell array
maxCStress = cell(1, numFiles);
o = 7; %The 08219 one that has a noise spike!
for i = 1:numFiles
    if zeroIndices{i} > 0
        % Find max of CStress after zeroIndex
        if i == o
            % This is a problem
            maxCStress{i} = 3.60852*1e-3;
            NormalCS{i} = CStress{i} / maxCStress{i}; % Normalize after zeroIndex
        else
            maxCStress{i} = max(CStress{i}(zeroIndices{i}:end));
            NormalCS{i} = CStress{i} / maxCStress{i}; % Normalize after zeroIndex
        end
    else
        NormalCS{i} = []; % Empty if zeroIndex is invalid
    end
end

%% Fix spike for 
% Define the spike range
  spikeRange = RelaxTime{o} >= 177.987 & RelaxTime{o} <= 186.366;
% 
% % Exclude spike data
  RelaxTime_clean = RelaxTime{o}(~spikeRange); % x-values without the spike
  NormalCS_clean = NormalCS{o}(~spikeRange);  % y-values without the spike
% 
% % Interpolate to fill the missing values
  RelaxTime_spike = RelaxTime{o}(spikeRange); % x-values in the spike range
  NormalCS_interpolated = interp1(RelaxTime_clean, NormalCS_clean, RelaxTime_spike, 'linear');
% 
% % Replace the spike in NormalCS with the interpolated values
  NormalCS{o}(spikeRange) = NormalCS_interpolated;

% Plot the cleaned data
 figure;
 plot(RelaxTime{o}, NormalCS{o}, 'b', 'LineWidth', 1.5);
 xlabel('Relaxation Time');
 ylabel('Normal Cauchy Stress');
 title('Filtered Data with Spike Replaced by Interpolation');
 grid on;

%% Plot all NormalCS{i} against RelaxTime{i} with RelaxTime >= 0
legendEntries = {'26 yrs M', '39 yrs F', '47 yrs M', '59 yrs M', '62 yrs F', '70 yrs M', '74 yrs F', '87 yrs M', '90 yrs F'}; % Adjust according to numFiles
figure;
hold on; % Allow multiple lines on the same plot

% Use the provided exact colors
colors = {
    '#cd001a',   % Teal
    '#1961ae',   % Blue
    '#f2cd00',   % Light blue
    '#61007d',   % Yellow
    '#ef6a00',   % Orange
    '#79c300',   % Purple
    '#5A5AFF',  %
    '#4B0082',  %
    '#D4A017'
};  

% Convert hex colors to RGB values for MATLAB
colors_rgb = cellfun(@(c) sscanf(c(2:end), '%2x%2x%2x', [1 3])/255, colors, 'UniformOutput', false);

% Darken the colors by scaling each RGB component by 0.7
colors_rgb_darker = cellfun(@(c) c * 0.7, colors_rgb, 'UniformOutput', false);

markers = {'.', '.', '.', '.', '.', '.', '.', '.', '.', '.'}; % Define a set of markers
markerSize = 5; % Increase marker size

for i = 1:numFiles
    % Check if RelaxTime{i} is valid (not empty)
    if ~isempty(RelaxTime{i})
        % Filter out values where RelaxTime < 0
        validIndices = RelaxTime{i} >= 0;
        
        % Cycle through markers if there are more curves than defined markers
        markerIndex = mod(i-1, length(markers)) + 1;
        
        % Plot only the data where RelaxTime >= 0
        plot(RelaxTime{i}(validIndices), NormalCS{i}(validIndices), ...
            'LineStyle', 'none', ...
            'Color', colors_rgb_darker{i}, ...
            'Marker', markers{markerIndex}, ...
            'MarkerSize', markerSize); % Adjusted marker size
    else
        % Skip the plot for this dataset if RelaxTime is empty
        warning('Skipping plot for dataset %d due to invalid RelaxTime.', i);
    end
end

hold off;

% Add labels, legend, and adjust fonts
xlabel('Relaxation Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Cauchy Stress', 'FontSize', 14, 'FontWeight', 'bold');
legend(legendEntries, 'Location', 'best', 'FontSize', 14); % Adjusted legend font size

% Set axis limits and hide specific ticks
ax = gca;
ax.FontSize = 14; % Adjust axis number font size
xlim([-50, 1300]); % Start x-axis from 0
ylim([-0.05, 1.1]);       % Set y-axis limits from 0 to 1.2

% Set the plot aspect ratio
pbaspect([1 1 1]); % Set the plot aspect ratio
grid off;


%% Plot all LastCycleStress{i} against LastCycleStretch{i}
%legendEntries = {'26 yrs M', '39 yrs F', '47 yrs M', '70 yrs M', '74 yrs F', '87 yrs M', '90 yrs F'}; % Adjust according to numFiles
figure;
hold on; % Allow multiple lines on the same plot

% Use the provided exact colors
% colors = {
%     '#cd001a',   % Teal
%     '#1961ae',   % Blue
%     '#f2cd00',   % Light blue
%     '#61007d',   % Yellow
%     '#ef6a00',   % Orange
%     '#79c300',    % Purple
%     '#cd001a'   % Teal
% };  

% Convert hex colors to RGB values for MATLAB
colors_rgb = cellfun(@(c) sscanf(c(2:end), '%2x%2x%2x', [1 3])/255, colors, 'UniformOutput', false);

% Darken the colors by scaling each RGB component by 0.7
colors_rgb_darker = cellfun(@(c) c * 0.7, colors_rgb, 'UniformOutput', false);

markers = {'*', '*', '*', '*', '*', '*', '*', '*', '*', '*'}; % Define a set of markers
markerSize = 8; % Increase marker size

downsampleFactor = 20; % Set downsampling factor (adjust as needed)

for i = 1:numFiles
    % Cycle through markers if there are more curves than defined markers
    markerIndex = mod(i-1, length(markers)) + 1;
    
    % Downsample the x and y data
    downsampledStretch = LastCycleStretch{i}(1:downsampleFactor:end); 
    downsampledStress = (LastCycleStress{i}(1:downsampleFactor:end) - LastCycleStress{i}(1)) * 1e3;
    
    % Plot the downsampled data
    plot(downsampledStretch, downsampledStress, ...
        'LineStyle', 'none', ...
        'Color', colors_rgb{i}, ...
        'Marker', markers{markerIndex}, ...
        'MarkerSize', markerSize); % Adjusted marker size
end

hold off;

% Add labels, legend, and adjust fonts
xlabel('Stretch', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Cauchy Stress (kPa)', 'FontSize', 14, 'FontWeight', 'bold');
legend(legendEntries, 'Location', 'best', 'FontSize', 14); % Adjusted legend font size

% Set axis limits and hide specific ticks
ax = gca;
ax.FontSize = 14; % Adjust axis number font size
xlim([0.99, 1.11]); % Start x-axis from 0.99
ylim([-1, 20]);   % Start y-axis from -1


% Remove the starting values (0.99 and -1) from the tick labels
%xticks(ax.XTick(ax.XTick > 0.99)); % Show only ticks greater than 0.99
%yticks(ax.YTick(ax.YTick > -1));   % Show only ticks greater than -1

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


