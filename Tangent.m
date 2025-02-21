clc;
clear;

% Define the input file paths
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

%%Tangent
T_Strain = cell(1,numFiles);
T_CS = cell(1,numFiles);
tangentStiffness = cell(1,numFiles);
for i = 1:numFiles
    % Check if RelaxTime{i} is valid (not empty)
    if ~isempty(RelaxTime{i})
        % Filter out values where RelaxTime < 0
        validIndices = RelaxTime{i} >= 0;
        
        % Plot only the data where RelaxTime >= 0
        T_Strain{i} = LastCycleStretch{i}-1;
        T_CS{i} = (LastCycleStress{i}-LastCycleStress{i}(1))*1e3;
    else
        % Skip the plot for this dataset if RelaxTime is empty
        warning('Skipping plot for dataset %d due to invalid RelaxTime.', i);
    end
end

% Preallocate an array to store the average tangent stiffness for each dataset
tangentStiffnessAtLast = zeros(numFiles, 1);

for i = 1:numFiles
    % Check if RelaxTime{i} is valid (not empty)
    if ~isempty(RelaxTime{i})
        % Filter out values where RelaxTime < 0
        validIndices = RelaxTime{i} >= 0;
        
        % Fit a polynomial to Cauchy Stress vs. Strain
        degree = 5; % Choose degree of polynomial (adjust as needed)
        p = polyfit(T_Strain{i}, T_CS{i}, degree); % Polynomial fit
        
        % Compute the derivative of the polynomial (tangent stiffness)
        dp = polyder(p); % Derivative of the polynomial
        tangentStiffness{i} = polyval(dp, T_Strain{i}); % Evaluate the derivative at strain values

        % Get the last value of the tangent stiffness (at the last strain value)
        tangentStiffnessAtLast(i) = max(tangentStiffness{i}); % Store the last tangent stiffness value

        % Plot the original data (Cauchy Stress vs. Strain)
        figure;
        hold on;
        plot(T_Strain{i}, T_CS{i}, '.', 'DisplayName', 'Original Data'); % Plot original data
        % Plot the polynomial fit
        strain_fit = linspace(min(T_Strain{i}), max(T_Strain{i}), 100); % Generate a smooth range for fitting curve
        stress_fit = polyval(p, strain_fit); % Evaluate the polynomial at these strain values
        plot(strain_fit, stress_fit, '-', 'DisplayName', 'Polynomial Fit'); % Plot fitted curve
        xlabel('Strain');
        ylabel('Cauchy Stress (kPa)');
        title(['Cauchy Stress vs. Strain and Fit for Dataset ' num2str(i)]);
        legend show;
        hold off;
        
        % Plot the tangent stiffness (derivative of the fit)
        age = [26, 39, 47, 59, 62, 64, 70, 72, 74, 87, 90];
        % Fit a line to the data
        p = polyfit(age, tangentStiffnessAtLast, 1); % Linear fit (1st-degree polynomial)
        fittedLine = polyval(p, age); % Evaluate the fitted line at the age values

        
    else
        % Skip the plot for this dataset if RelaxTime is empty
        warning('Skipping dataset %d due to invalid RelaxTime.', i);
    end
end

%% Plot the data and fitted line
        figure;
        plot(age, tangentStiffnessAtLast, '.', 'LineStyle', 'none', 'MarkerSize', 25, 'DisplayName', 'Experimental Data'); % Original data
        hold on;
        plot(age, fittedLine, ':r', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit'); % Fitted line
        hold off;
        grid off;

        % Add labels and title
        xlabel('Age (years)', 'FontSize', 14, 'FontWeight','bold');
        ylabel('Tangent Stiffness (kPa)', 'FontSize', 14, 'FontWeight','bold');
       
        % Set axis limits and hide specific ticks
        ax = gca;
        ax.FontSize = 14; % Adjust axis number font size
        xlim([20, 100]);
        ylim([50, 380]);
       
        % Add legend
        legend('show', 'Location', 'best', 'Fontsize', 12);
        box off;
        
        % Set the plot aspect ratio
        pbaspect([1 1 1]); % Set the plot aspect ratio


%% Create a bar chart of the average tangent stiffness for each dataset
% figure;
% b = bar(tangentStiffnessAtLast); % Create a bar chart
% 
% % Define a custom color for each bar (for example, use a range of colors)
% colors = lines(numFiles); % Generates a colormap with distinct colors for each bar
% 
% % Apply the colors to each bar individually
% b.FaceColor = 'flat'; % Allow individual bar color settings
% b.CData = colors; % Set the colors for each dataset
% 
% % Set x-axis labels
% set(gca, 'XTickLabel', {'Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4', 'Dataset 5', 'Dataset 6'});
% 
% % Label for x-axis and y-axis
% xlabel('Datasets');
% ylabel('Tangent Stiffness (kPa)');
% 
% % Title of the chart
% title('Tangent Stiffness (at Stretch = 1.1) for Each Dataset');

