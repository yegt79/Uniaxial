function [constantIndex, zeroIndex, secondZeroIndex] = findConstantIndexWithDisplay(data, threshold, tolerance, minDistance)
    % This function finds the first index where at least 'threshold' 
    % consecutive elements of the data are the same.
    % It also searches for:
    %   - The nearest value close to 1 before the constant index
    %   - The second closest value close to 1 going backward from the first found zeroIndex
    %     (which is at least 'minDistance' indices before the first zeroIndex)
    %
    % Input:
    %   data - Array of numerical values to search through
    %   threshold - Number of consecutive equal values to consider constant
    %   tolerance - Acceptable range around 1 for the near-one search
    %   minDistance - Minimum number of indices between the first and second zeroIndex
    %
    % Output:
    %   constantIndex - Index where constant values start (0 if not found)
    %   zeroIndex - Closest index to constantIndex with a value close to 1 (0 if not found)
    %   secondZeroIndex - Second closest value to 1 going backward from zeroIndex (0 if not found)

    % Ensure inputs are valid
    if ~isnumeric(threshold) || ~isnumeric(tolerance) || ~isnumeric(minDistance)
        error('Threshold, tolerance, and minDistance must be numeric values.');
    end

    constantIndex = 0;        % Initialize constantIndex to 0 (not found)
    zeroIndex = 0;            % Initialize zeroIndex to 0 (not found)
    secondZeroIndex = 0;      % Initialize secondZeroIndex to 0 (not found)

    % Find the constant index
    for i = 1:(length(data) - threshold)
        if all(data(i:i+threshold-1) == data(i))  % Check if values are the same
            constantIndex = i;
            break;
        end
    end

    % If a constant index is found, search upward for the nearest value close to 1
    if constantIndex > 0
        for j = constantIndex:-1:1
            if abs(data(j) - 1) <= tolerance  % Check if value is within tolerance of 1
                zeroIndex = j;
                break;
            end
        end

        % If zeroIndex is found, go backward to find the second nearest value close to 1
        if zeroIndex > 0
            for k = zeroIndex-1:-1:1
                if abs(data(k) - 1) <= tolerance  % Check if value is within tolerance of 1
                    if zeroIndex - k >= minDistance
                        secondZeroIndex = k;  % Store the second index
                        break;  % Stop searching once a valid index is found
                    end
                end
            end
       end

        % Display results
        disp(['The constant value starts at index ', num2str(constantIndex), ...
              ' with a value of ', num2str(data(constantIndex))]);
        if zeroIndex > 0
            disp(['The nearest value close to 1 before the constant value is at index ', ...
                  num2str(zeroIndex), ' with a value of ', num2str(data(zeroIndex))]);
        else
            disp('No value close to 1 was found before the constant value.');
        end
        if secondZeroIndex > 0
            disp(['The second value close to 1 backward from zeroIndex is at index ', ...
                  num2str(secondZeroIndex), ' with a value of ', num2str(data(secondZeroIndex))]);
        else
            disp('No second value close to 1 was found before the first zeroIndex with the required distance.');
        end
    else
        disp('No constant value found in the data for the given threshold.');
    end
end
