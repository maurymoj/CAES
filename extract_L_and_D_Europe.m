% strVector = {
%     "{'diameter_mm': 600.0, 'end_year': 2050, 'is_H_gas': 1.0, 'is_bothDirection': 1.0, 'lat_mean': 54.3708315, 'length_km': 101.11970692844721, 'long_mean': 9.5367365, 'max_cap_M_m3_per_d': 8.41096, 'max_pressure_bar': 63.51782227405727, 'num_compressor': 1.0, 'operator_name': None, 'pipe_class_LKD': None, 'start_year': 1996.0, 'waterDepth_m': None}"
% };

strVector = param;
% Preallocate arrays to store the extracted values
numEntries = length(strVector);
diameterValues = NaN(numEntries, 1);
lengthValues = NaN(numEntries, 1);

% Loop through each cell in the string vector
for i = 1:numEntries
    % Get the current string
    currentStr = strVector{i};
    
    % Extract diameter_mm value
    diameterMatch = regexp(currentStr, "'diameter_mm':\s*([0-9.]+)", 'tokens');
    if ~isempty(diameterMatch)
        diameterValues(i) = str2double(diameterMatch{1}{1});
    end
    
    % Extract length_km value
    lengthMatch = regexp(currentStr, "'length_km':\s*([0-9.]+)", 'tokens');
    if ~isempty(lengthMatch)
        lengthValues(i) = str2double(lengthMatch{1}{1});
    end
end

% Display results
% disp("Diameter values:");
% disp(diameterValues);
% disp("Length values:");
% disp(lengthValues);

%% Get values for specific country
Country = string(country_code);
(strcmp(Country,"['DE', 'DE']"));
D_DE =diameterValues(strcmp(Country,"['DE', 'DE']"));
L_DE =lengthValues(strcmp(Country,"['DE', 'DE']"));