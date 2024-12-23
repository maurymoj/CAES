% File path
filePath = 'GT AR 2023 Part H.csv'; % Part H - MILES OF TRANSMISSION PIPE BY NOMINAL PIPE SIZE (NPS)

% Read the CSV file
data = readtable(filePath);

% Step 1: Filter data for "Natural Gas" only
gas_data = data(strcmp(data{:, 'PARTA5COMMODITY'}, 'Natural Gas'), :);  % Assuming column G is for fluid type

% Step 2: Calculate total mileage for each pipe diameter category
% Note: Adjust the column names if necessary to match your data

mileage_less_than_4 = sum(gas_data.PARTHON4LESS, 'omitnan');
mileage_6_inch = sum(gas_data.PARTHON6, 'omitnan');
mileage_8_inch = sum(gas_data.PARTHON8, 'omitnan');
mileage_10_inch = sum(gas_data.PARTHON10, 'omitnan');
mileage_12_inch = sum(gas_data.PARTHON12, 'omitnan');
mileage_14_inch = sum(gas_data.PARTHON14, 'omitnan');
mileage_16_inch = sum(gas_data.PARTHON16, 'omitnan');
mileage_18_inch = sum(gas_data.PARTHON18, 'omitnan');
mileage_20_inch = sum(gas_data.PARTHON20, 'omitnan');
mileage_22_inch = sum(gas_data.PARTHON22, 'omitnan');
mileage_24_inch = sum(gas_data.PARTHON24, 'omitnan');
mileage_26_inch = sum(gas_data.PARTHON26, 'omitnan');
mileage_28_inch = sum(gas_data.PARTHON28, 'omitnan');
mileage_30_inch = sum(gas_data.PARTHON30, 'omitnan');
mileage_32_inch = sum(gas_data.PARTHON32, 'omitnan');
mileage_34_inch = sum(gas_data.PARTHON34, 'omitnan');
mileage_36_inch = sum(gas_data.PARTHON36, 'omitnan');
mileage_38_inch = sum(gas_data.PARTHON38, 'omitnan');
mileage_40_inch = sum(gas_data.PARTHON40, 'omitnan');
mileage_42_inch = sum(gas_data.PARTHON42, 'omitnan');
mileage_44_inch = sum(gas_data.PARTHON44, 'omitnan');
mileage_46_inch = sum(gas_data.PARTHON46, 'omitnan');
mileage_48_inch = sum(gas_data.PARTHON48, 'omitnan');
mileage_52_inch = sum(gas_data.PARTHON52, 'omitnan');
mileage_56_inch = sum(gas_data.PARTHON56, 'omitnan');
mileage_58_and_over = sum(gas_data.PARTHON58OVER, 'omitnan');

% Step 3: Summing up mileage into required diameter ranges
total_mileage = mileage_less_than_4 + mileage_6_inch + mileage_8_inch + mileage_10_inch + ...
                mileage_12_inch + mileage_14_inch + mileage_16_inch + mileage_18_inch + ...
                mileage_20_inch + mileage_22_inch + mileage_24_inch + mileage_26_inch + ...
                mileage_28_inch + mileage_30_inch + mileage_32_inch + mileage_34_inch + ...
                mileage_36_inch + mileage_38_inch + mileage_40_inch + mileage_42_inch + ...
                mileage_44_inch + mileage_46_inch + mileage_48_inch + ...
                mileage_52_inch + mileage_56_inch + mileage_58_and_over;

mileage_less_12 = mileage_less_than_4 + mileage_6_inch + mileage_8_inch + mileage_10_inch;
mileage_12_34 = mileage_12_inch + mileage_14_inch + mileage_16_inch + mileage_18_inch + ...
                mileage_20_inch + mileage_22_inch + mileage_24_inch + mileage_26_inch + ...
                mileage_28_inch + mileage_30_inch + mileage_32_inch + mileage_34_inch;
mileage_more_34 = mileage_36_inch + mileage_38_inch + mileage_40_inch + mileage_42_inch + ...
                  mileage_44_inch + mileage_46_inch + mileage_48_inch + ...
                  mileage_52_inch + mileage_56_inch + mileage_58_and_over;

% Step 4: Calculate percentages
percent_less_12 = (mileage_less_12 / total_mileage) * 100;
percent_12_34 = (mileage_12_34 / total_mileage) * 100;
percent_more_34 = (mileage_more_34 / total_mileage) * 100;

% Step 5: Display results in a table
resultTable = table({'< 12 inches'; '12-34 inches'; '> 34 inches'}, ...
                    [mileage_less_12; mileage_12_34; mileage_more_34], ...
                    [percent_less_12; percent_12_34; percent_more_34], ...
                    'VariableNames', {'Diameter_Range', 'Total_Mileage', 'Percentage'});

disp(resultTable);

%% Weighted average diameter
% Sum of mileage for each diameter category
mileage_values = [
    sum(gas_data.PARTHON4LESS, 'omitnan'), sum(gas_data.PARTHON6, 'omitnan'), ...
    sum(gas_data.PARTHON8, 'omitnan'), sum(gas_data.PARTHON10, 'omitnan'), ...
    sum(gas_data.PARTHON12, 'omitnan'), sum(gas_data.PARTHON14, 'omitnan'), ...
    sum(gas_data.PARTHON16, 'omitnan'), sum(gas_data.PARTHON18, 'omitnan'), ...
    sum(gas_data.PARTHON20, 'omitnan'), sum(gas_data.PARTHON22, 'omitnan'), ...
    sum(gas_data.PARTHON24, 'omitnan'), sum(gas_data.PARTHON26, 'omitnan'), ...
    sum(gas_data.PARTHON28, 'omitnan'), sum(gas_data.PARTHON30, 'omitnan'), ...
    sum(gas_data.PARTHON32, 'omitnan'), sum(gas_data.PARTHON34, 'omitnan'), ...
    sum(gas_data.PARTHON36, 'omitnan'), sum(gas_data.PARTHON38, 'omitnan'), ...
    sum(gas_data.PARTHON40, 'omitnan'), sum(gas_data.PARTHON42, 'omitnan'), ...
    sum(gas_data.PARTHON44, 'omitnan'), sum(gas_data.PARTHON46, 'omitnan'), ...
    sum(gas_data.PARTHON48, 'omitnan'), sum(gas_data.PARTHON52, 'omitnan'), ...
    sum(gas_data.PARTHON56, 'omitnan'), sum(gas_data.PARTHON58OVER, 'omitnan')
];

% Corresponding nominal diameters for each mileage category
diameter_values = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, ...
                   34, 36, 38, 40, 42, 44, 46, 48, 52, 56, 58];

% Step 1: Calculate the total mileage
total_mileage = sum(mileage_values);

% Step 2: Calculate the weighted sum of diameters
weighted_sum = sum(mileage_values .* diameter_values);

% Step 3: Calculate the average diameter
average_diameter = weighted_sum / total_mileage;

% Display the result
fprintf('The average diameter of the pipelines is %.2f inches.\n', average_diameter);

D_mm = diameter_values*25.4;
L_km = 1.6093*mileage_values;

% Step 1: Calculate the total km
total_km = sum(L_km);

% Step 2: Calculate the weighted sum of diameters
weighted_sum_mm = sum(L_km.*D_mm);

% Step 3: Calculate the average diameter
average_D_mm = weighted_sum_mm / total_km;

fprintf('The average diameter of the pipelines is %.2f mm.\n', average_D_mm);

