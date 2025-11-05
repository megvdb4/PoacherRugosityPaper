% Load and plot poacher DPIV data 
% OH Hawkins, Meg L Vandenberg 08.20.2023
% Load PIV text files for a range of frames, visualize average X velocity.
% Text files have the following components:
        % Column 1: X coordinates in pixels        
        % Column 1: X coordinates in pixels
        % Column 2: Y coordinates in pixels (needs to be flipped)
        % Column 3: U (Velocity in x direction in pixels/frame)
        % Column 4: V (Velocity in y direction in pixels/frame)
        % Column 5: Vector type
        % Column 6: Vorticity (1/frame)
        % Column 7: Magnitude (pixel/frame)
        % Column 8: Divergence (1/frame)
        % Column 9: DCEV (1) 
        % Column 10: Simple shear (1/frame)
        % Column 11: Simple strain (1/frame)
        % Column 12: Vector direction (degrees)
      

% Set the range of frames
start_frame = 1;
end_frame = 5000; % Replace with the actual number of frames

% Load data for the first frame to get the grid dimensions
filename_first = sprintf('PIVlab_%04d.txt', start_frame);
data_first = readmatrix(filename_first);
x_first = data_first(:, 1);
y_first = data_first(:, 2);

[X_mm, Y_mm] = meshgrid(unique(x_first), unique(y_first));

% Initialize the sum array
sumArray = zeros(size(X_mm)); % Initialize with the size of the meshgrid

% Loop through frames
for frame = start_frame:end_frame
    % Load data for the current frame
    filename = sprintf('PIVlab_%04d.txt', frame); % Update the file name pattern
    data = readmatrix(filename);
    
   % Extract data columns
            x = data(:, 1); % x positions in pixels
            y = data(:, 2); % y positions in pixels
            u = data(:, 3); % u velocity components
            v = data(:, 4); % v velocity components
            
            % Conversion factor (pixels to millimeters)
            pixel_to_mm = 0.1; % Adjust this value based on your specific conversion factor
            
            % Conversion factor (frames to seconds)
            fps = 500; % Since we sampled every second frame
            
            % Convert pixel positions to millimeters
            x_mm = x * pixel_to_mm;
            y_mm = y * pixel_to_mm;
            
            % Convert pixels to mm and frames to seconds to get velocity in mm/s
            u_frames = (u * pixel_to_mm) * fps; 
            v_frames = (v * pixel_to_mm) * fps; 
            
            % Create a meshgrid for x and y positions
            X_mm = x_mm;
            Y_mm = y_mm;
            
            % Create a meshgrid for x and y positions in millimeters
            [X_mm, Y_mm] = meshgrid(unique(x *0.1), unique(y *0.1));

            % Reshape velocity components to match the meshgrid
            U = reshape(u_frames, size(X_mm));
            V = reshape(v_frames, size(X_mm));
            
    
    % Add the current frame's data to the sum array
    sumArray = sumArray + U;
end

% Calculate the average velocity array
AvgArray = sumArray / (end_frame - start_frame + 1); % Divide by the number of frames

%% Create the figure of average velocity in the x direction
AverageXVel = figure;
    contourf(X_mm, -Y_mm, AvgArray, 100, 'LineStyle', 'none');
    h = colorbar; % Get the colorbar handle
    h.Label.String = 'X Velocity (mm/s)'; % Set colorbar label
    xlabel('X Position (mm)');
    ylabel('Y Position (mm)');
    title('Contour Plot of Velocity in X Direction');
    axis equal;

% Get starting coordinates in x and y direction 
start_x = 0; 
start_y = -380;


% Add streamlines
hold on;


streamline(X_mm, Y_mm, U, V, start_x, start_y);

% Plot contour lines on top of contour filled plot
contour(X_mm, -Y_mm, AvgArray, 25, 'LineColor', 'k'); % You can adjust the number of contour lines and line color

% Release the hold on the plot
hold off;


% Save the figure
saveas(AverageXVel, 'AverageXVelocity_contour.pdf');



%% Create transects and pull average velocities out for each transect

% Calculate the step size for transects
transect_step = 6.5; % in millimeters

% Calculate the number of transects
num_transects = floor(max(max(X_mm)) / transect_step);

% Initialize an empty matrix to store transect data
transects = zeros(size(AvgArray, 1), num_transects);

%Initialize an empty matrix to store mask data
mask = zeros(size(Y_mm, 1), num_transects);

% Loop through transects
for i = 1:num_transects
    % Calculate the X position of the transect
    transect_x = (i - 0.5) * transect_step;
    
    % Find the closest X position in the meshgrid
    [~, x_idx] = min(abs(X_mm(1, :) - transect_x));

    % Extract the transect data from the AvgArray
    transect_data = AvgArray(:, x_idx);
    
    % Store the transect data as a new column in the transects matrix
    transects(:, i) = transect_data;
end

%%
% Assume AvgArray is your matrix with NaN values
% Assume Y_mm is your matrix with Y positions

% Step 1: Extract NaN values from AvgArray and store their positions
nanPositions = isnan(transects);  % Logical matrix with 1s where NaNs are present
nanIndices = find(nanPositions);  % Linear indices of NaN values
[rowIndices, colIndices] = ind2sub(size(transects), nanIndices);  % Convert linear indices to row and column indices

% Step 2: Extract corresponding Y positions from Y_mm
nanYPositions = Y_mm(nanIndices);

% Step 3: Create a matrix containing NaN positions and corresponding Y positions
nanDataMatrix = [rowIndices, colIndices, nanYPositions];

csv_filename = 'NaNpositions_A.chiloensis.csv';
writematrix(nanDataMatrix, csv_filename);


%% Plot all transects on one graph

% Create a figure for the plot
figure;

% Loop through transects to plot each one
for i = 1:num_transects
    plot( transects(:, i),-Y_mm, 'LineWidth', 1.5);
    hold on;
end

% Customize the plot
xlabel('Velocity');
ylabel('y_{mm}');
title('Transect Velocities');

% Define legend entries (change these as needed)
legend_entries = arrayfun(@(x) sprintf('Transect %d', x), 1:num_transects, 'UniformOutput', false);
legend(legend_entries);
grid on;

%%
% Write the matrix to a CSV file
csv_filename_1 = 'Transects_A.chiloensis.csv';
writematrix(transects, csv_filename_1);

% Write the matrix to a CSV file
csv_filename = 'YHeights_A.chiloensis.csv';
writematrix(Y_mm, csv_filename);
