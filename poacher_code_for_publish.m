%% Poacher DPIV Data
% Load and plot poacher DPIV data 
% OH Hawkins, Meg L Vandenberg 08.20.2023
% Load PIV text files for a range of frames, visualize average X velocity.
% Text files from PIVlabhave the following components:
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
      
%CHECK YOUR TXT FILES UNITS - this is assuming you have calibrated in PIV
%lab
% Set the range of frames
start_frame = 1;
end_frame = 4999; % Replace with the actual number of frames depending on your export from PIVlab

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
            
            % convert m to mm
            x_mm = x * 1000;
            y_mm = y * 1000;

            u_frames = u * 1000;
            v_frames = v * 1000;
            
            % Create a meshgrid for x and y positions
            X_mm = x_mm;
            Y_mm = y_mm;
            
            % Create a meshgrid for x and y positions in millimeters 
            [X_mm, Y_mm] = meshgrid(unique(x *1000), unique(y *1000));

            % Reshape velocity components to match the meshgrid
            U = reshape(u_frames, size(X_mm));
            V = reshape(v_frames, size(X_mm));
            
    
    % Add the current frame's data to the sum array
    sumArray = sumArray + U;
end

% Calculate the average velocity array
AvgXArray = sumArray / (end_frame - start_frame + 1); % Divide by the number of frames

%% Create transects and pull average velocities out for each transect

% Calculate the step size for transects
transect_step = 2; % in millimeters - this is personal preferance

% Calculate the number of transects
num_transects = floor(max(max(X_mm)) / transect_step);

% Initialize an empty matrix to store transect data
transects = zeros(size(AvgXArray, 1), num_transects);

%Initialize an empty matrix to store mask data
mask = zeros(size(Y_mm, 1), num_transects);

% Loop through transects
for i = 1:num_transects
    % Calculate the X position of the transect
    transect_x = (i - 0.5) * transect_step;
    
    % Find the closest X position in the meshgrid
    [~, x_idx] = min(abs(X_mm(1, :) - transect_x));

    % Extract the transect data from the AvgArray
    transect_data = AvgXArray(:, x_idx);
    
    % Store the transect data as a new column in the transects matrix
    transects(:, i) = transect_data;
end

% Loop through transects
for i = 1:num_transects
    % Calculate the X position of the transect
    transect_x = (i - 0.5) * transect_step;

    % Find the closest X position in the meshgrid
    [~, x_idx] = min(abs(X_mm(1, :) - transect_x));

    % Extract the transect data from the AvgArray
    transect_data = AvgXArray(:, x_idx);

    % Store the transect data as a new column in the transects matrix
    transects(:, i) = transect_data;

    % Store the exact X position of each transect
    transect_positions(i) = X_mm(1, x_idx);
end

%% Plot velocity in the X direction

%Matlab plots upside down for some reason so flip your data first to make
%it right side up
AvgXArray([1:size(AvgXArray,1)],:) = AvgXArray([size(AvgXArray,1):-1:1],:);
% Create the figure of average velocity in the x direction
AverageXVel = figure;
set(gcf, 'Color', 'black'); % Set background color to black

contourf(X_mm, Y_mm, AvgXArray, 100, 'LineStyle', 'none');
colormap('jet'); % You can choose a different colormap if desired
h = colorbar;
h.Label.String = 'Y Velocity (mm/s)';
h.Label.Color = 'white';

 % Set colorbar label color to white
xlabel('X Position (mm)', 'Color', 'white');
ylabel('Y Position (mm)', 'Color', 'white');
title('Contour Plot of Velocity in X Direction', 'Color', 'white');
axis equal;
ax = gca;
ax.Color = 'black'; % Set axes background color to black
ax.XColor = 'white'; % Set x-axis label and tick color to white
ax.YColor = 'white';
ax.XTickLabelColor = 'white'; % Set x-axis tick label color to white
% Set y-axis label and tick color to white

% Get starting coordinates in x and y direction 
start_x = 0; 
start_y = -380;

% Release the hold on the plot
hold off;

% Save the figure
saveas(AverageXVel, 'AverageXVelocity_contour.pdf');

%% Find the BL velocities

%Using the values you fins from the velocity plot for you boundary layer

%make row into column
transect_positions = transect_positions'

% get points where veolcity is between 10-17 mm/s - change these consistent with what
% your boundary layer is 
[row, col] = find(AvgXArray > 0 & AvgXArray< 16);
Xrow = X_mm(1,:);
Yrow = Y_mm(:,1);
col2 = Xrow(col);
row2 = Yrow(row);

merged = [col2', row2];
G = findgroups(merged(:,1));
C = splitapply(@(m){m},merged,G);

maxX = [];
maxY = [];

for i = 1:length(C)
    dat = C{i};
    maxX = [maxX; dat(1,1)];
    maxY = [maxY, max(dat(:,2))];
end

c = polyfit(maxX,maxY,5);
ysmooth = polyval(c, col2);

%% Plotting Boundary Layer on top 

% Create the figure of average velocity in the x direction
AverageXVel = figure;
set(gcf, 'Color', 'black'); % Set background color to black

contourf(X_mm, Y_mm, AvgXArray, 97, 'LineStyle', 'none');
colormap('jet'); % You can choose a different colormap if desired
h = colorbar;
h.Label.String = 'Y Velocity (mm/s)';
h.Label.Color = 'white';

 % Set colorbar label color to white
xlabel('X Position (mm)', 'Color', 'white');
ylabel('Y Position (mm)', 'Color', 'white');
title('Contour Plot of Velocity in Y Direction', 'Color', 'white');
axis equal;
ax = gca;
ax.Color = 'black'; % Set axes background color to black
ax.XColor = 'white'; % Set x-axis label and tick color to white
ax.YColor = 'white';

% Get starting coordinates in x and y direction 
start_x = 0; 
start_y = -380;

% Add BL line and points
hold on;
plot(col2, row2, 'w.', 'LineWidth', 2)
plot(col2, ysmooth, 'w-', 'LineWidth', 2)
hold off;

%% Saving
% Save the figure
saveas(AverageXVel, 'AverageXVelocity_contour.pdf');