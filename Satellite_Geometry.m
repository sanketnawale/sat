function [Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats, Sat_IDs, Latitudes, Longitudes,Sat_To_NodeRM_Delay] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector)

    Num_Satellites = size(oev, 1); 
    num_steps = length(Time_Vector);
    num_nodes = size(Node_Coordinates, 1);

    % 🌍 Convert Ground Node Positions to Cartesian (ECEF)
    Node_Positions = zeros(num_nodes, 3);
    for n = 1:num_nodes
        lat = Node_Coordinates(n, 1);
        lon = Node_Coordinates(n, 2);
        Node_Positions(n, :) = [
            Earth_Radius * cosd(lat) * cosd(lon),
            Earth_Radius * cosd(lat) * sind(lon),
            Earth_Radius * sind(lat)
        ];
    end

    % ✅ Initialize Output Arrays
    Distances = zeros(num_nodes, Num_Satellites, num_steps);
    Elevation_Angles = zeros(num_nodes, Num_Satellites, num_steps);
    Ground_Distances = zeros(num_nodes, Num_Satellites, num_steps);
    Visibility = zeros(num_nodes, Num_Satellites, num_steps);
    Num_Visible_Sats = zeros(num_nodes, num_steps);
    Sat_IDs = cell(num_nodes, num_steps);  % Store satellite IDs for each node
    % ✅ Initialize Output Arrays
    Sat_To_NodeRM_Delay = zeros(Num_Satellites, num_steps);


    % ✅ New Outputs: Latitude & Longitude Storage
    Latitudes = zeros(Num_Satellites, num_steps);
    Longitudes = zeros(Num_Satellites, num_steps);

    % 📂 Open a text file to store results
    logFileName = 'Satellite_Geometry_Log.txt';
    logFile = fopen(logFileName, 'w');
    if logFile == -1
        error('Error: Unable to create or open the log file.');
    end
   
    

    % 🔄 Loop Over Time Steps
    for t = 1:num_steps
        fprintf('\n⏳ Time %.2f min:\n', Time_Vector(t) / 60);  % Convert sec to min
        fprintf(logFile, '\n⏳ Time %.2f min:\n', Time_Vector(t) / 60);

        

        for s = 1:Num_Satellites
            % ✅ Compute satellite position dynamically at time t
            mean_motion = 2 * pi / 5400;  % Assume 5400s orbital period
            true_anomaly = oev(s, 6) + mean_motion * Time_Vector(t);

            [Sat_X, Sat_Y, Sat_Z] = orbital_to_cartesian(...
                oev(s, 1), oev(s, 2), oev(s, 3), ...
                oev(s, 5), oev(s, 4), true_anomaly);
            Satellite_Position = [Sat_X, Sat_Y, Sat_Z];
            % ✅ Convert ECEF to Latitude & Longitude
            Latitudes(s, t) = atan2d(Sat_Z, sqrt(Sat_X^2 + Sat_Y^2));  % Latitude
            Longitudes(s, t) = atan2d(Sat_Y, Sat_X);  % Longitude
             Vector_Sat_To_NodeRM = Satellite_Position - Node_Positions(3, :);
                Distance_Sat_To_NodeRM(s, t) = norm(Vector_Sat_To_NodeRM);

    % ✅ Compute Signal Delay (Satellite → NodeRM)
        Speed_of_Light = 3e8;  % Speed of light in meters/second
        Sat_To_NodeRM_Delay(s, t) = Distance_Sat_To_NodeRM(s, t) * 1e3 / Speed_of_Light;  % Convert km to meters

            % 🛰️ Print Satellite Position
            fprintf('🛰️ Sat %d Pos -> X: %.2f km, Y: %.2f km, Z: %.2f km\n', ...
                    s, Sat_X / 1e3, Sat_Y / 1e3, Sat_Z / 1e3);
            fprintf(logFile, '🛰️ Sat %d Pos -> X: %.2f km, Y: %.2f km, Z: %.2f km\n', ...
                s, Sat_X / 1e3, Sat_Y / 1e3, Sat_Z / 1e3);
            

            for n = 1:num_nodes
                % ✅ Compute Distance Between Node and Satellite
                Vector = Satellite_Position - Node_Positions(n, :);
                Distances(n, s, t) = norm(Vector);

                % ✅ Compute Ground Distance
                Ground_Distances(n, s, t) = sqrt(Distances(n, s, t)^2 - H^2);
                % ✅ Compute Signal Delay in Seconds
                 Speed_of_Light = 3e8; % Speed of light in meters/second
                Signal_Delay(n, s, t) = Ground_Distances(n, s, t) * 1e3 / Speed_of_Light; % Convert km to meters
                 
                % ✅ Compute Distance Between Satellite and NodeRM
                Vector_NodeRM = Satellite_Position - Node_Positions(3, :); % NodeRM is Node 3
                Sat_To_NodeRM_Distance(s, t) = norm(Vector_NodeRM); % Distance in km

% ✅ Compute Signal Delay from Satellite to NodeRM
                Sat_To_NodeRM_Delay(s, t) = Sat_To_NodeRM_Distance(s, t) * 1e3 / Speed_of_Light; % Convert km to meters

% ✅ Print the computed delay for debugging
                fprintf('🛰️ Sat %d → NodeRM Delay: %.6f seconds\n', s, Sat_To_NodeRM_Delay(s, t));
                fprintf(logFile, '🛰️ Sat %d → NodeRM Delay: %.6f seconds\n', s, Sat_To_NodeRM_Delay(s, t));

                % ✅ Print the computed delay for debugging (optional)
                fprintf('📡 Node %d → Sat %d Delay: %.6f seconds\n', n, s, Signal_Delay(n, s, t));
                fprintf(logFile, '📡 Node %d → Sat %d Delay: %.6f seconds\n', ...
                    n, s, Signal_Delay(n, s, t));

                % ✅ Compute Elevation Angle
                Elevation_Angles(n, s, t) = asind(dot(Vector, Node_Positions(n, :)) / ...
                                               (norm(Vector) * norm(Node_Positions(n, :))));

                % ✅ Visibility Condition
                if Elevation_Angles(n, s, t) > 10  % Satellites above 10° elevation are visible
                    Visibility(n, s, t) = 1;
                    Sat_IDs{n, t} = [Sat_IDs{n, t}, s];  % Store satellite ID
                end
            end
        end

        % ✅ Count Visible Satellites per Node
        Num_Visible_Sats(:, t) = sum(Visibility(:, :, t), 2);
        

        % ✅ Print visibility results with satellite IDs
        fprintf('🔍 Rome sees %d satellites: %s\n', Num_Visible_Sats(1, t), mat2str(Sat_IDs{1, t}));
        fprintf('🔍 Milan sees %d satellites: %s\n', Num_Visible_Sats(2, t), mat2str(Sat_IDs{2, t}));
         fprintf(logFile, '🔍 Rome sees %d satellites: %s\n', Num_Visible_Sats(1, t), mat2str(Sat_IDs{1, t}));
        fprintf(logFile, '🔍 Milan sees %d satellites: %s\n', Num_Visible_Sats(2, t), mat2str(Sat_IDs{2, t}));
    
    end
end
