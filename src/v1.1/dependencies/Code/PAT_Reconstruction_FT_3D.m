function pa_img = PAT_Reconstruction_FT_3D(P)
%% FFT image reconstruction
%inputs:
% inputData---------> raw data which is a 3D matrixt (scans*samples*sensors)
% sensorsPositions -> position of the sensors
% scanArray --------> array that shows the selected scans
% sensorArray ------> array that shows the selected sensors
% xRange -----------> range of the imaging area in x direction
% yRange -----------> range of the imaging area in y direction

% samplingRate -----> sampling rate in the data acquisition
% speedOfSound -----> speed of sound in the imaging medium (homogeneous medium)
% method -----------> reconstruction method including: 1=Delay and sum(DAS)...
% 2=universal back projection(UBP),3=Fourier transform, 4-time reversal, ...
% 5-compressed sensing
%% input parameters
inputData       = P.inputData;
samplingRate    = P.fs;
speedOfSound    = P.V_M;
numSensors      = P.numSensors;
interpolation   = P.interpolation;
sensorPitch     = P.sensorPitch/1000; %[m]
delay           = P.delay; %micro second
delay           = delay*1e-6; %second
slicing         = P.slicing;


tempMat = [];
for i = 1: size(inputData,1)
    tempMat = cat(2, tempMat, reshape(inputData(i,:,:), size(inputData,2), []));
end
inputData = tempMat';

% sensorsPositions(:,[1 2]) = sensorsPositions(:,[2 1]);
medium.sound_speed = speedOfSound/1000;
dt = 1/samplingRate;
if interpolation
    arrayLength = sensorPitch*(numSensors-1); %[m]
    dx = medium.sound_speed/samplingRate;          %[m]
    kgrid = makeGrid(size(inputData,2), dx, round(arrayLength/dx), dx );
    cart_sensor_mask (1,:) = zeros(1,numSensors);
    cart_sensor_mask (2,:) = -arrayLength/2:sensorPitch:arrayLength/2;
    binary_sensor_mask = zeros(size(inputData,2), round(arrayLength/dx));
    binary_sensor_mask(1,:)=1;
    
    % interpolate data to remove the gaps and assign to sensor structure
    inputData = interpCartData(kgrid, inputData, cart_sensor_mask, binary_sensor_mask);
else
    dx = sensorPitch;
end

pa_img = kspaceLineRecon(inputData', dx, dt, medium.sound_speed, 'PosCond',true);
% pa_img = rot90(double(pa_img));
% create the computational grid
% x = x_size; % mm
% y = y_size; % mm
% z = z_size;
%     Nx = 800;               % number of grid points in the x (row)direction % 3140
%     Ny = 800;               % number of grid points in the y (column)direction (sensors)
% dx = medium.sound_speed/sampling_rate;   % grid point spacing in the x direction [m] total 0.05888m
% dy =  dx;      % grid point spacing in the y direction [m]
% dz =  dx;
% z = dx*length(Samples);
% Nx = ceil( x/(1000*dx))+1;               % number of grid points in the x (row)direction % 3140
% Ny = ceil( y/(1000*dy))+1;               % number of grid points in the y
% Nz = length(Samples);

% maskTemp = inputLoc'./1000;
% mask(1,:)=maskTemp(1,:);
% mask(2,:)=maskTemp(3,:);
% binary_sensor_mask = zeros(kgrid.Nx, kgrid.Ny);
% binary_sensor_mask(1, :) = 1;
%
% SensorDataInterpolated = interpCartData(kgrid, sensor_data', mask, binary_sensor_mask);

% define a binary line sensor
%                     sensor.mask = zeros(Nx, Nz);
%                     sensor.mask(1, :) = 1;
% p_xy1 = kspaceLineRecon(SensorDataInterpolated', dx, dt, medium.sound_speed, 'PosCond',true);
%                     [Nx_recon, Ny_recon] = size(p_xy);
%                     kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy);
%                     % resample p_xy to be the same size as source.p0
%                     p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));

if 0
    kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);
    sensor_data = squeeze(selectedData);
    sensor.mask = inputLoc'./1000;
    figure;
    plot3(inputLoc(:,1), inputLoc(:,2), inputLoc(:,3), '.');
    %     xlabel(['[' prefix 'm]']);
    %     ylabel(['[' prefix 'm]']);
    %     zlabel(['[' prefix 'm]']);
    axis equal;
    grid on;
    box on;
    
    % define a binary planar sensor
    % sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    % sensor.mask(1, :, :) = 1;
    binary_sensor_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    binary_sensor_mask(:, :, 1) = 1;
    
    %                     [grid_data, order_index, reorder_index] = cart2grid(kgrid, sensor.mask);
    
    % select suitable axis scaling factor
    %     [x_sc, scale, prefix] = scaleSI(max(sensor.mask(:)));  %#ok<ASGLU>
    
    % create the figure
    
    kgrid.t_array = 0:dt:dt*Nz;
    
    %                     % input arguments
    %                     input_args = {'PlotLayout', true, 'PlotPML', false, ...
    %                         'DataCast', 'single', 'CartInterp', 'nearest'};
    
    % run the simulation
    % sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    % sensor_data = sensor_data';
    %                     load sensorData_1
    % load sensorData_1
    kgrid.Nt = length(Samples);
    sensor_data_rs = reshape(sensor_data', 15, 15, kgrid.Nt);
    
    pa_img = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
        medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true);
    %                     PlotGUI_2(app,pa_img);
    
    SensorDataInterpolated = interpCartData(kgrid, sensor_data', sensor.mask, binary_sensor_mask);
    % reshape sensor data to y, z, t
    % sensor_data_rs = reshape(sensor_data', Ny, Nz, kgrid.Nt);
    sensor_data_rs = reshape(SensorDataInterpolated, Nx, Ny, kgrid.Nt);
    % reconstruct the initial pressure
    pa_img = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
        medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);
end
end