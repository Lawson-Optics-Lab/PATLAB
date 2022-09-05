function pa_img = PAT_Reconstruction_TR_3D(P)
%% Time reversal image reconstruction
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
sensorsPositions= P.inputLoc;
xRange          = P.x_range; %low to high
yRange          = P.y_range; %high to low
zRange          = P.z_range; %low to high
samplingRate    = P.fs;
speedOfSound    = P.V_M;

% directivity     = P.directivity;
v               = P.directivityAngle;
% dMethod         = P.directivityMethpod;
delay           = P.delay; %micro second
slicing         = P.slicing;
delay = delay*1e-6; %second

tempMat = [];
for i = 1: size(inputData,1)
    tempMat = cat(2, tempMat, reshape(inputData(i,:,:), size(inputData,2), []));
end
inputData = tempMat';
%% calculate imaging area & parameters
Npixel_x = size(xRange, 2);
Npixel_y = size(yRange, 2);
Npixel_z = size(zRange, 2);
try
    dx = abs(xRange(2)-xRange(1))/1000;
catch
    dx = 0;
end
try
    dy = abs(yRange(2)-yRange(1))/1000;
catch
    dy = 0;
end
try
    dz = abs(zRange(2)-zRange(1))/1000;
catch
    dz = 0;
end
% sensorsPositions(:,[1 2]) = sensorsPositions(:,[2 1]);
if size(sensorsPositions,1)==2 && strcmp(slicing, 'xyz')
    %     2D Reconstruction
    kgrid_recon = makeGrid(Npixel_x, dx, Npixel_y, dy);
    sensor.mask = sensorsPositions; % [m]
    sensor.time_reversal_boundary_data = inputData;
    Sc = size(inputData,2);
    dt = 1/samplingRate;
    kgrid_recon.t_array = ((1:Sc)-1)*dt;
    medium.sound_speed = speedOfSound/1000;
    medium.density = 1000;  % [kg/m^3]
    source.p0 = 0;
    pa_img = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor);
%     pa_img = rot90(double(pa_img));
else
    %     3D Reconstruction
    if size(sensorsPositions,2)==2; sensorsPositions = cat(2,sensorsPositions, squeeze(zeros(size(sensorsPositions,1), 1 ,size(sensorsPositions,3)))); end
    TP = []; %Transducers position
    if strcmp(slicing, 'xyz')
        TP = sensorsPositions;
    elseif strcmp(slicing, ' yzx')
        TP(:,1) = sensorsPositions(:,2);
        TP(:,2) = sensorsPositions(:,3);
        TP(:,3) = sensorsPositions(:,1);
    elseif strcmp(slicing, '  zxy')
        TP(:,1) = sensorsPositions(:,3);
        TP(:,2) = sensorsPositions(:,1);
        TP(:,3) = sensorsPositions(:,2);
    end
    
    kgrid_recon = makeGrid(Npixel_x, dx, Npixel_y, dy, Npixel_z, dz);
    sensor.mask = TP'./1000; % [m]
    sensor.time_reversal_boundary_data = inputData;
    Sc = size(inputData,2);
    dt = 1/samplingRate;
    kgrid_recon.t_array = ((1:Sc)-1)*dt;
    medium.sound_speed = speedOfSound/1000;
    medium.density = 1000;  % [kg/m^3]
    source.p0 = 0;
    pa_img = kspaceFirstOrder3D(kgrid_recon, medium, source, sensor);
    pa_img = double(pa_img);
end



    

% save(filename,'p0_recon');
% figure;
% imshow(p0_recon(300:720,300:720),[-1 1]);
% axis image;

if 0
    lambda = 1000*speedOfSound / samplingRate; % Wavelength[mm]
    %% Reconstruction
    Nsample = size(inputData, 2);
    Nstep   = length(sensorArray);
    f = waitbar(0,'1','Name','BackProjection Reconstruction');
    % setappdata(f,'canceling',0);
    
    iter = 0;
    for i = scanArray
        binary_sensor_mask= zeros(Npixel_x,Npixel_y);
        for i= 1:lendgth(sensorsPositions,2)
            binary_sensor_mask(round(sensorsPositions(i,1)/2)+Npixel_x/2,round(sensorsPositions(i,2)/2)+Npixel_y/2) = 1;
        end
    end
    
    % assign the time reversal data
    sensor.time_reversal_boundary_data = sensor_data;
    % run the time reversal reconstruction
    p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end
end