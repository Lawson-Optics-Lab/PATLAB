function pa_img = TimeReversal(inputData, sensorsPositions, scanArray, sensorArray,Npixel_x, Npixel_y,Npixel_z,x_size,y_size,z_size,samplingRate, speedOfSound, method )
%% Time Reversal image reconstruction
Nsample = size(inputData, 2);
Nstep   = length(sensorArray);
 chndata = inputData;
sensor.mask = sensorsPositions';
% total_step = Nstep; 
% angle_per_step = 2*pi/total_step; 
%calculate the negative delay caused by the system impulse response from the raw data file
x = x_size/1000;              % grid size [m]
y = y_size/1000;              % grid size [m]
z = z_size/1000;              % grid size [m]
Nx = Npixel_x;           % number of pixels in the x (column) direction
Ny = Npixel_y;
Nz = Npixel_z;           % number of pixels in the z (row) direction
dx = x/Nx;          % pixel width [m]
dy = y/Ny;          % pixel width [m]
dz = z/Nz;          % pixel height [m]
% 
SOS_full = ones(Ny,Nx)*speedOfSound/1000;
% % SOSmap = 50./(SOSimg/40+50/(waterSOS));
% % SOS_full(201:800,201:800)=SOSmap;
% % % load SOS_1020_15004.mat;                      %%%%%%%% speed of sound %%%%%%%%%%%%%
medium.sound_speed = SOS_full;
medium.density = 1000;    % [kg/m^3]
kgrid_recon = makeGrid(Ny, dy, Nx, dx);


% assign the time reversal data
% load rawleaf2.mat;
% load chndata30.mat
% chndata = mean(chndata_3D(:,:,40:41),3);
% denoise = 0;
% if denoise
%     chndata = Curvelet_denoising_US(chndata);
% else
%     [Xsize,Ysize] = size(chndata);
%     delay = ones(512,1);
%     delay(1:8:end) = 2;delay(2:8:end) = 2;delay(3:8:end) = 2;delay(4:8:end) = 2;
%     temp = zeros(Xsize,Ysize);
%     ret = zeros(Xsize,Ysize);
%     for ii = 1:512
%         if delay(ii) == 1
%             temp(1:Xsize,ii,:) = chndata(1:end,ii);
%         else
%             temp(1:(Xsize-1),ii,:) = chndata(2:end,ii);   
%         end
%     end
%     chndata = temp; clear temp;
% end
% chndata = fliplr(chndata);
% chndata1 = chndata;
% ini_angle = 225;
% shif = (ini_angle-180)/360*512;
% if shif >0
%     chndata1(:,1:(512-shif))= chndata(:,(shif+1):512);
%     chndata1(:,(512-shif+1):512)=chndata(:,1:shif);
% else
%     shif = -shif;
%     chndata1(:,(shif+1):512)= chndata(:,1:(512-shif));
%     chndata1(:,1:shif)=chndata(:,(512-shif+1):512);
% end
% chndata = chndata1;
% % break
% chndata = chndata';
% chndata = chndata(:,14:end);
Sc = size(chndata,2);
cart_sensor_mask = sensor.mask;
sensor_radius = 25e-3;      % [m]
sensor_angle = 4*pi/2;      % [rad]
sensor_radius_grid_points = round(sensor_radius/kgrid_recon.dx);
binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, kgrid_recon.Nx/2, kgrid_recon.Ny/2, sensor_radius_grid_points, sensor_angle);
sensor.mask = binary_sensor_mask;
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, chndata, cart_sensor_mask, binary_sensor_mask);
% sensor.time_reversal_boundary_data = chndata;%sensor_data;
% create the time array
% [kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
dt = 2.5e-8; %s
kgrid.t_array = ((1:Sc)-1)*2.5e-8;
% set the input options
input_args = {'Smooth', false, 'PMLInside', false,'DataCast', ...
    'gsingle','RecordMovie', false, 'PlotSim', false};
% attach the original time array
kgrid_recon.t_array = kgrid.t_array;
% reset the initial pressure
source.p0 = 0;
p0_recon = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor);
p0_recon = double(p0_recon);
% save(filename,'p0_recon');
figure;
imshow(p0_recon,[-1 1]);
axis image;
end
