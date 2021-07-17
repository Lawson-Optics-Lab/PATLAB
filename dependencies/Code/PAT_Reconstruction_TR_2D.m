function pa_img = PAT_Reconstruction_TR_2D(P)
%% Back projection image reconstruction
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
scanArray       = P.Scans;
sensorArray     = P.Sensors;
xRange          = P.x_range;
yRange          = P.y_range;
samplingRate    = P.fs;
speedOfSound    = P.V_M;
method          = P.recMethod;
matchedFilter   = 0;
coherentFactor  = P.coherentFactor;
minimumVariance = P.minimumVariance;
directivity     = P.directivity;
v               = P.directivityangle;
%% calculate imaging area & parameters
Npixel_x = size(xRange, 2);
Npixel_y = size(yRange, 2);
x_img = ones(Npixel_y,1)*xRange;
y_img = yRange'*ones(1,Npixel_x);
pa_img = zeros(Npixel_y, Npixel_x); % image buffer
lambda = 1000*speedOfSound / samplingRate; % Wavelength[mm]
dat =  459.6*exp(-0.1342*v)+30.51*exp(-0.02917*v);
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
delete(f)
end