function pa_img = PAT_Reconstruction_3D(P)
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
xRange          = P.x_range; %low to high
yRange          = P.y_range; %high to low
zRange          = P.z_range; %low to high
samplingRate    = P.fs;
speedOfSound    = P.V_M;
method          = P.recMethod;
coherentFactor  = P.coherentFactor;
DerGin = P.DerivativesGain;
% minimumVariance = P.minimumVariance;
dPow     = P.directivityPower;
v               = P.directivityAngle;
dEx     = P.directivityExample;
dMethod         = P.directivityMethpod;
delay           = P.delay; %micro second
slicing         = P.slicing;
delay = delay*1e-6; %second

% sensorArray     = P.Sensors;
%% calculate imaging area & parameters
tic
Npixel_x = size(xRange, 2);
Npixel_y = size(yRange, 2);
Npixel_z = size(zRange, 2);
if strcmp(slicing, 'xyz')
    x_img = ones(Npixel_y,1)*xRange;
    y_img = yRange'*ones(1,Npixel_x);
    z_img_temp = ones(Npixel_y,Npixel_x);
    z_img = zeros(Npixel_y,Npixel_x,Npixel_z);
    for i = 1:Npixel_z
        z_img(:,:,i) = z_img_temp*zRange(i);
    end
    pa_img = zeros(Npixel_y, Npixel_x, Npixel_z); % image buffer
elseif strcmp(slicing, ' yzx')
    x_img = ones(Npixel_z,1)*flip(yRange);
    y_img = flip(zRange)'*ones(1,Npixel_y);
    z_img_temp = ones(Npixel_z,Npixel_y);
    z_img = zeros(Npixel_z,Npixel_y,Npixel_x);
    for i = 1:Npixel_x
        z_img(:,:,i) = z_img_temp*xRange(i);
    end
    pa_img = zeros(Npixel_z, Npixel_y, Npixel_x); % image buffer
elseif strcmp(slicing, '  zxy')
    x_img = ones(Npixel_x,1)*zRange;
    y_img = flip(xRange)'*ones(1,Npixel_z);
    z_img_temp = ones(Npixel_x,Npixel_z);
    z_img = zeros(Npixel_x,Npixel_z,Npixel_y);
    yRange = flip(yRange);
    for i = 1:Npixel_y
        z_img(:,:,i) = z_img_temp*yRange(i);
    end
    pa_img = zeros(Npixel_x, Npixel_z, Npixel_y); % image buffer
end

dat =  459.6*exp(-0.1342*v)+30.51*exp(-0.02917*v);
%% Reconstruction
Nsample = size(inputData,2);
Nstep   = size(inputData,3);
f = waitbar(0,'1','Name','BackProjection Reconstruction');
% setappdata(f,'canceling',0);
iter = 0;

if size(sensorsPositions,2)==2; sensorsPositions = cat(2,sensorsPositions, squeeze(zeros(size(sensorsPositions,1), 1 ,size(sensorsPositions,3)))); end
% if dMethod~=0; total_angle_weight = zeros(size(pa_img)); end
total_angle_weight = zeros(size(pa_img));
if coherentFactor; coherent_square = zeros(size(pa_img)); end
% if minimumVariance;  temp_img_MV = zeros(size(pa_img)); end
pa_img0 = zeros(size(pa_img)); % image buffer

if method ==2 %FBP
    Sm = size(inputData,2);
    NFFT = 2^nextpow2(Sm);
    angles = linspace(1, 180, NFFT/2);
    n = length(angles);
    delta = angles(2) - angles(1);
    K = (-n:n-1) * delta;
    odd = rem(n,2);
    phi = zeros(2*n, 1);
    phi(2-odd:2:end) = -1 ./ (pi^2 * K(2-odd:2:end).^2);
    phi(n+1) = 1 / (4 * delta^2);
    phi = phi * delta^2;
    %         pa_data = conv2(pa_data, phi, "same");
elseif method ==3 %UBP
    fs = samplingRate;  %sampl freq
    Sm = size(inputData,2);
    NFFT = 2^nextpow2(Sm);
    fv = fs/2*linspace(0,1,NFFT/2+1);
    fv2 = -flip(fv);
    fv2(1) = [];
    fv2 (end)= [];
    fv3 = cat(2,fv,fv2);  %ramp filter for positive and negative freq components
    k = 2*pi*fv3/speedOfSound; %wave vector
    k = k';
    t = linspace(0,1/fs,Sm);
end
tic
for i = 1:size(inputData,1) %scans
    TP = sensorsPositions(:,:,i); %Transducers position
    if strcmp(slicing, 'xyz')
        x_receive = TP(:,1);
        y_receive = TP(:,2);
        z_receive = TP(:,3);
        if dMethod==1 %Center of the sensors locations for every scan
            x_center = mean(x_receive,'all');
            y_center = mean(y_receive,'all');
            z_center = mean(z_receive,'all');
        elseif dMethod==2 %Global focal point
            x_center = P.dir(1);
            y_center = P.dir(2);
            z_center = P.dir(3);
        elseif dMethod==3 %Global unit vector
            x_dir = P.dir(1);
            y_dir = P.dir(2);
            z_dir = P.dir(3);
        elseif dMethod==4 %Load focal points
            x_center = P.dir(:,1,:);
            y_center = P.dir(:,2,:);
            z_center = P.dir(:,3,:);
        elseif dMethod==5 %Load unit vectors
            x_dir = P.dir(:,1,:);
            y_dir = P.dir(:,2,:);
            z_dir = P.dir(:,3,:);
        end
    elseif strcmp(slicing, ' yzx')
        x_receive = TP(:,2);
        y_receive = TP(:,3);
        z_receive = TP(:,1);
        if dMethod==1 %Center of the sensors locations for every scan
            x_center = mean(x_receive,'all');
            y_center = mean(y_receive,'all');
            z_center = mean(z_receive,'all');
        elseif dMethod==2 %Global focal point
            x_center = P.dir(2);
            y_center = P.dir(3);
            z_center = P.dir(1);
        elseif dMethod==3 %Global unit vector
            x_dir = P.dir(2);
            y_dir = P.dir(3);
            z_dir = P.dir(1);
        elseif dMethod==4 %Load focal points
            x_center = P.dir(:,2,:);
            y_center = P.dir(:,3,:);
            z_center = P.dir(:,1,:);
        elseif dMethod==5 %Load unit vectors
            x_dir = P.dir(:,2,:);
            y_dir = P.dir(:,3,:);
            z_dir = P.dir(:,1,:);
        end
    elseif strcmp(slicing, '  zxy')
        x_receive = TP(:,3);
        y_receive = TP(:,1);
        z_receive = TP(:,2);
        if dMethod==1 %Center of the sensors locations for every scan
            x_center = mean(x_receive,'all');
            y_center = mean(y_receive,'all');
            z_center = mean(z_receive,'all');
        elseif dMethod==2 %Global focal point
            x_center = P.dir(3);
            y_center = P.dir(1);
            z_center = P.dir(2);
        elseif dMethod==3 %Global unit vector
            x_dir = P.dir(3);
            y_dir = P.dir(1);
            z_dir = P.dir(2);
        elseif dMethod==4 %Load focal points
            x_center = P.dir(:,3,:);
            y_center = P.dir(:,1,:);
            z_center = P.dir(:,2,:);
        elseif dMethod==5 %Load unit vectors
            x_dir = P.dir(:,3,:);
            y_dir = P.dir(:,1,:);
            z_dir = P.dir(:,2,:);
        end
    end
    
    pa_data = reshape(inputData(i,:,:),size(inputData,2), []);
    
    for iStep = 1:Nstep % sensors
        iter  = iter+1;
        waitbar(iter/(Nstep*size(inputData,1)),f,'Please Wait ...');
        pa_data_tmp = pa_data(:,iStep);
        if method ==2 %FBP
            pa_data_tmp = ifft(abs(fft(phi)).*fft(pa_data_tmp,NFFT));
            pa_data_tmp = real(pa_data_tmp(1:Sm));
        elseif method ==3 %UBP
            pf = real(ifft(1i*k.*fft(pa_data_tmp,NFFT))); %apply ramp filter
            pa_data_tmp = pa_data_tmp-DerGin.*speedOfSound*t'.*pf(1:Sm);
            %             pa_data_tmp = pa_data_tmp-speedOfSound*t'.*pf(1:Sm);
            %            pa_data_tmp =  pa_data_tmp-10.*diff([pa_data_tmp;0]);
        end
        pa_data_tmp = [pa_data_tmp;0];
        dx = x_img - x_receive(iStep);
        dy = y_img - y_receive(iStep);
        dz = z_img - z_receive(iStep);
        rr0_dis = sqrt(dx.^2+dy.^2); % in plane distance to the detector
        rr0_dis = repmat(rr0_dis,[1, 1, size(pa_img,3)]);
        rr2 = sqrt(rr0_dis.^2 + dz.^2);           % distance
        IDX = round((rr2/speedOfSound-delay)*samplingRate);
        inrange = (IDX >=50) & (IDX <= Nsample);
        IDX = (inrange).*IDX + (1-inrange).*(Nsample+1);
        if dMethod ~=0
            if dMethod==1 || dMethod==2 %Center of the sensors locations for every scan OR Global focal point
                a = sqrt((x_receive(iStep)-x_center)^2+(y_receive(iStep)-y_center)^2+(z_receive(iStep)-z_center)^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_center).^2+(y_img-y_center).^2+(z_img-z_center).^2);
                c = rr2;
                %                 d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
                d = sqrt(c.^2-((c.^2+a.^2-b.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            elseif dMethod==3 %Global unit vector
                a = sqrt((x_dir)^2+(y_dir)^2+(z_dir)^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_receive(iStep)-x_dir).^2+(y_img-y_receive(iStep)-y_dir).^2+(z_img-z_receive(iStep)-z_dir).^2);
                c = rr2;
                d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            elseif dMethod==4 %Load focal points
                a = sqrt((x_receive(iStep)-x_center(iStep,i))^2+(y_receive(iStep)-y_center(iStep,i))^2+(z_receive(iStep)-z_center(iStep,i))^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_center(iStep,i)).^2+(y_img-y_center(iStep,i)).^2+(z_img-z_center(iStep,i)).^2);
                c = rr2;
                %                 d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
                d = sqrt(c.^2-((c.^2+a.^2-b.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            elseif dMethod==5 %Load unit vectors
                a = sqrt((x_dir(iStep,i))^2+(y_dir(iStep,i))^2+(z_dir(iStep,i))^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_receive(iStep)-x_dir(iStep,i)).^2+(y_img-y_receive(iStep)-y_dir(iStep,i)).^2+(z_img-z_receive(iStep)-z_dir(iStep,i)).^2);
                c = rr2;
                d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            end
            angle_weight = real(cos(alpha).^dat)./abs(rr2).^dPow;
            if i== dEx(1) && iStep==dEx(2)
                weEx = (angle_weight-min(angle_weight(:)))./(max(angle_weight(:))-min(angle_weight(:)));
                if size(weEx,3)==1
                figure; imshow(weEx,[]); title('The sensitivity beam fot the selected sensor');
                end
            end
            total_angle_weight = total_angle_weight + angle_weight;
            temp_img = pa_data_tmp(IDX).*angle_weight;
        else
            angle_weight = 1./abs(rr2).^dPow;
            total_angle_weight = total_angle_weight + angle_weight;
            temp_img = pa_data_tmp(IDX).*angle_weight;
        end
        pa_img0 = pa_img0 + temp_img;
        if coherentFactor; coherent_square = coherent_square + temp_img.^2; end
    end
end
if coherentFactor
    coherent_factor = pa_img0.^2./coherent_square/iter;
    pa_img =  pa_img0.*coherent_factor;
else
    pa_img =  pa_img0;
end
% if dMethod~=0; pa_img = pa_img./total_angle_weight; end
pa_img = pa_img./total_angle_weight;
toc
delete(f)
end