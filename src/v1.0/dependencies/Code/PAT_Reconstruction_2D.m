function pa_img = PAT_Reconstruction_2D(P)
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
zRange          = P.z_range;
samplingRate    = P.fs;
speedOfSound    = P.V_M;
method          = P.recMethod;
matchedFilter   = 0;
coherentFactor  = P.coherentFactor;
minimumVariance = P.minimumVariance;
directivity     = P.directivity;
v               = P.directivityangle;
dMethod         = P.directivityMethpod;
%% 
if size(sensorsPositions,2)==2
    sensorsPositions(:,3,:)=zeros(size(sensorsPositions,1),1,size(sensorsPositions,3));
end
%% calculate imaging area & parameters
Npixel_x = size(xRange, 2);
Npixel_y = size(yRange, 2);
x_img = ones(Npixel_y,1)*xRange;
y_img = yRange'*ones(1,Npixel_x);
z_img = zRange*ones(Npixel_y,Npixel_x);
pa_img = zeros(Npixel_y, Npixel_x); % image buffer
lambda = 1000*speedOfSound / samplingRate; % Wavelength[mm]
dat =  459.6*exp(-0.1342*v)+30.51*exp(-0.02917*v);
%% Reconstruction
Nsample = size(inputData, 2);
Nstep   = length(sensorArray);
f = waitbar(0,'1','Name','BackProjection Reconstruction');
% setappdata(f,'canceling',0);
iter = 0;
if length(dMethod) ==3
    x_center = dMethod(1);
    y_center = dMethod(2);
    z_center = dMethod(3);
end
if directivity; total_angle_weight = zeros(Npixel_y, Npixel_x); end
for i = scanArray
    TP = sensorsPositions(sensorArray,:, i);%Transducers position
    x_receive = TP(:,1);
    y_receive = TP(:,2);
    z_receive = TP(:,3);
    if length(dMethod) ==1
        x_center = sum(x_receive,'all');
        y_center = sum(y_receive,'all');
        z_center = sum(z_receive,'all');
    end
    pa_data = reshape(inputData(i,:,sensorArray),size(inputData,2), []) ;
    if method ==2
        Sm = size(pa_data,1);
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
    elseif method ==3
        fs = samplingRate;  %sampl freq
            c = speedOfSound;
            Sm = size(pa_data,1);
            NFFT = 2^nextpow2(Sm);
            fv = fs/2*linspace(0,1,NFFT/2+1);
            fv2 = -flip(fv);
            fv2(1) = [];% delete(fv2,0);
            fv2 (end)= [];%np.delete(fv2,-1)
            fv3 = cat(2,fv,fv2);  %ramp filter for positive and negative freq components
            k = 2*pi*fv3/c; %wave vector
            k = k';
            t = linspace(0,1/fs,Sm);
    end
    
    pa_img0 = zeros(Npixel_y, Npixel_x); % image buffer
    
    if coherentFactor; coherent_square = zeros(Npixel_y, Npixel_x); end
    if minimumVariance;  temp_img_MV = zeros(Npixel_y, Npixel_x, Nstep); end
    if matchedFilter
        q=zeros(Npixel_x,Npixel_y);
        for ii=1:Npixel_x
            for jj=1:Npixel_y
                for kk=1:Npixel_z
                    p=zeros(1,Nstep);
                    for L=1:Ns
                        p(L)=sqrt((x_range_image(ii)-x_receive(i, L))^2+...
                            (y_range_image(jj)-y_receive(i, L))^2);
                    end
                end
                q(jj,ii, kk)= min(p');
            end
        end
    end
    delay = 0/samplingRate;
    tic
    for iStep = 1:Nstep
        iter  = iter+1;
        waitbar(iter/(Nstep*length(scanArray)),f,'Please Wait ...');
        %         pa_data_tmp=ifft(fft(pa_data(:,iStep)).*matchedfilter);
        %         pa_data_tmp = [pa_data; 0];
        pa_data_tmp = pa_data(:,iStep); 
        if method ==2
        %     pa_data_tmp= cumtrapz(pa_data_tmp')';
        %     t1 = diff(pa_data_tmp); % first derivative
        %     pa_data_tmp(2:end) = t1;
        %     pa_data_tmp= cumtrapz(pa_data_tmp')';
        %     pa_data_tmp = pa_data_tmp - pa_data_tmp2;
        pa_data_tmp = ifft(abs(fft(phi)).*fft(pa_data_tmp,NFFT));
        pa_data_tmp = real(pa_data_tmp(1:Sm));
        
        elseif method ==3
            pf = real(ifft(-1i*k.*fft(pa_data_tmp,NFFT))); %apply ramp filter
%             figure; ax1=subplot(311);plot(pa_data_tmp); ax2=subplot(312);plot(c*t'.*pf(1:Sm)); 
%             ax3=subplot(313);plot(pa_data_tmp-c*t'.*pf(1:Sm));linkaxes([ax1 ax2 ax3],'xy')
            pa_data_tmp = pa_data_tmp-c*t'.*pf(1:Sm);
        end
        dx = x_img - x_receive(iStep);
        dy = y_img - y_receive(iStep);
        dz = z_img - z_receive(iStep);
        dr = sqrt(dx.^2+dy.^2+dz.^2); % 
        IDX = round((dr/speedOfSound-delay)*samplingRate);
        inrange = (IDX >=50) & (IDX <= Nsample);
        IDX = (inrange).*IDX + (1-inrange).*(Nsample+1);
        if directivity
            if length(dMethod)==1 || length(dMethod)==3
                a = sqrt((x_receive(iStep)-x_center)^2+(y_receive(iStep)-y_center)^2+(z_receive(iStep)-z_center)^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_center).^2+(y_img-y_center).^2+(z_img-z_center).^2);
                c = dr;
%                 d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
d = sqrt(c.^2-((c.^2+a.^2-b.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            elseif length(dMethod) ==4
                a = sqrt((dMethod(1))^2+(dMethod(2))^2+(dMethod(3))^2);
                a = repmat(a,size(x_img));
                b = sqrt((x_img-x_receive(iStep)-dMethod(1)).^2+(y_img-y_receive(iStep)-dMethod(2)).^2+(z_img-z_receive(iStep)-dMethod(3)).^2);
                c = dr;
                d = sqrt(b.^2-((b.^2+a.^2-c.^2)./(2.*a)).^2);
                alpha = asin(d./c);
            end
            %             r0 = sqrt(x_receive(iStep)^2+y_receive(iStep)^2); %
            %    alpha2 = acos(min(abs((x_receive(iStep)*dx+y_receive(iStep)*dy+z_receiveI*dz)/r0_3D./rr2),0.999))+1.0e-5; % 3D angle
            %             alpha2 = acos(min(abs((x_receive(iStep)*dx + y_receive(iStep)*dy)/r0./dr), 0.999)) + 1.0e-5; % 3D angle
%             angle_weight = real(cos(alpha).^dat)./dr.^2;
            angle_weight = (real(cos(alpha).^dat)).^10;
            if iter ==1
                figure; imshow(angle_weight,[]); title('The sensitivity beam fot a single sensor');
            end
            total_angle_weight = total_angle_weight + angle_weight;
            temp_img = pa_data_tmp(IDX).*angle_weight;
        else
            temp_img = pa_data_tmp(IDX);
        end
        pa_img0 = pa_img0 + temp_img;
    end
    if directivity
%         pa_img = pa_img + pa_img0./total_angle_weight;
        pa_img = pa_img + pa_img0;
    else
        pa_img = pa_img + pa_img0;
    end
end   
if directivity
    pa_img = pa_img + pa_img0./total_angle_weight;
    %         pa_img = pa_img + pa_img0;
else
    pa_img = pa_img + pa_img0;
end
%         if method ==0
%             
%         elseif method ==1
%             
%         elseif method ==2
%             pa_img0 = pa_img0 + temp_img;
%             %             pa_img0 = pa_img0 + temp_img;
%             %             coherent_square = coherent_square + temp_img.^2;
%         elseif method ==3
%             warning('off')
%             temp_img_MV(:,:,:,iStep) = pa_data_tmp(IDX);
%         elseif method ==4
%             pa_img0 = pa_img0 + temp_img;
%             coherent_square = coherent_square + temp_img.^2;
%             temp_img_MV(:,:,:,iStep) = pa_data_tmp(IDX);
%         elseif method ==5
%             Ns = 77;                 % Number of Sensors
%             Ls = 60;                 % Lemgth of Sensor Array plane (Ls=<x)
%             pitch = Ls/Ns;           % distance between sensors in array [mm]
%             ws = 3/4 * pitch;
%             alpha = acosd(min(abs(q./rr),0.999))+1.0e-5;
%             df = (sin((pi*ws/lambda)*sind(alpha)))/(pi*ws/lambda*sind(alpha));
%             Df = (sin(pi*Ls/lambda*sind(alpha)))/(sin(pi*pitch/lambda*sind(alpha)));
%             total_D(:,:,i) = total_D(:,:,i) + df(:,:,iStep).*Df(:,:,iStep);
%             temp_img_2(:,:,iStep,i) = pa_data_tmp(idx).*df(:,:,iStep).*Df(:,:,iStep);
%             pa_img_2(:,:,i) = pa_img_2(:,:,i) + temp_img_2(:,:,iStep,i);
%         end      
%     if method == 0
%         pa_img = pa_img + pa_img0;
%     elseif method == 1
%         pa_img = pa_img + pa_img0./total_angle_weight;
%     elseif method == 2
%         pa_img = pa_img + pa_img0;
%         %         coherent_factor = pa_img0.^2./coherent_square/Nstep;
%         %         pa_img = pa_img + pa_img0.*coherent_factor;
%     elseif method == 3
%         L = 60;
%         cov_temp_all= zeros(L, L);
%         aa = ones(L, 1);
%         for ii=1:Npixel_y
%             for jj = 1:Npixel_x
%                 for kk = 1:Npixel_z
%                     for LStep_1 =1:Nstep-L+1
%                         temp_img_all =[];
%                         for LStep_2 = LStep_1:LStep_1+L-1
%                             temp_img_all = cat(1,temp_img_all ,temp_img_MV(ii,jj,kk,LStep_2));
%                         end
%                         mean_temp_all = (1/L)*sum(temp_img_all,1);
%                         temp_img_all_mines_mean = temp_img_all - repmat(mean_temp_all,L,1);
%                         cov_temp_all=cov_temp_all+temp_img_all_mines_mean*temp_img_all_mines_mean';
%                     end
%                     cov_temp_all=(1/(Nstep-L+1))*cov_temp_all;
%                     % calculation of Covariance Matrix
%                     W_MV= inv(cov_temp_all)*aa/(aa'*inv(cov_temp_all)*aa);
%                     pa_img_MV (ii, jj, kk)=(1/L)*W_MV'*temp_img_all;
%                 end
%             end
%         end
%         pa_img = pa_img +pa_img_MV;
%     elseif method == 4
%         coherent_factor = pa_img0.^2./coherent_square/Nstep;
%         L = 7;
%         cov_temp_all= zeros(L, L);
%         aa = ones(L, 1);
%         for ii=1:Npixel_y
%             for jj = 1:Npixel_x
%                 for kk = 1:Npixel_z
%                     for LStep_1 =1:Nstep-L+1
%                         temp_img_all =[];
%                         for LStep_2 = LStep_1:LStep_1+L-1
%                             temp_img_all = cat(1,temp_img_all ,temp_img_MV(ii,jj,kk,LStep_2));
%                         end
%                         mean_temp_all = (1/L)*sum(temp_img_all,1);
%                         temp_img_all_mines_mean = temp_img_all - repmat(mean_temp_all,L,1);
%                         cov_temp_all=cov_temp_all+temp_img_all_mines_mean*temp_img_all_mines_mean';
%                     end
%                     cov_temp_all=(1/(Nstep-L+1))*cov_temp_all;
%                     % calculation of Covariance Matrix
%                     W_MV= inv(cov_temp_all)*aa/(aa'*inv(cov_temp_all)*aa);
%                     pa_img_MV (ii, jj, kk)=(1/L)*W_MV'*temp_img_all;
%                 end
%             end
%         end
%         pa_img = pa_img +pa_img_MV.*coherent_factor;
%     elseif method ==5
%         pa_img_final(:,:,i) = pa_img_2(:,:,i)./total_D(:,:,i);
%     end
delete(f)
end