% Determine location of transducers relative to imaging volume
% using gaussian newton to solve the problem (http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm)
%
% original code taken from fileexchange -- updated by phil, apr15/15

%% Initialize Parameters
P.step = 1;                                 
P.xsize = 30;                              
P.ysize = 30;                             
P.zsize = 30;                              
P.d = P.xsize*P.step;                       
P.e = P.ysize*P.step;
P.f = P.zsize*P.step;
P.d1 = P.xsize;                             
P.e1 = P.ysize;
P.f1 = P.zsize; 
P.Fs = 50e6;           % sampling rate (Hz)
P.C = 1500*1000;       % speed of sound in water [mm/s]

N = 27000;  % number of gridpoints
M = 24;  % number of transducers

noisePow = 0;  % it means that the accuracy of distance measurement is 90 %
                % for instance the inaccuracy of a 1m measured distance
                % is around .1 meter.

networkSize = 200;  % cubic volume (mm3) that the transducers can wander
v2struct(P)

%% Grid point positions
if f1 < 1
    [c,r,p] = meshgrid(-e/2:step:(e/2),-d/2:step:(d/2),(-(f-1)/2:step:((f-1)/2)));
    c = c(1:d*(1/step),1:e*(1/step),1:f); 
    r = r(1:d*(1/step),1:e*(1/step),1:f); 
    p = p(1:d*(1/step),1:e*(1/step),1:f);    
else
    [c,r,p] = meshgrid(-e/2:step:(e/2),-d/2:step:(d/2),(-f/2:step:(f/2)));
    c = c(1:d*(1/step),1:e*(1/step),1:f*(1/step)); 
    r = r(1:d*(1/step),1:e*(1/step),1:f*(1/step)); 
    p = p(1:d*(1/step),1:e*(1/step),1:f*(1/step));  
end

p = p-((min(p(:))+max(p(:)))/2);
r = r-((min(r(:))+max(r(:)))/2);
c = c - ((min(c(:))+max(c(:)))/2);
position = cat(4,c,r,p); 
IOpos = reshape(position,[],3);

% % Plot
% figure, plot3(IOpos(:,1),IOpos(:,2),IOpos(:,3),'bo')
% % hold on, plot3(Array(:,1),Array(:,2),Array(:,3),'ro')
% % hold on, plot3(Array(center,1),Array(center,2),Array(center,3),'k+')
% axis equal, set(gcf, 'color','w')
% xlim([-networkSize networkSize])
% ylim([-networkSize networkSize])
% zlim([-networkSize networkSize])
           
%% Transducer Localization
plane = 1;
iter = 5;
numrot = 10;
dimension = 3;
numOfIteration = 25;

ArcArray = zeros(M,dimension,length(plane),dimension,iter,numrot,'single');

tic
for i0 = 1:numrot
    % Load IO
    fname = sprintf('IOwaterps_30x30x30_rot%d.mat',i0);
%     load(fname);
    load(fname,'D');

    % TOA
%     [~,D] = max(io,[],2); clear io
%     D = squeeze(D);
%     save(fname,'D','-append');

    % Delay
    delay = 3100; % timepoints (60 us)
    dist0 = (D+delay).*(C*1/Fs);

    % Initial guess 
%     mobileLocEst = networkSize*rand(M,3); %random locations
    mobileLocEst = zeros(M,3); % origin
    mobileLocEst(:,3) = 100; % 100 mm above origin

    % imgbuffer = zeros(d,e,M,'single');

    for dim = 1:dimension
        for m = 1:iter
            for k = plane
                for j = 1:M
                    distance = reshape(dist0(:,j),d,e,f);
                    distance = medfilt3(distance);
% 
%                     if dim == 1
%                         distance = distance(k,:,:);
%                     elseif dim == 2
%                         distance = distance(:,k,:);
%                     elseif dim == 3
%                         distance = distance(:,:,k);
%                     end
                    
                %     imgbuffer(:,:,j) = distance;
                    distance = reshape(distance,[],1);

                    % noisy measurements
                    distance = distance + distance.*noisePow./100.*(rand(N,1)-1/2);

                    for i = 1:numOfIteration
                        % Esimated distance
%                         temp = reshape(IOpos,d,e,f,[]);

%                         if dim == 1
%                             temp = reshape(temp(k,:,:,:),[],dimension);
%                         elseif dim == 2
%                             temp = reshape(temp(:,k,:,:),[],dimension);
%                         elseif dim == 3
%                             temp = reshape(temp(:,:,k,:),[],dimension);
%                         end
                        
                        temp = IOpos;
                        distanceEst = sqrt(sum((temp-repmat(mobileLocEst(j,:),N,1)).^2,2));

                        % computing the derivatives
                            % d0 = sqrt( (x-x0)^2 + (y-y0)^2 )
                            % derivatives -> d(d0)/dx = (x-x0)/d0
                            % derivatives -> d(d0)/dy = (y-y0)/d0

                        distanceDrv   = [(mobileLocEst(j,1)-temp(1:N,1))./distanceEst ... % x-coordinate
                                         (mobileLocEst(j,2)-temp(1:N,2))./distanceEst ... % y-coordinate
                                         (mobileLocEst(j,3)-temp(1:N,3))./distanceEst]; % z-coordinate
                                                                         
                        % delta 
                        [u,s,v] = svd(distanceDrv,'econ');
                        pseudo = v(:,:)*pinv(s(:,:))*u(:,:)';    
                        delta = -pseudo*(distanceEst-reshape(distance,[],1));   
                        
                        % constrain
                        delta(isnan(delta)==1) = 0;
                        delta(abs(delta)>networkSize) = 0;  
                        
                        % Updating the estimation
                        mobileLocEst(j,:) = mobileLocEst(j,:) + delta.';
                    end
                    ArcArray(j,:,k,dim,m,i0) = mobileLocEst(j,:);
                end
            end
        end
    end
    disp(i0)
end
toc

% Save
% ArcArray_filt = ArcArray;
% save('ArcArrayPos.mat','ArcArray_filt','-append')

%% Visualize array
M = 24; 
networkSize = 200;
for i = 1:10
    mobileLocEst = test(:,:,i)*-1;
    figure, plot3(IOpos(:,1),IOpos(:,2),IOpos(:,3),'bo')
    view(0,90)
    hold on, plot3(mobileLocEst(1:M,1),mobileLocEst(1:M,2),mobileLocEst(1:M,3),'ro')
    axis equal, set(gcf, 'color','w')
    xlim([-networkSize networkSize])
    ylim([-networkSize networkSize])
    zlim([-networkSize networkSize])
end

% figure, plot3(IOpos(:,1),IOpos(:,2),IOpos(:,3),'bo')
% hold on, plot3(mobileLocEst(1:M,1),mobileLocEst(1:M,2),mobileLocEst(1:M,3),'ro')
% axis equal, set(gcf, 'color','w')
% xlim([-networkSize networkSize])
% ylim([-networkSize networkSize])
% zlim([-networkSize networkSize])
% 
% % % Compute the Root Mean Squred Error
% % Err = sqrt(sum((mobileLocEst-mobileLoc).^2));
% % title(['Estimation error is ',num2str(Err),'meter'])


%% Average positions across rotations
load('ArcArrayPos.mat','Array'); % mm

Array = permute(Array,[1 3 2]);

% Re-orient
Array = Array.*-1;

N = size(Array,1); % Number of projections

% Plot scenario
Array = reshape(Array,[],3);
figure, plot3(Array(:,1),Array(:,2),Array(:,3),'ro')
axis equal, set(gcf, 'color','w')
Array = reshape(Array,N,[],3);

% Fit circle to find center
center = zeros(N,3);
for i = 1:N
    temp = squeeze(Array(i,:,1:2));
    Par = CircleFitByTaubin(temp);
    center(i,1:2) = Par(1:2) ;
    center(i,3) = mean(Array(i,:,3),2);
end

% Center of rotation
xyz = center;
r0=mean(xyz);
xyz=bsxfun(@minus,xyz,r0);
[~,~,V]=svd(xyz,0);
 
% center of rotation line --> x(t) = x0+u0*t
u = V(:,1);
x0 = [mean(center(:,1),1),mean(center(:,2),1),0]; 
deg = 18;

% Rotate
ArrayNew = zeros(N,10,3);
for i = 1:10    
    [XYZnew, R, t] = AxelRot(squeeze(Array(:,i,:))', deg*(i-1), u, x0);
    ArrayNew(:,i,:) = XYZnew';    
end

% Plot scenario
ArrayNew = reshape(ArrayNew,[],3);
figure, plot3(ArrayNew(:,1),ArrayNew(:,2),ArrayNew(:,3),'ro')
axis equal, set(gcf, 'color','w')
ArrayNew = reshape(ArrayNew,N,10,[]);

ArrayTrue = squeeze(mean(ArrayNew,2));
cent0 = [mean(center(:,1),1) mean(center(:,2),1)]; cent0(1,3) = 0;

% Center at origin
ArrayTrue = bsxfun(@plus,ArrayTrue,-1*cent0); 

% Plot scenario
figure, plot3(ArrayTrue(:,1),ArrayTrue(:,2),ArrayTrue(:,3),'ro')
axis equal, set(gcf, 'color','w')

% Rotate and plot
ArrayNew = zeros(N,10,3);
for i = 1:10    
    [ArrayTrue2, R, t] = AxelRot(ArrayTrue', -deg*(i-1), u, x0);
    ArrayNew(:,i,:) = ArrayTrue2';    
end
ArrayNew = reshape(ArrayNew,[],3);
figure, plot3(ArrayNew(:,1),ArrayNew(:,2),ArrayNew(:,3),'ro')
axis equal, set(gcf, 'color','w')
ArrayNew = reshape(ArrayNew,N,10,[]);

ArrayTrue = permute(ArrayNew,[1 3 2])*-1;

% Save ArrayTrue
% save('ArcArrayPos.mat','ArrayTrue','-append')
