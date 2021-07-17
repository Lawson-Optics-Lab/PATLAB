% Loads in experimental IO and re-orders it
%
% updated by phil, Mar9/15

function [IO,rawdata] = GenerateIO_daq128(P)

% Initialize
v2struct(P)

% Load IO
if IOtype == 0 % Experimental IO
    
    % Add IO file path
    addpath(IOfpath)
    
    calibration = zeros(length(vpositions),length(samples),length(trans),'single');
    tic
    k = 1;
    for i = vpositions
        [S,ERRMSG]=sprintf('%05d.pv3',i-1);    
        [~,DAQsettings,~]=importPAI128(S);
        RF=DAQ128settings2RF(DAQsettings);
        RF=(RF(trans,1:tpoints*256*numavg)')-8192;
        RF(RF>8000) = 8000;
        
        RF = reshape(RF,[],numavg,length(trans));
        RF = squeeze(mean(RF,2));
        
        calibration(k,:,:) = RF(samples,:);
        k = k+1;
    end
    rawdata.calibration = calibration;
    toc
    disp('Raw IO loaded in.')
    
    % Filter signals
    if IOfilt == 1
        tic
        for i = 1:length(trans)
            for j = 1:size(calibration,1)
                temp = calibration(j,:,i);
                calibration(j,:,i) = wvd_denoise(temp,P);
            end
        end
        toc
        disp('Filtered image data.')
    end
       
    % Re-order
    if unsnake == 1
        calibration = reshape(calibration,[],length(trans)*length(samples));
        [~,idx] = ReOrder(P);
        calibration = calibration(idx,:);
        calibration = reshape(calibration,[],length(trans)*length(samples));
        IO = calibration; clear calibration
        disp('Imaging operator reordered.')
    end
    
    % Rectify
    if rect == 1
       IO = abs(IO); 
    end
    
    % Sensitivity map per transducer
    if sensmap == 1
        SensMap = zeros(vpositions,length(trans),'single');
        for i = 1:length(trans)
            ct = sum(IO(:,:,i).*IO(:,:,i),2);   
            SensMap(:,i) = ct;
        end
        rawdata.SensMap = SensMap;
    end
    
    % TOA
    if TOA == 1
        M = squeeze(max(IO,[],2));

        for i = 1:length(vpositions)
            for j = 1:length(trans)
                a0 = IO(i,:,j);
                a0(a0<M(i,j)) = 0;        
                IO(i,:,j) = a0;
            end
        end
    end
    
    % Remove bad trans
    if badtranschannel > 0
        for b = 1:length(badtranschannel)           
            IO(:,:,trans(trans == badtranschannel(b))) = 0;
        end                
    end
    
    % Remove IOpath
    rmpath(IOfpath)
    
    IO = IO(:,P.samples,:);
    IO = reshape(IO,[],P.samples*length(trans));
    
    disp('Experimental IO loaded in.')        
    
else % Analytical IO
    
    if exist(P.IOfpath,'file') == 2           
        tic
        load(P.IOfpath)     
        toc
        disp('Analytical IO loaded in.')
    else
        IO = CreateAnalyticalIO(P);
        IO = IO(:,P.samples);
    end    
    rawdata = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-order IO
%     if unsnake == 1
%         tic
%         IO = zeros(length(vpositions),length(samples),length(trans),'single');
%         for i = 1:length(trans)
%             IOraw = calibration(:,:,i);
% 
%             IOraw = reshape(IOraw,xsize,ysize,zsize,[]);    
%             IOraw(:,2:2:end,:,:) = IOraw(end:-1:1,2:2:end,:,:);
%             
%             for k = 2:2:zsize
%                 for j = 1:size(IOraw,4)
%                     test = fliplr(squeeze(IOraw(:,:,k,j)));
%                     IOraw(:,:,k,j) = test;
%                 end
%             end
% 
%             IOraw = reshape(IOraw,xsize*ysize*zsize,[]);
% 
%             IO(:,:,i) = IOraw;   
%         end
%         toc
%         disp('Imaging operator reordered.')
%     end