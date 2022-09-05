function [ imdata ] = load_averaged_data( Averaging,scan_points_index,trans,samplesImg,numberoftransducers,AverageFluence,thresh,fluence_channel )
%LOAD_AVERAGED_DATA Load data from raw files, averaging if specified
%   Detailed explanation goes here

    if Averaging == 0
        parfor i = scan_points_index

            [S,~]=sprintf('%05d.pv3',i-1);    
            [~,DAQsettings,~]=importPAI128(S);
            RF=DAQ128settings2RF(DAQsettings);
            RF=(RF(trans,samplesImg)')-8192;
            %RF=abs(RF)
            %RF=-RF;
            RF(RF>8000) = 8000;
            %size(RF)
            RF = reshape(RF,[],1,length(trans));
            RF = squeeze(mean(RF,2));
            RF = detrend(RF,'constant');
            %threshold to ensure random noise doesn't enter
            RF((end-80):end,:)=0;
            RF(1,:)=0;  
            RF(abs(RF)<thresh) = 0;
            if AverageFluence ~= 1
                ShotPower = max(RF(samplesImg,fluence_channel));
                ShotRatio = AverageFluence/ShotPower;
                RFadjusted = RF.*ShotRatio;
                imdata(i,:,:) =RFadjusted(samplesImg,:);

            else
                imdata(i,:,:) =RF(samplesImg,:);
            end
        end
    else
        h = waitbar(0,'Please wait...');
        parfor j = scan_points_index
            
            k = ((j-1)*Averaging+1):((j-1)*Averaging+Averaging);
            imtempdata = zeros(Averaging,length(samplesImg),numberoftransducers,'single');
            waitbar(j / numel(scan_points_index))  
            for i = k
                [S,~]=sprintf('%05d.pv3',i-1);    
                [~,DAQsettings,~]=importPAI128(S);
                RF=DAQ128settings2RF(DAQsettings);
                RF=(RF(trans,samplesImg)')-8192;
                %RF=abs(RF)
                %RF=-RF;
                RF(RF>8000) = 8000;
                %size(RF)
                RF = reshape(RF,[],1,length(trans));
                RF = squeeze(mean(RF,2));
                RF = detrend(RF,'constant');
                %threshold to ensure random noise doesn't enter
                RF((end-80):end,:)=0;
                RF(1,:)=0;

                RF(abs(RF)<thresh) = 0;
                if AverageFluence ~= 1
                ShotPower = max(RF(samplesImg,fluence_channel));
                ShotRatio = AverageFluence/ShotPower;
                RFadjusted = RF.*ShotRatio;         

                imtempdata((i-((j-1)*Averaging)),:,:) = RFadjusted(samplesImg,:);
                else
                    imtempdata((i-((j-1)*Averaging)),:,:) = RF(samplesImg,:);
                end
            end
            
    %         imdata(j,:,:) = imtempdata(1,:,:);
            imdata(j,:,:) = mean(imtempdata,1);
        end
        delete(h);
    end

end

figure;
for q = 8
    for qd = 1:50
       plot (imdata(qd,:,q));
       hold on
    end
end
