%Written by Morteza Heydari Araghi

function [UserID, DAQsettings, tkdata]=importPAI128(filename)
     
  fid=fopen(filename,'r');
  
  version=fread(fid,1,'int16');
  UserID=fread(fid,1,'int16');
  
  DAQsettingssize=fread(fid,1,'int32');
  
%   % Pre-allocate
%   DAQsettings = cell(1,2);
  
  for ii=1:DAQsettingssize

      DAQsettings{ii}.NumPoints=fread(fid,1,'int32');
      DAQsettings{ii}.NumTriggers=fread(fid,1,'int32');  
      
      DAQsettings{ii}.NumChannels=fread(fid,1,'int32');         
      fread(fid,1,'int8');
      
      DAQsettings{ii}.TransferRate=fread(fid,1,'double')   ;
      DAQsettings{ii}.TimeStamp= fread(fid,1,'uint32');
      
      messagesize=fread(fid,1,'int32');
      message=fread(fid,messagesize,'uint8');
      
      datasize=fread(fid,1,'uint32')/2;
      
      DAQsettings{ii}.DataPointsCombined=fread(fid,datasize,'uint16');
      
      for jj=1:DAQsettings{ii}.NumChannels
         reshaped=reshape(DAQsettings{ii}.DataPointsCombined,datasize/DAQsettings{ii}.NumChannels,DAQsettings{ii}.NumChannels)';
         DAQsettings{ii}.DataPoints= swapbytes(uint16(reshaped(:,2:end)));              
      end
        
  end
  
  tkEnabled=fread(fid,1,'int8');
  
  if tkEnabled~=0
      tkdatasize=fread(fid,1,'int32');
      tkdata=fread(fid,6,'double');
  else
      tkdata=[];
  end
      

  fclose(fid);
  


  
  
end
