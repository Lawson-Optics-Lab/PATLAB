%Written by Morteza Heydari Araghi
%updated by Phil, Feb11/15

function RF=DAQ128settings2RF(DAQ128settings)

%    RF=zeros(numel(DAQ128settings)*DAQ128settings{1}.NumChannels,DAQ128settings{1}.NumPoints*256);
   RF=zeros(numel(DAQ128settings)*DAQ128settings{1}.NumChannels,DAQ128settings{1}.NumPoints*256*DAQ128settings{1}.NumTriggers,'single');

   for ii=1:numel(DAQ128settings)
       basen=((ii-1)*DAQ128settings{ii}.NumChannels);

       if DAQ128settings{ii}.NumTriggers > 1
           temp = DAQ128settings{ii}.DataPoints;
%            temp = reshape(temp,DAQ128settings{ii}.NumChannels,[],DAQ128settings{ii}.NumTriggers);
%            temp = mean(temp,3);
       else
           temp = DAQ128settings{ii}.DataPoints;
       end
       
%        RF( (basen+1):(basen+DAQ128settings{ii}.NumChannels) , : )= temp;
       RF( (basen+1):(basen+DAQ128settings{ii}.NumChannels) , : )= temp;
      
   end


end
