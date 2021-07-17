function [Transducer_Positions_Scanned] = transducer_position_scanning( Center,Transducer_Positions,Scan_Positions )
%TRANSDUCER_POSITION_SCANNING Summary of this function goes here
%   Detailed explanation goes here

delta_positions = zeros(length(Transducer_Positions),3);
for trans = 1:length(Transducer_Positions)
    x = Transducer_Positions(trans,1);
    y = Transducer_Positions(trans,2);
    z = Transducer_Positions(trans,3);
    xorig = Center(1);
    yorig = Center(2);
    zorig = Center(3);
    xdelta = x-xorig;
    ydelta = y-yorig;
    zdelta = z-zorig;
    delta_positions(trans,1) = xdelta;
    delta_positions(trans,2) = ydelta;
    delta_positions(trans,3) = zdelta;
end

Transducer_Positions_Scanned = zeros(length(Scan_Positions),length(Transducer_Positions),3);

for pos = 1:length(Scan_Positions)
    Transducer_Positions_Scanned(pos,:,:) = Scan_Positions(pos,:)+delta_positions;
end

figure;
for pos = 1:length(Scan_Positions)
    hold on
    scatter3(Transducer_Positions_Scanned(pos,:,1),Transducer_Positions_Scanned(pos,:,2),Transducer_Positions_Scanned(pos,:,3),'Filled','MarkerFaceColor',[0.643 0 0.427]);
end

end

