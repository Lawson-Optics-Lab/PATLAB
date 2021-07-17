figure;
for k=1:Number_Transducers
    figure;
    for j=1:Scan_Points
        plotname = "Transducer " + int2str(k)+", Scan Point "+int2str(j);
        plot(Image_Data(j,:,k),"DisplayName",plotname)
        hold on
    end
end

for L = 31:40
    figure;
    for k=1:Scan_Points
        plot(Image_Data(k,:,L))
        hold on
    end
end