function ret = get_total_RSS_SHL(filepath)
    temp =read_bf_file(filepath);
    error(nargchk(1,1,nargin));
    flag = cellfun(@isempty,temp);
    k=1;
    for i=1:size(temp,1)
        if flag(i,1)==0
            csi_trace(k,1)=temp(i,1);
            k=k+1;
        end
    end
    k=1;
    total_rss=0;
    for i=1:size(temp,1)
        if flag(i,1)==0         
            rssi_mag = 0 ;
            if csi_trace(k.1).rssi_a ~=0
                rssi_mag=dbinv((csi_trace(k.1).rssi_a)+rssi_mag;
            end
            if csi_trace(k.1).rssi_b ~=0
                rssi_mag=dbinv((csi_trace(k.1).rssi_b)+rssi_mag;
            end
            if csi_trace(k.1).rssi_c ~=0
                rssi_mag=dbinv((csi_trace(k.1).rssi_c)+rssi_mag;
            end
            rssi_mag = db(rssi_mag,'pow')-44-csi_trace(k.1).agc;
            k=k+1;
        end
        total_rss=total_rss+rssi_mag;
    end
    ret = total_rss/size(temp,1);
end
            
    