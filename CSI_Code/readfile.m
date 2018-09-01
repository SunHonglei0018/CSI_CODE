function csi_trace = readfile(filepath)
    % change the directory to matlab
% 	path('../linux-80211n-csitool-supplementary/matlab', path);
    % get things from filepath
	temp = read_bf_file(filepath);
    % use every elements of the temp as the parameter of the function
    % isempty
    flag=cellfun(@isempty,temp);
    k=1;
    for i=1:size(temp,1)
        if flag(i,1)==0
            csi_trace(k,1)=temp(i,1);
            k=k+1;
        end
    end
end
