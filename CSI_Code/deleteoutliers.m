function [b,idx,outliers] = deleteoutliers(a,alpha,rep);

if nargin == 1
	alpha = 0.05;
	rep = 0;
elseif nargin == 2
	rep = 0;
elseif nargin == 3
	if ~ismember(rep,[0 1])
		error('Please enter a 1 or a 0 for optional argument rep.')
	end
elseif nargin > 3
	error('Requires 1,2, or 3 input arguments.');
end

if isempty(alpha)
	alpha = 0.05;
end

b = a;
b(isinf(a)) = NaN;
%Delete outliers:
outlier = 1;
while outlier
	tmp = b(~isnan(b));
	meanval = mean(tmp);
	maxval = tmp(find(abs(tmp-mean(tmp))==max(abs(tmp-mean(tmp)))));
	maxval = maxval(1);
	sdval = std(tmp);
	tn = abs((maxval-meanval)/sdval);
	critval = zcritical(alpha,length(tmp));
	outlier = tn > critval;
	if outlier
		tmp = find(a == maxval);
		b(tmp) = NaN;
    end
end
if nargout >= 2
	idx = find(isnan(b));
end
if nargout > 2
	outliers = a(idx);
end
if ~rep
	%b=b(~isnan(b));
    b=b(~any(isnan(b), 2), :);
end
return
end
function zcrit = zcritical(alpha,n)
%ZCRIT = ZCRITICAL(ALPHA,N)
% Computes the critical z value for rejecting outliers (GRUBBS TEST)
tcrit = tinv(alpha/(2*n),n-2);
zcrit = (n-1)/sqrt(n)*(sqrt(tcrit^2/(n-2+tcrit^2)));
end