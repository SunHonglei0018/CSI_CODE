function [xymax,loc_max,xymin,smin] = extrema2(xy,t)


M = size(xy);
if length(M) ~= 2
 error('Entry must be a matrix.')
end
N = M(2);
M = M(1);

% Search peaks through columns:
[smaxcol,smincol] = extremos(xy);

% Search peaks through rows, on columns with extrema points:
im = unique([smaxcol(:,1);smincol(:,1)]); % Rows with column extrema
[smaxfil,sminfil] = extremos(xy(im,:).');

% Convertion from 2 to 1 index:
smaxcol = sub2ind([M,N],smaxcol(:,1),smaxcol(:,2));
smincol = sub2ind([M,N],smincol(:,1),smincol(:,2));
smaxfil = sub2ind([M,N],im(smaxfil(:,2)),smaxfil(:,1));
sminfil = sub2ind([M,N],im(sminfil(:,2)),sminfil(:,1));

% Peaks in rows and in columns:
smax = intersect(smaxcol,smaxfil);
smin = intersect(smincol,sminfil);

% Search peaks through diagonals?
% if nargin==1
 % Check peaks on down-up diagonal:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,1);

 % Check peaks on up-down diagonal:
 smax = intersect(smax,[M; (N*M-M); sextmax]);
 smin = intersect(smin,[M; (N*M-M); sextmin]);

 % Peaks on up-down diagonals:
 [iext,jext] = ind2sub([M,N],unique([smax;smin]));
 [sextmax,sextmin] = extremos_diag(iext,jext,xy,-1);

 % Peaks on columns, rows and diagonals:
 smax = intersect(smax,[1; N*M; sextmax]);
 smin = intersect(smin,[1; N*M; sextmin]);
% end

% Extrema points:
xymax = xy(smax);
xymin = xy(smin);

% Descending order:
[temp,inmax] = sort(-xymax); clear temp
xymax = xymax(inmax);
smax = smax(inmax);
%�ų�ĳЩС����
index = find(xymax > t);
xymax = xymax(index);
smax = smax(index);
[xymin,inmin] = sort(xymin);
smin = smin(inmin);
%ת��Ϊ��������
[row_max,col_max] = ind2sub(size(xy),smax);
loc_max = [row_max col_max];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [smax,smin] = extremos(matriz)
% Peaks through columns or rows.

smax = [];
smin = [];

for n = 1:length(matriz(1,:))
 [temp,imaxfil,temp,iminfil] = extrema(matriz(:,n)); clear temp
 if ~isempty(imaxfil)     % Maxima indexes
  imaxcol = repmat(n,length(imaxfil),1);
  smax = [smax; imaxfil imaxcol];
 end
 if ~isempty(iminfil)     % Minima indexes
  imincol = repmat(n,length(iminfil),1);
  smin = [smin; iminfil imincol];
 end
end
end

function [sextmax,sextmin] = extremos_diag(iext,jext,xy,A)
% Peaks through diagonals (down-up A=-1)

[M,N] = size(xy);
if A==-1
 iext = M-iext+1;
end
[iini,jini] = cruce(iext,jext,1,1);
[iini,jini] = ind2sub([M,N],unique(sub2ind([M,N],iini,jini)));
[ifin,jfin] = cruce(iini,jini,M,N);
sextmax = [];
sextmin = [];
for n = 1:length(iini)
 ises = iini(n):ifin(n);
 jses = jini(n):jfin(n);
 if A==-1
  ises = M-ises+1;
 end
 s = sub2ind([M,N],ises,jses);
 [temp,imax,temp,imin] = extrema(xy(s)); clear temp
 sextmax = [sextmax; s(imax)'];
 sextmin = [sextmin; s(imin)'];
end
end

function [i,j] = cruce(i0,j0,I,J)
% Indexes where the diagonal of the element io,jo crosses the left/superior
% (I=1,J=1) or right/inferior (I=M,J=N) side of an MxN matrix. 

arriba = 2*(I*J==1)-1;

si = (arriba*(j0-J) > arriba*(i0-I));
i = (I - (J+i0-j0)).*si + J+i0-j0;
j = (I+j0-i0-(J)).*si + J;

end
% Carlos Adri�n Vargas Aguilera. nubeobscura@hotmail.com