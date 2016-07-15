function [ filter ] = makeFilter( Radius, dim_x, dim_y,mode)
%makeFilter( Radius, dim_x, dim_y,mode) makes a Fourier Plane filter
%filter of dimensions (dim_x,dim_y) of a given Radius. See modes:
%   mode = 0 blocks central spot of a given Radius and lets everything else
%   through;
%   mode = 1 blocks everything except central spot of a given Radius;
filter=zeros(dim_x,dim_y);
Xc = dim_x /2;
Yc = dim_y / 2;

for i=1:dim_x
 for j=1:dim_y
     
        if (mode == 0 ) 
            if (((i-Xc)^2+(j-Yc)^2) >= Radius^2)
            filter(i,j)=1;
            end;
        elseif (mode == 1)
            if (((i-Xc)^2+(j-Yc)^2) <= Radius^2)
            filter(i,j)=1;
        end;

end;
 end;
end

