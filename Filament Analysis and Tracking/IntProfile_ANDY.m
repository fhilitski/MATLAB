function[IntPr]=IntProfile(contour, Img, w)
%Get Intensity Profile of filament

%Determine max dimensions of the image
s=size(Img);
xMax=s(2);
yMax=s(1);

%Filament length
L = length(contour);

for i = 1:L
    
    %Get i-th (X,Y) coordinates
    XY = contour(i,1:2);
    
    X = round(XY(1));
    if (X >= 1+w ) && (X <= xMax-w)
        X1=X-w;
        X2=X+w;
    else
        if X>=xMax-w
            X1=xMax;X2=X1;
            
        else
            X1=1;X2=1;
        end
    end
    
    Y=round(XY(2));
    if (Y >= 1+w) && (Y <= yMax-w)
        Y1=Y-w;
        Y2=Y+w;
    else
        if Y>=yMax-w
            Y1 = yMax;
            Y2 = yMax;
        else
            Y1=1;
            Y2=1;
        end
    end
       
    %Get Intensity Value
    subImg=Img(Y1:Y2,X1:X2);
    IntPr(i)=mean(mean(subImg)); %#ok<AGROW>
    
end
IntPr=IntPr-min(IntPr);
end