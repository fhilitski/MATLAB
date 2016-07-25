function[CL,coordsOut]=ContourLength(contour)
% Returns contour length of the filament CL from given contour. 
% Also, removes possible jumps in the filament shape, and returns smoother
% contour coordsOUT

 delta_x = contour(2:end,1)-contour(1:end-1,1);
 delta_y = contour(2:end,2)-contour(1:end-1,2);

 %find displacements between tracked points
 dS=sqrt(delta_x.^2+delta_y.^2);
 
 f=find(dS>4);
  if length(f)>0
     df = f(2:end) - f(1:end-1);
     tooth=find(df==1);
     contour(tooth,:)=[];
     
     dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);
     
     f=find(dS>4);
     
     if length(f)>0
         
         fb=find(f<=8);
         fb=f(max(fb));
         contour(1:fb,:)=[];
         dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);
         f=find(dS>4);
         if length(f)>0
             fe=find(f>=length(contour(:,1))-8);
             fe=f(min(fe));
             contour(fe:end,:)=[]; 
             dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);
         end
         f=find(dS>4);
         if length(f)>0
            contour(f+1,:)=[];
            dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);
         end
     end
     dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);
     f=find(dS>4);
     if length(f)==1
         L1=f(1);
         L2=length(dS)+1-L1;
         if L1<L2
             contour(1:f(1),:)=[];
         else
             contour(f(1):end,:)=[];
         end
        dS=sqrt((contour(2:end,1)-contour(1:end-1,1)).^2+(contour(2:end,2)-contour(1:end-1,2)).^2);

     end
 end
 
 CL(1)=0;
 for i=1:length(dS);
     CL(i+1)=CL(i)+dS(i);
 end;
 
 coordsOut=contour;

