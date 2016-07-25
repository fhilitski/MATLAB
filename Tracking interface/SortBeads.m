function trajectories = SortBeads( fBeadCoord )
%SortBeads sorts tracked data into separate columns.
%INPUTS fBeadsCoors,output of BeadTrack function
%  Detailed explanation goes here

traj=[];

%Nm = input('Enter Nmax: ');

for i=1:2
    f=find(fBeadCoord(:,4)==(i));
    if i>1
        lmax=length(traj(:,1));
        dL=length(f)-lmax;
        if dL>0
        f=f(1:end-dL);
        end
        if dL<0
            traj=traj(1:end+dL,:);
        end
        
    end


c=i*2-1;
traj(:,c:c+1)=fBeadCoord(f,1:2);


end

trajectories = traj;
end

