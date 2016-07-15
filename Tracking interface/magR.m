function [r]=magR(traj,t1,t2)
%magR does not get no simpler than this!!!

r=sqrt((traj(:,t1)-traj(:,t2)).^2+(traj(:,t1+1)-traj(:,t2+1)).^2);

end
