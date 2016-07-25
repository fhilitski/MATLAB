function[Vrod]=VrodTang(lb,b,d,ls,a)

%calculate the interaction energy/length between two cylnders with a center
%to center distance d, lb is bjerrum length, b is distance per charge on a
%filament, ls is the screening length, a is the filament radius

Vrod = zeros(1,length(d));

for i=1:length(d)
r = d(i);    
theta=0:0.001:2*pi;
rn=sqrt(a^2*sin(theta).^2+(r-a*cos(theta)).^2);
V = besselk(0,rn./ls);
Vrod(i)=sum(V)*.001;
Vrod(i)=4*pi*Vrod(i)*ls*lb/a/b^2/(2*pi)^2/besselk(1,a/ls)*1000;
end;


end