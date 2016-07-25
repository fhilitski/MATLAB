N_trials = 1;
stds = zeros(1,N_trials);
means = zeros(1,N_trials);

for j = 1:N_trials
%     
%     if j <= N_trials/2
%         set_sd = 1;
%     else
%         set_sd = 1.2;
%     end;
    set_sd = 1;
    n = 1000;
    a = random('normal',pi,set_sd,n,1); 
    %a = random('unif',0,2*pi,n,1);
    x = 0:0.1:2*pi;
    pdf = exp(-(x-pi).^2/2);
    pdf = ones(length(x),1);
    norm = sum(pdf)*0.1;
    p = pdf/norm;
    figure(100);
    plot(x,p,'.r');
    hold on;
    [h, xb] = hist(a);
    nhist = sum(h)*(xb(2)-xb(1));
    plot(xb,h/nhist,'ob');
    hold off;
      
   
    
   % polar(x,p,'or');
    figure(101);
    plot(a,'or');
    figure(102);
    hist(a);
    axis tight;
    stds(j) = std(a);
    means(j) = mean(a);
    
    figure(500);
    [t,r]=rose(a,10);
    %rose(a);
     %polar(xb,h,'r');
     %hold on;
    rose(a,10);
    figure(600);
    
    nhist = 1/2*sum(h.^2)*(abs(mean(dxb)));
    normp = h/sqrt(nhist);
    polar(xb,normp);
    sum(1/2 * normp.^2 * abs(mean(dxb)))
    
    
end;

% figure(202);
% plot(means,'ob');
% figure(203);
% plot(stds,'ok');



