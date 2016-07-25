function[U,R]=Umbrella(Rexp,Rcal,bin_size)
% Umbrella sampling based on distributions of separations for beads 
% in experimental and calibration runs

U=[];

%Make sure bins for each distribution start on integer pixels
R_start_Exp=min(Rexp)-mod(min(Rexp),1);
R_end_Exp=max(Rexp)-mod(min(Rexp),1)+1;

R_start_Cal=min(Rcal)-mod(min(Rcal),1);
R_end_Cal=max(Rcal)-mod(min(Rcal),1)+1;

% make bins for plotting distributions of bead separation
Rb_Exp=[R_start_Exp:bin_size:R_end_Exp];
Rb_Cal=[R_start_Cal:bin_size:R_end_Cal];

%Normalizing distributions and plotting them
W_Exp=hist(Rexp,Rb_Exp);
Pdf_Exp=W_Exp./sum(W_Exp);

W_Cal=hist(Rcal,Rb_Cal);
Pdf_Cal=W_Cal./sum(W_Cal);


figure(gcf+1);
plot(Rb_Cal,Pdf_Cal,'--oc',Rb_Exp,Pdf_Exp,'--or');
axis tight;
legend('Calibration','Bundle Present');
title('Probability distribution of bead position');


Ustart=input('Enter Start Position: ');
Uend=input('Enter End Position: ');

be=find(Rb_Exp<Ustart+bin_size/2 & Rb_Exp>Ustart-bin_size/2);
le=find(Rb_Exp<Uend+bin_size/2 & Rb_Exp>Uend-bin_size/2);
bc=find(Rb_Cal<Ustart+bin_size/2 & Rb_Cal>Ustart-bin_size/2);
lc=find(Rb_Cal<Uend+bin_size/2 & Rb_Cal>Uend-bin_size/2);

Pc=Pdf_Cal(bc:lc);
Wc=W_Cal(bc:lc);
Pe=Pdf_Exp(be:le);
We=W_Exp(be:le);
R=Rb_Exp(be:le);

U(:,1)=(log(Pc./Pe));
U(:,2)=Wc;
U(:,3)=We;

figure(gcf+1);
plot(R,U(:,1),'bs');
hold on;
xlabel('Bead separation, \mu m');
ylabel('Potential due to overlap, kT');
title('Bundling Energy VS Bead Separation');

% this function plots binding energy vs bead separation, not overlap length.
% overlap = total mts length - separation, 
% slope of the plot is \lambda in the potential
% lambda >= 0 

%this finds weights for our points and records them into U(i,4) for point i

for i=1:length(U(:,1))
    w1=U(i,2);
    w2=U(i,3);
    if w1>w2
        U(i,4)=w2;
    else
        U(i,4)=w1;
    end
end;  
% perform weighted linear fit    
    fitobject = fit(R',U(:,1),'poly1','Weight', U(:,4));
    
% find coefficients of the fit
    fit_coeffs = coeffvalues(fitobject);
    lambda = fit_coeffs(1);
    c = fit_coeffs(2);
    
   h_fit_plot = plot(fitobject,'--r');
   
   text('String',['\lambda = ', num2str(lambda),' k_BT/\mu m'],'Position',[Ustart+0.1 1 0]);
   %text(0,Ustart,['\lambda = ', num2str(lambda)]);
   
        
    


    



    

end
