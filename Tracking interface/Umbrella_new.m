function[U,R] = Umbrella_new(Rexp,Rcal,bin_size)
% Umbrella sampling based on distributions of separations for beads
% in experimental and calibration runs

U=[];

%Make sure bins for each distribution start on integer pixels
R_start_Exp=min(Rexp);%-mod(min(Rexp),1);
R_end_Exp = max(Rexp);%-mod(min(Rexp),1)+1;

R_start_Cal=min(Rcal);%-mod(min(Rcal),1);
R_end_Cal=max(Rcal);%-mod(min(Rcal),1)+1;

% make bins for plotting distributions of bead separation
Rb_Exp=[R_start_Exp:bin_size:R_end_Exp];
Rb_Cal=[R_start_Cal:bin_size:R_end_Cal];

%Normalizing distributions and plotting them
W_Exp=hist(Rexp,Rb_Exp);
Total = sum(W_Exp);
Pdf_Exp = W_Exp./ Total;

W_Cal=hist(Rcal,Rb_Cal);
Total_cal = sum(W_Cal);
Pdf_Cal = W_Cal./ Total_cal;


figure;
subplot(2,1,1);
plot(Rb_Cal,Pdf_Cal,'--oc',Rb_Exp,Pdf_Exp,'--or');
axis tight;
legend('Calibration','Bundle Present');
title('Probability distribution of bead position');
h_pdf_plot = gca;


%Also plot time series of calibration and experiment
subplot(2,1,2);
plot(Rexp,'.r');
hold on;
plot(Rcal,'.b');
hold off;
axis tight;
legend('Experiment','Calibraiton');
title('Time-series of bead separations');



%This part prompts user to manually enter range for umbrella sampling
start_string = 'Enter Start Position: ';
end_string = 'Enter End Position: ';
prompt = {start_string, end_string};
dialog_title = 'Umbrella Sampling Options';

user_output = inputdlg(prompt,dialog_title);
if (~isempty(user_output)) %only go ahead with sampling if Ok was pressed
    
    Ustart = str2double(user_output{1});
    Uend = str2double(user_output{2});
    
    %if Ustart and Uend are empty in the end, perform automatic search
    if (isnan(Ustart) || isnan(Uend))
     Threshold = 0.1;
     [Ustart, Uend,Index_s,Index_e] = get_umbrella_limits(Rb_Cal,Pdf_Cal,Rb_Exp,Pdf_Exp,Threshold);
    
    %if Ustart and Uend were given, find indeces Index_s and Index_n to
    %plot points on the graph;
    %start index is determined by ranges Rb_Cal 
    %end index is determined by Rb_Exp
    else
       Index_s = find(Rb_Cal <= Ustart,1,'last');
       Index_e = find(Rb_Exp >= Uend,1,'first');      
        
    end;
    
    
    %plot start- and end- points on the graph
    axes(h_pdf_plot);
    hold on;
    plot(h_pdf_plot, Ustart,Pdf_Cal(Index_s),'og','MarkerSize',15);
    plot(h_pdf_plot, Uend,Pdf_Exp(Index_e),'or','MarkerSize',15);
    hold off;
    
    be=find(Rb_Exp<Ustart+bin_size/2 & Rb_Exp>Ustart-bin_size/2);
    le=find(Rb_Exp<Uend+bin_size/2 & Rb_Exp>Uend-bin_size/2);
    bc=find(Rb_Cal<Ustart+bin_size/2 & Rb_Cal>Ustart-bin_size/2);
    lc=find(Rb_Cal<Uend+bin_size/2 & Rb_Cal>Uend-bin_size/2);
    
    Pc=Pdf_Cal(bc:lc);
    Wc=W_Cal(bc:lc);
    Pe=Pdf_Exp(be:le);
    We=W_Exp(be:le);
    R=Rb_Exp(be:le);
    
    
    %To eliminate possible index mismatch error
    %caused by binning
    %we truncate length of umbrella sampling 
    %and make dimensions for experiment and calibration same
    
    min_l= min(length(Pc),length(Pe));
    Pc = Pc(1:min_l);
    Pe = Pe(1:min_l);
    Wc = Wc(1:min_l);
    We = We(1:min_l);
    
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
    U(:,4) = min(U(:,2),U(:,3));
    
    % perform weighted linear fit
    fitobject = fit(R',U(:,1),'poly1','Weight', U(:,4));
    
    % find coefficients of the fit
    fit_coeffs = coeffvalues(fitobject);
    lambda = fit_coeffs(1);
    c = fit_coeffs(2);
    
    h_fit_plot = plot(fitobject,'--r');
    
    text('String',['\lambda = ', num2str(lambda),' k_BT/\mu m'],'Position',[Ustart+0.05 U(1:1) 0]);
    %text(0,Ustart,['\lambda = ', num2str(lambda)]);
else
    U=0;
    R=0;
end;











end
