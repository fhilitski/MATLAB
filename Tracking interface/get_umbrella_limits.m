function [StartPoint, EndPoint, start_index, end_index] = get_umbrella_limits (Range_Cal,Pdf_Cal,Range_Exp,Pdf_Exp,Threshold)
% Gets range for umbrella sampling from pdf and ranges for experment and
% calibration


%starting point determined by calibration
%rough estimate is
StartPoint = min(Range_Cal);

%endpoint is given by exp. run
%rough estimate is
EndPoint = max(Range_Exp);


%rough estimates are not enough, code errors out a lot, so a refinement
%follows. 04/05/2013
%starting point is determined by cal. run (with bundle);
max_cal_probability = max(Pdf_Cal);
startpoint_probability = Threshold * max_cal_probability;

i = 1;
while (Pdf_Cal(i) < startpoint_probability)
    StartPoint = Range_Cal(i);
    i = i + 1;
end;
start_index = i-1;

%end point is determined by exp. run (with bundle);
max_exp_probability = max(Pdf_Exp);
endpoint_probability = Threshold * max_exp_probability;

i = length(Range_Exp);
while (Pdf_Exp(i) < endpoint_probability)
    EndPoint = Range_Exp(i);
    i = i - 1;
end;
end_index = i+1;





