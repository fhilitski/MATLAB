function [ err_code, frames ] = make_animated_plot(t, x, style, title_str, x_string, y_string, block_size)
%MAKE_ANIMATED_PLOT makes a animated movie of x-vs-t plot
%MAKE_ANIMATED_PLOT(t, x, style, title_str, x_string, y_string)
%RETURN:
%err_code - the error code:
%   0 - normal operation;
%   1 - normal operation but video not saved (cancelled by the user);
%frames - structure of getframe results that can be played with
%movie(frames) command.
%INPUTS:
%t,x - the vectors for X vs T plot
%style - a string of LineSpec used for PLOT function
%title_str - plot title
%x_string - x-axis label
%y_string - y-axis label
%block_size - (optional) number of points that will be plotted in a single
%               frame; convenient if the data set is very large.
%(c) Feodor Hilitski 2015-2016 

proceed = true;
err_code = 0;

if nargin < 6
    warning('Not enough arguments!');
    proceed = false;
elseif nargin == 6
    block_size = 1;
end;

%frames is a returned movie frames
%get length of both x,t
l = length(x);
m = length(t);

if (l ~= m)
    msgbox('Unequal length of input vectors!','Plot Animator');
    proceed = false;
    err_code = -1;
    frames = struct;
    return;
else
    %pre-allocate memory for frames
    num_frames = ceil(l/block_size); %total number of frames after blocking
    rem_frames = mod(l,block_size); %remainder of frames not fitting in any block
    frames(num_frames) = struct('cdata',[],'colormap',[]);
    %frames(num_frames) = struct();
    
    %prepare the figure window
    h_fig = figure;
    plot(t, x, '.','Color',[77 190 238]./255);
    title(title_str,'FontSize',30,'FontName','Calibri');
    xlabel(x_string,'FontSize',25,'FontName','Calibri');
    ylabel(y_string,'FontSize',25,'FontName','Calibri');
    set(gca,'FontSize',25,'FontName','Calibri');
    hold on;
    
    %set axes limits
    xlim([min(t) max(t)]);
    ylim([min(x) max(x)]);
    
    block_start = 1;
    block_end = block_start + block_size - 1;
    frame_index = 1;
    while (block_start <= l)
       
        plot(t(block_start:block_end), x(block_start:block_end), style);
        frames(frame_index) = getframe(h_fig);
        block_start = block_end+1;
        block_end = block_start + block_size - 1;  
        disp(['Animating frames: ' num2str(block_start)]);
        if block_end > l
            block_end = l;
        end;
        frame_index = frame_index + 1;
    end; %frames for loop
    
    %close the helper figure
    close(h_fig);
    
    %display dialog to save the movie in a file
    dlg_title = 'Select file to save the movie...';
    fname = 'movie.avi';
    [fname,pathname,filter_index] = uiputfile('*.avi',dlg_title,fname);
    str_message = 'Done animating!';
    if (~strcmp('0',num2str(fname)))
        video_writer = VideoWriter([pathname fname],'Uncompressed AVI');
        open(video_writer);
        writeVideo(video_writer,frames);
        close(video_writer);
    else
        str_message = ['Video not saved! ' str_message];
        err_code = 1;
    end;
    msgbox(str_message,'Plot Animator');
    err_code = 0;
end %if

