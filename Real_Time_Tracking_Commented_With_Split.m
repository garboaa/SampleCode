%%%%%%%%%%%%%%                   Real- Time Tracking                  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     Edits to existing codes started on 3/22/17 (AG)    %%%%%%%%%%%%%%%%%
clear; close all;clc;
xlimits = [100 1000];
ylimits = [0 900];

winsize = 48;    % size of smoothing mask, # pixels
maxthresh = 0.6; % threshold for local maxima, relative to brightest pixel

pointsize = 6;   % size of plotted points on output display
linewidth = 2;   % width of plotted lines
fontsize=11;     % height of text on plots (pt)

% relative positions of reference points within ROI (between 0 and 1)
refpts = [1/3 1/2 2/3];

% form smoothing window
h = hann(winsize)*hann(winsize)';

%Edit made on 4/5/17 (AG)
%User Input will specify if data is from a video file (1) or if it is being
%collected in real-time (2)
res = menu('Are you using a video file or collecting live data?','Video File','Real-Time Tracking');

%Video File Tracking (AG)
if res == 1
    
    % Change speaker and prod variables according to which speaker and prod
    % iterations are being analyzed (AG)
    % IMPORTANT NOTE: Create folders named All variables, ar, Statistics and Displacement Figures in current directory  folder before running script
    speaker = 8;
    prod = {
        '_ar1',...
        '_ar2',...
        '_ar3',...
        '_ar4',...
        '_ar5',...
        '_ar6',...
        '_ar7',...
        '_ar8',...
        '_ar9',...
        '_ar10',...
        '_ar11',...
        '_ar12',...
        '_ar13',...
        '_ar14',...
        '_ar15',...
        '_ar16',...
        '_ar17',...
        '_ar18'};
    n = 1;
    tic
    for i = 1 : length(prod)
                
        file = ['SCI' num2str(speaker) cell2mat(prod(i)) '_r'];
        
        % if running from "Real-Time Tracking" folder on Box
        filenamebaseimages = ['../../Normals Data/Participant SCI' num2str(speaker) ...
                  '/SCI' num2str(speaker) '/ar/' file];

        % filenamebaseimages = (['ar/' file]);
        iframe = 2;
        Iorig = imread([filenamebaseimages num2str(iframe) '.jpg']);
        imaqhandle = []; % Empty vector for when variables are passed to the function (AG)
        figure(1);
        
        % Requires user to draw ROIs only once for all of the video files
        % (AG)
        if i == 1
            image(Iorig);
            axis equal; axis tight; axis off;
            xlim(xlimits); ylim(ylimits);
            title('Mark blade ROI');
            [tongueroi,roix,roiy] = roipoly;
            title('Mark dorsum ROI');
            [tongueroi2,roi2x,roi2y] = roipoly;
            title('Mark root ROI');
            [tongueroi3,roi3x,roi3y] = roipoly;
            title('');
            %%%%%%%%%%%%% Edit made on 5/1/17 (AG) %%%%%%%%%%%%%
            %             figure(1);
            %             subplot(2,2,1);
            %             image(Iorig);
            %             axis equal; axis tight; axis off;
            %             xlim(xlimits); ylim(ylimits);
            %             hold on;
            %             plot(roix,roiy,'b-',roi2x,roi2y,'g-',roi3x,roi3y,'r-',...
            %                 'LineWidth',linewidth);
            %             hold off;
            %             title('Dorsum/blade/root ROIs','FontSize',fontsize);
            
            %%%%%%%%%%%%% Edit made on 3/23/17 %%%%%%%%%%%%%
            %%%%%%% Calculates max and min x and y coordinates to limit
            %%%%%%% filtering (AG)
            xRefn=horzcat(roi2x,roi3x,roix);
            yRefn=horzcat(roi2y,roi3y,roiy);
            xPtsMax=round(max(max(xRefn)));
            xPtsMin=round(min(min(xRefn)));
            yPtsMax=round(max(max(yRefn)));
            yPtsMin=round(min(min(xRefn)));
            
            xRootMin = round(min(roi3x));
            xRootMax = round(max(roi3x));
            yRootMin = round(min(roi3y));
            yRootMax = round(max(roi3y));
            
            iframe = 1; % 1 = /a/ Frame; 2 = /r/ Frame (AG)
            
            %Uses FindPointsR Function and passes rhs variables from the
            %workspace to track the /a/ Frame (AG)
            [disp1, refdist] = FindPointsR(i,res,iframe,filenamebaseimages,imaqhandle,tongueroi,tongueroi2,tongueroi3,roix,roi3x,roiy,roi3y,xPtsMin,xPtsMax,yPtsMin,yPtsMax,xRootMin,xRootMax,yRootMin,yRootMax);
        end
        
        iframe = 2; % 1 = /a/ Frame; 2 = /r/ Frame (AG)
        
        %Uses FindPointsR Function and passes rhs variables from the
        %workspace to track the /r/ Frame (AG)
        [disp2, refdist] = FindPointsR(i,res,iframe,filenamebaseimages,imaqhandle,tongueroi,tongueroi2,tongueroi3,roix,roi3x,roiy,roi3y,xPtsMin,xPtsMax,yPtsMin,yPtsMax,xRootMin,xRootMax,yRootMin,yRootMax);
        
        disp = vertcat(disp1,disp2);
        
        eps_root_norm = -((disp(2,3) - disp(1,3))*(sqrt(2)))/refdist;
        eps_dorsum_norm = -(disp(2,2) - disp(1,2))/refdist;
        eps_blade_norm = -(disp(2,1) - disp(1,1))/refdist;
        Displacements_norm = [eps_root_norm eps_dorsum_norm eps_blade_norm];
        
        eps_root_nonorm = -((disp(2,3) - disp(1,3))*(sqrt(2)));
        eps_dorsum_nonorm = -(disp(2,2) - disp(1,2));
        eps_blade_nonorm = -(disp(2,1) - disp(1,1));
        Displacements_nonorm = [eps_root_nonorm eps_dorsum_nonorm eps_blade_nonorm];
        
        figure(1);
        %%%%%%%%%%%%% Edit made on 5/1/17 (AG) %%%%%%%%%%%%%
        %subplot(2,2,4);
        bar([eps_root_norm eps_dorsum_norm eps_blade_norm]);
        ylabel('Normalized displacement','FontSize',fontsize);
        ylim([-0.23 0.23]);
        set(gca,'XTickLabel',{'Root','Dorsum','Blade'},'FontSize',fontsize);
        
        %%%%%%%%%%%%% Edit made on 5/3/17 (AG) %%%%%%%%%%%%%
        %%%%%%%%%%%%% Creates a matrix for the root, dorsum, and blade displacements relative to time %%%%%%%%%%
        SaveDataRoot(i) = eps_root_norm;
        SaveDataDorsum(i) = eps_dorsum_norm;
        SaveDataBlade(i) = eps_blade_norm;
        Time(i) = toc;
        
        pause(.01)
    end
    %%%%%%%% Saves Displacement data and timestamps to an excel file (AG) %%%%%%%%%
    Time=Time';
    Root = SaveDataRoot';
    Dorsum = SaveDataDorsum';
    Blade = SaveDataBlade';
    T = table(Time,Root,Dorsum,Blade);
    writetable(T,'All Displacements Data.xls')
    
else
    filenamebaseimages = []; % Empty vector for when variables are passed to the function (AG)
    
    %%%%% Edits for image acquisition made on 4/3/17 (TDM) %%%%%%%%%
    % Initialize frame grabber, only needs to be done once
    imaqhandle=imaq.VideoDevice;
    
    % Grab a single frame; similar to frame grab from mp4 file
    Iorig=step(imaqhandle);
    % Fixes color scale (Black and White) (TDM)
    Iorig = hsv2rgb(0*Iorig(:,:,1),0*Iorig(:,:,1),Iorig(:,:,1));
    
    iframe = 2;
    i=1;
    
    % User will draw ROIs (AG)
    if i == 1
        image(Iorig);
        axis equal; axis tight; axis off;
        xlim(xlimits); ylim(ylimits);
        title('Mark blade ROI');
        [tongueroi,roix,roiy] = roipoly;
        title('Mark dorsum ROI');
        [tongueroi2,roi2x,roi2y] = roipoly;
        title('Mark root ROI');
        [tongueroi3,roi3x,roi3y] = roipoly;
        title('');
        %%%%%%%%%%%%% Edit made on 5/1/17 (AG) %%%%%%%%%%%%%
        %         figure(1);
        %         subplot(2,2,1);
        %         image(Iorig);
        %         axis equal; axis tight; axis off;
        %         xlim(xlimits); ylim(ylimits);
        %         hold on;
        %         plot(roix,roiy,'b-',roi2x,roi2y,'g-',roi3x,roi3y,'r-',...
        %             'LineWidth',linewidth);
        %         hold off;
        %         title('Dorsum/blade/root ROIs','FontSize',fontsize);
        
        %%%%%%%%%%%%% Edit made on 3/23/17 %%%%%%%%%%%%%
        %%%%%%% Calculates max and min x and y coordinates to limit
        %%%%%%% filtering (AG)
        xRefn=horzcat(roi2x,roi3x,roix);
        yRefn=horzcat(roi2y,roi3y,roiy);
        xPtsMax=round(max(max(xRefn)));
        xPtsMin=round(min(min(xRefn)));
        yPtsMax=round(max(max(yRefn)));
        yPtsMin=round(min(min(xRefn)));
        xRootMin = round(min(roi3x));
        xRootMax = round(max(roi3x));
        yRootMin = round(min(roi3y));
        yRootMax = round(max(roi3y));
        
        iframe = 1;
        [disp1, refdist] = FindPointsR(i,res,iframe,filenamebaseimages,imaqhandle,tongueroi,tongueroi2,tongueroi3,roix,roi3x,roiy,roi3y,xPtsMin,xPtsMax,yPtsMin,yPtsMax,xRootMin,xRootMax,yRootMin,yRootMax);
    end
    set(gcf,'currentchar','o')
    while get(gcf,'currentchar')~=' '
        tic
        iframe = 2;
        [disp2, refdist] = FindPointsR(i,res,iframe,filenamebaseimages,imaqhandle,tongueroi,tongueroi2,tongueroi3,roix,roi3x,roiy,roi3y,xPtsMin,xPtsMax,yPtsMin,yPtsMax,xRootMin,xRootMax,yRootMin,yRootMax);
        
        disp = vertcat(disp1,disp2);
        
        eps_root_norm = -((disp(2,3) - disp(1,3))*(sqrt(2)))/refdist;
        eps_dorsum_norm = -(disp(2,2) - disp(1,2))/refdist;
        eps_blade_norm = -(disp(2,1) - disp(1,1))/refdist;
        Displacements_norm = [eps_root_norm eps_dorsum_norm eps_blade_norm];
        
        eps_root_nonorm = -((disp(2,3) - disp(1,3))*(sqrt(2)));
        eps_dorsum_nonorm = -(disp(2,2) - disp(1,2));
        eps_blade_nonorm = -(disp(2,1) - disp(1,1));
        Displacements_nonorm = [eps_root_nonorm eps_dorsum_nonorm eps_blade_nonorm];
        
        figure(1);
        %%%%%%%%%%%%% Edit made on 5/1/17 (AG) %%%%%%%%%%%%%
        %subplot(2,2,4);
        bar([eps_root_norm eps_dorsum_norm eps_blade_norm]);
        ylabel('Normalized displacement','FontSize',fontsize);
        ylim([-0.23 0.23]);
        set(gca,'XTickLabel',{'Root','Dorsum','Blade'},'FontSize',fontsize);
        
        %%%%%%%%%%%%% Edit made on 5/3/17 (AG) %%%%%%%%%%%%%
        SaveDataRoot(i) = eps_root_norm;
        SaveDataDorsum(i) = eps_dorsum_norm;
        SaveDataBlade(i) = eps_blade_norm;
        
        pause(.001)
        toc
    end
    Root = SaveDataRoot';
    Dorsum = SaveDataDorsum';
    Blade = SaveDataBlade';
    T = table(Root,Dorsum,Blade);
    writetable(T,'All Displacements Data.xls')
    
end