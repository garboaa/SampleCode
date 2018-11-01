function [disp3, refdist] = FindPointsR(i,res,iframe,filenamebaseimages,imaqhandle,tongueroi,tongueroi2,tongueroi3,roix,roi3x,roiy,roi3y,xPtsMin,xPtsMax,yPtsMin,yPtsMax,xRootMin,xRootMax,yRootMin,yRootMax)
%%%%% Began editing function on 3/20/17  (AG) %%%%%
xlimits = [100,1000];
ylimits = [0, 900];

pointsize = 6;
linewidth = 2;
maxthresh = .6;
winsize = 48;    % size of smoothing mask, # pixels
fontsize=11;     % height of text on plots (pt)

% relative positions of reference points within ROI (between 0 and 1)
refpts = [1/3 1/2 2/3];
h = hann(winsize)*hann(winsize)';
if res == 1
    Iorig = imread([filenamebaseimages num2str(iframe) '.jpg']);
    I(:,:,iframe) = mean(double(Iorig),3);
    %%%%% Edit made on 4/6/17 (AG) %%%%%
    I(yPtsMin:yPtsMax,xPtsMin:xPtsMax,iframe) = conv2(I(yPtsMin:yPtsMax,xPtsMin:xPtsMax,iframe),h,'same');
else
    Iorig = step(imaqhandle);
    I(:,:,iframe) = Iorig(:,:,1);
    %%%%% Edit made on 4/6/17 (AG) %%%%%
    I(yPtsMin:yPtsMax,xPtsMin:xPtsMax,iframe) = conv2(I(yPtsMin:yPtsMax,xPtsMin:xPtsMax,iframe),h,'same');
end

% local maxima for measuring vertical displacement in ROI 1 (blade)
[Imax iymax] = max(I(:,:,iframe) .* tongueroi);
ixpeak = find(Imax > maxthresh*max(Imax));
iypeak = iymax(ixpeak);
plot(ixpeak,iypeak,'.b');

% local maxima for measuring vertical displacement in ROI 2 (dorsum)
[Imax iymax] = max(I(:,:,iframe) .* tongueroi2);
ixpeak2 = find(Imax > maxthresh*max(Imax));
iypeak2 = iymax(ixpeak2);
plot(ixpeak2,iypeak2,'.g');


%%% EDIT 2 %%%
% Create new variables for the number of rows and columns in the whole matrix of the frame of the ultrasound image (lines 102-103)
% Create a new variable for the matrix of extracted diagonals from the entire matrix of the frame of the ultrasound image (line 104)
r = size(I(:,:,iframe),1);
c = length(I(:,:,iframe));
I_diags = spdiags(I(:,:,iframe),-(r-1):c-1);

%%% EDIT 3 %%%
% Create new variables for the number of rows and columns of the region of interest (converted into matrix) for root drawn by user (lines 109-110)
% Create a new variable for the matrix of extracted diagonals from the drawn region of interest of the root (line 111)
r1 = size(tongueroi3(:,:),1);
c1 = length(tongueroi3(:,:));
t_diags = spdiags(tongueroi3(:,:),-(r1-1):c1-1);


%%% EDIT 4 %%%
% Finds the local maxima once the new region of interest matrix is applied to the new matrix of the frame of the ultrasound image (line 122)
% The values of each pixel remain the same in the diagonal image and original image. So, the maxima in the diagonal matrix is found, and then the locations of where those maxima are achieved in the original matrix is determined
[Imax iymax] = max(I_diags.*t_diags);
if max(Imax) ~= 0
    A = I(:,:,iframe).*tongueroi3(:,:);
    Imax2 = nonzeros(Imax);
    Imax2t = Imax2';
    M = zeros(length(Imax2t),2);
    k = 1;
    while (k<=length(Imax2t))
        while (Imax2t(k)> maxthresh*max(Imax2t))
            %%%%%%%%%%%%% Edit made on 4/19/17 (AG) %%%%%%%%%%%%%
            [x,y] = find(Imax2t(k) == A(yRootMin:yRootMax,xRootMin:xRootMax));
            M(k,1) = min(y)+yRootMin-1;
            M(k,2) = min(x)+xRootMin-1;
            break;
        end
        k = k + 1;
    end
    
    ixpeak3 = nonzeros(M(:,1))';
    iypeak3 = nonzeros(M(:,2))';
    plot(ixpeak3,iypeak3,'.r')
    
    hold off;
    
    % find reference points: ROI 1, vertical displacement
    nxpeak = length(ixpeak);
    refxindices = 1+ round(refpts*(nxpeak-1));
    ixpeakref = ixpeak(refxindices);
    iypeakref = iypeak(refxindices);
    
    % find reference points: ROI 2, vertical displacement
    nxpeak2 = length(ixpeak2);
    refxindices = 1+ round(refpts*(nxpeak2-1));
    ixpeak2ref = ixpeak2(refxindices);
    iypeak2ref = iypeak2(refxindices);
    
    % find reference points: ROI 3, diagonal displacement
    nypeak3 = length(iypeak3);
    nxpeak3 = length(ixpeak3);
    refxindices = 1+ round(refpts*(nxpeak3-1));
    refyindices = 1+ round(refpts*(nypeak3-1));
    ixpeak3ref = ixpeak3(refxindices);
    iypeak3ref = iypeak3(refyindices);
    
    % find reference difference for normalization: here, distance from
    %  center of blade ROI to center of root ROI
    xblade = mean(roix); yblade = mean(roiy);
    xroot = mean(roi3x); yroot = mean(roi3y);
    refdist = sqrt((yroot-yblade).^2 + (xroot - xblade)^2);
    
    %%% EDIT 5 %%%
    % pixel positions of each ROI for current frame
    % y reference point for blade, dorsum, and root
    disp3=[mean(iypeakref) mean(iypeak2ref) mean(iypeak3ref)];
    
    %%%%%%%%%%%%% Edit made on 5/1/17 (AG) %%%%%%%%%%%%%
    % if i == 1
    %     figure(1);
    %     subplot(2,2,1+iframe);
    %     image(Iorig);
    %     if res == 2
    %         imagesc(I(:,:,iframe),[0 255]); colormap(gray(256));
    %     end
    %     axis equal; axis tight; axis off;
    %     xlim(xlimits); ylim(ylimits);
    %     hold on;
    %     plot(ixpeakref,iypeakref,'ob','MarkerSize',pointsize,'LineWidth',linewidth);
    %     plot(ixpeak2ref,iypeak2ref,'og','MarkerSize',pointsize,'LineWidth',linewidth);
    %     plot(ixpeak3ref,iypeak3ref,'or','MarkerSize',pointsize,'LineWidth',linewidth);
    %     hold off;
    %     pause(.01);
    % end
end

