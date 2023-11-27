% USAGE 
%  Detect red and blue LEDs position in a video file and creates a 'led' file 
%
%    Process_DetectLED(videoFile,<options>)
%
%    videoFile      path to Basler video file, including the '.avi'
%                   extension or not
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'manualROI'   boolean: whether or not you want to manually adjust the
%                   area where LEDs are detected (default = 1)
%    =========================================================================
%
% DEPENDENCIES:
%
%   Computer vision toolbox
 
% Inspiration from Adrien Peyrache lab and from John D Long II
%  
%    =========================================================================
%How to run (Prawesh Dahal 2023) 
% Go to the desired folder with rat behavior video files
% detec_file = dir('*_test*avi');

% LEDpos = []; 

%for i = 1:length(detec_file) 
%    [~, fbasename, ~] = fileparts(detec_file(i).name); 
%    whl = Process_DetectLED_optimized(fbasename,2); 
%    LEDpos = [LEDpos; whl]; 
%    clear whl 
%end
% =========================================================================
 


function whl = Process_DetectLED_optimized(fbasename,factor,varargin)


manualROI = 1;

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'manualroi',
      manualROI = varargin{i+1};
      if ~isa(manualROI,'numeric')
        error('Incorrect value for property ''manualROI'' (type ''help Process_DetectLED'' for details).');
      end
   
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help Process_DetectLED'' for details).']);
  end
end


if strcmp(fbasename(end-3:end),'avi')
    fbasename = fbasename(end-3:end);
end
%%
file = [fbasename '.avi'];

if ~exist(file,'file')
    warning('No video file')
    keyboard
end

%%


% Creater readerobj of file
videoObj    = vision.VideoFileReader(file,'VideoOutputDataType','uint8');
videoSize   = (videoObj.info.VideoSize)./factor;
width       = videoSize(1);
height      = videoSize(2);
threshF     = 150;         % threshold on foreground pixel intensity
%%
% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);

% Define ROI (manual definition of the environment)
%try 
%First, we adapt the dynamical range of the pixel colors for manual
%selection of the ROI (i.e. the environment)
%load('firstFrame.mat');
%firstFrame = ans;


firstFrame  = step(videoObj);
firstFrame  = imresize(firstFrame,1/factor);


frame = imadjust(rgb2gray(firstFrame));
%%
if manualROI
    ok = 0;
    figure(1),clf

    while ~ok
        clf,imshow(frame)
        fprintf('Define your ROI. Click ''enter'' when finished\n')
        [x,y] = ginput;
        inArea = inpolygon(X(:),Y(:),x,y);
        inArea = double(reshape(inArea,[height width]));
        frame(~inArea) = 0;

        clf,imshow(frame)

        reply = input('OK with the result? Y/N [Y]:','s');
        if ~strcmp(reply,'N')
            ok = 1;
        end
    end
else
    load frame_crop.mat;
    inArea=imresize(inArea,1/factor);
end
%%

% Initialize background as a grayscale image of the first frame
bg_bw       = rgb2gray(firstFrame);

% Initialize foreground image and difference image
fg          = zeros(size(bg_bw));
fr_diff     = zeros(size(bg_bw));

% Initialize color mask
mask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');

% Initialize whl matrix
whl = [];
count = 0;

try
while ~isDone(videoObj)
%for i = Fint:readerobj.NumberOfFrames
    
    if count~=0
        backSp = repmat('\b',[1 length(num2str(count-1))]);
        fprintf(backSp)
    end
     fprintf('%i',count)

    if count~=0
        fr    = step(videoObj);
        fr  = imresize(fr,1/factor);

    else
        fr = firstFrame;
    end
    
    % convert frame to grayscale
    fr_bw = rgb2gray(fr);
    %%% Label Color Mask
    label         = repmat(logical(fr_bw>threshF) & inArea,[1 1 3]);
    mask(label)   = fr(label);
    mask(~label)  = 0;
    
    %%% Find centroid of remaining pixels %%%
    %Red
    bw_mask = squeeze(mask(:,:,1));
    [CC,nc] = bwlabel(bw_mask);
    
    if nc>0
%         pixels = regionprops(CC,'PixelList');
        pixels = regionprops(CC,'Area');
        
%         centroidSize = zeros(nc,1);
%         for ii=1:nc
%             p = pixels(ii).PixelList;
%             centroidSize(ii) = length(p);
%         end
%         [~,mxIx] = max(centroidSize);
        [~,mxIx] = max([pixels.Area]);
%         Rr   = round(mean(pixels(mxIx).PixelList,1));
        tmp  = zeros(size(CC),'uint8');
        tmp(CC == mxIx) = 1;
        chk  = regionprops(tmp,'centroid');
        Rr   = round(chk.Centroid);
    else
        Rr = [-1 -1];
    end
    
    %Blue
%     Br = [-1 -1];
    bw_mask = squeeze(mask(:,:,3));
    [CC,nc] = bwlabel(bw_mask);
%     pixels = regionprops(CC,'PixelList');
    
    if nc>0
%         pixels = regionprops(CC,'PixelList');
        pixels = regionprops(CC,'Area');
%         centroidSize = zeros(nc,1);
%         for ii=1:nc
%             p = pixels(ii).PixelList;
%             centroidSize(ii) = length(p);
%         end
%         [~,mxIx] = max(centroidSize);
        [~,mxIx] = max([pixels.Area]);
        tmp  = zeros(size(CC),'uint8');
        tmp(CC == mxIx) = 1;
        chk  = regionprops(tmp,'centroid');
%         Br   = round(mean(pixels(mxIx).PixelList,1));
        Br = round(chk.Centroid);
    else
         Br = [-1 -1];
    end
     
    % pre-allocate
     whl = [whl;[Rr(1),Rr(2),Br(1),Br(2)]];
    
    % End processing time
    
    % Display results every 1000 frame 
    if 0 %set it at 1 for debug
        if mod(count,1000)==0
%             figure(1),h1 = subplot(3,1,1);imshow(uint8(fr)), %title(sprintf('Frame %d, Time %f, Connected Components %d',i,t2,Nl))
%             subplot(3,1,2),imshow(uint8(floor(fr_diff)))
%             h2 = subplot(3,1,3);imshow(mask)
            %h3 = impoly(h1,[xj,yj]);
            %h4 = impoly(h2,[xi,yi]);
            %setColor(h3,'green')
            %setColor(h4,'yellow')
            %keyboard
            figure(2),
            plot(whl(:,1),whl(:,2),'r')
            hold on
            plot(whl(:,3),whl(:,4),'b')
            drawnow
         end
    end
    count = count+1;
end
catch
   keyboard
end
%write result file
dlmwrite([fbasename '.led'],whl.*factor,'\t')

fprintf('\n\n')
release(videoObj)