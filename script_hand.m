% inpath = 'E:\project\FluidFlow\Code\data\eccv2014\single_hand2';
inpath = '../';
outpath = '../outpath/Trial8short-v2/';
mkdir(outpath);

vidname = 'hand2_input';
cropRectList = [];
alphaList = 0.04;
limitFlow = [-2e-3,2e-3];
avgFlowStr=0.003;
ffBeta2 = 30;
ffBeta3 = 2e-2;
avgFFStrthList = [0.3];
threshFFStr = 0.1;
timeWin = 15;

% Initialization
mdFrame = 5;
tempThetaList = 2;
tcolormap = gray(256);
tVarColormap = jet(256);
quiver_step = 20;

%% Read video
% if exist(fullfile(inpath, [vidname, '.avi']), 'file')
%     vidReader = VideoReader(fullfile(inpath, [vidname, '.avi']));
% else
%     vidReader = VideoReader(fullfile(inpath, [vidname, '.mov']));
% end
% vid = read(vidReader);
% vidfirst10 = vid(1:10,1:5);
% vid = im2double(colorvid2gray(vid));


imlist = dir('/Users/pflomacpro/RLG/Airflow_OSX/Trial8short/*.pgm');
im = cell(length(imlist),1);
for i=1:length(imlist)
     im{i} = imread(strcat('/Users/pflomacpro/RLG/Airflow_OSX/Trial8short/',imlist(i).name));
end
vid = cell2mat(shiftdim(im,-3));
%vid = im2double(colorvid2gray(vid));

    
%% Filter the video
tw = tempThetaList*2;
f = fspecial('gaussian', [1,tw*2+1], tempThetaList);
vid = convn(vid, shiftdim(f,-1), 'same');
vid = vid(:,:,mdFrame+1:end-mdFrame);

%% optical flow
[flow,vA] = optflow_GMRF2(vid, 'alpha1', 1, 'alpha2', alphaList, 'verbose', 2, 'sorMethod', 0);
%for t=1:length(vA)
%    vA{t} = vA{t} * 1e6;
%end

%% Visualization
[~,~,~,nF] = size(flow);
writeVid(mat2colorvid(squeeze(flow(:,:,1,:)), 'limits', limitFlow, 'map', tcolormap), ...
    fullfile(outpath,'flow_x.avi'), 30, 'Motion JPEG AVI');
writeVid(mat2colorvid(squeeze(flow(:,:,2,:)), 'limits', limitFlow, 'map', tcolormap), ...
    fullfile(outpath,'flow_y.avi'), 30, 'Motion JPEG AVI');
motionColorVisualizeOverlay(flow, vid(:,:,1:nF).^0.2, [], ...
    fullfile(outpath,'flow_colorOverlap.avi'), 'color_scale', avgFlowStr);

%% Run Fluid Flow
[ffmean,ffvar] = fluidflow_GMRF3(flow, vA, 'beta2', ffBeta2, 'beta3', ffBeta3, ...
        'verbose', 1, 'startSigma', 8, 'endSigma', 1, 'nOuterIter', 4, 'nInnerIter', 1, 'timeWin', 15, 'outFrameValidOnly', true);    
 avgFFStr = avgFFStrthList;
 tballMag = 6;
 ffvardet = squeeze(ffvar(:,:,1,:).*ffvar(:,:,3,:) - ffvar(:,:,2,:).^2);
 ffvardet(ffvardet<0) = 1e20;
 ffvardet = sqrt(ffvardet);
 weightVid = (threshFFStr./ffvardet).^2; weightVid(weightVid>1)=1;
%%
 motionColorVisualizeOverlay(ffmean, repmat(permute(vid.^(1/3),[1,2,4,3]),[1,1,3,1]), weightVid, 'ff_colorOverlay.avi',...
     'color_scale', avgFFStr);
 motionColorVisualize(ffmean, 'ff_color.avi', 'color_scale', threshFFStr);
 motionQuiverVisualize(vid(:,:,1:size(ffmean,4))/2, ffmean, 'ff_quiver_small.avi', 'quiver_step', 20, ...
     'quiver_scale', avgFFStr, 'quiver_width', 2, 'vidRescale', 1, 'quiver_color', [1,0.6,0], ...
     'quiver_headsize', 2, 'abs_quiver_scale', 0.9);
