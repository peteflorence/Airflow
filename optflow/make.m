function make(mode)

if ~exist('mode', 'var')
    mode = 'release';
end

if strcmp(mode, 'release')
    mex -O -output ..\Coarse2FineTwoFrames Coarse2FineTwoFrames.cpp GaussianPyramid.cpp OpticalFlow.cpp
    %mex -O -output ..\Coarse2FineTwoFramesThreadsafe Coarse2FineTwoFramesThreadsafe.cpp GaussianPyramid.cpp OpticalFlow.cpp
elseif strcmp(mode, 'debug')
    mex -g -output ..\Coarse2FineTwoFrames Coarse2FineTwoFrames.cpp GaussianPyramid.cpp OpticalFlow.cpp
    %mex -g -output ..\Coarse2FineTwoFramesThreadsafe Coarse2FineTwoFramesThreadsafe.cpp GaussianPyramid.cpp OpticalFlow.cpp
else
    error('Incorrect mode!\n');
end

end

