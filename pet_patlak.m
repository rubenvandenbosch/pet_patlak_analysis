function pet_patlak(files,timingsfile,mask,outputdir)
% PET_PATLAK
% Patlak graphical analysis code for PET data
%
% DESCRIPTION
% The patlak analysis entails the following steps:
% 1. Read headers and image data of PET images and mask (ref ROI) image
% 2. Create input function from PET data in reference ROI
% 3. Voxel-wise Patlak analysis
%       - Integrate input function
%       - Calculate X of the Patlak plot function
%       - Per slice calculate Y of the Patlak function for each voxel
%       - Per slice calculate the slope (Ki) and intercept (Vd) of the 
%         Patlak function for each voxel
% 4. Save results of Ki and Vd voxel values in NIFTI images
% 
% INPUTS
% files         = PET images to analyze 
% timingsfile   = timing of PET scans 
% mask          = mask image of reference region
% outputdir     = output directory for patlak output
%
% OUTPUTS
% 1. NIFTI image of Ki map
% 2. NIFTI image of Vd map
% 
% DEPENDENCIES
% 1. Function: IntFunction = integration(Function,frames); 
%               - For input function calculation
% 2. Function: [s,int] = trend(A,t,dim)
%               - For slope calculation as a matrix operation
% 
% Both functions are included as subfunctions.
% -------------------------------------------------------------------------
% Ruben van den Bosch
% Donders Institute for Brain, Cognition and Behaviour
% Radboud University
% Nijmegen, The Netherlands
% January 2018
% 
% Based on the script GA_ParametricMapping_nifti_ABedit2015 from ...,
% University of Berkeley
%

% Specify time span for which to run the analysis (in minutes)
StartTime = 24;
StopTime = 89;

% Select PET images to analyze
nFile = length(files);

% Read in timing data
timings         = csvread(timingsfile,1,1);
frames.start    = timings(:,1)';
frames.end      = timings(:,2)';
frames.duration = timings(:,3)';

% Assert number of nii files to analyze matches number of time points
assert(numel(frames.start)==nFile, ...
    'Number of nii files to analyze does not match the number of timepoints (frames) specified');

% 1. Read headers and image data
% =========================================================================
% PET images
% -------------------------------------------------------------------------
% Get image dimensions of one of the images for pre-allocation purposes
V1 = spm_vol(files{1});

% Pre-allocate matrices for headers and image data
V = cell(1,nFile);
Y = nan([V1.dim,nFile]);

% Read PET headers and image data
for iFile = 1:nFile
    V{iFile} = spm_vol(files{iFile});
    Y(:,:,:,iFile) = spm_read_vols(V{iFile});
end

% Reference ROI
% -------------------------------------------------------------------------
VM = spm_vol(mask);

% Assert identical dimensions for ROI and data images
assert(all(VM.dim == V1.dim), 'Dimensions of mask image does not match with dimensions of data images');

% Read ROI mask data and make sure data type is uint
YM = spm_read_vols(VM);
YM = uint8(YM);

% 2. Calculate input function from data in reference ROI 
% =========================================================================
nTimePoints = nFile;

% Extract mean data value in mask area for each time point
for iTime = 1:nTimePoints
    timepoint_data = Y(:,:,:,iTime);
    input_function(iTime) = mean(timepoint_data(YM==1));
end

% Clear timepoint_data for clearing up memory
clear timepoint_data

% 3. Voxel-wise Patlak analysis
% =========================================================================
% Pre-allocate matrices for Ki and Vd maps
% -------------------------------------------------------------------------
Ki_map = zeros(V1.dim);
Vd_map = zeros(V1.dim);

% Determine frames to use in slope calculation
% -------------------------------------------------------------------------
time_minutes = frames.start/60;
frameidx = find((time_minutes >= StartTime) & (time_minutes <= StopTime));

if ~isempty(frameidx)
  start_frame = frameidx(1);
  end_frame = frameidx(end);
else
  error('Unable to run analysis: no frames within Start and End Times')
end

% Integrate input function (via subfunction)
% -------------------------------------------------------------------------
int_input_func = patlak_integration(input_function,frames);

% Calculate patlakX
% -------------------------------------------------------------------------
patlakX = int_input_func./input_function;

% Create matrix with replications of input function, size is
% sizeY(1)xsizeY(2)xnTimepoints. Matrix will be used to calculate patlakY
% values per slice.
% -------------------------------------------------------------------------
sizeY = size(Y);

input_func_mat = reshape(input_function,1,1,numel(input_function),1);
input_func_mat = repmat(input_func_mat,[sizeY(1),sizeY(2)],1);

% Voxel-wise calculation per slice: rate of net influx, Ki
% -------------------------------------------------------------------------
nSlices = sizeY(3);

for iSlice = 1:nSlices

    % After squeeze, dimension 3 is time
    sliceData = squeeze(Y(:,:,iSlice,:));

    % Calculate patlakY
    patlakY = sliceData./input_func_mat;

    % Calculate the slopes and intercepts
    [s,int] = trend(patlakY(:,:,start_frame:end_frame),patlakX(start_frame:end_frame),3);

    % Save Ki and Vd values of this slice
    Ki_map(:,:,iSlice) = s(:,:,:).*60;
    Vd_map(:,:,iSlice) = int(:,:,:);
end

% 4. Save results
% =========================================================================
fnameKi = fullfile(outputdir,'Ki_map.nii');
fnameVd = fullfile(outputdir,'Vd_map.nii');

% Create Ki and Vd volumes
% -------------------------------------------------------------------------
% IMPORTANT: set datatype to float! Otherwise the small Ki values get
% rounded to zero
V_Ki        = V1;
V_Ki.dt     = [16 0];
V_Ki.fname  = fnameKi;

V_Vd        = V1;
V_Ki.dt     = [16 0];
V_Vd.fname  = fnameVd;

% Write Ki and Vd data to NIFTI images
% -------------------------------------------------------------------------
Vki = spm_write_vol(V_Ki,Ki_map);
Vvd = spm_write_vol(V_Vd,Vd_map);

end

function IntFunction = patlak_integration(Function,frames)
% KRONOS INTEGRATION METHOD

MidFrame = frames.start + frames.duration/2;

IntFunction = MidFrame(1)*Function(1,1);

for i = 2:size(Function,2)
  IntFunction(i) = (MidFrame(i)-MidFrame(i-1))*Function(1,i) + IntFunction(i-1);
end
end

function [s,int] = trend(A,t,dim)
% [s,int] = trend(A) returns the linear trend of A by least squares. If A is a
% vector, s is the slope of the linear fit.  If A is an N-dimensional array, 
% s and int are (N-1)-dimensional matrices corresponding to the slope and 
% intercept along the dimension dim.  
% 
%
% SYNTAX: 
% s = trend(A) 
% s = trend(A,Fs) 
% s = trend(A,t) 
% s = trend(...,dim) 
% s = trend(A,[],dim)
% [s,int] = trend(...)
%
% 
% DESCRIPTION: 
% s = trend(A) returns the (N-1)-dimensional matrix s corresponding to the
% linear trend(s) along dimension 1 of A. Assumes data are evenly spaced
% along dimension 1. 
% 
% s = trend(A,Fs) declares sampling frequency Fs along trending dimension of A.     
% 
% s = trend(A,t) allows for unevenly-spaced data in the trending dimension
% with time vector t. length of t must equal the length of A along its
% trending dimension.  
% 
% s = trend(...,dim) returns the trend along dimension dim of A. 
% 
% s = trend(A,[],dim) assumes data are sampled at 1 Hz (or 1/(unit time) or
% 1/(unit space) or what-have-you). 
%
% [s,int] = trend(...) also returns the intercepts of the slope-intercept
% form. 
% 
% 
% *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
% Version history: 
% Version 1: April 25, 2014
% Version 2: April 28, 2014 now works for any N-dimensional inputs and
%            optionally returns intercepts.  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% Chad A. Greene
% Institute for Geophysics
% The University of Texas at Austin
% April 2014
% 
% Thanks to Matt J. for writing much of this code. 
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% 
% See also: polyfit. 

% Input checking: 
assert(~isscalar(A),'Input data A cannot be a scalar.'); 
sizeA = size(A); 
 
if exist('dim','var')
    assert(ndims(A)>=dim,'dim argument exceeds dimensions of A'); 
else
    dim = 1; 
end
dimLength = sizeA(dim); 

if exist('t','var')
    if isscalar(t)
        t = 0:1/t:(dimLength-1)/t; % for the case where t is declared as sampling frequency.
        if nargout == 2
            disp('Warning: independent variable vector starts at zero. This may affect intercept value.')
            disp('If you want full confidence in the intercept value, use time vector as input instead of scalar sampling frequency.') 
        end
    end
  
    if isempty(t)
        t = 0:dimLength-1; 
    end
    
else
    t = 0:dimLength-1;     
end
t = t(:); 

order=[dim, setdiff(1:ndims(A),dim)];

sizeA(dim)=1;
sizeA=sizeA(order);

data2D = reshape(permute(A,order),dimLength,[]);
coefficients = ([t ones(size(t))]\data2D); 
s = squeeze(ipermute(reshape(coefficients(1,:),sizeA),order));
if nargout==2
    int = squeeze(ipermute(reshape(coefficients(2,:),sizeA),order));
end

end

