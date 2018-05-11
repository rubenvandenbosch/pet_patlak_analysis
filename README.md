# Matlab code for Patlak graphical analysis of PET data

DESCRIPTION  
The patlak analysis entails the following steps:
1. Read headers and image data of PET images and mask (ref ROI) image
2. Create input function from PET data in reference ROI
3. Voxel-wise Patlak analysis
      - Integrate input function
      - Calculate X of the Patlak plot function
      - Per slice calculate Y of the Patlak function for each voxel
      - Per slice calculate the slope (Ki) and intercept (Vd) of the 
        Patlak function for each voxel
4. Save results of Ki and Vd voxel values in NIFTI images

INPUTS  
files         = PET images to analyze  
timingsfile   = timing of PET scans  
mask          = mask image of reference region  
outputdir     = output directory for patlak output  

OUTPUTS
1. NIFTI image of Ki map
2. NIFTI image of Vd map
 
DEPENDENCIES
1. Function: IntFunction = integration(Function,frames);  
               - For input function calculation
2. Function: [s,int] = trend(A,t,dim)  
               - For slope calculation as a matrix operation

Both functions are included as subfunctions.

-------------------------------------------------------------------------
Ruben van den Bosch  
Donders Institute for Brain, Cognition and Behaviour  
Radboud University  
Nijmegen, The Netherlands  
January 2018  

Based on the script GA_ParametricMapping_nifti_ABedit2015 from ...,  
University of Berkeley
