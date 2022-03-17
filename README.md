# Lymphocyte Detection Model (MATLAB-ML)
ML model in MATLAB to detect and differentiate lymphocyte from other cells (non-lymphocytes) on H&amp;E-stained digitized images.
## How to run

```MATLAB
clear; clc;

I=imread('data/case713_1001_501.png'); % original image
M=imread('data/case713_1001_501_mask.png'); % mask (nuclei segmentation)

%If mask is not available, use M=getWatershedMask(I,normalize,minRad,maxRad)
%(suggested params: normalize=true, minRad=4, & maxRad=10)
	
lympModel=load('data/lymp_svm_matlab_wsi.mat'); % trained models	

% centroids contains the position (x,y) of all detected nuclei
[nucCentroids,nucFeatures] = getNucLocalFeatures(I,M);
	
% isLympocyte indicates wheter or not each detected nucleus is a lymphocyte. (1 - lymphocyte, 0 - non-lymphocyte)
isLymphocyte = (predict(lympModel.model,nucFeatures(:,1:7)))==1;
	
% Optional: draws the image showing the class of each centroid (lymphocyte or non-lymphocyte)
drawNucleiCentroidsByClass(I,nucCentroids,isLymphocyte);
```
