function [ centroids,features,featureNames ] = getNucLocalFeatures( image, mask, logicalMask, orderZernikePol )
%GETNUCLOCALFEATURES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    logicalMask=true;
end

if nargin<4
    orderZernikePol=7;
end

image=normalizeStaining(image);
grayImg=rgb2gray(image);

%[h,~,~]=size(image);
%[mh,~,~]=size(mask);
%if h~=mh
%    image=imrotate(image,90);
%end

if logicalMask
    mask = bwareafilt(logical(mask(:,:,1)),[20 7000]);
    
    
    regionProperties = regionprops(mask,grayImg,'Centroid','Area',...
        'BoundingBox','Eccentricity','EquivDiameter','Image',...
        'MajorAxisLength','MaxIntensity','MeanIntensity','MinIntensity',...
        'MinorAxisLength','Orientation','PixelValues');
else
    if size(mask,3)>1
        mask=rgb2gray(mask);
    end
    vals=unique(mask)';
    
    regionProperties=[];
    for i=vals
        if i~=0
             regionProperties =  [regionProperties;regionprops(mask==i,...
                 grayImg,'Centroid','Area','BoundingBox','Eccentricity',...
                 'EquivDiameter','Image','MajorAxisLength','MaxIntensity',...
                 'MeanIntensity','MinIntensity','MinorAxisLength',...
                 'Orientation','PixelValues')];
        end
    end
end

centroids = cat(1, regionProperties.Centroid);


medRed=[];
medGreen=[];
medBlue=[];
entropyRed=[];
medIntensity=[];
entropyIntensity=[];
harFeat=[];
edgeMedIntensity=[];
zernFeat=[];

nucleiNum=size(regionProperties,1);
for i=1:nucleiNum
    nucleus=regionProperties(i);
    bbox = nucleus.BoundingBox;
    bbox = [round(bbox(1)) round(bbox(2)) (bbox(3) - 1) (bbox(4) - 1)];
    roi = image(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3), :);
    per = bwperim(nucleus.Image);
    
    gray=rgb2gray(roi);
    R=roi(:,:,1);
    G=roi(:,:,2);
    B=roi(:,:,3);
    
    R=R(nucleus.Image == 1);
    G=G(nucleus.Image == 1);
    B=B(nucleus.Image == 1);
    grayPix=gray(nucleus.Image == 1);
    perPix=gray(per==1);
    
    % Intensity features:
    medRed = [medRed;median(double(R))];
    medGreen = [medGreen;median(double(G))];
    medBlue = [medBlue;median(double(B))];
    medIntensity=[medIntensity;median(double(grayPix))];
    edgeMedIntensity=[edgeMedIntensity;median(double(perPix))];
    
    % Texture features:
    
    % Entropies
    entropyRed=[entropyRed; getNucEntropy(R)];
    entropyIntensity=[entropyIntensity;getNucEntropy(grayPix)];
    
    % Haralick features
    glcm = graycomatrix(gray);
    try
        h=haralickTextureFeatures(glcm,1:14);
    catch ex
        h=zeros(14,1);
    end
    harFeat=[harFeat;h'];
    
    % Zernike moments
    [w,h]=size(nucleus.Image);
    sqRoi=getSquareRoi(grayImg,nucleus.Centroid,w,h);
    zernFeat=[zernFeat;getZernikeMomentsImg(sqRoi,orderZernikePol)];
    
end

ratioMedRB=medRed./medBlue;
ratioMedRG=medRed./medGreen;

ratioAxes=[regionProperties.MajorAxisLength]./[regionProperties.MinorAxisLength];

features = horzcat([regionProperties.Area]',[regionProperties.Eccentricity]', ...
    ratioAxes',medRed,entropyRed,double([regionProperties.MinIntensity]'),...
    double([regionProperties.MaxIntensity]'),[regionProperties.EquivDiameter]',...
    [regionProperties.Orientation]',entropyIntensity,medIntensity,edgeMedIntensity,...
    ratioMedRB,ratioMedRG,harFeat,zernFeat...
    );

numZernikePol=size(zernFeat,2)/2;
zernNames={};
for i=1:numZernikePol
    num=num2str(i);
    zernNames=[zernNames, ['ZernPol' num '_A'], ['ZernPol' num '_Phi']];
end

harNames={'AngularSecondMoment','Contrast','Correlation',...
    'Variance','Homogeneity','SumAverage','SumVariance','SumEntropy',...
    'Entropy','DifferenceVariance','DifferenceEntropy','InfoMeasureCorrI',...
    'InfoMeasureCorrII','MaxCorrCoeff'};

featureNames=[{'Area','Eccentricity','RatioAxes','MedianRed','EntropyRed',...
    'MinIntensity','MaxIntensity','EquivDiameter','Orientation',...
    'EntropyIntensity','MedianIntensity','EdgeMeanIntensity','RatioMedianRedBlue',...
    'RatioMedianRedGreen'},harNames,zernNames];

end

