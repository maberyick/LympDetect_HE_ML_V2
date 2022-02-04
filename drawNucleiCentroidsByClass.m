function drawNucleiCentroidsByClass( image,centroids,class, markerSize )

 if nargin<4
     markerSize=5;    
 end

colors={'b','y','g','r'};

imshow(image);
hold on;

cl=unique(class);
numClass=length(cl);

for i=1:numClass
    cent=centroids(class==cl(i),:);
    plot(cent(:,1),cent(:,2),[colors{i} '*'],'MarkerSize',markerSize);
end

hold off;
end