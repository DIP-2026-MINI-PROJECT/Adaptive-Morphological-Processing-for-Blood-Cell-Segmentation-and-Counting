# Project Name
Adaptive Morphological Processing for Blood Cell Segmentation and Counting Under Variable Imaging Conditions
### Project Description:
#### Summary - Adaptive Morphological Processing for Blood Cell Segmentation and Counting Under Variable Imaging Conditions, focuses on automating the process of identifying and counting blood cells from microscopic images. We use image processing techniques such as filtering, thresholding, and morphological operations to improve image quality, separate individual cells, and obtain accurate cell counts

#### Course concepts used - 
DIP Concepts Used:
RGB colour space processing.
Gaussian filtering and noise reduction.
Otsu's automatic thresholding.
Morphological operations.
Connected component analysis.
Distance transform and watershed segmentation

   
#### Additional concepts used -
1. -
2. -
   
#### Dataset - 
Google Blood smear test images 

#### Novelty - 
As of now we do not have a major novelty but there is scope for future work as it can be made more adaptive

### Contributors:
Advika Dokania (PES1UG24EC020)
Khushi Shetty (PES1UG24EC103)
Lakshmi Sri R (PES1UG24EC110)
### Steps:
 Algorithm Steps
Image Acquisition & Preprocessing: Read the blood smear image and apply Gaussian filtering to reduce noise.
Cell Detection: Use automatic thresholding on the green channel to separate cells from the background.
Morphological Processing: Apply opening, closing and hole filling to remove noise and improve cell shapes.
Cell Size Estimation: Estimate typical cell area using region properties for adaptive processing.
WBC Detection: Identify dark purple/blue nuclei using intensity thresholding and colour differences.
RBC Detection: Detect remaining cells and separate overlapping cells using watershed segmentation.
Cell Counting & Visualization: Count RBCs and WBCs separately and display their boundaries on the original image.
MATLAB CODE
clc; clear; close all;

%% Read image
img = imread("3.jpeg");

if size(img,3) == 4
    img = img(:,:,1:3);
end

if size(img,3) ~= 3
    error('Please use an RGB blood smear image.');
end

img = im2uint8(img);

% Resize only if image is very large
scale = min(1,800/max(size(img,1),size(img,2)));
if scale < 1
    img = imresize(img,scale);
end

R = imgaussfilt(img(:,:,1),1);
G = imgaussfilt(img(:,:,2),1);
B = imgaussfilt(img(:,:,3),1);

%% 1. All cells
cells = G < 255*graythresh(G);

if mean(cells(:)) > 0.8
    cells = ~cells;
end

cells = imfill(imclose(cells,strel('disk',2)),'holes');
cells = imopen(cells,strel('disk',2));
cells = bwareaopen(cells,max(30,round(0.0002*numel(G))));

%% 2. Estimate cell size
s = regionprops(cells,'Area','Solidity');
a = [s.Area];

singles = a([s.Solidity] > 0.9);

if numel(singles) >= 3
    typ = median(singles);
else
    typ = median(a);
end

r = sqrt(typ/pi);

%% 3. WBC detection
tNuc = 255*graythresh(R(cells));
nuc = cells & (R < tNuc);

RB = double(R) - double(B);
refRB = median(RB(cells));

lab = bwlabel(nuc,8);
st = regionprops(lab,RB,'Area','MeanIntensity');

keep = find([st.Area] >= 0.3*typ & ...
            [st.MeanIntensity] < 0.5*refRB);

nuc = ismember(lab,keep);
nuc = imfill(nuc,'holes');

%% 4. Separate touching WBCs
CC = bwconncomp(nuc,8);

wbcLabels = zeros(size(nuc),'uint16');
wbcCount = 0;

for k = 1:CC.NumObjects

    areaNow = numel(CC.PixelIdxList{k});

    current = false(size(nuc));
    current(CC.PixelIdxList{k}) = true;

    if areaNow > 1.5*typ

        [yy,xx] = find(current);

        y1 = max(1,min(yy)-5);
        y2 = min(size(nuc,1),max(yy)+5);
        x1 = max(1,min(xx)-5);
        x2 = min(size(nuc,2),max(xx)+5);

        smallMask = current(y1:y2,x1:x2);

        D = bwdist(~smallMask);

        h = max(1,round(0.25*r));

        markers = imextendedmax(D,h);
        markers = markers & smallMask;

        markerCC = bwconncomp(markers,8);

        if markerCC.NumObjects >= 2

            L = watershed(imimposemin(-D, ...
                ~smallMask | markers));

            L(~smallMask) = 0;

            props = regionprops(L,'Area');

            valid = find([props.Area] >= 0.25*typ);

            L(~ismember(L,valid)) = 0;

            ids = unique(L(:));
            ids(ids == 0) = [];

            for j = 1:numel(ids)

                wbcCount = wbcCount + 1;

                temp = wbcLabels(y1:y2,x1:x2);
                temp(L == ids(j)) = wbcCount;
                wbcLabels(y1:y2,x1:x2) = temp;

            end

        else

            wbcCount = wbcCount + 1;
            wbcLabels(current) = wbcCount;

        end

    else

        wbcCount = wbcCount + 1;
        wbcLabels(current) = wbcCount;

    end
end

wbcMask = imdilate(wbcLabels > 0, ...
    strel('disk',max(1,round(0.5*r))));

%% 5. RBC detection
rbc = cells & ~wbcMask;

rbc = imopen(rbc,strel('disk',2));
rbc = bwareaopen(rbc,round(0.3*typ));

D = imgaussfilt(bwdist(~rbc),1);

markers = imregionalmax(imhmax(D,max(1,0.2*r))) & rbc;

L = watershed(imimposemin(-D,markers));

rbcMask = (L > 0) & rbc;
rbcMask = bwareaopen(rbcMask,round(0.3*typ));

[~,rbcCount] = bwlabel(rbcMask,4);

%% 6. Results
fprintf('RBC count: %d\n',rbcCount);
fprintf('WBC count: %d\n',wbcCount);

figure('Name','Blood Cell Detection');

imshow(img);
hold on;

visboundaries(rbcMask,'Color','g','LineWidth',0.5);
visboundaries(wbcLabels > 0,'Color','b','LineWidth',1.5);

title(sprintf('RBC: %d (green) | WBC: %d (blue)', ...
    rbcCount,wbcCount));

hold off;

### Outputs:
<img width="605" height="342" alt="image" src="https://github.com/user-attachments/assets/b9385c53-8774-4256-bb52-a0bfc5054f4c" />


### References:
A few IEEE papers 

https://ieeexplore.ieee.org/document/7950584 

https://ieeexplore.ieee.org/document/9244549 

https://ieeexplore.ieee.org/document/8899602 

 

An article  

https://www.tandfonline.com/doi/full/10.1080/23311916.2024.2304484#abstract 
   
### Limitations and Future Work:
As of now we do not have a major novelty but there is scope for future work as it can be made more adaptive
