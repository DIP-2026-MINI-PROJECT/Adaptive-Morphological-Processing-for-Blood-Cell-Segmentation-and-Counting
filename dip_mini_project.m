clc; 
clear all; 
close all;

%%\ Read image
img = imread("10.jpeg");
if size(img,3) == 4
    img = img(:,:,1:3); 
end
if size(img,3) ~= 3
    error('Please use an RGB blood smear image.'); 
end
img = im2uint8(img);

R = imgaussfilt(img(:,:,1),1);
G = imgaussfilt(img(:,:,2),1);
B = imgaussfilt(img(:,:,3),1);

%%\ 1. All cells: automatic threshold on green channel
cells = G < 255*graythresh(G);
if mean(cells(:)) > 0.8
    cells = ~cells; 
end        % safety if background is dark
cells = imfill(imclose(cells,strel('disk',2)),'holes');
cells = imopen(cells,strel('disk',2));
cells = bwareaopen(cells,max(30,round(0.0002*numel(G))));

%%\ 2. Estimate size of one cell from the image itself
s = regionprops(cells,'Area','Solidity');
a = [s.Area];
singles = a([s.Solidity] > 0.9);
if numel(singles) >= 3
    area = median(singles);
else
    area = median(a); 
end
r = sqrt(area/pi);

%%\ 3. WBCs: dark nuclei that are purple/blue (not pink)
nuc = cells & (R < 255*graythresh(R(cells)));

RB    = double(R) - double(B);          % pink cells: large, purple nuclei: small
refRB = median(RB(cells));
lab   = bwlabel(nuc);
st    = regionprops(lab,RB,'Area','MeanIntensity');
keep  = find([st.Area] >= 0.3*area & [st.MeanIntensity] < 0.5*refRB);
nuc   = ismember(lab,keep);
nuc   = imfill(imclose(nuc,strel('disk',max(1,round(0.2*r)))),'holes');

% split touching nuclei, then grow EACH nucleus separately (no merged mask)
nucLab   = splitTouching(nuc, max(1,0.15*r));
nucLab   = nucLab .* bwareaopen(nucLab>0, round(0.2*area));
wbcLab   = imdilate(nucLab, strel('disk',max(1,round(0.6*r))));
wbcMask  = wbcLab > 0;
wbcCount = numel(unique(nucLab(nucLab>0)));

%%\ 4. RBCs: remaining cells, touching ones split by watershed
rbc = imopen(cells & ~wbcMask, strel('disk',2));
rbc = bwareaopen(rbc, round(0.3*area));
rbcLab = splitTouching(rbc, max(1,0.2*r));
rbcMask = bwareaopen(rbcLab > 0, round(0.3*area));
[~,rbcCount] = bwlabel(rbcMask,4);

%%\ Results
fprintf('RBC count: %d\n',rbcCount);
fprintf('WBC count: %d\n',wbcCount);

figure('Name','Steps');
subplot(1,3,1); imshow(img);  title('Original');
subplot(1,3,2); imshow(cells); title('Cell mask');
subplot(1,3,3); imshow(nuc);  title('WBC nuclei');

figure('Name','Result');
imshow(img); hold on;
visboundaries(rbcMask,'Color','g','LineWidth',0.5);
for k = unique(wbcLab(wbcLab>0))'               % one outline per WBC
    visboundaries(wbcLab==k,'Color','b','LineWidth',1.5);
end
title(sprintf('RBC: %d (green) | WBC: %d (blue)',rbcCount,wbcCount));

%%\ Helper: split touching blobs with a distance-transform watershed
function L = splitTouching(mask,h)
D = imgaussfilt(bwdist(~mask),1);
markers = imregionalmax(imhmax(D,h)) & mask;
W = watershed(imimposemin(-D,markers));
L = bwlabel((W>0) & mask, 4);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
