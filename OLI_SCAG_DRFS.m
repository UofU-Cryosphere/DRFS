%% OLI-SCAG
clear; clc; close all

%Produces a delta vis geotiff (reduction in visible albedo due to light absorbing particles) for 
%Landsat-8 OLI scenes using grain size and snow fraction from OLI-SCAG.

%These files need to be in the path: multiband Landsat-8 reflectance file (.bip), 
% snow fraction file (.snow.tif), snow grainsize file (.grnsz.tif), and appropriate
% zenith angle look up table (.mat)

%Spatial reference data is found in basic scene metadata (.txt), the only
%thing that should need updating is the epsg (coordRefSysCode) code when moving 
%to a new scene in a different UTM zone

%Data Ranges: 0-100 Delta Vis, 245 Not Snow, 255 Off Grid/Not Run/No Data

%Load lookup table data
load lookup_table % variables generated from OLI.z30.LIB.csv

myDir=uigetdir; %choose parent folder
myFiles = dir(fullfile(myDir)); %all file names saved as struct
numfile = length(myFiles);
% 
for p=28:length(myFiles)

%insert file name here, excluding the file extension
%filename = 'LC08_L2SP_042034_20170616_20200903_02_T1';
filename = myFiles(p).name;

%get geotiff reference data to use when saving geotiff
[A,geodat] = readgeoraster([filename,'.snow.tif']);
[row col] = size(A);

%load OLI bands
OLI = multibandread([filename, '.bip'],[row col 6],'uint16',0,'bip','ieee-le');

%b1test = imread([filename,'.snow.tif']);
OLI_b1 = OLI(:,:,1);
OLI_b3 = OLI(:,:,3);
OLI_b4 = OLI(:,:,4);
OLI_b5 = OLI(:,:,5);

%load fractional snow cover
frac = imread([filename,'.snow.tif']);
frac = double(frac);

%load SCAG grain sizes
scag=imread([filename, '.grnsz.tif']);
scag = double(scag);
%generate find reflectance at each band in .tif based on lookup table



% finding the index for the grain size for use in the lookup table
[row col] = size(scag);
gs_ind = NaN(row,col); %pre allocate
for i = 1:row
    for j= 1:col
        if scag(i,j) == 65535 %set all nans to nan
            gs_ind(i,j) = NaN;
        elseif scag(i,j) == 0 %set all nans to nan
            gs_ind(i,j) = NaN;
        else
           gs_ind(i,j) = find(gs==scag(i,j)); 
        end
    end
end

%find lookup table reflectance
scag_refl = NaN(row,col,4); % preallocate 5 band reflectance from lookup table
for i = 1:row
    for j= 1:col
        for k=1:4 %4 bands
        if isnan(gs_ind(i,j))
        scag_refl(i,j,k)=NaN;    
        else
        scag_refl(i,j,k)=sim_refl(k,gs_ind(i,j))*10; %change order of magnitude to match OLI
        end
        end
    end
end

%Find reflectance difference at band 5 = wvl(4) = 0.87
diff = NaN(row,col);
for i = 1:row
    for j= 1:col
        if isnan(scag_refl(i,j,4))
            diff(i,j)=NaN;
        else
            diff(i,j) = scag_refl(i,j,4)-OLI_b5(i,j);
        end
    end
end

%add difference to OLI bands
OLI_diff=NaN(row,col,4);
for i = 1:row
    for j= 1:col
        if isnan(diff(i,j))
            OLI_diff(i,j)=NaN;
        else
            OLI_diff(i,j,1) = diff(i,j)+OLI_b1(i,j);
            OLI_diff(i,j,2) = diff(i,j)+OLI_b3(i,j);
            OLI_diff(i,j,3) = diff(i,j)+OLI_b4(i,j);
            OLI_diff(i,j,4) = diff(i,j)+OLI_b5(i,j);
        end
    end
end

%calculate scag - OLI for raditative forcing
radf=NaN(row,col,4);
for i = 1:row
    for j= 1:col
        for k =1:4 %4 bands
        if isnan(diff(i,j))
            radf(i,j,k)=NaN;
        else
            radf(i,j,k) = scag_refl(i,j,k) - OLI_diff(i,j,k);
        end
    end
    end
end

deltavis = mean(radf,3); %average bands
deltavis=deltavis/1000; %change to reflectance scale

%apply FSCA mask
%remove all pixels that are <90% FSCA
[row col] = size(frac);
mask=frac; %Before modifying frac use 0 and 255 values to place final mask
for i= 1:row
    for j=1:col
   if frac(i,j) < 90 || frac(i,j) > 100;
       frac(i,j) = NaN;
   end
    end
end

for i= 1:row
    for j=1:col
    if isnan(frac(i,j))
       deltavis(i,j) = 245;
   end
    end
end

%Use fractional snow cover to set no data values to 255 and no snow to 245
for i= 1:row
    for j=1:col
    if mask(i,j)==255
       deltavis(i,j) = 255;
    end
    if mask(i,j) == 0
        deltavis(i,j)=245;
   end
    end
end

%Any remaining NaN's to 255 (i.e. may contain snow, but no grnsz value)
for i= 1:row
    for j=1:col
    if isnan(deltavis(i,j))
       deltavis(i,j) = 255;
    end
    end
end


%plot for visualization
clim = [0 10];
figure(1)
imagesc(deltavis,clim);
figure(2)
histogram(deltavis)

%save as geotiff 
coordRefSysCode = 32642; %UTM 42N
% imwrite(uint8(deltavis),[filename '_deltavis.tiff'],'tiff')
geotiffwrite([filename '_deltavis.tif'],deltavis,geodat, 'CoordRefSysCode',coordRefSysCode);

%clear variables for next loop
clear A deltavis diff frac geodat mask OLI OLI_b1 OLI_b3 OLI_b4 OLI_b5 radf scag scag_refl row col

end
