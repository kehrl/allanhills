function [I,x,y,info] = L8read(foldername,varargin) 
% The L8read function reads full-resolution Level 1 Landsat 8 images and has the 
% ability to subset the image before loading.   
% 
% Developed by Chad Greene, University of Texas, January 2016.
% https://www.mathworks.com/examples/matlab/community/19932-l8read-documentation#6

%% Syntax
% 
%  I = L8read(foldername) 
%  I = L8read(...,lat,lon)
%  I = L8read(...,x,y)
%  I = L8read(...,buffer)
%  I = L8read(...,'band',band) 
%  I = L8read(...,'stretch') 
%  [I,x,y,info] = L8read(...)
%  L8read(...)
% 
%% Description
% 
% I = L8read(foldername) returns a full natural color RGB image. Loading an entire Level 1 Landsat 
% 8 image will take a few seconds. The foldername is the name of the folder containing band-level
% tiff images.  The foldername is what Earth Explorer calls the Entity ID. 
% 
% I = L8read(...,lat,lon) loads only a region large enough to encompass all geo coordinates given
% by lat,lon. Specifying |lat,lon| requires Matlab's Mapping Toolbox, but specifying |x,y| does not. 
% 
% I = L8read(...,x,y) loads only a region large enough to encompass all projected coordinates given
% by x,y. Coordinates x and y are in meters and projection information can be determined by 
% exploring the information of one individual band tiff files with  imfinfo(filename). 
% 
% I = L8read(...,buffer) specifies a buffer around any specified lat,lon or x,y in units of 
% meters.  To center a 20 km wide image on a single point, specify the single point as lat,lon or 
% x,y and set the buffer value to 10e3 meters.  The buffer value may be a two-element vector 
% to specify different buffer sizes in width and height, respectively. 
% 
% I = L8read(...,'band',band) specifies a Landsat 8 band to load.  band can be any number 1 through
% 11, 'qa' for quality, or 'rgb' for a natural color RGB image.  Default band is 'rgb' which 
% concatenates bands 4, 3, and 2, respectively to create an RGB image. 
% 
% I = L8read(...,'stretch') uses the Image Processing Toolbox functions imadjust and stretchlim 
% to stretch the limits of the image, thereby maximizing contrast. 
% 
% [I,x,y,info] = L8read(...) returns x and y coordinates in meters as well as meta data info. 
% 
% L8read(...) without any outputs displays the image with imagesc.
% 
%% Examples
% From Earth Explorer (http://earthexplorer.usgs.gov/ Earth Explorer) I downloaded path 27, row 39 from August 31, 2013, 
% Entity ID LC80270392013243LGN00.  When you download the Level 1 data, it will be about a gigabyte of 
% zipped data and when you unzip the folder the foldername will be the same as the Earth Explorer Entity ID. 
% For a quick look at the entire scene, make sure your current folder contains the LC80270392013243LGN00 and
% try this: 
% 
%   L8read('LC80270392013243LGN00') 
% 
%% 
% If you prefer a 20 km wide image centered on the Texas State Capitol building, specify the Capitol 
% coordinates and a 10 km buffer around the center coordinates.  I'll also use the 'stretch' option
% to give the image more contrast: 
% 
%   L8read('LC80270392013243LGN00',30.274722,-97.740556,10e3,'stretch'); 
% 
%% 
% If you prefer to get the image data into the Matlab workspace, simply include function outputs. 
% Here we load band 5 which covers the near-infrared range.  And instead of specifying a single centerpoint, 
% I'll say we have 7 locations we're interested in that surround the Austin airport.  We only want to load
% enough data to surround the airport, plus an additional 3000 meters on each side, 1500 meters on top, and 
% 1500 meters below the seven points.  The seven points surrounding the airport may be specified in geo coordinates
% as above, or in meters, UTM Zone 14: 
% 
%   Xpts = [626975;626650;626609;627910;629659;630106;629293];
%   Ypts = [3343378;3341833;3340043;3339068;3339108;3341467;3343378]; 
% 
%    [I,x,y] = L8read('LC80270392013243LGN00',Xpts,Ypts,[3000 1500],'band',5); 
%    
%   figure
%   imagesc(x,y,I)
%   colormap(gray)
%   axis xy equal
%   hold on
%   plot(Xpts,Ypts,'gp')
% 
%% Author Info
% The L8read function and supporting documentation were written by Chad A. Greene of the University 
% of Texas at Austin's Institute for Geophysics (UTIG), January 2016. L8read uses Aslak Grinsted's 
% geoimread function to accomplish regional loading.  
% 
% See also: geotiffread, imread 


%% Error checks: 

narginchk(1,inf) 
assert(~isnumeric(foldername)==1,'Input error: foldername must be a string.') 
assert(exist(foldername,'dir')==7,['Error: Folder ',foldername,' not found.'])

%% Set defaults: 

band = 'rgb'; 
stretchlimits = false; 

%% Input parsing: 

if nargin>1
    tmp = strcmpi(varargin,'band'); 
    if any(tmp)
        band = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
    end
    
    % Stretch limits of image? 
    tmp = strcmpi(varargin,'stretch'); 
    if any(tmp) 
        stretchlimits = true; 
        assert(license('test','image_toolbox')==1,'Error: L8read''s stretch option requires Matlab''s Image Processing Toolbox.')
        varargin = varargin(~tmp); 
    end
    
end

%% Load Data: 

switch lower(band) 
    case {1,2,3,4,5,6,7,8,9,10,11} 
        Ifile = dir(fullfile(foldername,['*_B',num2str(band),'.TIF'])); 
        [I,x,y,info] = geoimread(fullfile(foldername,Ifile.name),varargin{:}); 
        
    case 'qa'
        Ifile = dir(fullfile(foldername,'*_BQA.TIF'));
        [I,x,y,info] = geoimread(fullfile(foldername,Ifile.name),varargin{:}); 
        
    case 'rgb'
        
        % Red band: 
        Rfile = dir(fullfile(foldername,'*_B4.TIF'));
        [R,x,y,info] = geoimread(fullfile(foldername,Rfile.name),varargin{:}); 

        % Green band: 
        Gfile = dir(fullfile(foldername,'*_B3.TIF'));
        G = geoimread(fullfile(foldername,Gfile.name),varargin{:}); 

        % Blue band: 
        Bfile = dir(fullfile(foldername,'*_B2.TIF'));
        B = geoimread(fullfile(foldername,Bfile.name),varargin{:}); 

        % Concatenate: 
        I = cat(3,R,G,B); 
        
    otherwise 
        error('Unrecognized Landsat band.') 
end

%% Adjust image: 

if stretchlimits
    I = imadjust(I,stretchlim(I)); 
end

%% Show image: 

if nargout==0
    imagesc(x,y,I); 
    axis xy image 
    
    clear I
end


end















%% Aslak Grinsed's geomread as a subfunction: 


function [A,x,y,I] = geoimread(filename,varargin)
%GEOIMREAD reads a sub region of a geotiff or geojp2 image.  
% 
% 
%% Syntax
%
% A = geoimread(filename) 
% A = geoimread(filename,xlim,ylim)
% A = geoimread(filename,latlim,lonlim)
% A = geoimread(...,buffer)
% [A,x,y,I] = geoimread(...)
% geoimread(...) 
%
% 
%% Description
%
% A = geoimread(filename) returns the full image given by a filename. This 
% syntax is equivalent to A = geotiffread(filename). 
% 
% A = geoimread(filename,xlim,ylim) limits the region of the geotiff file to 
% the limits given by xlim and ylim, which are map units (usually meters) relative
% to the data projection. For example, if the geotiff is projected in Texas Centric 
% Mapping System/Lambert Conformal coordinates, xlim and ylim will have units of 
% meters relative to the origin (100 W, 18 N). xlim and ylim can be multimensional, 
% in which case the limits of the map will be taken as the limits of the limits of 
% the distribution of all points in xlim, ylim.  
% 
% A = geoimread(filename,latlim,lonlim) if no values in xlim, ylim exceed 
% normal values of latitudes and longitudes, geoimread assumes you've entered
% limits in geographic coordinates of degrees latitude and longitude. The first 
% input is latitude, the second input is longitude.  
%
% A = geoimread(...,buffer) adds a buffer in map units (usually meters or feet) to the 
% limits of the region of interest.  This may be useful if you want to load an image 
% surrounding scattered lat/lon data.  If you'd like an extra 2 kilometers of image around
% your data, enter 2000 as the buffer.  If buffer is a two-element vector, the first
% element is applied to the left and right extents of the image, and the second element 
% is applied to the top and bottom extents of the image.    
%
% [A,x,y,I] = geoimread(...) also returns pixel center coordinates (x,y) of the 
% output image and a geotiff info structure I. I is a useful input for projfwd and projinv. 
%
% geoimread(...) without any outputs shows the output image A without loading 
% any data into the workspace.  
%
% 
%% Examples: 
% 
% % Show a whole geotiff: 
% geoimread('boston.tif');
% 
% % Compare results from above to a subset geotiff: 
% mapx = [765884 766035 766963]; % units are feet
% mapy = [2959218 2957723 2958972]; 
% geoimread('boston.tif',mapx,mapy)
% 
% % Or if you have coordinates in lat/lon and you want a 500 foot buffer: 
% lat = [42.3675288 42.3634246 42.3668397];
% lon = [-71.0940009 -71.0934685 -71.0900125];
% 
% geoimread('boston.tif',lat,lon,500);
%
%% Author Info: 
% 
% (c) Aslak Grinsted 2014 (http://www.glaciology.net/)
%     & Chad A. Greene (http://www.chadagreene.com/)
% 
%% 
% 
% See also GEOTIFFREAD, GEOTIFFINFO, PIXCENTERS, and PROJFWD. 

%%  CHAD'S CHANGES:
% The following changes were made by Chad A. Greene (http://chadagreene.com/)
% of the University of Texas Institute for Geophysics (UTIG) on Sept 17, 2014: 
% 
% * More input checking and error messages. 
% 
% * Clarified syntax and description in header. 
%
% * The fileparts line now only writes the file extension because other outputs went unused.  
% 
% * If geographic coordinate limits are inferred, they are now ordered latlim,lonlim. <-- **FUNCTIONALITY CHANGE** 
% 
% * Limits xlim, ylim or latlim, lonlim can be scalar, vector, or matrix-- xlim is now taken 
%   as xlim = [min(xlim(:)) max(xlim(:))]. This will save a small step if you have some data points
%   given by x,y.  Now you can simply enter your data coordinates and geoimread will
%   figure out the limits.  
% 
% * Output variable I now has correct corner coordinates for the subset image instead of
%   the previous version which returned corner coordinates of the full original image.  
% 
% * A buffer can be specified around input limits. 
% 
% 
% Syntax before the changes: 
% geoimread('myimage.tif',[min(easting)-buffer_m max(easting)+buffer_m],...
%     [min(northing)-buffer_m max(northing)+buffer_m]);
% 
% Syntax after the changes:
% geoimread('myimage.tif',easting,northing,buffer_m)
% 
% 

%TODO: support downsampling (ReductionLevel parameter in imread)
%TODO: use map2pix and latlon2pix instead of projfwd and pixcenters. more robust if it is a rotational coordinate system.


%% Set defaults: 

usegeocoordinates = false; 
returnNONsubsetImage = true; 
buffer_x = 0; 
buffer_y = 0; 

%% Input checks: 

% Check for mapping toolbox: 
% assert(license('test','map_toolbox')==1,'geoimread requires Matlab''s Mapping Toolbox.')

% Check file type: 
assert(isnumeric(filename)==0,'Input filename must be a string.') 
[~,~,ext] = fileparts(filename);
switch upper(ext)
    case {'.JP2' '.JPEG2000' '.GEOJP2'}
        I = jp2tiffinfo(filename);
    case {'.TIF' '.TIFF' '.GTIF' '.GTIFF'}
        I = robustgeotiffinfo(filename);
    otherwise
        error('Unrecognized image file type. Must be tif, tiff, gtif, gtiff, jp2, jpeg2000, or geojp2.')
end


% Parse optional inputs: 
if nargin>1
    returnNONsubsetImage = false; 
    assert(nargin>2,'If you specify an xlim or latlim, you must specify a corresponding ylim or lonlim.')
    
    % Parse limits: 
    xlimOrLatlim = varargin{1}(:); 
    ylimOrLonlim = varargin{2}(:); 
    
    assert(isnumeric(xlimOrLatlim)==1,'Input xlim or latlim must be numeric.')
    assert(isnumeric(ylimOrLonlim)==1,'Input ylim or lonlim must be numeric.')
    
    % Assume geo coordinates if no input limits exceed normal lat/lon values: 
    if max(abs(xlimOrLatlim))<=90 && max(abs(ylimOrLonlim))<=360
        usegeocoordinates = true; 
    end
    
    % Parse buffer: 
    if nargin>3
        buffer_m = varargin{3}; 
        assert(isnumeric(buffer_m)==1,'Buffer value must be either a scalar or a two-element array.')
        assert(numel(buffer_m)<3,'Buffer value must be either a scalar or a two-element array.')
        buffer_x = buffer_m(1); 
        if isscalar(buffer_m)
            buffer_y = buffer_m(1);
        else
            buffer_y = buffer_m(2); 
        end
    end
    
    if nargin>4
        error('Too many inputs in geoimread.') 
    end    
    
end


%% Begin work: 

% Get pixel coordinates of full (non-subset) image: 
[x,y]=robustpixcenters(I);

% Set xlim and ylim depending on user inputs: 
if returnNONsubsetImage
    xlimOrLatlim = x(:); 
    ylimOrLonlim = y(:); 
end

if usegeocoordinates
    % lat/lon limits switch to x/y limits here: 
    if ~strcmp(I.ModelType,'ModelTypeGeographic')
        assert(license('test','map_toolbox')==1,'Mapping toolbox needed to project between lat/lon limits and x,y limits. Specify limits in x,y coordinates.')
        [xlimOrLatlim,ylimOrLonlim]=projfwd(I,xlimOrLatlim,ylimOrLonlim);
    end
end


xlim = [min(xlimOrLatlim)-buffer_x max(xlimOrLatlim)+buffer_x]; 
ylim = [min(ylimOrLonlim)-buffer_y max(ylimOrLonlim)+buffer_y]; 


% Rows and columns of pixels to read: 
rows=find((y>=ylim(1))&(y<=ylim(2)));
cols=find((x>=xlim(1))&(x<=xlim(2)));




%% Display messages if region of interest is partly or wholly outside the image region:

if xlim(1)<min(x)||xlim(2)>max(x)
    disp('geoimread limits extend beyond the available image output in the x direction.')
end

if ylim(1)<min(y)||ylim(2)>max(y)
    disp('geoimread limits extend beyond the available image output in the y direction.')
end

if isempty(rows)||isempty(cols)
    error('No image coordinates can be found inside the specified limits.')
end

%% Load region of interest: 
reductionlevel=0;
if reductionlevel==0
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
else
    %% Load region of interest:
    dpx=2^reductionlevel;
    rows=round(rows/dpx);cols=round(cols/dpx);
    
    rows=sort(rows([1 end]));
    cols=sort(cols([1 end]));
end
x=x(cols(1):cols(end));
y=y(rows(1):rows(end));

A=imread(filename,'PixelRegion',{rows cols});


%% Update info structure to more accurately reflect the new image: 

if nargout == 4
    I.FileSize = numel(A); 
    I.Height = size(A,1); 
    I.Width = size(A,2); 
    try
        I.TiePoints.WorldPoints.X = x(1);
        I.TiePoints.WorldPoints.Y = y(1);
        I.SpatialRef.RasterSize = [size(A,1),size(A,2)];
        I.RefMatrix(3,1) = x(1);
        I.RefMatrix(3,2) = y(1);
        I.BoundingBox = [min(x) min(y); max(x) max(y)];
        I.CornerCoords.X = [min(x) max(x) max(x) min(x)];
        I.CornerCoords.Y = [max(y) max(y) min(y) min(y)];
        %TODO: check whether GTRasterTypeGeoKey is RasterPixelIsArea or RasterPixelIsPoint
        I.CornerCoords.Row = .5 + [0 0 size(A,1) size(A,1)]; %TODO: is this .5 always true?  
        I.CornerCoords.Col = .5 + [0 size(A,2) size(A,2) 0];
        [I.CornerCoords.Lat,I.CornerCoords.Lon] = projinv(I,I.CornerCoords.X,I.CornerCoords.Y);
        I.GeoTIFFTags.ModelTiepointTag(4) = x(1);
        I.GeoTIFFTags.ModelTiepointTag(5) = y(1);
        I.SpatialRef.XLimWorld = [min(x),max(x)];
        I.SpatialRef.YLimWorld = [min(y),max(y)];
    catch,end
end

%% Clean up: 

if nargout==0
    imshow(A,'XData',x,'YData',y)
    axis xy
    clear A x y I
end

end



%--- BELOW is to make it more robust if no mapping toolbox ---

function I = robustgeotiffinfo(fname)
if license('test','map_toolbox')
    I=geotiffinfo(fname);
else
    I=imfinfo(fname);
%     %TODO: generate home-made refmatrix(?)....
%     if isfield(tags, 'ModelTransformationTag') && numel(tags.ModelTransformationTag) >= 8
%         geoimread does not work for rotated systems
%         
%     else %use ModelPixelScaleTag instead
%         dx =  I.ModelPixelScaleTag(1); dy = -I.ModelPixelScaleTag(2);
%         x0 = I.ModelTiepointTag(4) - dx * I.ModelTiepointTag(1);
%         y0 = I.ModelTiepointTag(5) - dy * I.ModelTiepointTag(2);
%         J = [dx 0; 0 dy];
%     end
%     I.RefMatrix=[flipud(J); x0-J(1,1)-J(1,2)];
end
end

function [x,y]=robustpixcenters(I)
if license('test','map_toolbox')
    [x,y]=pixcenters(I);
else
    %I have not read documentation... but this only works for rectilinear systems.
    assert(I.ModelPixelScaleTag(3)==0,'unexpected ModelPixelScaleTag format.');
    assert(all(I.ModelTiepointTag(1:3)==0),'unexpected ModelTiepointTag format.');
    x=((0:I.Width-1)-I.ModelTiepointTag(1))*I.ModelPixelScaleTag(1)+I.ModelTiepointTag(4);
    y=((0:I.Height-1)-I.ModelTiepointTag(2))*-I.ModelPixelScaleTag(2)+I.ModelTiepointTag(5);
end
end
