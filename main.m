
%% Reading the image
clc
clear all
close all

inImg = StdIP.readImg3D(0.5,4);
[nr,nc,nd] = size(inImg);
StdIP.write3DStack(inImg,'original');

%% Segmentation with TuFF: Vessel Enhancement 
horizontal_red = 8;
vert_red = 2;
mask = zeros(size(inImg));
mask(horizontal_red:nr-horizontal_red,horizontal_red:nc-horizontal_red,vert_red:nd-vert_red) = 1;
enhImg = TuFF.brightVesselEnhanceFrangi3D(inImg,3,1).*mask;
enhImg = enhImg/max(enhImg(:));
StdIP.write3DStack(enhImg,'enhanced');

%% Segmentation with TuFF: Initialization
phi0 = TuFF.initRegionOtsu(enhImg,0.9,1000);
StdIP.imshow3D(phi0>0,0.4);

%% Segmentation with TuFF: Curve Evolution

param_Tuff.phi0         = phi0;
param_Tuff.enhI         = enhImg;
param_Tuff.max_iter     = 300;
param_Tuff.error        = 1;
param_Tuff.cnvg_cnt     = 5;
param_Tuff.interval     = 30;
param_Tuff.disp_Img     = [];
param_Tuff.magnify      = 200;
param_Tuff.col          = 'y';
%-- parameters to tune
param_Tuff.dt           = 1.5;
param_Tuff.pull_back    = 0.0;  % 1/0:do(not)use deflation; essential for nice shape
param_Tuff.p            = 2;    % exponent in 1/1+()^p
param_Tuff.smooth_trm   = 0.01;
param_Tuff.edge_trm     = 0.2;
param_Tuff.attr_trm     = 0.8;
param_Tuff.edge_range   = 7;
param_Tuff.attr_range   = 7;

TuFF_obj = TuFF(inImg);  %-- object for TuFF class
phi = TuFF_obj.runEATuFF3D(param_Tuff);

%%
bwI = TuFF.extractObject(phi,500,0,1); % phi,area_open, preserve_largest, perform_close
StdIP.imshow3D(bwI,0.4);
StdIP.write3DStack(bwI,'segmented');

%% Create the Graph from skeleton

skel = StdIP.find3DSkeleton(bwI);
% skeleton_graph = StdIP.skel2Graph3D(skel);

%%
StdIP.saveAllswc(skel,2,1,6);

%% Delete the objects

delete(TuFF_obj);



