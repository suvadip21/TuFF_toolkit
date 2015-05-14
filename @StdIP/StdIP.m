classdef StdIP
    %STDIP Class of standard image processing functions
    %   Includes various general purpose IP routine
    %   Author: Suvadip Mukherjee (sm5vp@virginia.edu)
    
    properties
    end
    
    methods(Static)
        % ---------------------------------------------------------------
        % -------------------- 2D routines ------------------------------
        % ---------------------------------------------------------------
        
        function [rgbImg, grayImg, fname, pname] = readImg2D(res)
            [fname,pname] = uigetfile('*.*','Input Data');
            rgbImg = imread(strcat(pname,fname));
            dim = ndims(rgbImg);
            rgbImg = imresize(rgbImg,res);
            if(dim == 3)
                grayImg = rgb2gray(rgbImg);
            else
                grayImg = rgbImg;
            end
            grayImg = im2double(grayImg);
            
        end
        
        function h = imageOverlay(im1,im2,mag,col)
            %-- Overlay im2 over im1. Generally, im2 is binary
            figure; imshow(im1,'InitialMagnification',mag); hold on;
            R = col(1)*im2; G = col(2)*im2; B = col(3)*im2;
            rgbim2 = cat(3,R,G,B);
            h = imshow(rgbim2,'InitialMagnification',mag); hold off;
            alpha_data = im2;
            set(h,'AlphaData',alpha_data);
        end
        
        % ---------------------------------------------------------------
        % -------------------- 3D routines ------------------------------
        % ---------------------------------------------------------------
        
        skel = find3DSkeleton(bwI)
        
        function skeleton_graph = skel2Graph3D(skel)
            %-- Create a graph from a skeleton. No smoothing is performed
            skel_pts = find(skel > 0);
            N = length(skel_pts);
            G = zeros(N);
            parent = zeros(N,1);
            for ii = 1 : N
                for jj = 1:N
                    p1 = skel_pts(ii);
                    p2 = skel_pts(jj);
                    if p1 ~= p2
                        [p1_r,p1_c,p1_d] = ind2sub(size(skel),p1);
                        [p2_r,p2_c,p2_d] = ind2sub(size(skel),p2);
                        if (abs(p1_r-p2_r)<=1) && (abs(p1_c-p2_c)<=1) && (abs(p1_d-p2_d)<=1)
                            G(ii,jj) = 1;
                        end
                    end
                end
            end

%             [disc, pred,~] = graphtraverse(sparse(G),1,'Depth',Inf,'Directed',false);
            node_list   = [1:N]';
            for ii = 1 : N
                for jj = ii+1:N
                    if G(ii,jj)
                        parent(jj) = ii;
                    end
                end
            end
            parent(1) = -1;
            t = find(parent == 0);
            parent(t) = -1;
            skeleton_graph = struct('SparseGraph',[],...
                'NodeList',[],...
                'ParentList',[],...
                'NodePosition',[]);
            skeleton_graph.SparseGraph  = sparse(G);
            skeleton_graph.NodeList     = node_list;
            skeleton_graph.ParentList   = parent;
            skeleton_graph.NodePosition = skel_pts;
        end
        
        
        
        function Img = readImg3D(sliceReduce,stackSkip)
            dname = uigetdir('','Input the image folder');
            D = dir(strcat(dname,'\*.tif'));
            buffer = double(imread(strcat(dname,'\',D(1).name)));
            buffer = imresize(buffer,sliceReduce);
            [maxr,maxc] = size(buffer);
            Img = zeros(maxr,maxc,length(1:stackSkip:length(D)));
            count=0;
            for I=1:stackSkip:length(D),
                count=count+1;
                buffer = im2double(imread(strcat(dname,'\',D(I).name)));
                buffer = imresize(buffer,sliceReduce);
                Img(:,:,count) = buffer;
                
            end
            Img = (Img-min(Img(:)))/(max(Img(:))-min(Img(:)));
        end
        
        
        
        function [] = imshow3D(Img,isoval)
            theta1 = -150;
            theta2 = 50;
            [nr,nc,nd] = size(Img);
            lims = [1 nc 1 nr 1 nd];
            [x,y,z] = meshgrid(1:nc,1:nr,1:nd);
            data = Img ;
            cdata = smooth3(ones(size(data)),'box',7);
            p = patch(isosurface(x,y,z,data,isoval));
            isonormals(x,y,z,data,p)
            isocolors(x,y,z,cdata,p)
            p.FaceColor = 'interp';
            p.EdgeColor = 'none';
            view(theta1,theta2)
            daspect([1 1 1.5])
            axis(lims);
        end
        
        
        
        function write3DStack(Img,newDir)
            
            message = ['Write The Stacks (' newDir ')'];
            dname = uigetdir('./Results',message);
            [~,msg,~] = mkdir(dname,newDir);
            p_write = strcat(dname,'\',newDir,'\');
            if ~isempty(msg), % directory already exists
                delete(strcat(p_write,'*.tif')); % remove the exixting tiffs
            end
            fname = 'slice';
            for f = 1:1:size(Img,3)
                buffer = Img(:,:,f);
                imwrite(uint8(buffer.*255),[p_write fname num2str(f) '.tif']);
            end
        end
        
        
        function [] = saveASswc(skel,skeleton_graph)
            %--[node_id type x y z rad parent]
            [row,col,depth] = size(skel);
            node_id = skeleton_graph.NodeList;
            [y,x,z] = ind2sub(size(skel),skeleton_graph.NodePosition);
            y = mod(-round(y),row+1);
            type = 4*ones(length(x),1);
            rad = 1*ones(length(x),1);
            parent = skeleton_graph.ParentList;
            swc_mat = [node_id type x y z rad parent];
            [fname,pname]=uiputfile('*.swc','Save SWC File as');
            p_write = strcat(pname,fname);
            fid = fopen(p_write,'w');
            header='## n,type,x,y,z,radius,parent';
            fprintf(fid,'%s\n\n',header);
            for ii = 1 : length(x)
                data = swc_mat(ii,:);
                fprintf(fid,'%s\n',num2str(data));
            end
            fclose(fid);
        end
        
        
    end
    
end

