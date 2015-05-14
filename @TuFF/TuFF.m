classdef TuFF < handle
    %TUFF
    % Code for enhancement via Frangi
    % Code for enhancement via LDE
    % Code for segmentation using TuFF (IEEE TIP,14)
    % Code for segmentation using EATuFF (IEEE Trans. ?)
    
    properties
        Img;   %-- original signal, not the vesselness
    end
    
    methods
        
        function obj = TuFF(I)
            obj.Img = I;
        end
        
        segI   = runTuFF2D(obj,param_Tuff)
        segI   = runEATuFF2D(obj,param_Tuff)
        segI   = runEATuFF3D(obj,param_Tuff)
        segI   = runFastEATuFF2D(obj,param_Tuff)
        
    end % End methods
    
    
    methods(Static)
       
        phi0   = initRegionOtsu(Img,level,min_area)
        phi0   = initRegionManual(Img,type,magnify)

        enhImg = brightVesselEnhanceFrangi(I,scales, gamma)
        enhImg = darkVesselEnhanceFrangi(I,scales, gamma)
        enhImg = vesselEnhanceLDE(vess_param)
        
        function enhImg = brightVesselEnhanceFrangi3D(Img,scales,gamma)
            ImStack = [];
            for I=1:length(scales),
                s = scales(I);
                filterRad = 4*s;
                [x1,x2,x3] = ndgrid(-filterRad:1:filterRad, -filterRad:1:filterRad, -filterRad:1:filterRad);
                gaussFilt = exp((-x1.^2-x2.^2-x3.^2)/(2*s^2));
                gaussFilt = gaussFilt/sum(gaussFilt(:));
                filtImg = imfilter(Img,gaussFilt,'same');
                ImStack = [ImStack s^gamma*TuFF.vesselness3D(filtImg)];
                disp(strcat('scale = ',num2str(s),' in imEnhance computed '));
            end
            ImStack = ImStack';
            eImg = max(ImStack,[],1);
            enhImg = reshape(eImg,size(Img));
        end
        
        vslI = vesselness3D(Img)
        
        F = attrForce(bwI,range)
        F = attrForce3D(bwI,range)
        Ker = AM_VFK(dim,rad,type,gamma)
        
               
        % -- extract the neuron segment from the sdf
        function bwI = extractObject(phi,min_area,preserve_one_cc,do_close)
            
            xcut  = 1;
            ycut  = 1;
            zcut  = 1;
            initial_segment = phi >= 0;
            bwI = bwareaopen(initial_segment,min_area);
            maskImg = zeros(size(phi));
            [nrow,ncol,nslice] = size(phi);
            maskImg(xcut+1:nrow-xcut,ycut+1:ncol-ycut,zcut+1:nslice-zcut) = 1;
            bwI = bwI.*maskImg;
            
            if do_close
                sz   = 1;
%                 bwI = imopen(bwI,ones(sz,sz,sz));
                bwI = imclose(bwI,ones(sz,sz,sz));
            end
            if preserve_one_cc
                CC    = bwconncomp(bwI);      % preserve the largest component only
                numCC = CC.NumObjects;
                
                t = 0;
                for i = 1 : numCC
                    sz = length(CC.PixelIdxList{i});
                    if sz > t
                        t = sz ;
                    end
                end
                bwI = bwareaopen(bwI,t-2,26);
            end
            
        end
        
        
        
        function phi0 = computeSDF(bwI)
            % inside >= 0, outside < 0
            lsf = bwdist(bwI)-bwdist(1-bwI)+im2double(bwI)-.5;
            phi0 = -double(lsf); 
        end 
        
        
        
        function c = convergence(p_mask,n_mask,thresh,c)
            diff = p_mask - n_mask;
            n_diff = sum(abs(diff(:)));
            if n_diff < thresh
                c = c + 1;
            else
                c = 0;
            end
        end
        
        %-------------------------------------------------------------------------
        %               Check boundary condition
        function g = NeumannBoundCond(f)
            
            [nrow,ncol] = size(f);
            g = f;
            g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
            g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
            g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
        end
        
        function g = NeumannBoundCond3D(f)
            
            [nrow,ncol,depth] = size(f);
            g = f;
            g([1 nrow],[1 ncol],[1 depth]) = f([2 nrow-1],[2 ncol-1],[2 depth-1]);
            g(2:nrow-1,[1 ncol],[1 depth]) = f(2:nrow-1,[2 ncol-1],[2 depth-1]);
            g([1 nrow],2:ncol-1,[1 depth]) = f([2 nrow-1],2:ncol-1,[2 depth-1]);
            g([1 nrow],[1 ncol], 2:depth-1) = f([2 nrow-1],[2 ncol-1], 2:depth-1);
        end
        
        
        %-------------------------------------------------------------------------
        %          Reinitialize LSF by Sussman reinitialization method   
        
        function [D] = SussmanReinitLS(D,dt)
     
            %D  : level set function
            a = D - shiftR(D); % backward
            b = shiftL(D) - D; % forward
            c = D - shiftD(D); % backward
            d = shiftU(D) - D; % forward
            
            a_p = a;  a_n = a; % a+ and a-
            b_p = b;  b_n = b;
            c_p = c;  c_n = c;
            d_p = d;  d_n = d;
            
            a_p(a < 0) = 0;
            a_n(a > 0) = 0;
            b_p(b < 0) = 0;
            b_n(b > 0) = 0;
            c_p(c < 0) = 0;
            c_n(c > 0) = 0;
            d_p(d < 0) = 0;
            d_n(d > 0) = 0;
            
            dD = zeros(size(D));
            D_neg_ind = find(D < 0);
            D_pos_ind = find(D > 0);
            dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
            dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
            
            D = D - dt .* TuFF.sussman_sign(D) .* dD;
        end
        
        
        function [D] = SussmanReinitLS3D(D,dt)
     
            %D  : level set function
            [nr,nc,nz] = size(D);
            a = D - cat(2, D(:,1,:), D(:,1:nc-1,:)); % backward
            b = cat(2, D(:,2:nc,:), D(:,nc,:)) - D;  % forward
            c = D - cat(1, D(1,:,:), D(1:nr-1,:,:)); % down
            d = cat(1, D(2:nr,:,:), D(nr,:,:)) - D;  % down
            e = D - cat(3, D(:,:,1), D(:,:,1:nz-1));   % away- z dirn
            f = cat(3, D(:,:,2:nz), D(:,:,nz)) - D;  % towards - z dirn
            
            a_p = a;  a_n = a; % a+ and a-
            b_p = b;  b_n = b;
            c_p = c;  c_n = c;
            d_p = d;  d_n = d;
            e_p = e;  e_n = e;
            f_p = f;  f_n = f;
            
            a_p(a < 0) = 0;
            a_n(a > 0) = 0;
            b_p(b < 0) = 0;
            b_n(b > 0) = 0;
            c_p(c < 0) = 0;
            c_n(c > 0) = 0;
            d_p(d < 0) = 0;
            d_n(d > 0) = 0;
            e_p(e < 0) = 0;
            e_n(e > 0) = 0;
            f_p(f < 0) = 0;
            f_n(f > 0) = 0;
            
            dD = zeros(size(D));
            D_neg_ind = find(D < 0);
            D_pos_ind = find(D > 0);
            dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                            + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2) ...
                            + max(e_p(D_pos_ind).^2, f_n(D_pos_ind).^2)) - 1;
            dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                            + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)...
                            + max(e_n(D_neg_ind).^2, f_p(D_neg_ind).^2)) - 1;
            
            D = D - dt .* TuFF.sussman_sign(D) .* dD;
        end
        
        
        function S = sussman_sign(D)
            S = D ./ sqrt(D.^2 + 1);
        end
        
        %-------------------------------------------------------------------------
        %               Compute curvature of a function u(lsf)
        
        function k = compute_divergence(u)
            % Computes div(\nabla u|\nabla u|)
            [ux,uy] = gradient(u);
            normDu = sqrt(ux.^2+uy.^2+1e-10);	
            Nx = ux./normDu;
            Ny = uy./normDu;
            nxx = gradient(Nx);
            [~,nyy] = gradient(Ny);
            k = nxx+nyy;                        
        end
   
        
        function k = compute_divergence3D(u)
            % Computes div(\nabla u|\nabla u|)
            [ux,uy,uz] = gradient(u);
            normDu = sqrt(ux.^2+uy.^2+uz.^2+1e-10);	
            Nx = ux./normDu;
            Ny = uy./normDu;
            Nz = uz./normDu;
            nxx = gradient(Nx);
            [~,nyy,~] = gradient(Ny);
            [~,~,nzz] = gradient(Nz);
            k = nxx+nyy+nzz;                        
        end
        
        
        %-------------------------------------------------------------------------
        %                   Display the evolving curve
        
        function showCurveAndPhi(phi,magnify,Img,cl)
            
            imshow(Img,[],'InitialMagnification',magnify);
            hold on;
            [c,h] = contour(phi,[0 0],cl,'Linewidth',3); hold off;
            test = isequal(size(c,2),0);
            while (test==false)
                s = c(2,1);
                if ( s == (size(c,2)-1) )
                    t = c;
                    hold on; plot(t(1,2:end)',t(2,2:end)',cl,'Linewidth',3);
                    test = true;
                else
                    t = c(:,2:s+1);
                    hold on; plot(t(1,1:end)',t(2,1:end)',cl,'Linewidth',3);
                    c = c(:,s+2:end);
                end
            end
        end
        
    end % End static methods
    
end % End class TuFF


%-------------------------------------------------------------------------
%                       Utility functions
%-------------------------------------------------------------------------

function shift = shiftD(M)
shift = shiftR(M')';
end
%-------------------------------------------------------------------------

function shift = shiftL(M)
shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];
end
%-------------------------------------------------------------------------

function shift = shiftR(M)
shift = [ M(:,1) M(:,1:size(M,2)-1) ];
end
%-------------------------------------------------------------------------

function shift = shiftU(M)
shift = shiftL(M')';
end
%--------



