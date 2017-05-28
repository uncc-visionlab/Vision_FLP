classdef RGBDImage < handle
    properties (SetAccess = private)
        % RGBD camera parameters
        width                     % width of the image
        height                    % height of the image
        numPixels                 % number of pixels in the image
        region_of_interest        % top-left and bottom-right corners define
                                  % a rectangular subset of the image for processing
        MAXRANGE                  % maximum range (in m.) of the sensor
        % The distortion parameters, size depending on the distortion model.
        % For "plumb_bob", the 5 parameters are: (k1, k2, t1, t2, k3).
        D                         % vector of distortion coefficients
        % Intrinsic camera matrix for the raw (distorted) images.
        %     [fx  0 cx]
        % K = [ 0 fy cy]
        %     [ 0  0  1]
        K
        % Rectification matrix (stereo cameras only)
        % A rotation matrix aligning the camera coordinate system to the ideal
        % stereo image plane so that epipolar lines in both stereo images are
        % parallel.
        R
        % Projection/camera matrix
        %     [fx'  0  cx' Tx]
        % P = [ 0  fy' cy' Ty]
        %     [ 0   0   1   0]
        % By convention, this matrix specifies the intrinsic (camera) matrix
        %  of the processed (rectified) image. That is, the left 3x3 portion
        %  is the normal camera intrinsic matrix for the rectified image.
        % It projects 3D points in the camera coordinate frame to 2D pixel
        %  coordinates using the focal lengths (fx', fy') and principal point
        %  (cx', cy') - these may differ from the values in K.
        % For monocular cameras, Tx = Ty = 0. Normally, monocular cameras will
        %  also have R = the identity and P[1:3,1:3] = K.
        P
        f                          % focal length (in pixels) in (x,y)
        inv_f
        c                          % principal/center point of image
        
        % image of angle weights from RGBD camera calibration parameters
        tan_theta_x
        tan_theta_y
        
        % integral images from RGBD camera calibration parameters
        iImg_tan_theta_x
        iImg_tan_theta_y
        iImg_tan_theta_x_sq
        iImg_tan_theta_y_sq
        iImg_tan_theta_xy
        % for RGBD range fitting and depends on measured data
        iImg_ix
        iImg_iy
        iImg_iz
        iImg_iz_sq
        
        % for XYZ point cloud fitting and depends on measured data
        I                         % sensed depth image
        iImg_numInvalidPts        % integral image counting invalid points
        X
        Y
        Z
        iImg_x
        iImg_y
        iImg_z
        
        % grid decomposition variables
        blocksize                 % size of a block in (x,y) pixels
        x_numfaces                % num of blocks in x-direction
        y_numfaces                % num of blocks in y-direction
        numfaces                  % num of blocks
        faceXPts                  % x coordinates of block corners
        faceYPts                  % y coordinates of block corners
    end
    
    methods (Static)
        function demo()
            % Pre-computed values and constants
            width=640;
            height=480;
            D=[0.0, 0.0, 0.0, 0.0, 0.0];
            K= [567.8377685546875, 0.0, 319.5, 0.0, 567.8377685546875, 239.5, 0.0, 0.0, 1.0];
            R= [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            P= [567.8377685546875, 0.0, 319.5, 0.0, 0.0, 567.8377685546875, 239.5, 0.0, 0.0, 0.0, 1.0, 0.0];
            K=reshape(K,3,3);
            R=reshape(R,3,3);
            P=reshape(P,4,3)';
            rgbdImg = RGBDImage();
            rgbdImg.setCameraParameters(width,height,D,K,R,P);
            rgbdImg.initializeEstimationMatrices();
            
            blocksize=70;
            rgbdImg.imposeGridDecomposition(blocksize);
            
            % read an image from the disk / sensor
            fid=fopen('image0.npy');
            val=fread(fid,rgbdImg.numPixels,'single');
            I=reshape(val,width,height)';
            
            rgbdImg.setImage(I);
            
            figure(1), hold on, imshow(I,[0 rgbdImg.MAXRANGE]);
            colormap('jet');
            
            MAXLEVEL=4;
            ERROR_THRESHOLD=1.0e-3/8;
            %ERROR_THRESHOLD=1.0e-3/15;
            
            planeCoeffs = rgbdImg.getPlanes(MAXLEVEL, ERROR_THRESHOLD);
            
        end
    end % Static methods
    
    methods
        function obj = RGBDImage()
            obj.MAXRANGE = 6;
        end
        
        function setCameraParameters(obj, width, height, D, K, R, P)
            obj.width = width;
            obj.height = height;
            obj.numPixels=width*height;
            obj.D = D;
            obj.K = K;
            obj.R = R;
            obj.P = P;
            f=diag(P);            
            obj.f=f(1:2);
            obj.inv_f=1./f;
            obj.c=P(1:2,3);
        end
        
        % selects a rectangular subset of the image to process
        % Specified as the 2x2 matrix [[x1,y1]; [x2,y2]] where 
        % [x1,y1] = top-left-corner
        % [x2,y2] = bottom-right-corner
        function setRegionOfInterest(obj, region_of_interest)
            obj.region_of_interest = region_of_interest;
        end
        
        function initializeEstimationMatrices(obj)
            [x,y] = meshgrid(0:obj.width-1,0:obj.height-1);
            obj.tan_theta_x = (x-obj.c(1))/obj.f(1);
            obj.tan_theta_y = (y-obj.c(2))/obj.f(2);
            obj.iImg_tan_theta_x = myIntegralImage(obj.tan_theta_x);
            obj.iImg_tan_theta_y = myIntegralImage(obj.tan_theta_y);
            obj.iImg_tan_theta_x_sq = myIntegralImage(obj.tan_theta_x.^2);
            obj.iImg_tan_theta_y_sq = myIntegralImage(obj.tan_theta_y.^2);
            obj.iImg_tan_theta_xy = myIntegralImage(obj.tan_theta_x.*obj.tan_theta_y);
        end
        
        function p2 = project(obj, p3)
            p2 = Point2d(0,0);
            p2.x = p3.x / (p3.z * obj.inv_f(1)) + obj.c(1);
            p2.y = p3.y / (p3.z * obj.inv_f(2)) + obj.c(2);
        end
        
        function lineseg2 = getVisibleLineSegment2d(obj, line3) 
            lineseg2 = LineSegment2d(Point2d(0,0),Point2d(0,0),0,0);
            ip3 = Point3d(0,0,0);
            %            if (abs(line3.v.x) > abs(line3.v.y)) {
            lambda = -line3.p0.x / line3.v.x;
            ip3 = line3.setPoint(ip3, lambda); % intersect w/ x=0 plane
            %            } else {
            %                lambda = -line3.p0.y / line3.v.y;
            %                line3.setPoint(ip3, lambda); // intersect w/ y=0 plane
            %            }
            lineseg2.p0 = obj.project(ip3); % project pt at intersection into 2D image
            lineseg2.v.x = line3.v.x;
            lineseg2.v.y = line3.v.y;

            ip = Point2d(0,0);
            lineseg2.l_start = -lineseg2.p0.x / lineseg2.v.x; % guess line passes through left side of image
            ip = lineseg2.setPoint(ip, lineseg2.l_start);
            if (ip.y < 0 || ip.y > obj.height)  % line does not pass through left side of image
                lineseg2.l_start = -lineseg2.p0.y / lineseg2.v.y; % guess line passes through top of image
                ip = lineseg2.setPoint(ip, lineseg2.l_start);
                if (ip.x < 0 || ip.x > obj.width)  % line does not pass through top side of image
                    lineseg2.l_start = (obj.height - 1 - lineseg2.p0.y) / lineseg2.v.y; % get right intersection point
                    lineseg2.l_end = (obj.width - 1 - lineseg2.p0.x) / lineseg2.v.x; % get bottom intersection point
                    return; % line through bottom and right edges
                end
                lineseg2.l_end = (obj.height - 1 - lineseg2.p0.y) / lineseg2.v.y; % guess line passes through bottom side of image
                ip = lineseg2.setPoint(ip, lineseg2.l_end);
                if (ip.x < 0 || ip.x > obj.width) % line does not pass through bottom side of image
                    lineseg2.l_end = (obj.width - 1 - lineseg2.p0.x) / lineseg2.v.x;
                    return; % line through top and right edges
                else 
                    % line through top and bottom edges (already computed))
                    return;
                end
            else  % guess line passes through right side of image
                lineseg2.l_end = (obj.width - 1 - lineseg2.p0.x) / lineseg2.v.x;
                ip = lineseg2.setPoint(ip, lineseg2.l_end);
                if (ip.y < 0 || ip.y > obj.height) % wrong
                    lineseg2.l_end = -lineseg2.p0.y / lineseg2.v.y; % guess top edge
                    ip = lineseg2.setPoint(ip, lineseg2.l_end);
                    if (ip.x < 0 || ip.x > obj.width) % wrong
                        lineseg2.l_end = (obj.height - 1 - lineseg2.p0.y) / lineseg2.v.y; % has to be bottom edge
                    end
                end
            end
            %std::cout << "Visible 2D line segment = " << lineseg2 << " lineseg2.start " << lineseg2.start << std::endl;
            %std::cout << "Visible 3D line segment = " << line3 
            %        << " start " << p0 << " start 2d " << project(p0)
            %        << " end " << p1  << " end 2d " << project(p1)
            %        << std::endl;
        end
        
        function setImage(obj, I)
            %I(I==0)=1000;
            %[dx,dy]=gradient(I);
            %gradMag = sqrt(dx.^2+dy.^2);
            %gradMag(gradMag>10)=0;
            %imshow(log10(gradMag+1),[]);
            %colormap(jet);
            %hold on, imagesc(I*50); colormap(jet);
            obj.I = I;
            obj.I(obj.I<0.1)=NaN;
            obj.I(obj.I>obj.MAXRANGE)=NaN;
            obj.I=imgaussfilt(obj.I,1);
            [obj.iImg_z, obj.iImg_numInvalidPts] = myIntegralImage(obj.I,1);
            iX=obj.tan_theta_x./obj.I;
            iY=obj.tan_theta_y./obj.I;
            %I(isnan(I))=0;
            obj.X=obj.tan_theta_x.*obj.I; % obj.X has NaN values
            obj.Y=obj.tan_theta_y.*obj.I; % obj.Y has NaN values
            obj.Z=obj.I;
            obj.iImg_x = myIntegralImage(obj.X);
            obj.iImg_y = myIntegralImage(obj.Y);
            %obj.iImg_z = myIntegralImage(obj.Z);
            obj.iImg_ix = myIntegralImage(iX);
            obj.iImg_iy = myIntegralImage(iY);
            obj.iImg_iz = myIntegralImage(1./obj.I);
            obj.iImg_iz_sq = myIntegralImage(1./(obj.I.^2));
            % S=[reshape(X,numpixels,1) reshape(Y,numpixels,1) reshape(Z, numpixels,1)];
        end
        
        function imposeGridDecomposition(obj, blocksize)
            aratio=obj.height/obj.width;  % aspect ratio
            obj.blocksize = blocksize;
            if (prod(size(obj.region_of_interest))~=4)
                margin=round(min(1*[aratio*blocksize,blocksize],[aratio*60,60]));
                obj.region_of_interest=[margin(2) margin(1); ...
                    obj.width-margin(2) obj.height-margin(1)];
            end
            [iy_cpt,ix_cpt]=meshgrid( ...
                obj.region_of_interest(1,2):blocksize:obj.region_of_interest(2,2), ...
                obj.region_of_interest(1,1):blocksize:obj.region_of_interest(2,1));
            [x_res, y_res] = size(ix_cpt);
            num_cpt = x_res*y_res;
            ix_cpt = reshape(ix_cpt,num_cpt,1);
            iy_cpt = reshape(iy_cpt,num_cpt,1);
            %I_idx = sub2ind(size(obj.I), iy_cpt, ix_cpt);
            obj.x_numfaces = x_res-1;
            obj.y_numfaces = y_res-1;
            obj.numfaces=obj.x_numfaces*obj.y_numfaces;
            fcpt_idxs=zeros(obj.numfaces,5);
            obj.faceXPts=zeros(obj.numfaces,5);
            obj.faceYPts=zeros(obj.numfaces,5);
            fidx=1;
            for yface=1:obj.y_numfaces
                for xface=1:obj.x_numfaces
                    idx = xface + (yface-1)*obj.x_numfaces;
                    fcpt_idxs(idx,:) = yface - 1 + ...
                        [idx idx+1 idx+obj.x_numfaces+2 idx+obj.x_numfaces+1 idx];
                    obj.faceXPts(fidx,:) = ix_cpt(fcpt_idxs(idx,:));
                    obj.faceYPts(fidx,:) = iy_cpt(fcpt_idxs(idx,:));
                    fidx = fidx + 1;
                end
            end
        end
        
        function subDivideFace(obj, faceIdx)
            xVals=obj.faceXPts(faceIdx,:);
            yVals=obj.faceYPts(faceIdx,:);
            midpt=zeros(2,4);
            for edgeIdx=1:4
                midpt(:,edgeIdx)=round(0.5*[(xVals(edgeIdx)+xVals(edgeIdx+1)); ...
                    (yVals(edgeIdx)+yVals(edgeIdx+1))]);
            end
            centerpt=round(mean([xVals(1:4)' yVals(1:4)']))';
            newFace=[[xVals(1); yVals(1)] ...
                midpt(:,1) centerpt midpt(:,4) ...
                [xVals(1); yVals(1)]];
            %hold on, plot(newFace(1,:),newFace(2,:),'y');
            obj.faceXPts=[obj.faceXPts; newFace(1,:)];
            obj.faceYPts=[obj.faceYPts; newFace(2,:)];
            newFace=[midpt(:,1) [xVals(2); yVals(2)] midpt(:,2) centerpt midpt(:,1)];
            obj.faceXPts=[obj.faceXPts; newFace(1,:)];
            obj.faceYPts=[obj.faceYPts; newFace(2,:)];
            newFace=[midpt(:,4) centerpt midpt(:,3) [xVals(4); yVals(4)] midpt(:,4)];
            obj.faceXPts=[obj.faceXPts; newFace(1,:)];
            obj.faceYPts=[obj.faceYPts; newFace(2,:)];
            newFace=[centerpt midpt(:,2) [xVals(3); yVals(3)] midpt(:,3) centerpt];
            obj.faceXPts=[obj.faceXPts; newFace(1,:)];
            obj.faceYPts=[obj.faceYPts; newFace(2,:)];
        end
        
        function [planeCoeffs,areas,locations, errors] = getPlanes(obj, MAXLEVEL, ERROR_THRESHOLD)
            USE_INTEGRAL_IMAGES = false;
            FAST_FIT = USE_INTEGRAL_IMAGES;
            faceIdx = 1;
            level = 1;
            %ptsPerEdge=4;
            nfaces(level)=size(obj.faceXPts,1);
            processedFaces = 0;
            planeCoeffs = [];
            areas = [];
            locations = [];
            % continue processing blocks until done
            while(nfaces(level) > 0 && level <= MAXLEVEL)
                %loc_blocksize=obj.blocksize/2^(level-1);
                %pathskip=max(1,round(obj.blocksize/ptsPerEdge));
                lastFaceIdx=sum(nfaces(1:level));
                % processes all the blocks for a given resolution/level
                while(faceIdx <= lastFaceIdx)
                    if (mod(faceIdx-processedFaces-1,nfaces(level)/10)==0)
                        fprintf(1,'Processing face %d of %d faces in level %d.\n',...
                            faceIdx-processedFaces,nfaces(level),level);
                    end
                    xVals=obj.faceXPts(faceIdx,:);
                    yVals=obj.faceYPts(faceIdx,:);
                    %sidelength=floor(obj.blocksize/pathskip)+1;
                    rectPts = [xVals', yVals'];
                    rectPts(1,:)=rectPts(1,:)-1;
                    rectPts(2,2)=rectPts(2,2)-1;
                    rectPts(4,1)=rectPts(4,1)-1;
                    numInvalidPts = getSum(obj.iImg_numInvalidPts, rectPts);
                    numPts = (rectPts(3,1)-rectPts(1,1))*(rectPts(3,2)-rectPts(1,2));
                    if (numInvalidPts == 0)
                        % fit and analyze error
                        if (USE_INTEGRAL_IMAGES)
                            if (FAST_FIT==1)
                                [coeffs, error] = obj.integralImageFit(rectPts);
                            end
                        else
                            if (FAST_FIT==1)
                                [coeffs, error] = FittingLibrary.fastFitPlane( ...
                                    reshape(obj.tan_theta_x(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                    reshape(obj.tan_theta_y(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                    reshape(obj.Z(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1));
                            else
                                [coeffs, error] = FittingLibrary.fitPlane( ...
                                    [reshape(obj.X(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                    reshape(obj.Y(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                    reshape(obj.Z(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1)]);
                                % [coeffs, error] = FittingLibrary.fitQuadric( ...
                                % [reshape(obj.X(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                % reshape(obj.Y(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1), ...
                                % reshape(obj.Z(yVals(1):yVals(3),xVals(1):xVals(3)),numPts,1)]);
                                % pause;
                            end
                        end
                        % average depth of points in the block
                        %avgZ = getSum(obj.iImg_z, rectPts)/numPts;
                        if (coeffs(3) < 0) % inward (away from camera) orientation
                            coeffs = -coeffs;
                        end
                        if (error < ERROR_THRESHOLD*level)
                            % draw the block in the image (green)
                            figure(1), hold on, plot(xVals, yVals,'g');
                            planeCoeffs(faceIdx,:) = coeffs;
                            areas(faceIdx,1) = numPts;
                            locations(faceIdx,:) = ...
                                [(rectPts(3,1)+rectPts(1,1)), ...
                                (rectPts(3,2)+rectPts(1,2))];
                            errors(faceIdx,1) = error;
                        else
                            planeCoeffs(faceIdx,:) = zeros(1,4);
                            areas(faceIdx,1) = numPts;
                            locations(faceIdx,:) = ...
                                [(rectPts(3,1)+rectPts(1,1)), ...
                                (rectPts(3,2)+rectPts(1,2))];
                            errors(faceIdx,1) = error;
                            if (level < MAXLEVEL)
                                % draw the block in the image (blue)
                                figure(1), hold on, plot(xVals, yVals,'b');
                                obj.subDivideFace(faceIdx);       % subdivide
                            else
                                % draw the block in the image (red)
                                hold on, figure(1), plot(xVals, yVals,'r');
                            end
                        end
                    else
                        planesCoeffs(faceIdx,:) = zeros(1,4);
                        areas(faceIdx,1) = numPts;
                        locations(faceIdx,:) = ...
                            [(rectPts(3,1)+rectPts(1,1)), ...
                            (rectPts(3,2)+rectPts(1,2))];
                        errors(faceIdx,1) = Inf;
                        if (numInvalidPts < 0.5*numPts && level < MAXLEVEL)
                            % draw the block in the image (blue)
                            hold on, figure(1), plot(xVals, yVals,'b');
                            obj.subDivideFace(faceIdx);        % subdivide
                        else
                            % draw the block in the image (red)
                            hold on, figure(1), plot(xVals, yVals,'r');
                        end
                    end
                    faceIdx = faceIdx + 1;
                end
                level=level+1;
                processedFaces=processedFaces+nfaces(level-1);
                nfaces(level)=size(obj.faceXPts,1)-sum(nfaces(1:level-1));
            end
            fprintf(1,'Processed %d faces\n',processedFaces);
            %[classifications,meanVals]=kmeans([planeCoeffs(:,1:2) planeCoeffs(:,4)],21);
        end
        
        function [coeffs, error] = integralImageFit(obj, rectPts)
            numPts = (rectPts(3,1)-rectPts(1,1))*(rectPts(3,2)-rectPts(1,2));
            sVals=zeros(4,4);
            sVals(2,1) = getSum(obj.iImg_tan_theta_xy, rectPts);
            sVals(3,1) = getSum(obj.iImg_tan_theta_x, rectPts);
            sVals(3,2) = getSum(obj.iImg_tan_theta_y, rectPts);
            sVals(4,1) = getSum(obj.iImg_ix, rectPts);
            sVals(4,2) = getSum(obj.iImg_iy, rectPts);
            sVals(4,3) = getSum(obj.iImg_iz, rectPts);
            sVals = sVals + sVals';
            sVals(1,1) = getSum(obj.iImg_tan_theta_x_sq, rectPts);
            sVals(2,2) = getSum(obj.iImg_tan_theta_y_sq, rectPts);
            sVals(3,3) = numPts;
            sVals(4,4) = getSum(obj.iImg_iz_sq, rectPts);
            %   sVals
            [evec,eval]=eig(sVals);
            coeffs = evec(:,1);  % vector associated with smallest eigenvalue
            normf = norm(coeffs(1:3));
            coeffs = coeffs'./normf;
            error = sqrt(eval(1,1))/(normf*numPts);
            
            % sVals2 = zeros(3,3);
            % sVals2(2,1) = sVals(2,1);
            % sVals2(3,1) = sVals(3,1);
            % sVals2(3,2) = sVals(3,2);
            % sVals2 = sVals2 + sVals2';
            % sVals2(1,1) = sVals(1,1);
            % sVals2(2,2) = sVals(2,2);
            % sVals2(3,3) = numPts;
            % sbVals(1,1) = sVals(4,1);
            % sbVals(2,1) = sVals(4,2);
            % sbVals(3,1) = sVals(4,3);
            % coeffs2 = inv(sVals2)*sbVals;
            % coeffs2 = [coeffs2; -1];
            % normf = norm(coeffs2(1:3));
            % coeffs2 = coeffs2'./normf
            % coeffs=coeffs2;
        end
    end
end