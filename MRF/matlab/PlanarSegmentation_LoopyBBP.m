classdef PlanarSegmentation_LoopyBBP < LoopyBBP
    properties (SetAccess = protected)
        planeCoeffs
        %pts3D
        X
        Y
        Z
        horz_bonds
        vert_bonds
        BETA
        region_of_interest
        plane_cost_matrix
        USE_BONDS
        BOND_WEIGHT
    end
    methods(Static)
        function demo()
            clear all;
            close all;
            clc;
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
            % read an image from the disk / sensor
            fid=fopen('image0.npy');
            val=fread(fid,rgbdImg.numPixels,'single');
            I=reshape(val,width,height)';
            
            % using NYU dataset
            %nyu_index = 165; 
            %nyu_index = 845; 
            %nyu = matfile('nyu_depth_v2_labeled.mat');
            %I = nyu.rawDepths(:,:,nyu_index);
            
            rgbdImg.setImage(I);
            
            figure(1), hold on, imshow(I,[0 rgbdImg.MAXRANGE]);
            colormap('jet');
            
            blocksize=50;
            ROI=[50, 50; 590, 450];
            rgbdImg.setRegionOfInterest(ROI);
            rgbdImg.imposeGridDecomposition(blocksize);
            
            MAXLEVEL=1;
            ERROR_THRESHOLD=1.0e-3;
            tic
            % detect planes
            [planeCoeffs, areas1, locations, errors] = rgbdImg.getPlanes(MAXLEVEL, ERROR_THRESHOLD);
            planeCoeffs=planeCoeffs(errors<(ERROR_THRESHOLD/4),:);
            planeCoeffs2=planeCoeffs(any(planeCoeffs,2),:);
            errors2=errors(errors <(ERROR_THRESHOLD/4));
            errors2=errors2(any(planeCoeffs,2));
            assert(size(planeCoeffs2,1) == size(errors2,1));
            
            fprintf('%i planes detected\n', size(planeCoeffs2,1));
            % cluster planes
            BETA = 0.4; % weight d1-d2 relative to angle difference cos(theta)
            LOC_THRESH=0.2; % threshold magnitude of difference for merge
            % loops over the 4 plane parameters [a b c d]
            sortedPlanes = planeCoeffs2;
            for param=4:-1:1
                %for param=1:4
                % sort the planes by the parameter having index param
                [~, order] = sort(sortedPlanes(:,param));
                sortedPlanes = sortedPlanes(order,:); % planes sorted by index param
                sortedErrors = errors2(order);
                planeIdx=1;                           % start at plane with index 1                
                while (planeIdx < size(sortedPlanes,1))
                    cplaneIdx=planeIdx+1; % compare with planes at index 2->end of list
                    errorVal=0;
                    while (cplaneIdx <= size(sortedPlanes,1) && errorVal < LOC_THRESH)
                        errorVal = 1-abs(sortedPlanes(planeIdx,1:3)*sortedPlanes(cplaneIdx,1:3)') + ...
                            BETA*abs(sortedPlanes(planeIdx,4)-sortedPlanes(cplaneIdx,4));
                        if (errorVal < LOC_THRESH)
                            if(sortedErrors(planeIdx) < sortedErrors(cplaneIdx))
                               sortedPlanes(cplaneIdx,:)=[];
                            else
                               sortedPlanes(planeIdx,:) = sortedPlanes(cplaneIdx,:);
                               sortedErrors(planeIdx) = sortedErrors(cplaneIdx);
                               sortedPlanes(cplaneIdx,:)=[];                               
                            end
                        else
                            cplaneIdx = cplaneIdx + 1;
                        end
                    end
                    planeIdx = planeIdx + 1;
                end
            end
            numPlanes = size(sortedPlanes,1);
            fprintf('%i planes after clustering\n', numPlanes);

            planeList=cell(numPlanes,1);
            for planeIdx=1:numPlanes
                planeList{planeIdx} = Plane3d(sortedPlanes(planeIdx,1), ...
                    sortedPlanes(planeIdx,2), sortedPlanes(planeIdx,3), ...
                    sortedPlanes(planeIdx,4));
            end
            
            for planeA_Idx=1:numPlanes
                planeA = planeList{planeA_Idx};
                for planeB_Idx=planeA_Idx+1:numPlanes
                    planeB = planeList{planeB_Idx};
                    line3 = planeA.intersect(planeB);
                    lineseg2 = rgbdImg.getVisibleLineSegment2d(line3);
                    ptA = lineseg2.getPoint(lineseg2.l_start);
                    ptB = lineseg2.getPoint(lineseg2.l_end);
                    ptA.x = ptA.x + 1;
                    ptA.y = ptA.y + 1;
                    ptB.x = ptB.x + 1;
                    ptB.y = ptB.y + 1;
                    %dihedral_angle = (180/pi)*acos(planeA.x*planeB.x + ...
                    %    planeA.y*planeB.y + planeA.z*planeB.z)
                    if (ptA.x >= 1 && ptA.x <= 640 && ...
                        ptA.y >= 1 && ptA.y <= 480 && ...
                        ptB.x >= 1 && ptB.x <= 640 && ...
                        ptB.y >= 1 && ptB.y <= 480)
                        figure(1), plot([ptA.x ptB.x],[ptA.y ptB.y],'w');
                    end
                end
            end

            if (true)
                fileID = fopen('dimage.dat','w');
                fprintf(fileID,'%f %f\n', size(K));
                fprintf(fileID,'%f ', K);
                fprintf(fileID,'\n%f %f\n', size(sortedPlanes));
                fprintf(fileID,'%f ', sortedPlanes');
                fprintf(fileID,'\n%f %f\n', size(I));
                fprintf(fileID,'%f ', I');
                %formatSpec = '%s %d %2.1f %s\n';
                %[nrows,ncols] = size(C);
                %for row = 1:nrows
                %    fprintf(fileID,formatSpec,);
                %end
                fclose(fileID);
                %save dimage.txt -ascii K R P I sortedPlanes -ascii
            end
            
            % run multi-scale BBP to segment points to the N plane classes
            % - data cost = weighted sum of orientation and location errors
            % - truncated data cost = max depth error+20%
            % - smoothness cost = Ising model
            % - truncated discontinuity cost = depth discontinuity cost            
            imgSeg = PlanarSegmentation_LoopyBBP();

            VALUES = size(sortedPlanes,1);
            planedist = zeros(VALUES,VALUES);
            for i=1:VALUES
                for j=(i+1):VALUES
                    errorVal = 1 - (sortedPlanes(i,1:3)*sortedPlanes(j,1:3)') + ...
                        BETA*abs(sortedPlanes(i,4)-sortedPlanes(j,4));
                    planedist(i,j)=errorVal;
                end
            end
            planedist = planedist + planedist';
            
            IMAGE_TO_3PLANES=1;
            if (IMAGE_TO_3PLANES==1)
                BETA = 0.4;
                LOC_THRESH=0.4;
                ITER = 5;           % number of BP iterations at each scale
                LEVELS = 1;         % number of scales
                LAMBDA = .5;       % weighting of data cost
                imgSeg.setParameters(ITER, LEVELS, VALUES, LAMBDA);
                imgSeg.DISC_K = 2;         % truncation of discontinuity cost
                imgSeg.DATA_K = 0.5;         % truncation of data cost
                imgSeg.INF = 1e10;          % large cost
                imgSeg.USE_BONDS = false;
                imgSeg.BOND_WEIGHT = 1/50; % weighting of the bond cost
            end
            imgSeg.setRegionOfInterest(ROI);
            imgSeg.setData(rgbdImg.X, rgbdImg.Y, rgbdImg.Z, sortedPlanes);
            imgSeg.plane_cost_matrix = planedist;
            seg = imgSeg.plane_ms(rgbdImg);
            toc
            figure(2), imshow(seg,[]);
            
            %tlc = [rgbdImg.faceXPts(1,1), ...
            %    rgbdImg.faceYPts(1,1)];
            %brc = [rgbdImg.faceXPts(rgbdImg.numfaces,3), ...
            %    rgbdImg.faceYPts(rgbdImg.numfaces,3)];
            %nheight=brc(2)-tlc(2);
            %nwidth=brc(1)-tlc(1);
            
            % N plane classes, create 1 background class
            seg = seg + 1;
            [nheight,nwidth]=size(seg);
            imgSeg2 = zeros(nheight,nwidth);
            for y=1:nheight
                for x=1:nwidth
                    myplane = sortedPlanes(seg(y,x),:);
                    error = abs([imgSeg.X(y,x), imgSeg.Y(y,x), imgSeg.Z(y,x), 1]*myplane');
                    if (error > 0.15)
                        imgSeg2(y,x) = 0;
                    else
                        imgSeg2(y,x)=seg(y,x);
                    end
                end
            end
            %figure(3), imshow(imgSeg2,[]);
            figure(3), hold off, imshow(imgSeg.Z,[]);
            rgb=label2rgb(imgSeg2);
            figure(3), hold on;
            figure(3), h=imshow(rgb);
            set(h,'AlphaData',0.1);
            aaa=1;
        end
        
        function demo2()
            img=zeros(10,10);
            [rows,cols]=size(img);
            imgSeg = PlanarSegmentation_LoopyBBP();
            neighbors = imgSeg.getManhattanNeighbors(1,5,3,rows,cols);
            for i=1:length(neighbors)
                y = ceil(neighbors(i)/cols);
                x = neighbors(i)-(y-1)*cols;
                img(y,x)=1;
            end
            img
            aaa=1;
        end
    end
    
    methods
        function obj = PlanarSegmentation_LoopyBBP(varargin)
            obj@LoopyBBP();          % call superclass constructor
        end
        
        
        % selects a rectangular subset of the image to process
        % Specified as the 2x2 matrix [[x1,y1]; [x2,y2]] where
        % [x1,y1] = top-left-corner
        % [x2,y2] = bottom-right-corner
        function setRegionOfInterest(obj, region_of_interest)
            obj.region_of_interest = region_of_interest;
        end
        
        function setData(obj, X, Y, Z, planeCoeffs, subsample)
            obj.planeCoeffs = planeCoeffs;
            [height,width]=size(X);
            %margin=35;
            if (prod(size(obj.region_of_interest))~=4)
                aratio=obj.height/obj.width;  % aspect ratio
                margin=round(min(1*[aratio*blocksize,blocksize],[aratio*60,60]));
                obj.region_of_interest=[margin(1) margin(2); ...
                    obj.height-margin(1) obj.width-margin(2)];
            end
            if (exist('subsample')==0)
                subsample=3;
            end
            
            ysubsamp=obj.region_of_interest(1,2):subsample:obj.region_of_interest(2,2);
            xsubsamp=obj.region_of_interest(1,1):subsample:obj.region_of_interest(2,1);
            %xsubsamp=margin:skip:(width-margin);
            %ysubsamp=margin:skip:(height-margin);
            
            % set NaNs to zero
%             Znonans = Z;
%             Znonans(isnan(Znonans)) = 0;
            
            % delta X in horizontal direction
            obj.horz_bonds(:,:,1) = [(X(ysubsamp, xsubsamp(1:end-1) + subsample) - ...
                X(ysubsamp, xsubsamp(1:end-1))), zeros(size(ysubsamp, 2), 1)];
            % delta Z in horizontal direction
            obj.horz_bonds(:,:,2) = [(Z(ysubsamp, xsubsamp(1:end-1) + subsample) - ...
                Z(ysubsamp, xsubsamp(1:end-1))), zeros(size(ysubsamp, 2), 1)];
            % delta Y in vertical direction
            obj.vert_bonds(:,:,1) = [(Y(ysubsamp(1:end-1) + subsample, xsubsamp) - ...
                Y(ysubsamp(1:end-1), xsubsamp)); zeros(1, size(xsubsamp, 2))];
            % delta Z in vertical direction
            obj.vert_bonds(:,:,2) = [(Z(ysubsamp(1:end-1) + subsample, xsubsamp) - ...
                Z(ysubsamp(1:end-1), xsubsamp)); zeros(1, size(xsubsamp, 2))];
            
            obj.X = X(ysubsamp,xsubsamp);
            obj.Y = Y(ysubsamp,xsubsamp);
            obj.Z = Z(ysubsamp,xsubsamp);
            %obj.pts3D = cat(3, X, Y, Z);
        end
        
        function neighbors = getManhattanNeighbors(obj, yv, xv, dist, rows, cols)
            neighbors = [];
            level=0;
            for y=yv-dist:yv+dist
                for x=(xv-level):(xv+level)
                    if (y >= 1 && y <= rows && x >= 1 && x <= cols)
                        offset = (y-1)*cols+x;
                        neighbors = [neighbors, offset];
                    end
                end
                if (y < yv)
                    level=level+1;
                else
                    level=level-1;
                end
            end
        end
        
        % multiscale belief propagation for planar segmentation
        function Iout = plane_ms(obj, Irgbd)
            %obj.ITER = 3;           % number of BP iterations at each scale
            %obj.LEVELS = 1;         % number of scales
            %obj.DISC_K = 10;         % truncation of discontinuity cost
            %obj.DATA_K = 0.5;         % truncation of data cost
            %obj.LAMBDA = 100;       % weighting of data cost
            %obj.INF = 20000;          % large cost
            %obj.VALUES = size(obj.planeCoeffs,1);         % number of possible plane values
            
            obj.DT_L1 = true;
            
            % data costs
            data{1} = obj.comp_data_planeseg(Irgbd);
            Iout = obj.loopybbp(data);
            aa=1;
        end
        
        % computation of data costs (plane segmentation)
        function [data] = comp_data_planeseg( obj, img)
            [height, width] = size(obj.X);
            data.values = zeros( height, width, obj.VALUES);
            data.width = width;
            data.height = height;
            
            for y=1:height
                for x=1:width
                    for value=1:obj.VALUES
                        val = abs( [obj.X(y,x), obj.Y(y,x), obj.Z(y,x), 1]*obj.planeCoeffs(value,:)'); % sq intensity diff
                        data.values(y,x,value) = obj.LAMBDA * min(val, obj.DATA_K);
                    end
                end
            end
        end
        
        % belief propagation using checkerboard update scheme
        function [u, d, l, r] = bp_cb( obj, u, d, l, r, data, iter)
            width = data.width;
            height = data.height;
            for t=0:(obj.ITER-1)
                fprintf('iter %d of %d\n', t+1, iter);
                for y=2:(height-1)
                    % alternate message passing to even/odd indexed columns
                    for x=mod(y+t,2)+2:2:(width-1)
                        u.values(y,x,:) = obj.msg(u.values(y-1,x,:), l.values(y,x-1,:), ...
                            r.values(y,x+1,:), data.values(y,x,:), obj.vert_bonds(y,x,:), 2);
                        d.values(y,x,:) = obj.msg(d.values(y+1,x,:), l.values(y,x-1,:), ...
                            r.values(y,x+1,:), data.values(y,x,:), -obj.vert_bonds(y-1,x,:), 2);
                        r.values(y,x,:) = obj.msg(u.values(y-1,x,:), d.values(y+1,x,:), ...
                            r.values(y,x+1,:), data.values(y,x,:), -obj.horz_bonds(y,x-1,:), 1);
                        l.values(y,x,:) = obj.msg(u.values(y-1,x,:), d.values(y+1,x,:), ...
                            l.values(y,x-1,:), data.values(y,x,:), obj.horz_bonds(y,x,:), 1);
                        %if (x==round(width/2) && y==round(height/2))
                        %    figure(2), plot(reshape(u.values(y,x,:),1,obj.VALUES))
                        %    pause(0.5);
                        %end
                    end
                end
            end
        end
        
        % overrides LoopyBBP.msg()
        function outgoingmsg = msg(obj, neighbor1, neighbor2, neighbor3, data, bonds, pass_dir)
            num_labels = obj.VALUES;
            outgoingmsg = zeros(1, num_labels);
            
            if obj.USE_BONDS == true
                dxy = bonds(1);
                dz = bonds(2);
                bond_weight = obj.BOND_WEIGHT;
            end
            
            for i = 1:num_labels
                min_val = obj.INF;
                
                for j = 1:num_labels
                    planej = obj.planeCoeffs(j,:);
                    
                    p = 0;
                    p = p + data(j);
                    p = p + neighbor1(j) + neighbor2(j) + neighbor3(j);
                    p = p + obj.smoothnessCost(i, j);
                    
                    if obj.USE_BONDS == true
                    
                        if pass_dir == 1 % horizontal
                            bond_cost = bond_weight*abs(dz/dxy + planej(1)/planej(3));
                        else % vertical
                            bond_cost = bond_weight*abs(dz/dxy + planej(2)/planej(3));
                        end
                        if isnan(bond_cost) || isinf(bond_cost)
                            bond_cost = 0;
                        end

                        bond_truc = 5*bond_weight;
                        p = p + min(bond_cost, bond_truc);
                        
                    end
                    
                    min_val = min(min_val, p);
                end
                
                outgoingmsg(i) = min_val;
                
            end
            
            % truncate costs to min(outgoingmsg) + DISC_K
            for label = 1:num_labels
                outgoingmsg(label) = min(outgoingmsg(label), min(outgoingmsg) + obj.DISC_K);
            end
            
            % normalization
            outgoingmsg = outgoingmsg - mean(outgoingmsg);

            % alternative normalization
%             outgoingmsg = outgoingmsg - log(sum(exp(outgoingmsg)))
            
            % alternative normalization
%             outgoingmsg = outgoingmsg - min(outgoingmsg)
        end
        
        function cost = smoothnessCost(obj, i, j)
            cost = obj.plane_cost_matrix(i, j);
        end
        
    end
end
