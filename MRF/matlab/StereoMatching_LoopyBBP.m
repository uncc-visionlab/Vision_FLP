classdef StereoMatching_LoopyBBP < LoopyBBP
    properties (SetAccess = protected)
        SIGMA
    end
    methods(Static)
        function demo()
            imgLeft=imread('tsukuba1.pgm');
            imgRight=imread('tsukuba2.pgm');
            stereo = StereoMatching_LoopyBBP();
            Iout = stereo.stereo_ms(double(imgLeft),double(imgRight));
            imshow(Iout,[0,16])
        end
    end
    
    methods
        function obj = StereoMatching_LoopyBBP(varargin)
            obj@LoopyBBP();          % call superclass constructor
        end
                
        % multiscale belief propagation for stereo matching
        function Iout = stereo_ms(obj, Ileft, Iright)
            obj.ITER = 5;            % number of BP iterations at each scale
            obj.LEVELS = 5;          % number of scales
            obj.DISC_K = 1.7;        % truncation of discontinuity cost
            obj.DATA_K = 15;         % truncation of data cost
            obj.LAMBDA = 0.07;       % weighting of data cost
            obj.INF = 1e20;          % large cost
            obj.VALUES = 16;         % number of possible graylevel values
            
            obj.SIGMA = 0.7;         % amount to smooth the input images
            obj.DT_L1 = true;
            
            % data costs
            data{1} = obj.comp_data_stereo(Ileft,Iright);
            Iout = obj.loopybbp(data);
        end
        
        % computation of data costs (stereo matching)
        function data = comp_data_stereo( obj, img1, img2)
            [height, width] = size(img1);
            data.values = zeros(width, height, obj.VALUES);
            data.width = width;
            data.height = height;
            
            if (obj.SIGMA >= 0.1)
                sm1 = imgaussfilt(img1, obj.SIGMA);
                sm2 = imgaussfilt(img2, obj.SIGMA);
            else
                sm1 = img1;
                sm2 = img2;
            end
            
            for y = 1:height
                for x = obj.VALUES:width
                    for value = 1:obj.VALUES
                        % abs instensity diff for disparity d=(value-1)
                        val = abs(sm1(y,x)-sm2(y,x-(value-1)));
                        data.values(y,x,value) = obj.LAMBDA * min(val, obj.DATA_K);
                    end
                end
            end
        end
    end
end