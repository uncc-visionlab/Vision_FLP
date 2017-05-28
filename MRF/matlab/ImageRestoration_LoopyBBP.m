classdef ImageRestoration_LoopyBBP < LoopyBBP
    properties (SetAccess = protected)
    end
    methods(Static)
        function demo()
            I=imread('penguin.pgm');
            J = imnoise(I,'gaussian',0,.01);
            subplot(2,3,1), imshow(I,[0 255]);
            subplot(2,3,2), imshow(J,[0 255]);
            %J = imnoise(I,'localvar',V)
            %J = imnoise(I,'localvar',image_intensity,var)
            %J = imnoise(I,'poisson')
            %J = imnoise(I,'salt & pepper',d)
            %J = imnoise(I,'speckle',v)
            imrestore = ImageRestoration_LoopyBBP();
            imrestore.setParameters(5, 5, 256, 0.05);
            Iout = imrestore.restore_ms(double(J));
            subplot(2,3,3), imshow(Iout,[0 255]);
            
            [rows,cols]=size(I);
            compVal=8;
            I2=round(I/compVal);
            I2=imresize(I2,0.25);
            [rows,cols]=size(I2);
            J2 = imnoise(I2,'gaussian',0,.0001);
            %levels = log2(max(rows,cols))-3;
            %values = 256/compVal;
            subplot(2,3,4), imshow(I2,[0 256/compVal]);
            subplot(2,3,5), imshow(J2,[0 256/compVal]);
            imrestore.setParameters(3, 2, 32, 0.5);
            imrestore.DATA_K = 100;
            imrestore.DISC_K = 10;
            Iout = imrestore.restore_ms(double(J2));
            subplot(2,3,6), imshow(Iout,[0 255/compVal]);
        end
    end
    
    methods
        function obj = ImageRestoration_LoopyBBP(varargin)
            obj@LoopyBBP();          % call superclass constructor
            obj.ITER = 5;            % number of BP iterations at each scale
            obj.LEVELS = 5;          % number of scales
            obj.DISC_K = 200;        % truncation of discontinuity cost
            obj.DATA_K = 10000;      % truncation of data cost
            obj.LAMBDA = 0.05;       % weighting of data cost
            obj.INF = 1e10;          % large cost
            obj.VALUES = 256;        % number of possible graylevel values
        end
        
        % multiscale belief propagation for image restoration
        function Iout = restore_ms(obj, Iin)            
            obj.DT_L2 = true;            
            % data costs
            data{1} = obj.comp_data_restore(Iin);
            Iout = obj.loopybbp(data);
        end
        
        % computation of data costs (image restoration)
        function data = comp_data_restore( obj, img)
            [height, width] = size(img);
            data.values = zeros( height, width, obj.VALUES);
            data.width = width;
            data.height = height;
            for y=1:height
                for x=1:width
                    for value=1:obj.VALUES
                        val = obj.square( img(y,x) - value - 1); % sq intensity diff
                        data.values(y,x,value) = obj.LAMBDA * min(val, obj.DATA_K);
                    end
                end
            end
        end
    end
end