classdef LoopyBBP < handle
    % 
    % Code follows closely C++ code provided by Pedro Felzenszwalb
    % from his article:
    %
    % Efficient Belief Propagation, 2006
    % Felzenszwalb, Huttenlocher
    %
    % Andrew Willis
    % UNC Charlotte
    % January 15, 2017
    %
    % Copyright (C) 2006 Pedro Felzenszwalb
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program; if not, write to the Free Software
    % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
    %
    properties (SetAccess = protected)
        ITER
        LEVELS
        DISC_K
        DATA_K
        LAMBDA
        INF
        VALUES
        DT_L1
        DT_L2
        debug
    end
    
    methods (Static)
        function xsq = square(x)
            xsq=x*x;
        end
    end
    
    methods
        function obj = LoopyBBP(varargin)
            obj.ITER = 0;            % number of BP iterations at each scale
            obj.LEVELS = 0;          % number of scales
            obj.DISC_K = 0;          % truncation of discontinuity cost
            obj.DATA_K = 0;          % truncation of data cost
            obj.LAMBDA = 0;          % weighting of data cost
            obj.INF = 0;             % large cost
            obj.VALUES = 0;          % number of possible graylevel values
        end
        
        function setParameters( obj, ITER, LEVELS, VALUES, LAMBDA)
            obj.ITER = ITER;         % number of BP iterations at each scale
            obj.LEVELS = LEVELS;     % number of scales
            obj.DISC_K = 200;        % truncation of discontinuity cost
            obj.DATA_K = 10000;      % truncation of data cost
            obj.LAMBDA = LAMBDA;       % weighting of data cost
            obj.INF = 1e10;          % large cost
            obj.VALUES = VALUES;     % number of possible graylevel values
        end
        
        function Iout = loopybbp(obj, data)
            % data pyramid
            for i=2:obj.LEVELS
                old_width = data{i-1}.width;
                old_height = data{i-1}.height;
                new_width = ceil(old_width/2);
                new_height = ceil(old_height/2);
                assert(new_width >= 1);
                assert(new_height >= 1);
                data{i}.values = zeros(new_height, new_width, obj.VALUES);
                data{i}.width = new_width;
                data{i}.height = new_height;
                % sum (column,row) neighbors to create 1/2 resolution image
                for y=1:old_height
                    yhalf=ceil(y/2);
                    for x=1:old_width
                        xhalf=ceil(x/2);
                        for value=1:obj.VALUES
                            data{i}.values(yhalf,xhalf,value) = ...
                                data{i}.values(yhalf,xhalf,value) + ...
                                data{i-1}.values(y,x,value);
                        end
                    end
                end
            end
            
            % run bp from coarse to fine
            for i=obj.LEVELS:-1:1
                width = data{i}.width;
                height = data{i}.height;
                % allocate & init memory for messages
                if (i == obj.LEVELS)
                    % in the coarsest level messages are initialized to zero
                    u{i}.values = zeros(height, width, obj.VALUES);
                    d{i}.values = zeros(height, width, obj.VALUES);
                    l{i}.values = zeros(height, width, obj.VALUES);
                    r{i}.values = zeros(height, width, obj.VALUES);
                else
                    % initialize messages from values of previous level
                    u{i}.values = zeros(height, width, obj.VALUES);
                    d{i}.values = zeros(height, width, obj.VALUES);
                    l{i}.values = zeros(height, width, obj.VALUES);
                    r{i}.values = zeros(height, width, obj.VALUES);
                    
                    for y=1:height
                        yhalf=ceil(y/2);
                        for x=1:width
                            xhalf=ceil(x/2);
                            for value=1:obj.VALUES
                                u{i}.values(y,x,value) = u{i+1}.values(yhalf,xhalf,value);
                                d{i}.values(y,x,value) = d{i+1}.values(yhalf,xhalf,value);
                                l{i}.values(y,x,value) = l{i+1}.values(yhalf,xhalf,value);
                                r{i}.values(y,x,value) = r{i+1}.values(yhalf,xhalf,value);
                            end
                        end
                    end
                    % delete old messages and data
                    u{i+1}=[];
                    d{i+1}=[];
                    l{i+1}=[];
                    r{i+1}=[];
                    data{i+1}.values=[];
                end
                
                % BP
                [u{i}, d{i}, l{i}, r{i}] = obj.bp_cb( u{i}, d{i}, l{i}, r{i}, data{i}, obj.ITER);
            end
            Iout = obj.output( u{1}, d{1}, l{1}, r{1}, data{1});
            clear u;
            clear d;
            clear l;
            clear r;
            clear data;
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
                            r.values(y,x+1,:), data.values(y,x,:));
                        d.values(y,x,:) = obj.msg(d.values(y+1,x,:), l.values(y,x-1,:), ...
                            r.values(y,x+1,:), data.values(y,x,:));
                        r.values(y,x,:) = obj.msg(u.values(y-1,x,:), d.values(y+1,x,:), ...
                            r.values(y,x+1,:), data.values(y,x,:));
                        l.values(y,x,:) = obj.msg(u.values(y-1,x,:), d.values(y+1,x,:), ...
                            l.values(y,x-1,:), data.values(y,x,:));
                        %if (x==round(width/2) && y==round(height/2))
                        %    figure(2), plot(reshape(u.values(y,x,:),1,obj.VALUES))
                        %    pause(0.5);
                        %end
                    end
                end
            end
        end
        
        % compute message
        function dst = msg(obj, s1, s2, s3, s4) % s1,s2,s3,s4,dst have VALUES elements
            dst = zeros(1,obj.VALUES);
            % aggregate/sum messages to node and find min label value
            minimum = obj.INF;
            for value = 1:obj.VALUES
                dst(value) = s1(value) + s2(value) + s3(value) + s4(value);
                if (dst(value) < minimum)
                    minimum = dst(value);
                end
            end
            
            % dt
            if (obj.DT_L2)
                tmp = obj.dt_l2(dst);
            elseif (obj.DT_L1)
                tmp = obj.dt_l1(dst);
            end
            
            % truncate costs to min(dst)+DISC_K and store in destination vector
            minimum = minimum + obj.DISC_K;
            for value=1:obj.VALUES
                dst(value) = min(tmp(value), minimum);
            end
            % normalize: computes mean of function range, shift it to 0
            val = 0;
            for value = 1:obj.VALUES
                val = val + dst(value);
            end
            val = val / obj.VALUES;            
            for value = 1:obj.VALUES
                dst(value) = dst(value) - val;
            end
            clear tmp;
        end
        
        % dt of 1d function having linear (absolute difference) error behaviors
        % use for L1 objective functions
        % Distance Transform between label values (x-axis) and costs
        % (y-axis) for a site in the MRF for a linear (L1) cost function
        function f = dt_l1(obj, f)
            %d = zeros(1,obj.VALUES);
            for q=2:obj.VALUES
                prev = f(q-1) + 1;
                if (prev < f(q))
                    f(q) = prev;
                end
            end
            for q=(obj.VALUES-1):-1:1
                prev = f(q+1) + 1;
                if (prev < f(q))
                    f(q) = prev;
                end
            end
        end
        
        % dt of 1d function having quadric (parabolic) error behaviors
        % use for L2 objective functions
        % Distance Transform between label values (x-axis) and costs
        % (y-axis) for a site in the MRF for a quadratic (L2) cost function
        function d = dt_l2( obj, f)
            d = zeros(1,obj.VALUES);
            v = zeros(1,obj.VALUES);
            z = zeros(1,obj.VALUES+1);
            
            k = 1;
            v(1) = 0;
            z(1) = -obj.INF;
            z(2) = +obj.INF;
            
            for q=1:(obj.VALUES-1)
                s  = ((f(q+1)+obj.square(q))-(f(v(k)+1)+obj.square(v(k)))) / (2*(q-v(k)));
                while (s <= z(k))
                    k = k-1;
                    s  = ((f(q+1)+obj.square(q))-(f(v(k)+1)+obj.square(v(k)))) / (2*(q-v(k)));
                end
                k=k+1;
                v(k) = q;
                z(k) = s;
                z(k+1) = obj.INF;
            end
            k = 1;
            for q=0:(obj.VALUES-1)
                while (z(k+1) < q)
                    k=k+1;
                end
                d(q+1) = obj.square(q-v(k)) + f(v(k)+1);
            end
            %clear v;
            %clear z;
        end        
        
        % generate output from current messages
        function Iout = output(obj, u, d, l, r, data)
            width = data.width;
            height = data.height;
            Iout = zeros(height,width);
            
            for y = 2:(height-1)
                for x = 2:(width-1)
                    % keep track of best value for current pixel
                    best = 0;
                    best_val = obj.INF;
                    for value = 1:obj.VALUES
                        val = u.values(y-1,x,value) + ...
                            d.values(y+1,x,value) + ...
                            l.values(y,x-1,value) + ...
                            r.values(y,x+1,value) + ...
                            data.values(y,x,value);
                        if (val < best_val)
                            best_val = val;
                            best = value-1;
                        end
                    end
                    Iout(y,x)=best;
                end
            end
        end        
    end
end