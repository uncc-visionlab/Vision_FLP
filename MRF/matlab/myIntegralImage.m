function [iImg,num_inValidImg] = myIntegralImage(I, computeInvalidImg)

if(exist('computeInvalidImg','var')==0)
    computeInvalidImg=0;
end
[height,width]=size(I);
iImg=zeros(height,width);
for y=1:height
    for x=1:width,
        if (~isfinite(I(y,x)) || I(y,x)==0)
            integral_pixel = 0.0;
        else
            integral_pixel = I(y,x);
        end
        left_value = 0; top_left_value = 0; top_value = 0;
        if (x > 1)
            left_value = iImg(y,x-1);
            if (y > 1)
                top_left_value = iImg(y-1,x-1);
            end
        end
        if (y > 1)
            top_value = iImg(y-1,x);
        end
        iImg(y,x) = integral_pixel + left_value + top_value - top_left_value;
    end
end

if (computeInvalidImg)
    num_inValidImg=zeros(height,width);
    for y=1:height
        for x=1:width,
            if (~isfinite(I(y,x)) || I(y,x)==0)
                invalid_point = 1;
            else
                invalid_point = 0;
            end
            left_invalid_points = 0; top_left_invalid_points = 0; top_invalid_points = 0;
            if (x > 1)
                left_invalid_points = num_inValidImg(y,x-1);
                if (y > 1)
                    top_left_invalid_points = num_inValidImg(y-1,x-1);
                end
            end
            if (y > 1)
                top_invalid_points = num_inValidImg(y-1,x);
            end
            num_inValidImg(y,x) = invalid_point + left_invalid_points + top_invalid_points - top_left_invalid_points;
        end
    end
end

        