function iVal = getSum(iImg, rectPts)

iVal = iImg(rectPts(1,2), rectPts(1,1)) + ...
    iImg(rectPts(3,2), rectPts(3,1)) - ...
    iImg(rectPts(2,2), rectPts(2,1)) - ...
    iImg(rectPts(4,2), rectPts(4,1));
    