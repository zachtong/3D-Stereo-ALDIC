% This function is used to change 32bits depth image to 24bits
function ChangeImageIf32bits(ImageInfo)

if imfinfo(ImageInfo{1,end}).BitDepth == 32
    for i = 1:size(ImageInfo,2)
        img_rgb = imread(ImageInfo{1,i});
        imwrite(img_rgb(:,:,1:3),[ImageInfo{2,i},'\',ImageInfo{1,i}]);
    end
end

end

