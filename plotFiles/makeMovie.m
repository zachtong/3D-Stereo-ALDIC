files = dir('*.tiff'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

v = VideoWriter('video_oht_cfrp_RawDICImages.mp4');
v.FrameRate = 5;
open(v);
figure,
for tempk = [ 1 : length(im) ]
    clf; 
    %imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow(imread(im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    %waitbar(tempk/length(files));
    
end
close(v);