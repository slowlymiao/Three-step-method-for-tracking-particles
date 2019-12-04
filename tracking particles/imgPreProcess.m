
function [Img1]=imgPreProcess(Img0,se)
Img1=imtophat(Img0,se);

    h=imhist(Img1);
    Pmax=percentile2i(h,0.98);
    Img1=imadjust(Img1,[0, Pmax],[0, 1]);

end