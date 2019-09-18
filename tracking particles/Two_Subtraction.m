%two subtraction%%%%%%%%%%%%%%%%
function [cent,fc]=Two_Subtraction(Image1,Image2,rad,pts)

    I0=Image2-Image1;
    for i=1:2
        if i==1
            mark=I0<0;
        else
            mark=I0>0;
        end
        I1=abs(I0.*mark);
        se=fspecial('gaussian',5,1);
        I2=filter2(se,I1,'same');
        THRE=graythresh(I2);
        if THRE<(mean(Image1(:))*0.3+THRE)/2
            THRE=(mean(Image1(:))*0.3+THRE)/2;
        end
        I3=im2bw(I2,THRE);
        I4=imerode(I3,strel('disk',4));
        I5=imdilate(I4,strel('disk',4));

        [L,N]=bwlabel(I5,4);
        bar=zeros(N,3);
        for n=1:N
            [r,c]=find(L==n);
            marker=(L==n);
            bar(n,2)=round(mean(r));
            bar(n,1)=round(mean(c));
            bar(n,3)=ceil(sqrt(sum(marker(:))));
        end
        bar(bar(:,3)<pts,:)=[];%
        if i==1
            cent1=bar;
            cent1(:,4)=1;
            cent1(:,5)=0;
            clear I1 I2 I3 L N bar
        else
            cent2=bar;
            cent2(:,4)=2;
            cent2(:,5)=0;
        end
    end

m1=cent1';
m2=cent2';
cent1=[cent1,zeros(size(cent1,1),2)];
cent2=[cent2,zeros(size(cent2,1),2)];
[rp1,rp2]=Data_Association(m1,m2,2*rad);

if isempty(rp2)==0 && sum(rp2)~=0
    rtemp=tabulate(rp1);
    ind=find(rtemp(:,2)>1);
    ind2=[];
    if isempty(ind)==0
        for i=1:length(ind)
            ind1=find(rp1==rtemp(ind(i),1));
            ind2=[ind2,ind1(2)];
        end
        rp1(ind2)=[];
        rp2(ind2)=[];
    end
    for ii=1:length(rp1)
    if cent1(rp1(ii),2)>cent2(rp2(ii),2) 
        cent1(rp1(ii),1:2)=round(0.75*cent1(rp1(ii),1:2)+0.25*cent2(rp2(ii),1:2));
        cent1(rp1(ii),5)=1;
        cent1(rp1(ii),6:7)=round(0.25*cent1(rp1(ii),1:2)+0.75*cent2(rp2(ii),1:2));
        cent2(rp2(ii),:)=0;
    else
        cent2(rp2(ii),1:2)=round(0.75*cent2(rp2(ii),1:2)+0.25*cent1(rp1(ii),1:2));
        cent2(rp2(ii),5)=1;
         cent2(rp2(ii),6:7)=round(0.25*cent2(rp2(ii),1:2)+0.75*cent1(rp1(ii),1:2));
         cent1(rp1(ii),:)=0;
    end
    end
    cent=[cent1;cent2];
    cent(cent(:,1)==0,:)=[];
else
    cent=[cent1;cent2];
end

aa=find(cent(:,5)==0 & cent(:,3)<1.5*pts);
if isempty(aa)==0
cent(aa,:)=[];
end
bb=find(cent(:,4)==2 & cent(:,5)==0);
if isempty(bb)==0
cent(bb,:)=[];
end
fc=cent(:,4:5);
cent(:,4:7)=[];





