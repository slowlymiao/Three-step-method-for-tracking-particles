%-

function [centn]=S_centn(x,y,Image1,rad,t)

[BL,BR,BU,BD,~]=Cro_boundary([x,y],[t,t],[round(t*1.2),t],size(Image1,1),size(Image1,2),1,rad);
if BU<rad
    centn=[];
else
A=Image1(BU:BD,BL:BR);
oc=[x,y]-[BL,BU]+1;%原始质心
g=histeq(A,256);
level=graythresh(g);
g1=im2bw(g,level);
r_maxl_minl=2;%初始值
c=5; i=0;
while r_maxl_minl>1.2 && i<=3
    
    g2=imerode(g1,strel('disk',c-i));
    g2=imdilate(g2,strel('disk',c-i-1));
    s=regionprops(g2,'Centroid','MajorAxisLength','MinorAxisLength','EquivDiameter');
    B=bwboundaries(g2);
    if isempty(s)==0
        temp=struct2cell(s);
        cent1=cell2mat(temp(1,:));
        cent=reshape(cent1,2,length(cent1)/2);
        maxl=cell2mat(temp(2,:));
        minl=cell2mat(temp(3,:));
        d0=cell2mat(temp(4,:));
        clear temp cent1 s
      
        dist=sqrt((cent(1,:)-oc(1)).^2+(cent(2,:)-oc(2)).^2);
        dist(dist>rad+5)=1000;
        bizhi=dist./d0;
        [~,ind]=min(dist);
        if length(ind)>1
            ll=maxl./minl;
            ii=find(ll(ind)==min(ll(ind)));
            ind1=ind(ii);
            clear ind
            ind=ind1;
            clear ind1
        end
        centn=cent(:,ind)';
        centn(:,3)=d0(ind);
        r_maxl_minl=maxl(ind)/minl(ind);
        i=i+1;
    else
        r_maxl_minl=1.1;
        centn(1:2)=oc;
        centn(3)=15;
    end
    
end
centn(:,1:2)=centn(:,1:2)+[BL,BU]-1;
end


