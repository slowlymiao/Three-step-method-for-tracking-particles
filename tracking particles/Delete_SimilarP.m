
function [cent,fc1]=Delete_SimilarP(cent1,rad,fc)
cent=[];
a=[];b=[];

 for j=1:length(cent1(:,1))
     dist=sqrt((cent1(j,1)-cent1(:,1)).^2+(cent1(j,2)-cent1(:,2)).^2);
     a=find(dist<=rad);%
    if length(a)>1
        b=find(fc(a,1)~=0 & fc(a,2)~=0);
        if isempty(b)==0
            c=setdiff(a,a(b(1)));
            cent1(c,3)=1000;
        else
            cent1(a(2:end),3)=1000;
        end      
    end
 end
 a1= cent1(:,3)==1000;
 cent1(a1,:)=[];
 fc(a1,:)=[];

cent=cent1;fc1=fc;