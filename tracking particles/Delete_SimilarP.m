%----delete the similar paticle-------------------------------
function [a,b]=Delete_SimilarP(m1,m2,rad,mode)
a=[];b=[];
dist=1000*ones(size(m1,2),size(m2,2));
for j=1:size(m1,2)
    
    if mode==1 % auto_association
        dist(j,j+1:end)=sqrt((m1(1,j)-m2(1,j+1:end)).^2+(m1(2,j)-m2(2,j+1:end)).^2);
    elseif mode==2
        dist(j,:)=sqrt((m1(1,j)-m2(1,:)).^2+(m1(2,j)-m2(2,:)).^2);
    end
end
if mode==1
    [a,b]=find(dist<=rad);
elseif mode==2
    for jj=1:size(dist,1)
        aa=find(dist(jj,:)<=rad & dist(jj,:)==min(dist(jj,:)));
        if isempty(aa)==0
            b(jj)=aa(1);
            a(jj)=jj;
            
        else
            b(jj)=0;
        end
    end
    a(a==0)=[];
    b(b==0)=[];
    
end



