

function [rp1,rp2]=Data_Association(m1,m2,rad)

c=[];
if isempty(m1)==0 && isempty(m2)==0 
   for ii=1:length(m1(1,:))
    vs(:,ii)=sqrt((m1(1,ii)-m2(1,:)).^2+(m1(2,ii)-m2(2,:)).^2);
end
for jj=1:size(vs,1)
    a=find(vs(jj,:)<=rad & vs(jj,:)==min(vs(jj,:)));
    if isempty(a)==0
        b(jj)=a(1);
        c(jj)=jj;
        if length(a)>1
            aa=a(2:end);
            bb=jj;
        end
        else
        b(jj)=0;
    end
end  
b(b==0)=[];
c(c==0)=[];
rp1=b;rp2=c;
else
    rp1=0;rp2=0;
end
