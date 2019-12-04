
function [c11]=Associate_BC(c,d,th1,th2,mu,ss)
c11=c(:,1:4);
D=bsxfun(@minus, c(:,3), c(:,4)');
D=D./ss-1;
D(D<-th1)=1000;
[D1,I]=sort(D);
D1(6:end,:)=[];
I(6:end,:)=[];
D1=D1';
I=I';
m=1:length(I);
clear D
for i=1:5
   atemp=find(I(:,i)==m');
   DIS(:,i)=sqrt((d(I(:,i),1)-d(:,3)).^2+(d(I(:,i),2)-d(:,4)).^2).*(-d(I(:,i),2)+d(:,4))./abs(d(I(:,i),2)-d(:,4)) ;
   DIS(atemp,i)=nan;
end

[ix,iy]=find(DIS./abs(DIS).*D1./abs(D1)>0 & abs(DIS)<abs(th2*mu*0.8));
ind=sub2ind(size(DIS),ix,iy);
ctemp=tabulate(ix);
cind=ctemp(ctemp(:,2)>1,1);
[cinds,ii]=setdiff(ix,cind);
c11(cinds,5)=c11(I(ind(ii)),1);
ii1=setdiff(1:length(ix),ii);
 [ix1, ind1]=unique(ix(ii1), 'first');
c11(ix1,5)=c11(I(ind(ii1(ind1))));
c11(cinds,6)=D1(ind(ii));
c11(cinds,7)=DIS(ind(ii));
c11(ix1,6)=D1(ind(ii1(ind1)));
c11(ix1,7)=DIS(ind(ii1(ind1)));
clear DIS D1 m I ix iy ind ctemp cind cinds ii ii1 ind1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ctable=tabulate(c11(:,5));
atemp=find(ctable(:,2)>1 & ctable(:,1)~=0);
if isempty(atemp)==0
    for i=1:length(atemp)
        btemp=find(c11(:,5)==ctable(atemp(i),1));
        ind=find(c11(btemp,6)==0);
        if length(ind)>1
            ctemp=abs(c11(btemp,7));
        else  
            ctemp=abs(c11(btemp,6).*c11(btemp,7));
        end
         tt=find(ctemp>min(ctemp));
        c11(btemp(tt),5:7)=0;  
        clear btemp ctemp tt
    end
end
clear ctable atemp

c11(:,8)=0;
atemp=find(c11(:,5)>0);
c11(atemp,8)=c11(atemp,1);
ctemp=[c11(atemp,1),c11(atemp,5)];
[atemp2,~]=find(ctemp(:,1)==ctemp(:,2).');
ctemp(atemp2,:)=[];
cind=setdiff(atemp,atemp2);
ic=1;in=1;
while isempty(ic)==0 
[atemp1,~]=find(c11(:,1)==ctemp(:,in+1).' );
cc=find(ctemp(:,in+1)>0);
cind(cc,in+1)=atemp1;
ctemp1=c11(atemp1,5);
ctemp1(isnan(ctemp1))=0;
 ic=find(ctemp1>0);
 iic=find((ctemp1(ic)-ctemp(cc(ic),in+1))<0);
 ctemp1(ic(iic))=0;
 if sum(ctemp1)>0
ctemp(cc,in+2)=ctemp1;
 else
     ic=[];
 end
 c11(atemp1,8)=ctemp(cc,1);
in=in+1;
end
clear atemp1 atemp2 ctemp1 ic iic 



