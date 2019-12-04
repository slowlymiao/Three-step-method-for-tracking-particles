%%%%repair the broken chain%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nn1,p1]=Repair_BC(c,nn,n,p,d,ss)
nn1=nn;
p2=zeros(size(p));
flag=nn(:,8);
nn2=[];
nf=[];
ind=[];
fp=[];
rad=mean(nn(:,7));
atemp=find(c(:,5)>0);


    for i=1:length(atemp)
            tt=find(c(:,1)==c(atemp(i),5));
            ind1=find(flag==c(atemp(i),1));
            ind2=find(flag==c(tt,1));
            [nf22,p22,nn22]=Par_Repair_BC(nn,p,atemp(i),c,d,n,ind1,ind2,tt,ss,rad);
            p2=p2+p22;
            nf=[nf;nf22];
            nn2=[nn2;nn22];
            ind=[ind;ind2];
            fp=[fp;c(atemp(i),8)*ones(length(ind2),1)];
    end
   
    p1=p+p2;
    nn1(nf,:)=[];
    nn1=[nn1;nn2];
    nn1(ind,8)=fp;
    
