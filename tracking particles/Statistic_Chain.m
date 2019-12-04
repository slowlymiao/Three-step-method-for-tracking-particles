

function [c,d,av]=Statistic_Chain(nn1)

 ntable=tabulate(nn1(:,8));
 ntable(ntable(:,2)==0,:)=[];
 [B, I1]=unique(nn1(:,8), 'first');
 [~,I2]=unique(nn1(:,8),'last');
c(:,1)=B;
 c(:,2)=ntable(:,2);
 c(:,3)=nn1(I1,1);
 c(:,4)=nn1(I2,1);
 c(:,5)=nn1(I2,3)+nn1(I2,5)-nn1(I1,3);
 c(:,6)=nn1(I2,4)+nn1(I2,6)-nn1(I1,4);
 d(:,1:2)=nn1(I1,3:4);
 d(:,3:4)=nn1(I2,3:4)+nn1(I2,5:6);
 av(:,1:2)=nn1(I1,5:6);
 av(:,3:4)=nn1(I2,5:6);
 
 nc=(c(:,4)-c(:,3))/2+1;
 ind=find((nc-c(:,2))~=0);
 
