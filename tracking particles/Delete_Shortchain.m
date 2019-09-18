%%%删除过短的颗粒链及相应的p与nn1的值
function [c1,d1,av1,p1,nn1]=Delete_Shortchain(c,d,av,n,p,nn,th1,th2)
c1=c;
d1=d;
av1=av;
p1=p;
nn1=nn;
atemp=find(abs(c(:,6))<th1 | c(:,2)<=th2); %找到需要删除的链
if isempty(atemp)==0
[ind1,~]=find(nn1(:,8)==c(atemp,1).');%需要删除的nn1的索引
ptable=tabulate(nn1(ind1,1));
ptable(ptable(:,2)==0,:)=[];
[I1,~]=find(n==ptable(:,1).');%需要修正的p的索引

c1(atemp,:)=[];
d1(atemp,:)=[];
av1(atemp,:)=[];
nn1(ind1,:)=[];
p1(I1)=p1(I1)-ptable(:,2);
end

