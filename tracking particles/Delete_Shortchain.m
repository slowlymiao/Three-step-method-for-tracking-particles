%%%ɾ�����̵Ŀ���������Ӧ��p��nn1��ֵ
function [c1,d1,av1,p1,nn1]=Delete_Shortchain(c,d,av,n,p,nn,th1,th2)
c1=c;
d1=d;
av1=av;
p1=p;
nn1=nn;
atemp=find(abs(c(:,6))<th1 | c(:,2)<=th2); %�ҵ���Ҫɾ������
if isempty(atemp)==0
[ind1,~]=find(nn1(:,8)==c(atemp,1).');%��Ҫɾ����nn1������
ptable=tabulate(nn1(ind1,1));
ptable(ptable(:,2)==0,:)=[];
[I1,~]=find(n==ptable(:,1).');%��Ҫ������p������

c1(atemp,:)=[];
d1(atemp,:)=[];
av1(atemp,:)=[];
nn1(ind1,:)=[];
p1(I1)=p1(I1)-ptable(:,2);
end

