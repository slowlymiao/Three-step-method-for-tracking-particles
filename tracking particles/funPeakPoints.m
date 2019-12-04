function [XY,m] = funPeakPoints(E,Th)
a=size(E,1);b=size(E,2);
Ex=diff(E); 
i=a-1;
while i>0
    Ex(i+1,:)=Ex(i,:);
    i=i-1;
end   
E=E';           
Ey=diff(E);
i=b-1;
while i>0
    Ey(i+1,:)=Ey(i,:);
    i=i-1;
end  
Ey=Ey';E=E';
clear i j;

m=0; 
for i=1:a-1
    for j=1:b-1
        if Ex(i,j)*Ex(i+1,j)<0 && Ey(i,j)*Ey(i,j+1)<0     
            m=m+1;
            xy(m,1)=i;   
            xy(m,2)=j;
            xy(m,3)=E(i,j);
        end
    end
end
if m==0
    XY=0;
    m=0;
else
    ind=find(xy(:,3)>Th);
    m=length(ind);
    XY=xy(ind,:);
end