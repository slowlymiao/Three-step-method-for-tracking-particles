
function [BL,BR,BU,BD,cent]=Cro_boundary(cent,radz,radx,m,n,mode,rad)

BL=cent(1)-radz(1);
BR=cent(1)+radz(2);
BU=cent(2)-radx(1);
BD=cent(2)+radx(2);

if BL<1 && mode==1
    BL=1;
    BR=BR-(cent(1)-radz(1))+1;
    cent(1)=cent(1)+(cent(1)-radz(1))-1;
elseif BL<rad+1 && mode==2
    BL=rad+1;
    BR=BR-(cent(1)-radz(1))+1+rad;
end
if BR>n && mode==1
    BR=n;
    BL=BL-((cent(1)+radz(2))-n);
    cent(1)=cent(1)-((cent(1)+radz(2))-n);
elseif BR>n-rad && mode==2
    BR=n-rad;
    BL=BL-((cent(1)+radz(2))-n)-rad;
end
if BU<1 && mode==1
    BU=1;
    BD=BD-(cent(2)-radx(1))+1;
    cent(2)=cent(2)+(cent(2)-radx(1))-1;
elseif BU<rad+1 && mode==2
    BU=rad+1;
    BD=BD-(cent(2)-radx(1))+rad+1;
end
if BD>m && mode==1
    BD=m;
    BU=BU-(cent(2)+radx(2)-m);
    cent(2)=cent(2)-(cent(2)+radx(2)-m);
elseif BD>m-radz(2) && mode==2
    BD=m-rad;
    BU=BU-(cent(2)+radx(2)-m)-rad;
  
end
BL=round(BL);
BR=round(BR);
BU=round(BU);
BD=round(BD);
