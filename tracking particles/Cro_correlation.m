
function [s,dis,fp]=Cro_correlation(s1,Image1,Image2,fp,th,QC,RC,rad)

dis=zeros(5,size(s1,2));%displacement

for i=1:size(s1,2)
    r=round(s1(3,i)/2*0.8);
    r(r>rad)=round(r*0.8);
    if r<rad
        r=rad;
    end
    [BL,BR,BU,BD,s1(1:2,i)]=Cro_boundary(s1(1:2,i),[r,r],[r,r],size(Image1,1),size(Image1,2),1,r);
    A=Image1(BU:BD,BL:BR);
    B=Image2(BU:BD,BL:BR);
    COR=corr2(A,B);
    if (COR>0.9 && s1(4,i)==0) || COR>0.95
        s1(:,i)=0;
    else
        [L,R,U,D,~]=Cro_boundary(s1(1:2,i),[r,r],[ceil(4*r),0],size(Image1,1),size(Image1,2),2,r);
        ii=1;
        for nn=L:R
            jj=1;
            for mm=U:D
                B1=Image2(mm-r:mm+r,nn-r:nn+r);
                COR2(ii,jj)=corr2(A,B1);
                jj=jj+1;
            end
            ii=ii+1;
        end
        a=s1(1,i)-L+1;
        b=s1(2,i)-U+1;
        [XY,m]=funPeakPoints(COR2,th);%
        if m==0
            s1(:,i)=0;
        else
            if s1(4,i)==0
                XYmax=max(XY(:,3));
                [a1,b1]=find(COR2==XYmax);
                if length(a1)>1
                    a1(2:end)=[];
                    b1(2:end)=[];
                end
                if a1==1 || a1==size(COR2,1)
                    z=a1;
                else
                    z=three_point_interplot_Gauss(a1,COR2(a1-1,b1),COR2(a1,b1),COR2(a1+1,b1),a);
                end
                if b1==1 || b1==size(COR2,2)
                    x=b1;
                else
                    x=three_point_interplot_Gauss(b1,COR2(a1,b1-1),COR2(a1,b1),COR2(a1,b1+1),b);
                end
                if x<=0 && sqrt(x^2+z^2)>=1 && (XYmax-COR)/abs(COR)>0.05
                    dis(1,i)=z;
                    dis(2,i)=x;
                else
                    s1(:,i)=0;
                end
            else
                dis1(1)=s1(1,i)+fp(i,1);
                if s1(4,i)==1
                    dis1(2)=s1(2,i)+fp(i,2);
                else
                    dis1(2)=s1(2,i)+fp(i,2)+0.5*(fp(i,2)-fp(i,3));
                end
                
                z1=s1(1,i)+XY(:,1)-a;
                x1=s1(2,i)+XY(:,2)-b;
                dist=sqrt((dis1(1)-z1).^2+(dis1(2)-x1).^2);
                clear z1 x1
                [dmin,d]=min(dist./max(dist)./XY(:,3));
                XYmax=XY(d,3);
                [a1,b1]=find(COR2==XYmax);
                if length(a1)>1
                    a1(2:end)=[];
                    b1(2:end)=[];
                end
                
                if a1==1 || a1==size(COR2,1)
                    z=a1;
                else
                    z=three_point_interplot_Gauss(a1,COR2(a1-1,b1),COR2(a1,b1),COR2(a1+1,b1),a);
                end
                if b1==1 || b1==size(COR2,2)
                    x=b1;
                else
                    x=three_point_interplot_Gauss(b1,COR2(a1,b1-1),COR2(a1,b1),COR2(a1,b1+1),b);
                end
                while x>0 && m>1
                    XY(d,:)=[];
                    dist(d)=[];
                    [dmin,d]=max(XY(:,3));
                    XYmax=XY(d,3);
                    [a1,b1]=find(COR2==XYmax);
                    if length(a1)>1
                        a1(2:end)=[];
                        b1(2:end)=[];
                    end
                    
                    if a1==1 || a1==size(COR2,1)
                        z=a1;
                    else
                        z=three_point_interplot_Gauss(a1,COR2(a1-1,b1),COR2(a1,b1),COR2(a1+1,b1),a);
                    end
                    if b1==1 || b1==size(COR2,2)
                        x=b1;
                    else
                        x=three_point_interplot_Gauss(b1,COR2(a1,b1-1),COR2(a1,b1),COR2(a1,b1+1),b);
                    end
                    m=m-1;
                end
                if x<0 && sqrt(x^2+z^2)>=1 && dmin*XYmax<1.2*rad
                    kg=fp(i,4:5)./(fp(i,4:5)+RC);
                    p=(1-kg).*fp(i,4:5);
                    fp(i,4:5)=p+QC;
                    dis2(1)=dis1(1)+kg(1)*(s1(1,i)+z-dis1(1));
                    dis2(2)=dis1(2)+kg(2)*(s1(2,i)+x-dis1(2));
                    dis(1,i)=dis2(1)-s1(1,i);
                    dis(2,i)=dis2(2)-s1(2,i);
                else
                    if s1(4,i)>0
                        dis(1,i)=dis1(1)-s1(1,i);
                        dis(2,i)=dis1(2)-s1(2,i);
                    else
                        s1(:,i)=0;
                    end
                end
            end
            dis(3,i)=fp(i,2);
        end
        clear COR2
    end
    
end
dis(4:5,:)=fp(:,4:5)';
s=s1(:,s1(1,:)~=0)';
btemp=find(s1(1,:)==0);
dis(:,btemp)=[];
dis=dis';
fp(btemp,:)=[];

