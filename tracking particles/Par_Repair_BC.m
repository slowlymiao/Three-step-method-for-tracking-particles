%%--------------------------------------
function [nf,p2,nn2]=Par_Repair_BC(nn,p,atemp,c,d,n,ind1,ind2,tt,ss,rad)
nf=[];
nn2=[];
p2=zeros(size(p));
 if c(atemp,6)>0 
            ptemp1=c(atemp,4)+ss:ss:c(tt,3)-ss;
            xyb(:,1:2)=nn(ind1(end-1:end),3:4);%��һ���ӵ������֡λ��
            xya(:,1:2)=nn(ind2(1:3),3:4);%��һ���ӵ��ʼ����֡��λ��
            if d(atemp,4)>d(tt,2)
                x(1)=d(atemp,3);%����ĵ�һ֡
                y(1)=d(atemp,4);
                
                if c(atemp,6)>1 %%����֡����������ֵ
                    xn=[1:3,3+c(atemp,6):5+c(atemp,6)];%��������Ա���
                    xni=4:2+c(atemp,6);%��Ҫ��ֵ���Ա���
                    xs=[xyb(:,1);x(1);xya(:,1)];%������
                    ys=[xyb(:,2);y(1);xya(:,2)];%������
                    xi=interp1(xn,xs,xni,'spline');%������ֵ
                    yi=interp1(xn,ys,xni,'spline');%������ֵ
                    xi=[x(1),xi];
                    yi=[y(1),yi];
                    disx=[diff(xi),xya(1,1)-xi(end)];
                    disy=[diff(yi),xya(1,2)-yi(end)];
                else
                    xi=x(1);
                    yi=y(1);
                    disx=xya(1,1)-xi;
                    disy=xya(1,2)-yi;
                end
               
                clear xyb xya
            else %%�Ȳ�ֵ��У��------------------------------------------
                xn=[1:2,3+c(atemp,6):5+c(atemp,6)];
                xni=3:2+c(atemp,6);
                xs=[xyb(:,1);xya(:,1)];
                ys=[xyb(:,2);xya(:,2)];
                xi=interp1(xn,xs,xni,'spline');
                yi=interp1(xn,ys,xni,'spline');
                disx=[diff(xi),xya(1,1)-xi(end)];
                disy=[diff(yi),xya(1,2)-yi(end)];
                
                %----����У��---------------------------------------------
                nxi=flipud([nn(ind1,3);xi']);%��ҪУ����������
                nyi=flipud([nn(ind1,4);yi']);
                pxi(1)=nn(ind2(1),3)-nn(ind2(1),5);%Ԥ��ֵ
                pyi(1)=nn(ind2(1),4)-nn(ind2(2),6);
                ndisx=flipud([nn(ind1,5);disx']);
                ndisy=flipud([nn(ind1,6);disy']);
                
                disp=sqrt((pxi(1)-nxi(1))^2+(pyi(1)-nyi(1))^2);
                j=1;
                while disp>8 && j<=length(nxi)
                    nxi(j)=pxi(j)+0.1*(nxi(j)-pxi(j)); %У�����λ��
                    nyi(j)=pyi(j)+0.1*(nyi(j)-pyi(j));
                    if j==1
                        ndisx(1)=nn(ind2(1),3)-nxi(1);
                        ndisy(1)=nn(ind2(1),4)-nyi(1);
                    else
                        ndisx(j)=nxi(j-1)-nxi(j);%У�����λ��
                        ndisy(j)=nyi(j-1)-nyi(j);
                    end
                    if ndisy(j)>0
                        ndisy(j)=-1;
                    end
                   if j==length(nxi)
                       j=j+1;
                   else
                    pxi(j+1)=nxi(j)-ndisx(j);
                    pyi(j+1)=nyi(j)-ndisy(j);
                    j=j+1;
                    disp=sqrt((pxi(j)-nxi(j))^2+(pyi(j)-nyi(j))^2);
                   end
                end
                %----�Բ�ֵ�������¸�ֵ---------------------------------
               
                xi=flipud(nxi(1:length(xi)));
                yi=flipud(nyi(1:length(yi)));
                disx=flipud(ndisx(1:length(disx)));
                disy=flipud(ndisy(1:length(disy)));
                %----��ind1�������¸�ֵ(ɾ�������)---------------------------------
                nn1(:,1:2)=nn(ind1,1:2);
                nn1(:,3)=flipud(nxi(length(xi)+1:end));
                nn1(:,4)=flipud(nyi(length(xi)+1:end));
                nn1(:,5)=flipud(ndisx(length(xi)+1:end));
                nn1(:,6)=flipud(ndisy(length(xi)+1:end));
                nn1(:,7)=nn(ind1,7);
                nn1(:,8)=c(atemp,8);
                nf=ind1';
            end
             %---��ֵ---------------------------------------------------
                  n1=c(atemp,4)+ss:ss:c(tt,3)-ss;
                  m=zeros(length(xi),1);
                  r=rad*ones(size(m));
                  ff=c(atemp,8)*ones(size(m));
                  nntemp=[n1',m,xi',yi',disx',disy',r,ff];
                  if d(atemp,4)>d(tt,2)
                  nn2=nntemp;%��ֵ�ĵ�
                  else
                      nn2=[nn1;nntemp];
                  end
                  
                  [ip,~]=find(n==n1'.');%��Ҫ�仯��p2��λ��
                  p2(ip)=p2(ip)+1;
            clear ind1 ind2 ptemp1  
            
        elseif c(atemp,6)<0 %��Ҫ��ȥ��֡
            ptemp2=c(tt,4)+ss+ss:ss:c(atemp,4)-ss;
            [ip,~]=find(n==ptemp2'.');%��Ҫ�仯��p2��λ��
            p2(ip)=p2(ip)-1;
            ftemp=ind2(1:length(ptemp2));
            nf=ftemp';
        end