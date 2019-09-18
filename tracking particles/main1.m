%%%%%%%%%three step tracking method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
%%%%%%%%%%%%%%%%%first two step with Kalman filter%%%%%%%%%%%%%%%%%%%
%-------------pre-processing parameters----------------------------------
format long g
th=0.4;
rad=8;
rangec=1.5*rad;

se=strel('disk',30);
se1=strel('disk',4);
se2=strel('disk',3);
%--------------kalman filter parameters-------------------------------------------
Q=[0.3,3];
R=[0.2,2];
P=[0.3,3];
%--------------file information----------------------------------------------------
StartFile=1;%第一帧结果已知
EndFile=22;
FileStep=1;
IPath='C:\Users\PIVGroup\Desktop\image\';
maxn=0;

for nFile=StartFile:FileStep:EndFile-1
    
    if nFile>StartFile && points(nFile-FileStep)~=0
        fcentt=centi(:,1:2)+dis(:,1:2);
        fcentt(:,3:5)=centi(:,3:5);
        fdis=dis;
        clear centt dis centi centt
    end
    %----------Pre-processing----------------------------------
    filename1=[IPath,num2str(nFile),'.mat'];
    filename2=[IPath,num2str(nFile+FileStep),'.mat'];
    load(filename1); Image1=image;
    load(filename2); Image2=image; clear image
    Image1=im2double(imgPreProcess(Image1,se));%预处理
    Image2=im2double(imgPreProcess(Image2,se));%预处理
    
    
    %---------Two subtraction------------------------------------------
    [cent1,fc]=Two_Subtraction(Image1,Image2,rad,5);
    cent2=[];
    if isempty(cent1)==0
        [cent2,fc1]=FindCent(Image1,cent1,fc,rad,rangec);
        if isempty(cent2)==0
            [cent2,~]=Delete_SimilarP(cent2,ceil(1.0*rad),fc1);
        end
    end
    clear fc fc1 cent1
    %-------initialization---------------------------------------------------------
    if isempty(cent2)==0
        fp=zeros(length(cent2(:,1)),5);
        fp(:,4)=P(1);
        fp(:,5)=P(2);
        cent2(:,4:5)=0;
        %------------Associate with previous frame-------------------------------------
        if nFile>StartFile && points(nFile-FileStep)~=0
            m1=fcentt';
            m2=cent2(:,1:3)';
            [rp1,rp2]=Data_Association(m1,m2,2*rad);
            if isempty(rp1)==0 %%存在关联值
                atemp=tabulate(rp1);
                btemp=find(atemp(:,2)>1);
                ctemp=setdiff(rp1,btemp);
                if isempty(btemp)==0
                    for ii=1:length(btemp)
                        temp=find(rp1==atemp(btemp(ii)));
                        dd=sqrt((cent2(rp2(temp),1)-fcentt(rp1(temp),1)).^2+(cent2(rp2(temp),2)-fcentt(rp1(temp),2)).^2);
                        tt=find(dd==min(dd));
                        fcentt(rp1(temp(1)),1:3)=round(0.5*cent2(rp2(temp(tt)),1:3)+0.5*fcentt(rp1(temp(1)),1:3));
                    end
                end
                if isempty(ctemp)==0
                    for ii=1:length(ctemp)
                        temp=find(rp1==atemp(ctemp(ii)));
                        fcentt(rp1(temp),1:3)=round(0.5*cent2(rp2(temp),1:3)+0.5*fcentt(rp1(temp),1:3));
                        clear temp dd
                        %                         end
                    end
                end
                clear atemp btemp ctemp
                cent2(rp2,:)=[];
                fp(rp2,:)=[];
                
            end
            cent=[round(fcentt);round(cent2)];
            fp=[fp1;fp];
        else
            cent=round(cent2);
        end
        clear cent2
    else
        if nFile>StartFile && points(nFile-FileStep)~=0
            cent=round(fcentt);
            fp=fp1;
        else
            cent=[];
            fp=[];
        end
    end
    
    points(nFile)=0;
    centi=zeros(1,5);
    dis=zeros(1,3);
    centt=zeros(1,2);
    
    
    if isempty(cent)==0
        %-------------cross-correlation------------------------------------------
        [s1,dis,fp1]=Cro_correlation(cent',Image1,Image2,fp,th,Q,R,rad);%centi为初始点位置标志，dis为位移；采用原始图片
        
        
        if sum(sum(s1))>0
            btemp=find(s1(1,:)==0);
            centi=s1(:,s1(1,:)~=0)';
            atemp=find(centi(:,5)==0);
            centi(atemp,5)=maxn+1:maxn+length(atemp);
            centi(:,4)=centi(:,4)+1;
            dis(:,btemp)=[];
            dis=dis';
            fp1(btemp,:)=[];
            fp1(:,3)=dis(:,3);
            
            [~,a,b]=Auto_Association(centt',rad);
            if isempty(a)==0
                m1=[];
                for aa=1:length(a)
                    if centi(a(aa),4)>centi(b(aa),4)
                        m1=[m1,b];
                    elseif centi(a(aa),4)<centi(b(aa),4)
                        m1=[m1,a];
                    else
                        m1=[m1,b];
                    end
                end
                
                centi(m1,:)=[];
                dis(m1,:)=[];
                fp1(m1,:)=[];
                centt(m1,:)=[];
            end
            points(nFile)=size(dis,1);
            maxn=maxn+length(atemp);
            clear atemp btemp a b cent
        end
    end
    
    Result1(nFile,1)=nFile;
    Result1(nFile,2)=points(nFile);%number of active particles
    Result2(sum(points(1:nFile))-points(nFile)+1:sum(points(1:nFile)),1)=nFile;
    Result2(sum(points(1:nFile))-points(nFile)+1:sum(points(1:nFile)),2)=1:points(nFile);
    Result2(sum(points(1:nFile))-points(nFile)+1:sum(points(1:nFile)),3:8)=[centi(:,1:2),dis(:,1:2),centi(:,3),centi(:,5)];

end
 %%%%%%%%%step3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th1=2;
th2=10;
ss=1;
%----statistic chain----------------------------------------------------------------
[c1,d1,av1]=Statistic_Chain(Result2);
rad1=round(mean(Result2(:,7))*2);

%----------------------------------------------------

mu=mean(c1(:,6)./c1(:,2));%平均位移
[c11]=Associate_BC(c1,d1,th1,th2,mu,ss);%链接
if isempty(find(c11(:,5)>0))==0
    [nn3,p3]=Repair_BC(c11,Result2,Result1(:,1),Result1(:,2),d1,ss);
    [c3,~,~]=Statistic_Chain(nn3);
    [c4,~,~,p4,nn4]=Delete_Shortchain(c3,d1,av1,n1,p3,nn3,rad1,th);
end







