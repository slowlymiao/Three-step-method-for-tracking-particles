%%%%%re-find center----------------------------------
function [cent1,fc1]=FindCent(Image1,cent,fc,rad,t)

cent1=[];
fc1=[];

j=1;
for i=1:length(cent(:,1))
    
    [centn]=S_centn(cent(i,1),cent(i,2),Image1,rad,t,fc(i,:));%t=20ΪѰ�ҷ�Χ

    if size(centn,2)>1
        dist=sqrt((cent(i,1)-centn(:,1)).^2+(cent(i,2)-centn(:,2)).^2);
        a=find(dist==min(dist));
       
        if isempty(a)==0 && centn(a(1),3)~=0
            cent1(j,1:3)=centn(a(1),1:3);
            fc1(j,:)=fc(i,:);
            
            j=j+1;
        end
    elseif isempty(centn)==0
        cent1(j,1:3)=centn(1,1:3);
        fc1(j,:)=fc(i,:);
        
        j=j+1;
    end
    
    
end


