%%%%%re-find center----------------------------------
function [cent1]=FindCent(Image1,cent,rad,t)

cent1=[];
j=1;
for i=1:length(cent(:,1)) 
    [centn]=S_centn(cent(i,1),cent(i,2),Image1,rad,t);
    if size(centn,2)>1
        dist=sqrt((cent(i,1)-centn(:,1)).^2+(cent(i,2)-centn(:,2)).^2);
        a=find(dist==min(dist));  
        if isempty(a)==0 && centn(a(1),3)~=0
            cent1(j,1:3)=centn(a(1),1:3);
            j=j+1;
        end
    elseif isempty(centn)==0
        cent1(j,1:3)=centn(1,1:3);
        j=j+1;
    end
    
    
end


