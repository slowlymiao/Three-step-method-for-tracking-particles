
function [m1,a,b]=Auto_Association(m,rad)
m1=m;
vs=1000*ones(length(m),length(m));
if isempty(m)==0
    for i=1:length(m(1,:))-1
        vs(i+1:end,i)=sqrt((m(1,i)-m(1,i+1:end)).^2+(m(2,i)-m(2,i+1:end)).^2);
    end
    [a,b]=find(vs<rad);
    m1(:,a)=[];
end
    