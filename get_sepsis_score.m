function [score, label] = get_sepsis_score(x,model)

xname={'HR','O2Sat','Temp','SBP','MAP','DBP','Resp','EtCO2','BaseExcess','HCO3',...
    'FiO2','pH','PaCO2','SaO2','AST','BUN','Alkalinephos','Calcium','Chloride',...
    'Creatinine','Bilirubin_direct','Glucose','Lactate','Magnesium','Phosphate',...
    'Potassium','Bilirubin_total','TroponinI','Hct','Hgb','PTT','WBC','Fibrinogen',...
    'Platelets','Age','Gender','Unit1','Unit2','HospAdmTime','ICULOS'}';    

tmax=cat(1,model.tmax);
t=x(:,end);
tt=t(end);
if tt>max(tmax)
    label=1;
    score=20;
    return
end
tmin=cat(1,model.tmin);
if tt<1,tt=1;end
k=find(tt>=tmin,1,'last');
tmodel=model(k);

nx=length(xname);

[x,xmeas]=lastsample(x,t);

%Demographic data
for c=35:39
    if isnan(x(c))
        x(c)=0;
    end
    xmeas(c)=0;
end
maxmeas=336;
meas=tmodel.meas;
for i=1:length(meas)
    if meas(i)>=maxmeas,continue,end
    if xmeas(i)>=meas(i)
        x(i)=NaN;
    end
end

res=tmodel.res;
for i=1:length(res)
    if isnan(x(i)),continue,end    
    if isnan(res(i)),continue,end
    x(i)=round(res(i)*x(i))/res(i);
end

binary=tmodel.binary;
C=tmodel.C;
v1=tmodel.v1;
v2=tmodel.v2;
b0=tmodel.b(1);
b=tmodel.b(2:end);
z=zeros(1,nx);
for i=1:nx
    xx=x(i);    
    if binary(i)
        if xx==1
            j=find(C==i);
            if length(j)==1
                z(i)=b(j);                
            end
        end
        continue
    end
    if isnan(xx)
        j=find(C==i&isnan(v1));
    else
        j=find(C==i&~isnan(v1));
        x1=v1(j);
        x2=v2(j);
        k=find(xx>=x1&xx<x2);
        if length(k)==1
            j=j(k);
        end
    end
    if length(j)==1
        z(i)=b(j);
    end
end

score=b0+sum(z);
label=score>=0;
    
end

function [x,xmeas]=lastsample(data,t)

maxtime=28*24;
[nr,nx]=size(data);
t2=t(nr);
x=data(nr,:);
xmeas=maxtime*ones(1,nx);
for i=1:nx
    if ~isnan(x(i))
        xmeas(i)=0;
        continue
    end
    xx=data(:,i);
    k=find(~isnan(xx),1,'last');
    if isempty(k),continue,end
    t1=t(k);
    x(i)=xx(k);
    xmeas(i)=t2-t1;
end

end
