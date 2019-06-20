function [score, label] = get_sepsis_score(x,model)

xname={'HR','O2Sat','Temp','SBP','MAP','DBP','Resp','EtCO2','BaseExcess','HCO3',...
    'FiO2','pH','PaCO2','SaO2','AST','BUN','Alkalinephos','Calcium','Chloride',...
    'Creatinine','Bilirubin_direct','Glucose','Lactate','Magnesium','Phosphate',...
    'Potassium','Bilirubin_total','TroponinI','Hct','Hgb','PTT','WBC','Fibrinogen',...
    'Platelets','Age','Gender','Unit1','Unit2','HospAdmTime','ICULOS'}';    

if x(end,end)>=60
    label=1;
    score=14524/82011;
    return
end

nx=length(xname);

%Preprocess demographic data
c=strmatch('Age',xname);
x(x(:,c)==100,1)=90;
for c=36:40
    x(isnan(x(:,c)),c)=0;
end

% [x,xmeas]=samplehold(x,t);
[x,xmeas]=lastsample(x,x(:,end));
vmeas=model.vmeas;
for i=1:length(vmeas)    
    if xmeas(i)>=vmeas(i)
        x(i)=NaN;
    end
end
res=model.res;
for i=1:length(res)
    if isnan(x(i)),continue,end    
    if isnan(res(i)),continue,end
    x(i)=round(res(i)*x(i))/res(i);
end
gc=model.gc;
v1=model.v1;
v2=model.v2;
gz=model.gz;
z=zeros(1,nx);
for i=1:(nx-1)
    xx=x(i);
    if isnan(xx)
        j=find(gc==i&isnan(v1));
    else
        j=find(gc==i&~isnan(v1));
        x1=v1(j);        
        k=find(xx>=x1,1,'last');
        if length(k)==1
            j=j(k);
        end
    end
    if length(j)==1
        z(i)=gz(j);
    end
end

binary=(36:38);
z(binary)=x(binary);

tz=model.tz;
t=x(end);
t=round(t);
if t<1,t=1;end
if t>length(tz)
    t=length(tz);
end
z(end)=tz(t);

z(isnan(z))=0;
mod=model.mod;
tmod=model.tmod;
b=model.b;
X=[1 z(mod) t*z(tmod)];

xb=X*b;
s=round(10*xb);
T=model.T;
label=s>=T;
pb=model.pb;
pxb=pb(1)+pb(2)*xb;
score=exp(pxb)/(1+exp(pxb));
    
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
