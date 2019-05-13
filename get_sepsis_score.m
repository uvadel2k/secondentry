function [score, label] = get_sepsis_score(x,model)

xname={'HR','O2Sat','Temp','SBP','MAP','DBP','Resp','EtCO2','BaseExcess','HCO3',...
    'FiO2','pH','PaCO2','SaO2','AST','BUN','Alkalinephos','Calcium','Chloride',...
    'Creatinine','Bilirubin_direct','Glucose','Lactate','Magnesium','Phosphate',...
    'Potassium','Bilirubin_total','TroponinI','Hct','Hgb','PTT','WBC','Fibrinogen',...
    'Platelets','Age','Gender','Unit1','Unit2','HospAdmTime','ICULOS'}';    

nx=length(xname);

%Preprocess demographic data
c=strmatch('Age',xname);
x(x(:,c)==100,1)=90;
c=strmatch('Unit1',xname);
x(isnan(x(:,c)),c)=0;
c=strmatch('Unit2',xname);
x(isnan(x(:,c)),c)=0;
c=strmatch('HospAdmTime',xname);
x(isnan(x(:,c)),c)=0;

mod=model.mod;
t=x(:,end);

% [x,xmeas]=samplehold(x,t);
[x,xmeas]=lastsample(x,t);
z=steprisk(model,x,xmeas);
z(isnan(z))=0;
b=model.b;
score=glmval(b,z(end,mod),'logit');
T=model.Tmax;
label=score>T;
    
end

function [x,xmeas]=samplehold(x,t)

[nr,nx]=size(x);
maxtime=28*24;
xmeas=maxtime*ones(nr,nx);
nt=length(t);
x1=x(1,:);
t1=NaN*ones(1,nx);
j=find(~isnan(x1));
t1(j)=t(1);
xmeas(1,j)=0;
for k=2:nt
    x2=x(k,:);
    j=find(~isnan(x2));
    xmeas(k,j)=0;
    x1(j)=x2(j);
    t1(j)=t(k);        
    t2=t(k)*ones(1,nx);        
    j=find(isnan(x2)&~isnan(x1));        
    x2(j)=x1(j);
    x(k,:)=x2;
    xmeas(k,j)=t2(j)-t1(j);        
end

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

function z=steprisk(model,x,xmeas)

[nr,nx]=size(x);
z=zeros(nr,nx);

if nargin<3
    xmeas=zeros(nr,nx);
end

dc=model.dc;
v1=model.v1;
v2=model.v2;
dx=model.dx;
maxmeas=model.maxmeas;

nd=length(dc);

for i=1:nd
    c=dc(i);
    sub=xmeas(:,c)<=maxmeas;
    if ~isnan(v1(i))
        sub=sub&x(:,c)>=v1(i);
        if v2(i)<Inf
            sub=sub&x(:,c)<v2(i);
        end
    else
        sub=xmeas(:,c)>maxmeas;
    end
    z(sub,c)=dx(i);
end

end

function p=evalmodel(model,x)

b=model.b;
mod=model.mod;
p=glmval(b,x(:,mod),'logit');

end

function r=logit(p)
%r=logit(p)

r=log(p./(1-p));

end
