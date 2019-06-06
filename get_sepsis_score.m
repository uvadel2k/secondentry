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
for c=36:40
    x(isnan(x(:,c)),c)=0;
end
binary=false(nx,1);
binary(36:38)=1;

t=x(:,end);

% [x,xmeas]=samplehold(x,t);
[x,xmeas]=lastsample(x,t);
z=steprisk(model,x,xmeas,binary);
z(isnan(z))=0;
disp(z)
t1=model.t1;
nt=length(t1);
k=find(t(end)>=t1,1,'last');
if isempty(k),k=nt;end
b=model.b(:,k);
T=model.T(k);
mod=model.mod;
xb=[1 z(mod)]*b;
label=xb>T;
score=exp(xb)/(1+exp(xb));
%score=glmval(b,z(end,mod),'logit');
%T=model.Tmax;
%label=score>T;
    
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

function z=steprisk(model,x,xmeas,binary)

[nr,nx]=size(x);
z=zeros(nr,nx);
z(:,binary)=x(:,binary);
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
    if binary(c),continue,end
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

