function [score, label] = get_sepsis_score(data, model)

xname={'HR','O2Sat','Temp','SBP','MAP','DBP','Resp','EtCO2','BaseExcess','HCO3',...
    'FiO2','pH','PaCO2','SaO2','AST','BUN','Alkalinephos','Calcium','Chloride',...
    'Creatinine','Bilirubin_direct','Glucose','Lactate','Magnesium','Phosphate',...
    'Potassium','Bilirubin_total','TroponinI','Hct','Hgb','PTT','WBC','Fibrinogen',...
    'Platelets','Age','Gender','Unit1','Unit2','HospAdmTime','ICULOS'}';    

t=data(:,end);
[x,xmeas]=samplehold(data,t);
z=steprisk(model,x,xmeas);
p=evalmodel(model,z);

T=model.T;
score=p(end);
label=score>T;
    
end

function [x,xmeas]=samplehold(x,t)

[nr,nx]=size(x);
xmeas=-ones(nr,nx);
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

[nr,nx]=size(x);
t2=t(nt);
x=data(nt,:);
for i=1:nx
    if ~isnan(x(i)),continue,end
    xx=data(:,i);
    k=find(~isnan(xx),1,'last');
    if isempty(k),continue,end
    x(i)=xx(k);
end

end

function z=steprisk(model,x,xmeas)

[nr,nx]=size(x);
z=NaN*ones(nr,nx);

if nargin<3
    xmeas=zeros(nr,nx);
end

pc=model.pc;
v1=model.v1;
v2=model.v2;
dp=model.dp;
uy=model.uy;

zp=logit(dp)-logit(uy);

nd=length(pc);

for i=1:nd
    j=pc(i);
    if ~isnan(v1(i))
        sub=xmeas(:,j)>=0&x(:,j)>=v1(i);
        if v2(i)<Inf
            sub=sub&x(:,j)<v2(i);
        end
    else
        sub=xmeas(:,j)<0;
    end
    z(sub,j)=zp(i);
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
