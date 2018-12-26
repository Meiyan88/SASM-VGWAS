% perform GSIS procedure

function [pp, W]=globalWald(SNP,YW,PX,deta,Flag_Z)


E3=eye(size(PX,1))-PX;
E4=YW*(repmat(deta,[1,size(YW,1)]).*YW');
if Flag_Z
    c=size(SNP,2)/2;%number of SNPs
    t1=E3*E4*E3;
    j=1;
    W=zeros(c,1);
    for i=1:c
        t2=SNP(:,j:j+1)'*E3*SNP(:,j:j+1);
        t3=SNP(:,j:j+1)'*t1*SNP(:,j:j+1);
        W(i)=trace(t2\t3);
        j=j+2;
    end
else
    Q=SNP'*E3;
    t1=sum(Q.^2,2);
    t1=t1.^(-1);
    [U,D]=eig(E4);
    D = real(D);
    D(D(:)<0)=0;
    t4=U*(D^(1/2))*U';
    t4=real(t4);
    t2=sum((Q*t4).^2,2);%C*1
    W=t1.*t2;
end
W=W./length(deta);

k1=mean(W);
k2=var(W);
k3=mean((W-k1).^3);
a=k3/(4*k2);
b=k1-2*k2^2/k3;
d=8*k2^3/k3^2;
pv=1-cdf('chi2',(W-b)/a,d);
pp=-log10(pv);
