function [ p,px ] = Analytical(SNR,K,M,N,method)
sum1=zeros(K,M);
sum2=zeros(M,M);
sum3=zeros(M,K);
sum4=zeros(K,K);
a=randi([0 1],2*K,1);
cnt=zeros(16);
cntx=zeros(4);
a=a*2-1;
x=modulate(a)';
H=normrnd(0,1,M,K)+j*normrnd(0,1,M,K);
Pt=(sum(abs(x).^2))/K;
sigma=Pt*10^(-SNR/10);
n=normrnd(0,sqrt(sigma),M,1) + j*normrnd(0,sqrt(sigma),M,1);
r=sqrt(Pt)*H*x+n;
y=quantizer(r,M);
%%求A的方法
switch(method)
    case 1  %LS
        for i=1:N
            b=randi([0 1],2*K,1);
            b=b*2-1;
            xn(:,i)=modulate(b)';
            rn=sqrt(Pt)*H*xn(:,i)+n;
            yn(:,i)=quantizer(rn,M);
            sum1=sum1+xn(:,i)*yn(:,i)';
            sum2=sum2+yn(:,i)*yn(:,i)';
        end
        A=sum1*inv(sum2);
    case 2  %MRC
        for i=1:N
            b=randi([0 1],2*K,1);
            b=b*2-1;
            xn(:,i)=modulate(b)';
            rn=sqrt(Pt)*H*xn(:,i)+n;
            yn(:,i)=quantizer(rn,M);
            sum3=sum3+sqrt(Pt)*yn(:,i)*xn(:,i)';
            sum4=sum4+Pt*xn(:,i)*xn(:,i)';
        end
        HE=sum3*inv(sum4);
        for ii=1:K
            A(:,ii)=HE(:,ii)/norm(HE(:,ii),2)^2;
        end
       
    case 3  %ZF
        for i=1:N
            b=randi([0 1],2*K,1);
            b=b*2-1;
            xn(:,i)=modulate(b)';
            rn=sqrt(Pt)*H*xn(:,i)+n;
            yn(:,i)=quantizer(rn,M);
            sum3=sum3+sqrt(Pt)*yn(:,i)*xn(:,i)';
            sum4=sum4+Pt*xn(:,i)*xn(:,i)';
        end
        HE=sum3*inv(sum4);
        A=inv(HE'*HE)*HE';
        A=A';
       
    case 4  %MRC,full CSI
        for ii=1:K
            A(:,ii)=H(:,ii)/norm(H(:,ii),2)^2;
        end
       
    case 5  %ZF,full CSI
        A=inv(H'*H)*H';
        A=A';
end
xse=A'*y;
xhe=demodulate(xse)';
for i=1:K
    if(sign(real(x(i)))==1&&sign(imag(x(i)))==1)
        cntx(1)=cntx(1)+1;
    elseif(sign(real(x(i)))==1&&sign(imag(x(i)))==-1)
        cntx(2)=cntx(2)+1;
    elseif(sign(real(x(i)))==-1&&sign(imag(x(i)))==1)
        cntx(3)=cntx(3)+1;
    elseif(sign(real(x(i)))==-1&&sign(imag(x(i)))==-1)
        cntx(4)=cntx(4)+1;
    end
end
for iii=1:4
    px(iii)=cntx(iii)/K;  %p(xk)
end

uk=zeros(1,4);
sk=zeros(1,4);
for k=1:4
    u=zeros(1,M);
    s2=zeros(1,M);
    sc=[1+j,1-j,-1-j,-1+j];
    si=zeros(1,M);
    tempu=0;
    temps=0;
    for i=1:M
        temph=0;
         for jj=1:K
            temph=temph+abs(H(i,jj))*abs(H(i,jj));
        end
        temph=temph-abs(H(i,k))*abs(H(i,k));
        si(i)=Pt*temph;
        de=sqrt(sigma+si(i));
        if (k==1)
            a=0.5*erfc((-sqrt(Pt)*real(H(i,k)*x(k)))/de);
            b=0.5*erfc((-sqrt(Pt)*imag(H(i,k)*x(k)))/de);
        elseif (k==2)
            a=0.5*erfc((-sqrt(Pt)*real(H(i,k)*x(k)))/de);
            b=1-0.5*erfc((-sqrt(Pt)*imag(H(i,k)*x(k)))/de);
        elseif (k==3)
            a=1-0.5*erfc((-sqrt(Pt)*real(H(i,k)*x(k)))/de);
            b=0.5*erfc((-sqrt(Pt)*imag(H(i,k)*x(k)))/de);
        elseif (k==4)
            a=1-0.5*erfc((-sqrt(Pt)*real(H(i,k)*x(k)))/de);
            b=1-0.5*erfc((-sqrt(Pt)*imag(H(i,k)*x(k)))/de);
        end       
        pi1=a*b;
        pi2=a*(1-b);
        pi3=(1-a)*(1-b);
        pi4=(1-a)*b;
        u(i)=H(i,k)'/norm(H(:,k),2)^2*(pi1*sc(1)+pi2*sc(2)+pi3*sc(3)+pi4*sc(4));
        s2(i)=abs(H(i,k))^2/norm(H(:,k),2)^4*(2-(abs(pi1*sc(1)+pi2*sc(2)+pi3*sc(3)+pi4*sc(4)))^2);
        tempu=tempu+u(i);
        temps=temps+s2(i);
    end
    uk(k)=tempu;
    sk(k)=temps;
    syms v;
    f1=1/(sqrt(2*pi*sk(k)))*exp(-(v-real(uk(k)))^2/(2*sk(k)));  %公式12的实部的概率密度函数
    f2=1/(sqrt(2*pi*sk(k)))*exp(-(v-imag(uk(k)))^2/(2*sk(k)));  %公式12的虚部的概率密度函数
    p(4*(k-1)+1)=double(int(f1,v,0,inf)*int(f2,v,0,inf));  %p(xk|xk,H)
    p(4*(k-1)+2)=double(int(f1,v,0,inf)*int(f2,v,-inf,0));
    p(4*(k-1)+3)=double(int(f1,v,-inf,0)*int(f2,v,0,inf));
    p(4*(k-1)+4)=double(int(f1,v,-inf,0)*int(f2,v,-inf,0));
end

end

