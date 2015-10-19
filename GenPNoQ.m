function [ p,px ] = GenPNoQ(SNR,K,M,N,method)
%K=20;
%M=32;
%N=40;
iteration = 100;
tempP=zeros(16);
tempPx=zeros(4);
for iter = 1:iteration
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
    y=r;  %没有量化
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
             A=A';
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
            if(xhe(i*2-1)==1&&xhe(i*2)==1)
                cnt(1)=cnt(1)+1;
            elseif(xhe(i*2-1)==1&&xhe(i*2)==-1)
                cnt(2)=cnt(2)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==1)
                cnt(3)=cnt(3)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==-1)
                cnt(4)=cnt(4)+1;
            end
        elseif(sign(real(x(i)))==1&&sign(imag(x(i)))==-1)
            cntx(2)=cntx(2)+1;
            if(xhe(i*2-1)==1&&xhe(i*2)==1)
                cnt(5)=cnt(5)+1;
            elseif(xhe(i*2-1)==1&&xhe(i*2)==-1)
                cnt(6)=cnt(6)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==1)
                cnt(7)=cnt(7)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==-1)
                cnt(8)=cnt(8)+1;
            end
        elseif(sign(real(x(i)))==-1&&sign(imag(x(i)))==1)
            cntx(3)=cntx(3)+1;
            if(xhe(i*2-1)==1&&xhe(i*2)==1)
                cnt(9)=cnt(9)+1;
            elseif(xhe(i*2-1)==1&&xhe(i*2)==-1)
                cnt(10)=cnt(10)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==1)
                cnt(11)=cnt(11)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==-1)
                cnt(12)=cnt(12)+1;
            end
        elseif(sign(real(x(i)))==-1&&sign(imag(x(i)))==-1)
            cntx(4)=cntx(4)+1;
            if(xhe(i*2-1)==1&&xhe(i*2)==1)
                cnt(13)=cnt(13)+1;
            elseif(xhe(i*2-1)==1&&xhe(i*2)==-1)
                cnt(14)=cnt(14)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==1)
                cnt(15)=cnt(15)+1;
            elseif(xhe(i*2-1)==-1&&xhe(i*2)==-1)
                cnt(16)=cnt(16)+1;
            end
        end
    end
    %cntx
    %cnt
    for ii=1:16
        tempP(ii)=tempP(ii)+cnt(ii)/cntx(fix(ii/4.1)+1);  %p(xk|xk,H)
    end
    for iii=1:4
        tempPx(iii)=tempPx(iii)+cntx(iii)/K;  %p(xk)
    end
end
    for ii=1:16
        p(ii)=tempP(ii)/iteration;
    end
    for iii=1:4
        px(iii)=tempPx(iii)/iteration;
    end
end

