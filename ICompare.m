clear; close all; clc;
K=20;
M=400;
N=2000;
MI=zeros(2,7);
m=1;
t=10;  %公式10求均值的循环次数
method=4;

s=1;
for SNR=-20:5:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px]=GenP(SNR,K,M,N,method);
        for k=1:4
            for kk=1:4
                pxk(k)=pxk(k)+p(4*(kk-1)+k)*px(kk);
            end
        end
        for j=1:16
            if(p(j)>0)
                I=I+p(j)*px(fix(j/4.1)+1)*log(p(j)/pxk(j-4*jj))/log(2);
            end
            if(mod(j,4)==0)
                jj=jj+1;
            end
        end
        temp=temp+I;
    end
    MI(m,s)=temp/t;
    s=s+1;
end
m=m+1;

s=1;
for SNR=-20:5:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px]=Analytical(SNR,K,M,N,method);
        for k=1:4
            for kk=1:4
                pxk(k)=pxk(k)+p(4*(kk-1)+k)*px(kk);
            end
        end
        for j=1:16
            if(p(j)>0)
                I=I+p(j)*px(fix(j/4.1)+1)*log(p(j)/pxk(j-4*jj))/log(2);
            end
            if(mod(j,4)==0)
                jj=jj+1;
            end
        end
        temp=temp+I;
    end
    MI(m,s)=temp/t;
    s=s+1;
end

%画图程序
hf = figure;
set( hf, 'color', 'white');
SNR=-20:5:10;
plot( SNR, MI(1,:), '-rs','LineWidth',1.5);
hold on;
plot( SNR, MI(2,:), '-bo','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
grid on;
legend('numerical','analytical')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}Mutual information per user (bit/channel use)');