clear; close all; clc;
K=20;
M=400;
N=1000;
SER=zeros(4,10);
m=1;
t=10;  %公式23求均值的循环次数
s=1;
for SNR=-20:1:-11
    temp=0;
    for h=1:t
    [p,px]=GenPNoQ(SNR,K,M,N,5);
    temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
    end
    SER(m,s)=temp/t;
    s=s+1;
end
m=m+1;

s=1;
for SNR=-20:1:-11
    temp=0;
    for h=1:t
    [p,px]=GenP(SNR,K,M,N,5);
    temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
    end
    SER(m,s)=temp/t;
    s=s+1;
end
m=m+1;

s=1;
for SNR=-20:1:-11
    temp=0;
    for h=1:t
    [p,px]=GenPNoQ(SNR,K,M,N,3);
    temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
    end
    SER(m,s)=temp/t;
    s=s+1;
end
m=m+1;

s=1;
for SNR=-20:1:-11
    temp=0;
    for h=1:t
    [p,px]=GenP(SNR,K,M,N,3);
    temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
    end
    SER(m,s)=temp/t;
    s=s+1;
end


%画图程序
hf = figure;
set( hf, 'color', 'white');
SNR=-20:1:-11;
semilogy( SNR, SER(1,:), '-rs','LineWidth',1.5);
hold on;
semilogy( SNR, SER(2,:), '--rs','LineWidth',1.5);
hold on;
semilogy( SNR, SER(3,:), '-bo','LineWidth',1.5);
hold on;
semilogy( SNR, SER(4,:), '--bo','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
grid on;
legend('Full CSI,no quantizer','Full CSI,quantizer','Channel estimation,no quantizer','Channel estimation,quantizer')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}SER per user');