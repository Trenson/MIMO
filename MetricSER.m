clear; close all; clc;
K=20;
M=400;
N=100;
t=10; %公式23求均值的循环次数
SER=zeros(8,11);
m=1;  %表示不同的线
for method=4:5
    s=1;
    for SNR=-20:2:0
        temp=0;
        for h=1:t
            [p,px]=GenP(SNR,K,M,N,method);
            temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
        end
        SER(m,s)=temp/t;
        s=s+1;
    end
    m=m+1;
end

%N=50K
for method=2:3
    s=1;
    for SNR=-20:2:0
        temp=0;
        for h=1:t
            [p,px]=GenP(SNR,K,M,10*N,method);
            temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
        end
        SER(m,s)=temp/t;
        s=s+1;
    end
    m=m+1;
end

%N=10K
for method=2:3
    s=1;
    for SNR=-20:2:0
        temp=0;
        for h=1:t
            [p,px]=GenP(SNR,K,M,2*N,method);
            temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
        end
        SER(m,s)=temp/t;
        s=s+1;
    end
    m=m+1;
end

%N=5K
for method=2:3
    s=1;
    for SNR=-20:2:0
        temp=0;
        for h=1:t
            [p,px]=GenP(SNR,K,M,N,method);
            temp=temp+px(1)*(p(2)+p(3)+p(4))+px(2)*(p(5)+p(7)+p(8))+px(3)*(p(9)+p(10)+p(12))+px(4)*(p(13)+p(14)+p(15));
        end
        SER(m,s)=temp/t;
        s=s+1;
    end
    m=m+1;
end

%画图程序
hf = figure;
set( hf, 'color', 'white');
SNR=-20:2:0;
semilogy( SNR, SER(1,:), '-b^','LineWidth',1.5);
hold on;
semilogy( SNR, SER(2,:), '-r^','LineWidth',1.5);
hold on;
semilogy( SNR, SER(3,:), '-bo','LineWidth',1.5);
hold on;
semilogy( SNR, SER(4,:), '-rs','LineWidth',1.5);
hold on;
semilogy( SNR, SER(5,:), '-.bo','LineWidth',1.5);
hold on;
semilogy( SNR, SER(6,:), '-.rs','LineWidth',1.5);
hold on;
semilogy( SNR, SER(7,:), '--bo','LineWidth',1.5);
hold on;
semilogy( SNR, SER(8,:), '--rs','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
grid on;
legend('MRC,full CSI','ZF,full CSI','MRC,channel estimation,N=50K','ZF,channel estimation,N=50K','MRC,channel estimation,N=10K','ZF,channel estimation,N=10K','MRC,channel estimation,N=5K','ZF,channel estimation,N=5K')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}SER per user');