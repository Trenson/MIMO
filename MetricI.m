clear; close all; clc;
K=20;
M=20;
N=100;
MI=zeros(8,19);
m=1;
t=10;  %公式10求均值的循环次数
for method=4:5
    s=1;
for SNR=-20:5/3:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px,A,y ]=GenP(SNR,K,M,N,method);
        for k=1:4
            for kk=1:4
                pxk(k)=pxk(k)+p(4*(kk-1)+k)*px(kk);
            end
        end
        for j=1:16
            if(p(j)>0&&pxk(j-4*jj)>0)
                I=I+p(j)*px(fix(j/4.1)+1)*log(p(j)/pxk(j-4*jj))/log(2);  %pxk(j-4*jj)可能为0
            end
            if(mod(j,4)==0)
                jj=jj+1;
            end
        end
        temp=temp+I;
    end
    MI(m,s)=temp/t;  % m表示图2的某条线
    s=s+1;  % s表示横轴的第几个值
end
    m=m+1;
end

for method=1:3
    s=1;
for SNR=-20:5/3:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px,A,y ]=GenP(SNR,K,M,N,method);
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
end
%设置N=50M
for method=1:3
    s=1;
for SNR=-20:5/3:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px,A,y ]=GenP(SNR,K,M,N*10,method);
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
end

MI = MI*2; % 08.06加
save MI MI;
%画图程序
hf = figure;
set( hf, 'color', 'white');
SNR=-20:5/3:10;
plot( SNR, MI(1,:), '-bo','LineWidth',1.5);
hold on;
plot( SNR, MI(2,:), '-rs','LineWidth',1.5);
hold on;
plot( SNR, MI(3,:), '-.k*','LineWidth',1.5);
hold on;
plot( SNR, MI(4,:), '-.bo','LineWidth',1.5);
hold on;
plot( SNR, MI(5,:), '-.rs','LineWidth',1.5);
hold on;
plot( SNR, MI(6,:), '--k*','LineWidth',1.5);
hold on;
plot( SNR, MI(7,:), '--bo','LineWidth',1.5);
hold on;
plot( SNR, MI(8,:), '--rs','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
grid on;
legend('MRC,full CSI','ZF,full CSI','LS,N=5M','MRC,channel estimation,N=5M','ZF,channel estimation,N=5M','LS,N=50M','MRC,channel estimation,N=50M','ZF,channel estimation,N=50M')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}Mutual information per user (bit/channel use)');