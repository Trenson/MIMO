clear; close all; clc;
K=20;
M=400;
N=100;
MI=zeros(8,19);
m=1;
t=1;  %公式10求均值的循环次数
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
                I=I+p(j)*px(fix(j/4.1)+1)*log(p(j)/pxk(j-4*jj))/log(2);
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

%N=50K
for method=2:3
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

%N=10K
for method=2:3
    s=1;
for SNR=-20:5/3:10
    temp=0;
    for h=1:t
        I=0;
        jj=0;
        pxk=zeros(1,4);
        [p,px,A,y ]=GenP(SNR,K,M,N*2,method);
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

%N=5K
for method=2:3
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

%MI = MI*2; % 08.06加
%画图程序
hf = figure;
set( hf, 'color', 'white');
SNR=-20:5/3:10;
plot( SNR, MI(1,:), '-bo','LineWidth',1.5);
hold on;
plot( SNR, MI(2,:), '-rs','LineWidth',1.5);
hold on;
plot( SNR, MI(3,:), '--bo','LineWidth',1.5);
hold on;
plot( SNR, MI(4,:), '--rs','LineWidth',1.5);
hold on;
plot( SNR, MI(5,:), ':bo','LineWidth',1.5);
hold on;
plot( SNR, MI(6,:), ':rs','LineWidth',1.5);
hold on;
plot( SNR, MI(7,:), '-.bo','LineWidth',1.5);
hold on;
plot( SNR, MI(8,:), '-.rs','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
grid on;
legend('MRC,full CSI','ZF,full CSI','MRC,channel estimation,N=50K','ZF,channel estimation,N=50K','MRC,channel estimation,N=10K','ZF,channel estimation,N=10K','MRC,channel estimation,N=5K','ZF,channel estimation,N=5K')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}Mutual information per user (bit/channel use)');