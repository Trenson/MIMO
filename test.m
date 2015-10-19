clear; close all; clc;
K=20;
M=20;
N=100;
MI=zeros(1,2);
m=1;
method=2;
s=1;
t=10;  %公式10求均值的循环次数
SNR=-5;
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
                I=I+p(j)*px(fix(j/4.1)+1)*log(p(j)/pxk(j-4*jj))/log(2)  %pxk(j-4*jj)可能为0
            end
            if(mod(j,4)==0)
                jj=jj+1;
            end
        end
        temp=temp+I;
    end
    MI(s)=temp/t  % m表示图2的某条线
    s=s+1  % s表示横轴的第几个值