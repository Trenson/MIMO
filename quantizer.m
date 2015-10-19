function [ q ] = quantizer( v,M )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i=1:M
    if (real(v(i))>0&&imag(v(i))>0)
        q(i)=1+j;
    elseif (real(v(i))>0&&imag(v(i))<0)
        q(i)=1-j;
    elseif (real(v(i))<0&&imag(v(i))>0)
        q(i)=-1+j;
    elseif (real(v(i))<0&&imag(v(i))<0)
        q(i)=-1-j;
    elseif (real(v(i))==0&&imag(v(i))>0)
        q(i)=j;
    elseif (real(v(i))==0&&imag(v(i))<0)
        q(i)=-j;
    elseif (real(v(i))==0&&imag(v(i))==0)
        q(i)=0;
    elseif (real(v(i))>0&&imag(v(i))==0)
        q(i)=1;
    elseif (real(v(i))<0&&imag(v(i))==0)
        q(i)=-1;
end
q=q';
end

