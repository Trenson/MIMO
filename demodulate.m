function [d]=demodulate(bit)
for n=1:2:length(bit)*2-1
    p=real(bit((n+1)/2));
    imp=imag(bit((n+1)/2));
    if   p>=0
        if imp>=0
        d(n)=1; 
        d(n+1)=1; 
        elseif imp<0
        d(n)=1; 
        d(n+1)=-1; 
        end
    elseif p<0
        if imp>=0
        d(n)=-1; 
        d(n+1)=1; 
        elseif imp<0
        d(n)=-1; 
        d(n+1)=-1;
        end                                       
    end
end
end