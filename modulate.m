function [d]=modulate(bit)
for n=1:length(bit)/2
    imp=bit(2*n);
    p=bit(2*n-1);
    if   p==1
        if imp==1
        d(n)=1/sqrt(2)+j*1/sqrt(2); 
        elseif imp==-1
        d(n)=1/sqrt(2)-j*1/sqrt(2); 
        end
    elseif p==-1
        if imp==1
        d(n)=-1/sqrt(2)+j*1/sqrt(2); 
        elseif imp==-1
        d(n)=-1/sqrt(2)-j*1/sqrt(2); 
        end                                       
    end
end