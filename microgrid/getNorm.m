function [ ret ] = getNorm( e,sigma )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

ret=zeros(5,1);
idx=zeros(4,1);
normidx=[-1.28,-0.52,0.52,1.28];

for i=1:4
    idx(i)=normidx(i)*sigma+e;
end

ret(1)=-sigma/sqrt(2*pi)*exp(-(idx(1)-e)^2/(2*sigma^2))+e*0.1;
ret(2)=-sigma/sqrt(2*pi)*exp(-(idx(2)-e)^2/(2*sigma^2))+sigma/sqrt(2*pi)*exp(-(idx(1)-e)^2/(2*sigma^2))+e*0.2;
ret(3)=-sigma/sqrt(2*pi)*exp(-(idx(3)-e)^2/(2*sigma^2))+sigma/sqrt(2*pi)*exp(-(idx(2)-e)^2/(2*sigma^2))+e*0.4;
ret(4)=-sigma/sqrt(2*pi)*exp(-(idx(4)-e)^2/(2*sigma^2))+sigma/sqrt(2*pi)*exp(-(idx(3)-e)^2/(2*sigma^2))+e*0.2;
ret(5)=sigma/sqrt(2*pi)*exp(-(idx(4)-e)^2/(2*sigma^2))+e*0.1;

ret(1)=ret(1)/0.1;
ret(2)=ret(2)/0.2;
ret(3)=ret(3)/0.4;
ret(4)=ret(4)/0.2;
ret(5)=ret(5)/0.1;
end

