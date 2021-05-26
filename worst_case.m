%MIT License
%Copyright (c) [2021] [Shuang Qin]
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
function x = worst_case(X,U,t,bias)
% X is the sensor coordinates vector, e.g., [-20,-20;20,-20;20,20;-20,20]
% for 4 sensors.
% t is the measurement of each sensor
% bias is pre-set
% U is Z or PZ in proposition 2
xr = [];
txr = [];
for i =1:size(U, 1)
    l = find(U(i,:) ==  1);
    m = find(U(i,:) == -1);
    if t(i) >= 0
        xr = [xr,[m;l]];
        txr = [txr,t(i)];
    else
        xr = [xr,[l;m]]; 
        txr = [txr,-t(i)];
    end
end
xc = unique(xr);
cvx_begin sdp quiet
variable y(2+length(xc)+1)
variable Y(2+length(xc)+1,2+length(xc)+1)  hermitian 
variable lam(length(xr))
minimize (sum(lam))
for m =1:length(xr) 
    j = find(xc==xr(1,m));
    r = txr(m);
    si = X(xr(2,m)+1,:)';
    sj = X(xr(1,m)+1,:)';
    c1 = bias^2+2*bias*r+r^2+norm(sj)^2-norm(si)^2;
    c2 = -bias^2+2*bias*r-r^2-norm(sj)^2+norm(si)^2; 
    subject to
    B1=zeros(2+length(xc)+1,2+length(xc)+1);%(18b)
    B1(1:2,1:2)=4*(si-sj)*(si-sj)';
    B1(1:2,2+j)=4*(bias+r)*(si-sj);
    B1(1:2,end)=2*c1*(si-sj);
    B1(2+j,end)=2*c1*(bias+r);
    B1(2+j,1:2)=4*(bias+r)*(si-sj)';
    B1(end,1:2)=2*c1*(si-sj)';
    B1(end,2+j)=2*c1*(bias+r);
    B1(2+j,2+j)=4*(bias+r)^2;
    B1(end,end)=c1^2;
    trace(B1*Y) <= lam(m);
    B2=zeros(2+length(xc)+1,2+length(xc)+1);%(18c)
    B2(1:2,1:2)=4*(si-sj)*(si-sj)';
    B2(1:2,2+j)=-4*(bias-r)*(si-sj);
    B2(1:2,end)=-2*c2*(si-sj);
    B2(2+j,end)=2*c2*(bias-r);
    B2(2+j,1:2)=-4*(bias-r)*(si-sj)';
    B2(end,1:2)=-2*c2*(si-sj)';
    B2(end,2+j)=2*c2*(bias-r);
    B2(2+j,2+j)=4*(bias-r)^2;
    B2(end,end)=c2^2;
    trace(B2*Y) <= lam(m);
    B=zeros(2+length(xc)+1,2+length(xc)+1);%(18a)
    B(1:2,1:2)=eye(2);
    B(1:2,end)=-sj;
    B(end,1:2)=-sj';
    B(end,end)=sj'*sj;
    B(2+j,2+j)=-1;
    trace(B*Y)==0;
    norm(y(1:2)-sj) <= y(2+j);
end
[Y,y;y',1] >= 0;
y(end) == 1;
Y(end,end) == 1;
cvx_end
x = y(1:2);
end