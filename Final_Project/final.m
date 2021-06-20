%% Task1. Implement Metropolis-Hastings
clear all;
rng(0);
 
p = @(x) exp(-x)/(1+exp(-x))^2*(x>=0);
n = 10000; 
s = zeros(n,1); %Markov Chain
 
s(1) = 1; % Initialize MC
% Metropolis acceptane ratio 'a', either A(x',x)=1 or A(x,x')=1 and either
% way, the condition is satisfied
 
for i = 2 : n
    xprime = s(i-1)+1*randn; % Generate candidate
    alpha = min(1,p(xprime)/p(s(i-1)));
    if rand < alpha % f(x')/f(x) = p(x')=p(x) because f(x), normal distribution is proportional to target
        s(i) = xprime; % accept
    else
        s(i) = s(i-1); % reject
    end
end
del = 0.02; smin = min(s); smax = max(s);
bincenters = smin:del:smax;
prb_est = zeros(1,length(bincenters));
for i=1:length(bincenters)
    prb_est(i) = nnz(s > bincenters(i)-del/2 & s <= bincenters(i)+del/2)...
    /n/del;
end
plot(bincenters, prb_est, 'g'); hold on;
% Comparison with the estimated PDF with true PDF
figure(1)
x=linspace(smin,smax,1000);
plot(x,2.*exp(-x)./(1+exp(-x)).^2,'r'); 
legend('Metropolis-Hastings', 'True PDF'); hold off;

%% Task2. estimate and plot the PDF of X[n] for n = 1,3,10,30,100
s2 = zeros(5,n); 
for i = 1:5
    s2(i,1) = 4;
end
for i = 1:n
    for j = 2:101
        xprime = s2(j-1)+randn;
        if rand < p(xprime)/p(s2(j-1)) 
                s2(j) = xprime;
        else
                s2(j) = s2(j-1);
        end
    end
    s2(1,i) = s2(2);  
    s2(2,i) = s2(4);
    s2(3,i) = s2(11); 
    s2(4,i) = s2(31);
    s2(5,i) = s2(101);
end
 
del = 0.02;
bincenters = 0:del:100;
prb_est = zeros(5,length(bincenters));
 
for i = 1 : 5
    for j= 1 : length(bincenters)
        prb_est(i,j) = nnz(s2(i,:)>bincenters(j)-del/2 & s2(i,:)<bincenters(j)+del/2)/n/del;
    end
end
figure(2)
x=linspace(0,5,100);
plot(x, 2.*exp(-x)./(1+exp(-x)).^2, bincenters, prb_est(1,:), ...
    bincenters,prb_est(2,:), bincenters, prb_est(3,:), ...
    bincenters, prb_est(4,:), bincenters, prb_est(5,:))
axis([0,10,0,1])
  
  
%% Task3. Repeat Task2 for the other conditional distributions g:
s2 = zeros(5,n); 
for i = 1:5
    s2(i,1) = 0;
end
tic
for i = 1:n
    for j = 2:101
        xprime = s2(j-1)+10*randn;
        if rand < p(xprime)/p(s2(j-1)) 
                s2(j) = xprime;
        else
                s2(j) = s2(j-1);
        end
    end
    s2(1,i) = s2(2);  
    s2(2,i) = s2(4);
    s2(3,i) = s2(11); 
    s2(4,i) = s2(31);
    s2(5,i) = s2(101);
end
toc
del = 0.02;
bincenters = 0:del:100;
prb_est = zeros(5,length(bincenters));
 
for i = 1 : 5
    for j= 1 : length(bincenters)
        prb_est(i,j) = nnz(s2(i,:)>bincenters(j)-del/2 & s2(i,:)<bincenters(j)+del/2)/n/del;
    end
end
figure(3)
x=linspace(0,5,100);
plot(x, 2.*exp(-x)./(1+exp(-x)).^2, bincenters, prb_est(1,:), ...
    bincenters,prb_est(2,:), bincenters, prb_est(3,:), ...
    bincenters, prb_est(4,:), bincenters, prb_est(5,:))
axis([0,10,0,1])
 
 
s2 = zeros(5,n); 
for i = 1:5
    s2(i,1) = 0;
end
tic
for i = 1:n
    for j = 2:101
        xprime = s2(j-1)+0.1*randn;
        if rand < p(xprime)/p(s2(j-1)) 
                s2(j) = xprime;
        else
                s2(j) = s2(j-1);
        end
    end
    s2(1,i) = s2(2);  
    s2(2,i) = s2(4);
    s2(3,i) = s2(11); 
    s2(4,i) = s2(31);
    s2(5,i) = s2(101);
end
toc
del = 0.02;
bincenters = 0:del:100;
prb_est = zeros(5,length(bincenters));
 
for i = 1 : 5
    for j= 1 : length(bincenters)
        prb_est(i,j) = nnz(s2(i,:)>bincenters(j)-del/2 & s2(i,:)<bincenters(j)+del/2)/n/del;
    end
end
figure(4)
x=linspace(0,5,100);
plot(x, 2.*exp(-x)./(1+exp(-x)).^2, bincenters, prb_est(1,:), ...
    bincenters,prb_est(2,:), bincenters, prb_est(3,:), ...
    bincenters, prb_est(4,:), bincenters, prb_est(5,:))
axis([0,10,0,1])
 
 
s2 = zeros(5,n); 
for i = 1:5
    s2(i,1) = 0;
end
tic
for i = 1:n
    for j = 2:101
        xprime = s2(j-1)+(rand-0.5);
        if rand < p(xprime)/p(s2(j-1)) 
                s2(j) = xprime;
        else
                s2(j) = s2(j-1);
        end
    end
    s2(1,i) = s2(2);  
    s2(2,i) = s2(4);
    s2(3,i) = s2(11); 
    s2(4,i) = s2(31);
    s2(5,i) = s2(101);
end
toc
del = 0.02;
bincenters = 0:del:100;
prb_est = zeros(5,length(bincenters));
 
for i = 1 : 5
    for j= 1 : length(bincenters)
        prb_est(i,j) = nnz(s2(i,:)>bincenters(j)-del/2 & s2(i,:)<bincenters(j)+del/2)/n/del;
    end
end
figure(5)
x=linspace(0,5,100);
plot(x, 2.*exp(-x)./(1+exp(-x)).^2, bincenters, prb_est(1,:), ...
    bincenters,prb_est(2,:), bincenters, prb_est(3,:), ...
    bincenters, prb_est(4,:), bincenters, prb_est(5,:))
axis([0,10,0,1])

%% Task4. Choose an image
img = imread('BangWool.jpg'); 
figure(6), image(img);


clear all;
rng(0);
%% Task5. 
% The color of palette consists of 'k' colors, each of which consists of 3 numbers 
%representing the red, green, and blue channels.
% (i.e.) The palette consists of 3k random variables denoted P_ij
%where i=1,...,k and j=1,2,3.
% Each of palette numbers should be between 0 and 255
 
% The color index for each pixel: Z_ij for the pixel position (i,j),
% where i=1,...,m and j=1,...,n assuming m*n is resolution of the image
% Each of the index numbers Z_ij={1,2,...,k} indicates which color from the
% color palette this pixel takes on
 
p = 255*rand(256,3);
z = randi([1 4],300,300);
X = zeros(300,300,3);
for i=1:300
    for j=1:300
        X(i,j,1)=p(z(i,j),1);    
        X(i,j,2)=p(z(i,j),2);    
        X(i,j,3)=p(z(i,j),3);    
    end
end
test=uint8(X);
imwrite(test,'test.jpg'); 
figure(7), image(test);
title('Mosaic')
r=0.5;
c_=1; % 'c' does mean nothing here since we care about the ratio
pz = @(x) c_*exp(r*x);
 
 
 
for a=1:100
    count = zeros(300,300); % Direc Delta Counting
    i=1; j=1;
    if (z(i,j)==z(i,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j))
        count(i,j) = count(i,j) + 1;
    end
    i=300; j=1;
    if (z(i,j)==z(i,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j))
        count(i,j) = count(i,j) + 1;
    end
    i=1; j=300;
    if (z(i,j)==z(i,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j))
        count(i,j) = count(i,j) + 1;
    end      
    i=300; j=300;
    if (z(i,j)==z(i,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j))
        count(i,j) = count(i,j) + 1;
    end
    
    i=1;j=2;
    while(j<300)
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        j = j+1;
    end
    i=300;j=2;
    while(j<300)
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        j = j+1;
    end
    
    i=2;j=300;    
    while(i<300)
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        i = i+1;
    end       
    i=2;j=1;    
    while(i<300)
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        i = i+1;
    end      
    for i=2:299
        for j=2:299
            if z(i,j)==z(i+1,j-1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i+1,j)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i+1,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j-1)
                count(i,j)=count(i,j)+1;
            end            
            if z(i,j)==z(i,j-1)
                count(i,j)=count(i,j)+1;
            end            
        end
    end
    i = randi(300);
    j = randi(300);
    for k = 1 : 600000
        a=randi(3);
        b=randi(3);
        if a == 1
            xprime = i+1;
        elseif a == 2
            xprime = i-1;
        elseif a ==3
            xprime = i;
        end
        if b == 1
            yprime = j+1;
        elseif b == 2
            yprime = j-1;
        elseif b == 3
            yprime = j;
        end
        if (xprime>1 && xprime<300) && (yprime>1 && yprime<300)
            alpha = min(1, pz(count(xprime,yprime))/pz(count(i,j)));
                    
            if rand < alpha
                z(i,j) = z(xprime, yprime);
                i = xprime;
                j = yprime;
            else
                z(xprime,yprime) = z(i,j);
                i = randi(300);
                j = randi(300);
            end
        end
    end
end
X_2=zeros(300,300,3);
for i=1:300
    for j=1:300
        X_2(i,j,1)=p(z(i,j),1);
        X_2(i,j,2)=p(z(i,j),2);    
        X_2(i,j,3)=p(z(i,j),3);   
    end
end
 
test2=uint8(X_2);
imwrite(test2,'test2.jpg'); 
figure(8), image(test2);
title('Sampled')

  
clear all;
rng(0);
z = randi([1 8],300,300); 
p=[234 220 85 ; 132 163 207 ; 171 161 149 ; 248 250 117  ; ...
    158 193 225 ; 73 101 138 ; 131 157 194 ; 172 255 0];
 
 
C = zeros(300,300,3);
for i=1:300
    for j=1:300
        C(i,j,1) = p(z(i,j),1);
        C(i,j,2) = p(z(i,j),2);
        C(i,j,3) = p(z(i,j),3);
    end
end
init = uint8(C);
figure(1); image(init);
 
A = imread('sample.jpg');
image(A); title('Original')
A = double(A);
 
 
% Definition of Pz(z) and Ppz(p,z)
c = 1; r = 1; s = 0.5;
Pz = @(x) c*exp(r*x);
Ppz = @(x,y) c*(255*rand)*Pz(x)*exp((-s)*y);
 
 
S=zeros(300,300,3);
 
for k=1:100
 
     
    % Direc-Delta Function
    count = zeros(300,300);
    i=1; j=1;
    if (z(i,j)==z(i,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j))
        count(i,j) = count(i,j) + 1;
    end
    i=300; j=1;
    if (z(i,j)==z(i,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j+1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j))
        count(i,j) = count(i,j) + 1;
    end
    i=1; j=300;
    if (z(i,j)==z(i,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i+1,j))
        count(i,j) = count(i,j) + 1;
    end      
    i=300; j=300;
    if (z(i,j)==z(i,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j-1))
        count(i,j) = count(i,j) + 1;
    elseif (z(i,j)==z(i-1,j))
        count(i,j) = count(i,j) + 1;
    end
    
    i=1;j=2;
    while(j<300)
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        j = j+1;
    end
    i=300;j=2;
    while(j<300)
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        j = j+1;
    end
    
    i=2;j=300;    
    while(i<300)
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j-1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        i = i+1;
    end       
    i=2;j=1;    
    while(i<300)
        if z(i,j)==z(i-1,j)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i-1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j+1)
            count(i,j)=count(i,j)+1;
        end
        if z(i,j)==z(i+1,j)
            count(i,j)=count(i,j)+1;
        end
        i = i+1;
    end      
    for i=2:299
        for j=2:299
            if z(i,j)==z(i+1,j-1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i+1,j)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i+1,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j+1)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j)
                count(i,j)=count(i,j)+1;
            end
            if z(i,j)==z(i-1,j-1)
                count(i,j)=count(i,j)+1;
            end            
            if z(i,j)==z(i,j-1)
                count(i,j)=count(i,j)+1;
            end            
        end
    end
    % abs(Cijl-Aijl)
    S = zeros(300,300);
    for i=1:300
        for j=1:300
            for l=1:3
                S(i,j) = S(i,j)+abs(C(i,j,l)-A(i,j,l));
            end
        end
    end
    
    %% Metropolis-Hastings
    i = randi(300);
    j = randi(300);
    for m=1:50000
        a = randi(3);
        b = randi(3);
        if a == 1
            iprime = i+1;
        elseif a == 2
            iprime = i-1;
        elseif a == 3
            iprime = i;
        end
        if b == 1
            jprime = j+1;
        elseif b == 2
            jprime = j-1;
        elseif b == 3
            jprime = j;
        end
        
        if (iprime>1 && iprime<300) && (jprime>1 && jprime<300)
            alpha = min(1,Ppz(count(iprime,jprime),S(iprime,jprime))/...
                Ppz(count(i,j),S(i,j)));
            if rand < alpha
                S(i,j) = S(iprime,jprime);
                z(i,j) = z(iprime,jprime);
                C(i,j,:)=C(iprime,jprime,:);
                i=iprime;
                j=jprime;
            else
                S(iprime,jprime)=S(i,j);
                z(iprime,jprime)=z(i,j);
                C(iprime,jprime,:)=C(i,j,:);
                i=randi(300);
                j=randi(300);
            end
        end
    end
end
result = uint8(C);
figure(2); image(result); title('Sampled')
  
