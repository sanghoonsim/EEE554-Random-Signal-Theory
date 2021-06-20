clear all;
rng(0);
% Task 1
syms x   % definite integrate of g(x)
x_1 = 0; % given initial value
ya = (exp(-x))/((1+exp(-x))^2); % g(x)
integral = int(ya); % int(y,x1,Inf) = 0.5, thus c = 2
c = 2; % since fs = c*g(x) = 1
fs_function = 2*(exp(-x))/((1+exp(-x))^2);
fs = @(x) 2*(exp(-x))/((1+exp(-x))^2); % PDF
% Task 2
Fs = @(x) (1-exp(-x))/(1+exp(-x)); % CDF
 
% Task 3, 4
nsamples = 1000000;      
s = zeros(nsamples,1); 
u = rand(nsamples,1); 
tol = 10^(-12);
tic
for i = 1 : nsamples
    xmin = x_1;
    xmax = 20;
    while true
        x = (xmax+xmin)/2;
        if xmax-xmin < tol
            break;
        elseif Fs(x) > u(i)
            xmax = x;
        else
            xmin = x;
        end
    end
    s(i) = x;
end
toc
 
del = 0.02; xmin = min(s); xmax = max(s);
bincenters = xmin:del:xmax;
pdf_est = zeros(1,length(bincenters));
for i=1:length(bincenters)
    pdf_est(i) = nnz(s > bincenters(i)-del/2 & s ...
    <= bincenters(i)+del/2)/nsamples/del;
end
plot(bincenters, pdf_est, 'g'); hold on;
% Comparison with the estimated PDF with true PDF
figure(1)
x=linspace(xmin,xmax,100);
plot(x,2.*exp(-x)./(1+exp(-x)).^2,'r'); 
legend('Estimated PDF', 'True PDF'); hold off;
 
 
% Task 5
anum = 4;   % the number of A
xa = zeros(anum,1);   % an array of storing x-values
Amax = 50;  % Initial value of max(A) 
Amin = 0;   % Initial value of min(A)
while true   % Outer bisection 
    Area = (Amax+Amin)/2;
    for i = 2 : anum   % Inner bisection
        xmin = 0;   % Initial value of min(x)
        xmax = 10;   % Initial value of max(x)
        while true   
            xa(i) = (xmax+xmin)/2;
            if (xmax-xmin)<tol
                break;
            % comparing the area of rectangle with guessed A
            elseif xa(i)*(fs(xa(i-1))-fs(xa(i))) > Area
                xmax = xa(i);
            else
                xmin = xa(i);
            end
        end
    end
    if (Amax-Amin)<tol
        break;
    % Comparing the nth area with guessed A
    elseif xa(anum)*fs(xa(anum))+1-Fs(xa(anum)) < Area
        Amax = Area;
    else
        Amin = Area;
    end
end
 
ya = zeros(anum,1);
for i = 1 : anum
    ya(i) = fs(xa(i));
end
va = zeros(anum,1);
for i = 1 : anum-1
    va(i) = xa(i+1)*(ya(i)-ya(i+1));
end
va(anum) = xa(anum)*(ya(anum))+1-Fs(xa(anum));
 
% Task 6
% Area of Rectangular in base / Area
p = xa(anum)*ya(anum)/va(anum);   
 
 
% Task 7
t = 10000;
ta = zeros(t,1);
c_prime = (exp(261/80)+1); % c' when there are 4 areas
lam = 0.5;
M = 1;
g = 1-exp(-lam*xa(anum));
for i=1:t
        while true
            tx = -1/lam*log((exp(-lam*xa(anum)))*rand);                   
            ty = M*(lam*exp(-lam*tx))*rand;
            if ty < 2*exp(-tx)/(1+exp(-tx))^2
                ta(i) = tx;
            break;
            end
        end
end            
del = 0.05; tmin = min(ta); tmax = max(ta);
bincenters = tmin:del:tmax;
pdf_est = zeros(1,length(bincenters));
for i=1:length(bincenters)
    pdf_est(i) = nnz(ta > bincenters(i)-del/2 & ta ...
    <= bincenters(i)+del/2)/t/del;
end
figure(2)
plot(bincenters, pdf_est, 'g'); hold on;
x=linspace(tmin,tmax,100);
plot(x,c_prime.*exp(-x)./(1+exp(-x)).^2,'r'); hold off;
 
% Task 8, 9, 10
z = 1000000;
zsamples = zeros(z,1);
tracking_a = 0;
tracking_b = 0;
tracking_c = 0;
tracking_d = 0;
tracking_e = 0;
tic
for i = 1 : z
    if rand <= (1-p)/anum % do first, because 
        while true
            zx = -1/lam*log((exp(-lam*xa(anum))*rand));                   
            zy = M*lam*exp(-lam*zx)*rand;
            if zy < 2*exp(-zx)/(1+exp(-zx))^2
                zsamples(i) = zx;
                tracking_e = tracking_e + 1;
                break;
            end
        end
    elseif rand <= p/anum
        zsamples(i) = rand*(xa(anum));
        tracking_d = tracking_d + 1;
    else
        ai = randi([1, anum-1]);
        zx = rand*(xa(ai+1));
        if ai == 1
           zx = rand*(xa(2));
           zy = rand*(ya(1)-ya(2))+ya(2);
           if zy < fs(zx)
               zsamples(i) = zx;
               tracking_b = tracking_b + 1;
           else
              tracking_c = tracking_c + 1;
           end
        elseif zx >= xa(ai)
           zx = rand*(xa(ai+1)-xa(ai))+xa(ai);
           zy = rand*(ya(ai)-ya(ai+1))+(ya(ai+1));
           if zy < fs(zx)
              zsamples(i) = zx;
              tracking_b = tracking_b + 1;
           else
              tracking_c = tracking_c + 1;
           end
        else
            zsamples(i) = zx;
            tracking_a = tracking_a + 1;
        end
    end   
end
toc
zsamples(zsamples==0)=nan;
tracking_sum = tracking_a + tracking_b + tracking_c + ...
tracking_d + tracking_e;
del = 0.01; zmin = min(zsamples); zmax = max(zsamples);
bincenters = zmin:del:zmax;
pdf_est = zeros(1,length(bincenters));

for i=1:length(bincenters)
    pdf_est(i) = nnz(zsamples > bincenters(i)-del/2 & ...
    zsamples <= bincenters(i)+del/2)/ ...
    (tracking_sum-tracking_c)/del;
end
figure(3)
plot(bincenters, pdf_est, 'b'); hold on;
x=linspace(zmin,zmax,100);
plot(x,2.*exp(-x)./(1+exp(-x)).^2,'r'); 
legend("Ziggurat Sampling", "Ture PDF"); hold off;
