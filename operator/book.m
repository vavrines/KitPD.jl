%  Madenci et al., (2019) Peridynamic Differential Operator for Numerical Analysis 
%  MATLAB code for the N-th order derivative of a function with M dimensions.
%
clc;
clear;
close all;
echo off;
global icount;
global p;
global ndiv;
global coord;
global dx;
global M;
global N;
%
% Input parameters
%
% M is currently limited to 10.
% If M>10 increase the size of arrays num[] and num1[] everywhere.
M = str2num(input('Specify the number of dimensions= ','s'));
N = str2num(input('Specify the maximum polynomial order in TSE= ','s'));
%
fprintf ('Specify the domain length for each dimension\n')
length = zeros(M,1);
for i=1:M
    fprintf('x%d ',i);
    length(i,1) = str2num(input('= ' ,'s'));
end
%
fprintf ('Specify the number of intervals for each dimension\n')
ndiv = zeros(M,1);
for i=1:M  
    fprintf('ndivx%d ',i);
    ndiv(i,1) = str2num(input('= ' ,'s'));
end
%
fprintf ('Specify the order of derivatives for each dimension\n')
porder = zeros(M,1);
for i=1:M
    fprintf('p%d ',i);
    porder(i,1) = str2num(input('= ' ,'s'));
end
%
% Compute total number of PD points
totnode= prod(ndiv);
% Initialize coordinate array
coord=zeros(totnode,M);
% Compute the number of TS terms based on N
icount=0;
loop1(1,N);
nsize=icount;       % nsize = size of PDDO (number of TS terms)
%
p=zeros(nsize,M);   % 2d array storing polynomial orders of TS terms
num=zeros(10,1);
icount=0;
loop2(1,N,num);     % Assign polynomial orders in each TS term.
bb = zeros(nsize,1);
for ii=1:nsize
    bb(ii)=1;
    for mm=1:M
        bb(ii) = bb(ii)*factorial(p(ii,mm));
    end
end
%
dmag = 0;
dEntity = 1;
dx=zeros(M,1);      % Initialize interval size array
for ii=1:M
    dx(ii) = length(ii)/ndiv(ii);  % Compute interval size
    delta(ii) = dx(ii)*(N+1);      % Compute horizon size in each dimension
    dmag = dmag + delta(ii)*delta(ii); % Compute max. size of horizon
    dEntity = dEntity*dx(ii);      % Compute jacobian of domain integration
end
dmag = sqrt(dmag);
num=zeros(10,1);
icount=0;

% Generate peridynamic grid and store the coordinates to coord array.
loop3(1,num);
%
% fvec (input): Function values at peridynamic points
fvec = zeros(totnode,1);
for k=1:totnode
    % Specify function
    % f(x1,x2,x3,x4) = x1^2 + x2^2 + x3^2 + x4^2 -- Example
    fvec(k) = coord(k,1)^2 + coord(k,2)^2 + coord(k,3)^2 + coord(k,4)^2;
end
%
% Initialize family array of a material point.
% Increase the size of this array if family members exceed 10000.
nodefam = zeros(10000,1);
% Array storing abs(xsi(m),m=1,M)
idist=zeros(M,1);
nmax = 0;
% Array storing xsi(1)^p1*xsi(2)^p2*...*xsi(M)^pM for each TS term
pvec = zeros(nmax,1) ;
% Array storing weight for each TS term
weight = zeros(nmax,1);
fileID = fopen('PDDO.out','w');
% Array storing output of PDDO
dfvec=zeros(totnode,1);
for k=1:totnode
    if int32(k)/int32(1000)*int32(1000) == int32(k)
        fprintf("k = %d\n",k)
    end
    %
    % Generate family of material point k
    numfam =1;
    nodefam(1) = k;
    for j = 1:totnode
        if j == k
            continue;
        end
        for mm=1:M
            idist(mm) = abs(coord(j,mm) - coord(k,mm));
        end
        inside = true;
        for mm=1:M
            if idist(mm) > delta(mm)
                inside=false;
                break;
            end
        end
        if inside
            numfam = numfam + 1;
            nodefam(numfam) = j;
        end
    end
    % fprintf("k = %d , numfam = %d\n", k, numfam)
    % End of generation of family of material point k
    %
    if(numfam>nmax)
        nmax = numfam;
        fprintf("numfam = %d\n",numfam);
    end
    %
    % Compute shape matrix Amat
    Amat = zeros(nsize,nsize);
    bvec = zeros(nsize,1);
    for kk=1:numfam
        j = nodefam(kk);
        xsimag = 0;
        for mm=1:M
            xsi(mm) = coord(j,mm) - coord(k,mm);
            xsimag = xsimag + xsi(mm)*xsi(mm);
        end
        xsimag = sqrt(xsimag);
        for ii=1:nsize
            pvec(ii) = 1.0;
            for mm=1:M
                pvec(ii) = pvec(ii)*xsi(mm)^p(ii,mm);
            end
            weight(ii) = exp(-4*(xsimag/dmag)^2);
            %weight[ii] = 1.0;
        end
        for ii=1:nsize
            for jj=1:nsize
                Amat(ii,jj) = Amat(ii,jj) + ...
                    weight(ii)*pvec(ii)*pvec(jj)*dEntity;
            end
        end
    end
    % End of computation of shape matrix Amat
    %
    % Invert shape matrix and store the inversion to AmatInv
    AmatInv = inv(Amat);
    %
    % Compute rhs vector for computing PD function for requested porder
    for ii=1:nsize
        imatch = true;
        for mm=1:M
            if p(ii,mm)~= porder(mm)
                imatch = false;
                break;
            end
        end
        if imatch
            bvec(ii) = bb(ii);
            break;
        end
    end
    %
    % Compute the coefficients for compution PD function associated with
    % porder
    avec = AmatInv*bvec;
    %
    %
    % Compute derivative based on specified porder
    dfval = 0;
    for kk=1:numfam
        j = nodefam(kk);
        ff = fvec(j);
        xsimag = 0;
        for mm=1:M
            xsi(mm) = coord(j,mm) - coord(k,mm);
            xsimag = xsimag + xsi(mm)*xsi(mm);
        end
        xsimag = sqrt(xsimag);
        %
        % Compute PD function
        gfun = 0;
        for ii=1:nsize
            pvec(ii) = 1;
            for mm=1:M
                pvec(ii) = pvec(ii)*xsi(mm)^p(ii,mm);
            end
            weight(ii) = exp(-4*(xsimag/dmag)^2);
            gfun = gfun + avec(ii)*weight(ii)*pvec(ii);
        end
        % Compute the derivative from PD integration
        dfval = dfval + fvec(j)*gfun*dEntity;
    end
    dfvec(k) = dfval;
    for mm=1:M
        fprintf(fileID,"%f , ",coord(k,mm));
    end
    fprintf(fileID,"%f , %f\n",fvec(k),dfvec(k));
end
%
% Function to compute factorial of an integer number
function f = factorial1(n)
if(n<=1)
    f=1;
else
    f = n*factorial1(n-1);
end
end
%
% Function to count the total number of TS terms
function loop1(m, n)
global M;
global icount;
if(m>0) && (m<=M)
    for i=0:n
        loop1(m+1,n-i);
        if m==M
            icount=icount+1;
        end
    end
end
end
%
% Function to identify and assign the powers of TS terms
function loop2(m, n, num)
global icount;
global p;
global M;
num1=zeros(10,1);
for ii=1:10
    num1(ii)=num(ii);
end
if(m>0) && (m<=M)
    for i=0:n
        num1(m) = i;
        loop2(m+1,n-i,num1);
        if m==M
            icount=icount+1;
            for m=1:M
                p(icount,m) = num1(m);
            end
        end
    end
end
end
%
% Function to generate the grid of the multidimensional domain
function loop3(m, num)
global icount;
global ndiv;
global coord;
global M;
global dx;
num1 = zeros(10,1);
for ii=1:10
    num1(ii)=num(ii);
end
if(m>0 && m<=M)
    for i=1:ndiv(m)
        num1(m)=i;
        loop3(m+1,num1);
        if(m==M)
            icount=icount+1;
            for ii=1:M
                coord(icount,ii) = dx(ii)/2.0+(num1(ii)-1)*dx(ii);
            end
        end
    end
end

end
