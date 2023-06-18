% three-DIMENSIONAL MLS APPROXIMATION
%by amir ashkan mehrabi
%this code made for three dimension relevanted data as input and one output
clc
clear all
%INITIAL SETTING 
m=3;
scale=3;
type='GAUSS';
% SETTING THE COORDINATES OF POINTS
random=1:105;
Count=1:105;
count=zeros(105,104);
for j=1:105
count(j,1:104)=Count(~(random(j)==Count));
end
X=csvread('dTA.csv',0,2,[0 2 0 106]);
Y=csvread('dTA.csv',1,2,[1 2 1 106]);
Z=csvread('dTA.csv',2,2,[2 2 2 106]);
rain_mean1=csvread('dTA.csv',3,2,[3 2 3 106]);
xI=X(random);
yI=Y(random);
zI=Z(random);
rainresultI=rain_mean1(random);
nnodes = length(xI);
Distance =40000*ones(1,nnodes);
result=zeros(1,nnodes);
for j=1:105
x=X(count(j,:));
y=Y(count(j,:));
z=Z(count(j,:));
rain_mean=rain_mean1(count(j,:));
npoints = length(x);
% SET UP NODAL COORDINATES (EVALUATION POINTS)

% DETERMINE RADIUS OF SUPPORT OF EVERY NODE
% Evaluate MLS shape function at all evaluation points x
[F,ARcond] = MLS3DShape(m, nnodes, xI,yI,zI, npoints, x,y,z,rain_mean, Distance, type, 3.0 ,j);
while (abs(F-rainresultI(j))/rainresultI(j))>=0.03
    Distance(j)=Distance(j)+2000;
 [F,ARcond] = MLS3DShape(m, nnodes, xI,yI,zI, npoints, x,y,z,rain_mean, Distance, type, 3.0 ,j);
 if ARcond>=5e-19
     break
 end
end

result(j)=F;
end

% Curve fitting. y = peaks(x,y)

fid1 = fopen('fun.dat','w');

fprintf(fid1,'%10s%10s%10s%10s%10s\n', ' ', 'Exact','Approximate');

for j = 1 : nnodes
   fprintf(fid1,'%10.4f', j);
    fprintf(fid1,'%10.4f%10.4f%10.4f\n', rainresultI(j), result(j));
end
fclose(fid1);
%r2 ceoficient
Mean=mean(rainresultI);
A1=zeros(105,1);
A2=zeros(105,1);
for i=1:105
 A1(i)=(rainresultI(i)-result(i))^2;   
 A2(i)= (rainresultI(i)-Mean)^2;  
end
R2=1-(sum(A1)/sum(A2))
 


