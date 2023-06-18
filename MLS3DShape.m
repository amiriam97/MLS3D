function [F,ARcond] = MLS3DShape(m, nnodes, xI,yI,zI, npoints, x,y,z,rain_mean, Distance, type, para,j) 
% SHAPE FUNCTION OF 3D MLS APPROXIMATION
% INPUT PARAMETERS
%    m - Total number of basis functions (1: Constant basis;  3: Linear
%    basis, 6:  basis
%    nnodes  - Total number of nodes used to construct MLS approximation
%    npoints - Total number of points whose MLS shape function to be evaluated
%    xI,yI(nnodes) - Coordinates of nodes used to construct MLS approximation
%    xi,yi(npoints) - Coordinates of points whose MLS shape function to be evaluated
%    Distance(nnodes) - Radius of support of nodes
%    type - Type of weight function
%    para  - Weight function parameter

wI= zeros (1, npoints);  % Weight funciton
xII = zeros (1, npoints);
yII = zeros (1, npoints);
zII = zeros (1, npoints);
% INITIALIZE SHAPE FUNCTION MATRICES
counter=zeros(1,npoints);
% LOOP OVER ALL EVALUATION POINTS TO CALCULATE VALUE OF SHAPE FUNCTION Fi(X)

for i = 1:npoints
r= sqrt((x(i)-xI(j)).^2+(y(i)-yI(j)).^2+(z(i)-zI(j)).^2)/ Distance(j);
if r<=1
 counter(i)=i;   
end
end    


for i =find(counter)
    r= sqrt((x(i)-xI(j)).^2+(y(i)-yI(j)).^2+(z(i)-zI(j)).^2)/ Distance(j);
     [wI(i)] =Weight3D(type, para,r);
      xII(1,i)=x(i);
       yII(1,i)=y(i);
        zII(1,i)=z(i);
end
mm=length(find(counter));
xII=xII(find(xII));
yII=yII(find(yII));
zII=zII(find(zII));
w=diag(wI(find(wI))); 
 if (m == 1)  % Shepard function
  
      p = ones(length(find(counter)), 1); 
      B    = p' * w;
   elseif (m == 3)     
      p = [ones(mm,1), xII',yII',zII']; 
      B    = p' * w;
      elseif (m == 6)     
      p = [ones(mm,1), xII',yII',zII',[xII.*xII]',[xII.*yII]',[yII.*yII]',[xII.*zII]',[zII.*zII]']; 
      B    = p' * w;
 end
	
A=B*p;	

ARcond = rcond(A);


if m==1
 pp=1;   
elseif m==3
 pp=[1,xI(j),yI(j),zI(j)];
else
 pp=[1,xI(j),yI(j),zI(j),xI(j).*xI(j),xI(j).*yI(j),yI(j).*yI(j),xI(j).*zI(j),zI(j).*zI(j)];
end
 AInv = inv(A);
 
 u=rain_mean(find(counter));
  
 PHI = pp* AInv*B; % shape function
 F=PHI*u'; 
  