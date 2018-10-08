\\Title: 	A Picard-Fuchs Algorithm for Hyperelliptic Level Curves
\\Author: 	Bradley Klee
\\License: 	GPLv3

\\Input: Polynomial potential function V(x) 
\\Output: Matrix encoding of Hermite reduction 
{MatrixG(poly,deg)=
my(Gx,GL,GR);
Gz=matrix(2*deg-1,2*deg-1,j,k,if(j==k&&k<deg,z,0));
GL=matrix(2*deg-1,2*deg-1,j,k,if(0<j-k&&j-k<=deg&&k<deg,
	-2*polcoef(poly,j-k),0));
GR=matrix(2*deg-1,2*deg-1,j,k,if(k>=deg&&j-(k-deg)>0&&j-(k-deg)<=deg,
	-(j-(k-deg))*polcoef(poly,j-(k-deg)),0));
Gz+GL+GR};

\\Input: Polynomial potential function V(x)
\\Output: One-form reduction matrices
{ReductionData(poly,deg)=
my(ID,InvG,MU,MV,dMV);
ID=matrix(2*deg-1,2*deg-1,j,k,if(j==k,1,0));
InvG=matinverseimage(MatrixG(poly,deg),ID);
MU=matrix(deg-1,deg-1,j,k,InvG[j,k]);
MV=matrix(deg,deg-1,j,k,InvG[j+deg-1,k]);
dMV=matrix(deg-1,deg-1,j,k,j*MV[j+1,k]);
[MU,dMV,MV]};

\\Input: One-form and One-form reduction matrices
\\Output: Certified, hermite-reduced one form 
{HermiteReduce(dtform,cert,dat,n)=if(n==0,[dtform,cert],
HermiteReduce((dat[1]+1/(2*n-1)*dat[2])*dtform,
cert+1/((2*n-1)*y^(2*n-1))*dat[3]*dtform,
dat,n-1))};

\\n!!
{DoubleFactorial(n) = prod(i=0, (n-1)\2, n - 2*i )};

\\Input: Parity Boolean, one-form reduction matrices
\\Output: Raw Picard-Fuchs data
{PeriodBasis(sym,deg,dat)=
my(dt,bound,reduction);
dt=vector(deg-1,j,if(j==1,1,0))~;
bound=2*if(sym==0,floor((deg-1)/2),floor(deg/4));
decomposition=concat([[dt,vector(deg,j,0)~]],
vector(bound,j,HermiteReduce(((-1/2)^j)*DoubleFactorial(2*j-1)*dt,0,dat,j)));
[vector(bound+1,j,decomposition[j][1]),
vector(bound+1,j,vector(deg,j,x^(j-1))*decomposition[j][2]),
concat([1],vector(bound,j,((-1/(2*y^2))^(j)*DoubleFactorial(2*j-1))))]};

\\Input: Polynomial Potential Function V(x)
\\Output: Processed Picard-Fuchs data
{HyperellipticPicardFuchs(poly)=
my(deg,sym,dat,decomposition,PFCs);
deg=poldegree(poly);
sym=if(subst(poly,x,-x)-poly==0,1,0);
dat = ReductionData(poly,deg);
decomposition = PeriodBasis(sym,deg,dat);
PFCs=lindep(decomposition[1]);
[(1/2)*y^2+poly,decomposition[3]*dt,PFCs,decomposition[2]*PFCs]};

\\Input: Processed Picard-Fuchs data
\\Output: Boolean check value 0
{CheckCertificate(PFData)=
y^(2*length(PFData[2])-2)*(
deriv(PFData[4],x)*deriv(PFData[1],y)
+deriv(PFData[4],y)*deriv(-PFData[1],x)
+PFData[2]*PFData[3]/dt)%(2*PFData[1]-z)
};

