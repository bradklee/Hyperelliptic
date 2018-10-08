\\genus-degree formula: g = floor((d-1)/2) .
\\for d=3,4,5,6... ; g=1,1,2,2...  
v=[];
for(deg=3,8,{
RandCoeffs = vector(deg,j,random(10));
Monomials = vector(deg,j,x^j)~;
PFData = HyperellipticPicardFuchs(RandCoeffs*Monomials+x^3+x^deg);
CheckZeros = [a==0|a<-PFData[3]]*vector(length(PFData[3]),j,1)~;
v=concat(v,[[deg,(length(PFData[3])-1)/2,CheckZeros,CheckCertificate(PFData)]]);
})
v

\\apparent genus-degree formula: g' = floor(d/4) .
\\for d=4,6,8,10... ; g'=1,1,2,2...  
v=[];
for(deg=2,7,{
RandCoeffs = vector(deg,j,random(10));
Monomials = vector(deg,j,x^(2*j))~;
PFData = HyperellipticPicardFuchs(RandCoeffs*Monomials+x^(2*deg));
CheckZeros = [a==0|a<-PFData[3]]*vector(length(PFData[3]),j,1)~;
v=concat(v,[[2*deg,(length(PFData[3])-1)/2,CheckZeros,CheckCertificate(PFData)]]);
})
v
