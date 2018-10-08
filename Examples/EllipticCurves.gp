\\ Hypergeometric Quartic (Cf. OEIS: A000897). 
\\ 2F1 parameters: (a,b,c) = (1/4,3/4,1).
PFData = HyperellipticPicardFuchs((1/2)*(x^2-(1/4)*x^4));
[PFData[3],CheckCertificate(PFData)]

\\ Hypergeometric Cubics (Cf. OEIS: A113424). 
\\ 2F1 parameters: (a,b,c) = (1/6,5/6,1).
PFData = HyperellipticPicardFuchs((1/2)*(x^2-c*x^3));
[-vector(3,j,9/4*(polcoef(PFData[3][j],2,c)*(4/27)+polcoef(PFData[3][j],0,c)))~,CheckCertificate(PFData)]

\\ two-parameter elliptic stratum 
PFData = HyperellipticPicardFuchs((1/2)*(x^2+a*x^3+b*x^4));
[PFData[3],CheckCertificate(PFData)]
subst(subst(-PFData[3],a,0),b,-1/4)/(1-2*z+z^2)
pfdat = subst(-PFData[3],b,0)/a^4;
vector(3,j,9/4*(polcoef(pfdat[j],2,a)*4/27+polcoef(pfdat[j],0,a)))~
