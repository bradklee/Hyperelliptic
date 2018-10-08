\\ Parity symmetric sextic (Cf. OEIS: A113424). 
\\ 2F1 parameters: (a,b,c) = (1/6,5/6,1).
PFData = HyperellipticPicardFuchs((1/2)*(x^2-(4/27)*x^6));
[-PFData[3],CheckCertificate(PFData)]

\\ Asymmetric sextic. 
PFData = HyperellipticPicardFuchs((1/2)*(x^2-x^3-(4/27)*x^6));
[-PFData[3],CheckCertificate(PFData)]
length(PFData[3])
