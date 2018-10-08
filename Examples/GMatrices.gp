\\Example 1: Print a few G matrices. 
\\The G Matrix is always invertible. 
coeffs = [a,b,c,d,e,f];
vars = [x,x^2,x^3,x^4,x^5,x^6]~;
MatrixG(coeffs[1..3]*vars[1..3],3)
MatrixG(coeffs[1..4]*vars[1..4],4)
MatrixG(coeffs[1..5]*vars[1..5],5)
MatrixG(coeffs[1..6]*vars[1..6],6)
