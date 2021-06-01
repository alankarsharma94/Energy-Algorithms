%This is Rebonato method. The most general methodology to create a valid 
% correlation matrix. This method uses spectral decomposition. 
% Useful to decompose non positive determine matrix C as C=U*U'
function Correl = Rebonato(C)

% Find eigesystem S and enginvalues lambdas
[S,lambda]= eig(C);

% check eigenvectors for negative values
% if found replace them with zeros
for i=1:length(lambda)
  for j=1:length(lambda)  
      if ( lambda(i,j) < 0 )
        lambda(i,j) = 0;
      end;  
  end;
end;

% Define non-zero elements of the matrix T with respect to the eigensystem
% S by ti = Sum(S(i)^2*lambda(i))
sqS = S.^2;
N=length(sqS);
T = zeros(N,1);
for i=1:N
    T(:,1) = T(:,1)+ sqS(:,i)*lambda(i,i);
end;    

%Build diagonal matrix  of square root of 1/T
TT1 = diag(1./sqrt(T'));

% Multiply the vector S with their associated "corrected" eigenvalues lambda and
% normilize the row vectors to unit length TT1  
U = TT1*S*sqrt(lambda);
U1 = fliplr(U);

Correl = U1*U1';
%Correl = adj(U1);
