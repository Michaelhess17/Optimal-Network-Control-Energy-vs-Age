function [ x, u ] = optimalControlContinuous( A, B, rho, x0, xT, T )
nStep = 1000;
p = size(A,1);

% redefine the variables and reduce the equation to stardard linear
% equations
I = eye(p);
A_tilde = [A, -1.0/(2*rho)*(B*B');-2*I, -A'];
b = [zeros(p);I]*2*xT;
b_tilde = A_tilde\b;
b1 = b_tilde(1:p); b2 = b_tilde(p+1:end);

% compute the constant vector c = [c1;c2]
c1 = x0+b1;
E = expm(-A_tilde*T);

% E = [E11, E12; E21, E22];
E11 = E(1:p,1:p); E12 = E(1:p,p+1:2*p);
E21 = E(p+1:2*p,1:p); E22 = E(p+1:2*p, p+1:2*p);

% calculate PT, c2
PT = E12\(c1- E11*xT- E11*b1-E12*b2);
c2 = E21*xT + E22*PT + E21*b1 + E22*b2;
P0 = c2 - b2;
c = [c1;c2];

% calculate x,P, u
xstar = zeros(2*p, nStep);

tt = expm(1/nStep*A_tilde);
ct = c;
for i = 1:nStep
    ct = tt*ct;
    xstar(:,i) = ct - b_tilde;
end
xstar0 = [x0;P0];
xstar = [xstar0, xstar];
x = xstar(1:p,:);
P = xstar(p+1:2*p,:);
u = -B'*P/(2*rho);
end