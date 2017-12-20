function onedimension(k,s)
% Electron in crystal - 1 dimension
% Authors: Kelly C. Offield and Ed Bueler
% University of Alaska Fairbanks
% Updated 18 December 2017

% The quantum mechanical problem of the motion of an electron through
% a lattice of identical atoms. This program finds the matrix H (the 
% Hamiltonian operator) so that the eigen values (=allowed energies)
% may be found from the time-independent Schroedinger equation
% HPhi=EPhi. Now, H = -((hbar)^2)/(2mass)*(del^2) + V, where hbar is
% the quantum constant, mass is the mass of the electron, del is the
% delta operator, and V is the potential energy function that
% determines the motion of the particle. This program approximates the
% second derivative (del^2)Phi as the difference (1/dx^2)*
% [Phi(x+dx)-2Phi(x)+Phi(x-dx)] and this formula is used recursively
% in the first matrix calculated in this code, 'A'. The 2nd matrix
% calculated in this code, 'B' represents the 2nd term in the
% Hamiltonian, or V. if we let c = (hbar^2)/(2*mass), then the working
% equation for this program is (1/c)*H = (1/d^2)*[-Phi(n-1)+2*Phi(n)-
% Phi(n+1)]+(1/c)V(n)Phi(n) = [(1/d^2)*A + B]Phi.

if nargin < 1
    k = 40;             % number of atoms
end
if nargin < 2
    s = 50;             % number subintervals between atoms
end

% physical parameters and function
hbar = 1.05457e-34;     % planck's constant [Joule-seconds]
mass = 9.10938e-31;     % electron mass [kg]
c = (hbar^2)/(2*mass);
m = 3e-10;			% distance between atoms [meters]; 10^-10 m = 1 angstrom
                    % compare observed values at https://en.wikipedia.org/wiki/Lattice_constant

%%%%%%%%%%%%% HAS MAJOR FIXES
%V0 = 10e30*hbar;		% The depth of the potential wells;
%                        % = 1.05457e-4 Joules?
V0 = 1.54e-18;		% The depth of the potential wells; see below
V = @(x) -V0*cos(x*2*pi/m) - V0; 	% The potential energy function

% coulomb potential would be:
%   V_C(r) = k_e q Q / r
% where k_e = 8.987e9 N m^2 / C^2    is coulombs constant
%         q = -1.60217662e-19 C      is charge of electron
% and we can assume Q = -q (roughly); then
%   V_C(r) = -2.31e-28 / r    Joules (if r in m)
% so scale of V should be comparable to change in V_C(r) if r changes from
% (e.g.) crystal spacing to one-third of crystal spacing:
%   V_C(3e-10) - V_C(1e-11) = 1.54e-18
%%%%%%%%%%%%%

% numerical parameters
dx = m/s;			% spacing between grid points
L = k*m;			% total length of crystal
M = round(L/dx);    % size of MxM matrix
x = 0:dx:L-dx;  		% the grid on the interval [0,L)
% xx = mod(x,m);	 	% due to periodicity of V(x)...

fprintf('using %d atoms in crystal with spacing %.3e m\n',k,m)
fprintf('using %d x %d matrix to approximate the Hamiltonian\n',M,M)

figure(1), plot(x,V(x))
xlabel x, ylabel('energy (Joules)')

%%%%%%%%%%%%%%%%% ADDITIONAL PLOTTING
xxx = -L:dx:2*L-dx;
Vtrue = zeros(size(xxx));
for j = -k:2*k
    VC = 2.31e-28 ./ abs(xxx - j * m);   % coulomb potential for nucleus at x = j m
    Vtrue = Vtrue - VC;
end
hold on, plot(xxx,Vtrue,'r'), hold off
axis([0, L, -10*V0, V0])
%%%%%%%%%%%%%%%%%


%%%%% AS CONSTRUCTED, THIS MATRIX IS THE NEGATIVE OF dx^2 TIMES THE SECOND DERIVATIVES
% This first matrix 'A' will represent the first term in the Hamiltonian, or the second spatial derivative term.
	% Define the size of the matrix 'A'
		A = zeros(M,M);
	% Boundary Conditions + something else???
		A(1,1) = 2;  A(M,M)   = 2;
		A(1,2) = -1; A(M,M-1) = -1;
		A(1,M) = -1; A(M,1)   = -1;
	% Loop through A for difference rule for 2nd derivative
		for j = 2:M-1
		    A(j,j-1) = -1;
		    A(j,j)   = 2;
		    A(j,j+1) = -1;
		end

% This second matrix 'B' will represent the second term in the
% Hamiltonian, or the potential, V.
	% Loop through B for potential values
    	% Assign the size of matrix 'B'
		B = zeros(M,M);
		for j = 1:M
		 		B(j,j) = V(x(j));
		end

% The Hamiltonian for the time-independent Schroedinger equation, H*phi = E*phi is
% H = -c*(1/(dx^2))*A + B; where E will 
% be the eigenvalues. Plot these values and compare to V(max) and to theory. For the current configuration (23 Nov 2017),
% V(max) = 0.

%H = -c*(1/(dx^2))*A + B;
H = c*(1/(dx^2))*A + B;

E = eig(H);


fprintf('there are %d positive energies (scattering)\n',sum(E>0))
fprintf('there are %d negative energies (bound)\n',sum(E<=0))

figure(2)
plot(zeros(size(E)),E,'ko')
axis([-1, 1, -4*V0, 2*V0])
ylabel('energy (Joules)')
grid on

%stairs(x,E)
%scatter(x,E)
%xx = xx(1:length(xx)-1);
%stairs(xx,E)
%hold on
%V_max = max(V(x));  V_max(1:length(xx)) = V_max;
%V_max = max(V(x)) * ones(size(x));
%plot(x,V_max)
%hold off
