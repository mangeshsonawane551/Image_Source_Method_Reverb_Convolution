%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%Author: Mangesh Sonawane 
%Program Details: An implementation of Image source method for generation
% shoebox-type-room impulse response for the given room dimension using 
%'Vectorize method and without for loops. 
% Each wall has different absorption coefficient.The output generated has
% the dimension mentioned and changes as the parameter changes for room
% dimension.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;

%Sample rate
Fs = 44100;

% Speed of sound in air m/s
Cair = 343;

%-------------------------------------------------------------------------%
                            %Room dimensions
%-------------------------------------------------------------------------%                            

%DIMENSION 01
% Dimension in X direction (length) in meters
Lx = 11.8;
if isstring(Lx)==1
    error('Please enter numeric value for Length of room');
end
% Dimension in Y direction (height) in meters
Ly = 14.5;
if isstring(Ly)==1
    error('Please enter numeric value for height of room');
end
% Dimension in Z direction (width) in meters
Lz = 18.3;
if isstring(Lz)==1
    error('Please enter numeric value for width of room');
end


%Area of wall, sidewall and ceiling/floor 
A1 = Lx * Ly;
A2 = Lx * Lz;
A3 = Ly * Lz;

%Volume of the room
V = Lx * Ly * Lz;

%-------------------------------------------------------------------------%
            %Aborption coefficient and reflection coefficient
%-------------------------------------------------------------------------%
%'alpha'is the absorption coefficient and 'R' is the reflection coefficient

% Coefficients for wall 1
alpha1 = 0.01;   
if (alpha1>1 || alpha1<-1)==1
   error('Please enter value within range')
end
R1 = sqrt(1-alpha1);

% coefficients for wall 2
alpha2 = 0.01;  
if (alpha2>1 || alpha2<-1)==1
  error('Please enter value within range')
end
R2 = sqrt(1-alpha2);

% coefficients for floor
alpha3 = 0.02;   
if (alpha3>1 || alpha3<-1)==1
   error('Please enter value within range')
end
R3 = sqrt(1-alpha3);

% coefficients for ceiling
alpha4 = 0.20; 
if (alpha4>1 || alpha4<-1)==1
  error('Please enter value within range')
end
R4 = sqrt(1-alpha4);

%coefficients for side wall 1
alpha5 = 0.06;  
if (alpha5>1 || alpha5<-1)==1
   error('Please enter value within range')
end
R5 = sqrt(1-alpha5);

%coefficient for side wall 2
alpha6  = 0.04; 
if (alpha6>1 || alpha6<-1)==1
   error('Please enter value within range')
end
R6 = sqrt(1-alpha6);
%-------------------------------------------------------------------------%
                        % Calculation of T60
%-------------------------------------------------------------------------%
% Calculating denominator of T60 equation
DenominatorT60 = Cair *(alpha1*A1 + alpha2*A1 + alpha3*A2 + alpha4*A2 + ...
    alpha5*A3 + alpha6*A3 );

%Calculation of T60
T60 = (12*log(10)*V)/(DenominatorT60);

%Initialising vector for impulse respone output
impulse_resp = zeros(ceil(T60*Fs),1);

% N is the maximum integer image source order in X,Y,Z direction
Nx = Cair*T60/Lx;
Ny = Cair*T60/Ly;
Nz = Cair*T60/Lz;


%position of source where x=p, y=q, z=r directions i.e. (p, q,r)
p = 1;
q = 1;
r = Lz/sqrt(2);


%position of listener where x=a, y=b, z=c directions i.e. (a,b,c)
a = Lx/sqrt(2);
b = Ly/sqrt(2);
c = Lz/sqrt(2); 

%-------------------------------------------------------------------------%
         % Calculation of impulse response without using for loops
%-------------------------------------------------------------------------%
%Start stop watch timer
tic

%Creates 3D full grid vector
[d,e,f] = ndgrid(-Nx:Nx , -Ny:Ny, -Nz:Nz);

Ad = (d*Lx) + (mod(d,2)*Lx) +(-1.^d)*p - a;
Be = (e*Ly) + (mod(e,2)*Ly) +(-1.^e)*q - b;
Cf = (f*Lz) + (mod(f,2)*Lz) +(-1.^f)*r - c;
 %Calculation of distance between the listener and 
 %virtual source in a virtual room using Pythagoras theorem
L = sqrt((Ad.^2) + (Be.^2) + (Cf.^2));
                    
%Calculation of time of arrival of wave 
t = L/Cair;                    
%Calculation of number of collision between the walls
%and sound waves
W = abs(d) + abs(e) +abs(f);                    
%Calculation of total reflection coefficient
R = (R1.^W).*(R2.^W).*(R3.^W).*(R4.^W).*(R5.^W).*(R6.^W);
%Calculation of magnitude of each impulse
g = (R) ./ L;                    
%Calculation of sample bin 
bin = round(t*Fs);
%Impulse response and its magnitude
impulse_resp(bin) = g;
%Stop stop watch
toc
  
%-------------------------------------------------------------------------%
        % output impulse response
%-------------------------------------------------------------------------%

%Play impulse response
soundsc(impulse_resp, Fs);

plot(impulse_resp);






