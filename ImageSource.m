%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Program Details: An implementation of Image source method for generation
% shoebox-type-room impulse response for the given room dimension. Each 
% wall has different absorption coefficient.The output generated has
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
Lx = 2.9;
if isstring(Lx)==1
    error('Please enter numeric value for Length of room');
end
% Dimension in Y direction (height) in meters
Ly = 3.5;
if isstring(Ly)==1
    error('Please enter numeric value for height of room');
end
% Dimension in Z direction (width) in meters
Lz = 2.3;
if isstring(Lz)==1
    error('Please enter numeric value for width of room');
end


%%DIMENSION 02 (UNCOMMENT TO USE THESE VALUES)
% % Dimension in X direction (length) in meters
% Lx = 10.8;
% if isstring(Lx)==1
%    error('Please enter numeric value for Length of room');
% end
% % Dimension in Y direction (height) in meters
% Ly = 13.5;
% if isstring(Ly)==1
%    error('Please enter numeric value for height of room');
% end
% % Dimension in Z direction (width) in meters
% Lz = 20.3;
% if isstring(Lz)==1
%    error('Please enter numeric value for width of room');
% end


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

%----------------------------------------%
% WET 1
%----------------------------------------%
% % Coefficients for wall 1
% desc='WET';
% alpha1 = 0.01; 
% if (alpha1>1 || alpha1<-1)==1
%     error('Please enter value within range')
% end
% R1 = sqrt(1-alpha1);
% 
% % coefficients for wall 2
% alpha2 = 0.01;    
% if (alpha2>1 || alpha2<-1)==1
%     error('Please enter value within range')
% end
% R2 = sqrt(1-alpha2);
% 
% % coefficients for floor
% alpha3 = 0.02;  
% if (alpha3>1 || alpha3<-1)==1
%     error('Please enter value within range')
% end
% R3 = sqrt(1-alpha3);
% 
% % coefficients for ceiling
% alpha4 = 0.20;  
% if (alpha4>1 || alpha4<-1)==1
%     error('Please enter value within range')
% end
% R4 = sqrt(1-alpha4);
% 
% %coefficients for side wall 1
% alpha5 = 0.06;    
% if (alpha5>1 || alpha5<-1)==1
%     error('Please enter value within range')
% end
% R5 = sqrt(1-alpha5);
% 
% %coefficient for side wall 2
% alpha6  = 0.04;
% if (alpha6>1 || alpha6<-1)==1
%     error('Please enter value within range')
% end
% R6 = sqrt(1-alpha6);

%----------------------------------------%
% WET 2 (UNCOMMENT TO USE THESE VALUES)
%----------------------------------------%
% Coefficients for wall 1
desc = 'WET';
alpha1 = 0.03;     
if (alpha1>1 || alpha1<-1)==1
   error('Please enter value within range')
end
R1 = sqrt(1-alpha1);

% coefficients for wall 2
alpha2 = 0.03;  
if (alpha2>1 || alpha2<-1)==1
  error('Please enter value within range')
end
R2 = sqrt(1-alpha2);

% coefficients for floor
alpha3 = 0.21;   
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
alpha5 = 0.03; 
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

%----------------------------------------%
% DRY 1 (UNCOMMENT TO USE THESE VALUES)
%----------------------------------------%

% % Coefficients for wall 1
%desc = 'DRY';
% alpha1 = 0.55;    
%if (alpha1>1 || alpha1<-1)==1
%    error('Please enter value within range')
%end
% R1 = sqrt(1-alpha1);
% 
% % coefficients for wall 2
% alpha2 = 0.55;  
%if (alpha2>1 || alpha2<-1)==1
%   error('Please enter value within range')
%end
% R2 = sqrt(1-alpha2);
% 
% % coefficients for floor
% alpha3 = 0.26; 
%if (alpha3>1 || alpha3<-1)==1
%    error('Please enter value within range')
%end
% R3 = sqrt(1-alpha3);
% 
% % coefficients for ceiling
% alpha4 = 0.17; 
%if (alpha4>1 || alpha4<-1)==1
%   error('Please enter value within range')
%end
% R4 = sqrt(1-alpha4);
% 
% %coefficients for side wall 1
% alpha5 = 0.16;   
%if (alpha5>1 || alpha5<-1)==1
%    error('Please enter value within range')
%end
% R5 = sqrt(1-alpha5);
% 
% %coefficient for side wall 2
% alpha6  = 0.15;   
%if (alpha6>1 || alpha6<-1)==1
%    error('Please enter value within range')
%end
% R6 = sqrt(1-alpha6);

%----------------------------------------%
% DRY 2 (UNCOMMENT TO USE THESE VALUES)
%----------------------------------------%
% 
% % Coefficients for wall 1
% desc = 'DRY';
% alpha1 = 0.59;    
% if (alpha1>1 || alpha1<-1)==1
%    error('Please enter value within range')
% end
% R1 = sqrt(1-alpha1);
% 
% % coefficients for wall 2
% alpha2 = 0.55;   
% if (alpha2>1 || alpha2<-1)==1
%   error('Please enter value within range')
% end
% R2 = sqrt(1-alpha2);
% 
% % coefficients for floor
% alpha3 = 0.30;  
% if (alpha3>1 || alpha3<-1)==1
%    error('Please enter value within range')
% end
% R3 = sqrt(1-alpha3);
% 
% % coefficients for ceiling
% alpha4 = 0.37; 
% if (alpha4>1 || alpha4<-1)==1
%   error('Please enter value within range')
% end
% R4 = sqrt(1-alpha4);
% 
% %coefficients for side wall 1
% alpha5 = 0.25;  
% if (alpha5>1 || alpha5<-1)==1
%    error('Please enter value within range')
% end
% R5 = sqrt(1-alpha5);
% 
% %coefficient for side wall 2
% alpha6  = 0.28; 
% if (alpha6>1 || alpha6<-1)==1
%    error('Please enter value within range')
% end
% R6 = sqrt(1-alpha6);

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
            % Calculation of impulse response
%-------------------------------------------------------------------------%
%Start stop watch timer
tic
for d = -Nx:Nx
    if mod(d,2)~=0
        %when d is odd
        Ad = (d + 1) * (Lx) - p - a;
    else
        % when d is even
        Ad = d*Lx + p - a;
    end
    
        for e = -Ny:Ny
            %when e is odd
          if mod(e,2)~=0
              Be = (e + 1) * (Ly) - q - b;
          else
              % when e is even
              Be = e*Ly + q - b;
          end  
          
            for f = -Nz:Nz
                if mod(f,2)~=0
                    %when f is odd
                     Cf = (f + 1) * (Lz) - r - c;
                else
                    %when f is even
                      Cf = f*Lz + r - c;
                end
                    %Calculation of distance between the listener and 
                    %virtual source in a virtual room using Pythagoras
                    %theorem
                    L = sqrt((Ad.^2) + (Be.^2) + (Cf.^2));
                    
                    %Calculation of time of arrival of wave
                    t = L/Cair;
                    
                    %Calculation of number of collision between the walls
                    %and sound waves
                    W = abs(d) + abs(e) +abs(f);
                    
                    %Calculation of total reflection coefficient
                    R = (R1.^W)*(R2.^W)*(R3.^W)*(R4.^W)*(R5.^W)*(R6.^W);
                    
                    %Calculation of magnitude of each impulse
                    g = (R) / L;
                    
                    %Calculation of sample bin 
                    bin = round(t*Fs);
                    
                    %Impulse response and its magnitude
                    impulse_resp(bin) = g;
            end
        end
end

%Stop stop watch
toc
  
%-------------------------------------------------------------------------%
        % Naming convention for generated output impulse response
%-------------------------------------------------------------------------%
% Conversion of dimensions of room into string
dimension_1 = num2str(Lx);
dimension_2 = num2str(Ly);
dimension_3 = num2str(Lz);
filename = sprintf('IR_%sX%sX%s_%s_S1889125_Sonawane.wav',...
                    dimension_1,dimension_2,dimension_3,desc);
                
% Save generated impulse response as .wav file
audiowrite(filename,impulse_resp,Fs);

%Play impulse response
soundsc(impulse_resp, Fs);

%Plot impulse response
plot(impulse_resp);






