classdef ImageSource_ConvReverb_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        FloorEditField               matlab.ui.control.NumericEditField
        CeilingEditFieldLabel        matlab.ui.control.Label
        CeilingEditField             matlab.ui.control.NumericEditField
        Sidewall1EditFieldLabel      matlab.ui.control.Label
        Sidewall1EditField           matlab.ui.control.NumericEditField
        Sidewall2EditFieldLabel      matlab.ui.control.Label
        Sidewall2EditField           matlab.ui.control.NumericEditField
        Wall1EditFieldLabel          matlab.ui.control.Label
        Wall1EditField               matlab.ui.control.NumericEditField
        Wall2EditFieldLabel          matlab.ui.control.Label
        Wall2EditField               matlab.ui.control.NumericEditField
        LengthEditFieldLabel         matlab.ui.control.Label
        LengthEditField              matlab.ui.control.NumericEditField
        WidthEditFieldLabel          matlab.ui.control.Label
        WidthEditField               matlab.ui.control.NumericEditField
        HeightEditFieldLabel         matlab.ui.control.Label
        HeightEditField              matlab.ui.control.NumericEditField
        PlayImpulseResponseButton    matlab.ui.control.Button
        PlayConvolveSoundButton      matlab.ui.control.Button
        ImageSourceMethodReverberationLabel  matlab.ui.control.Label
        FloorLabel                   matlab.ui.control.Label
        AbsorptionCoefficientsLabel  matlab.ui.control.Label
    end

    methods (Access = private)

        % Value changed function: FloorEditField
        function FloorEditFieldValueChanged(app, event)
            value = app.FloorEditField.Value;
            
        end

        % Value changed function: CeilingEditField
        function CeilingEditFieldValueChanged(app, event)
            value = app.CeilingEditField.Value;
            
        end

        % Value changed function: Sidewall1EditField
        function Sidewall1EditFieldValueChanged(app, event)
            value = app.Sidewall1EditField.Value;
            
        end

        % Value changed function: Sidewall2EditField
        function Sidewall2EditFieldValueChanged(app, event)
            value = app.Sidewall2EditField.Value;
            
        end

        % Value changed function: Wall1EditField
        function Wall1EditFieldValueChanged(app, event)
            value = app.Wall1EditField.Value;
            
        end

        % Value changed function: Wall2EditField
        function Wall2EditFieldValueChanged(app, event)
            value = app.Wall2EditField.Value;
            
        end

        % Value changed function: LengthEditField
        function LengthEditFieldValueChanged(app, event)
            value = app.LengthEditField.Value;
            
        end

        % Value changed function: WidthEditField
        function WidthEditFieldValueChanged(app, event)
            value = app.WidthEditField.Value;
            
        end

        % Value changed function: HeightEditField
        function HeightEditFieldValueChanged(app, event)
            value = app.HeightEditField.Value;
            
        end

        % Button pushed function: PlayImpulseResponseButton
        function PlayImpulseResponseButtonPushed(app, event)
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Program Details: An implementation of Image source method for generation
% shoebox-type-room impulse response for the given room dimension. Each 
% wall has different absorption coefficient.The output generated has
% the dimension mentioned and changes as the parameter changes for room
% dimension.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%Sample rate
Fs = 44100;

% Speed of sound in air m/s
Cair = 343;

%-------------------------------------------------------------------------%
                            %Room dimensions
%-------------------------------------------------------------------------%                            


% 
%DIMENSION 04 (UNCOMMENT TO USE THESE VALUES)
% Dimension in X direction (length) in meters
Lx = app.LengthEditField.Value;
if isstring(Lx)==1
    error('Please enter numeric value for length');
end
% Dimension in Y direction (height) in meters
Ly = app.HeightEditField.Value;
if isstring(Ly) ==1
    error('Please enter numeric value for height');
end
% Dimension in Z direction (width) in meters
Lz = app.WidthEditField.Value;
if isstring(Lz)==1
    error('Please enter numeric valur for width ');
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


% Coefficients for wall 1
alpha1 = app.Wall1EditField.Value;    
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R1 = sqrt(1-alpha1);

% coefficients for wall 2
alpha2 = app.Wall2EditField.Value;
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R2 = sqrt(1-alpha2);

% coefficients for floor
alpha3 = app.FloorEditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R3 = sqrt(1-alpha3);

% coefficients for ceiling
alpha4 = app.CeilingEditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R4 = sqrt(1-alpha4);

%coefficients for side wall 1
alpha5 = app.Sidewall1EditField.Value;  
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R5 = sqrt(1-alpha5);

%coefficient for side wall 2
alpha6  = app.Sidewall2EditField.Value;  
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
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
     % output impulse response
%-------------------------------------------------------------------------%

                
%Play impulse response
soundsc(impulse_resp, Fs);









        end

        % Button pushed function: PlayConvolveSoundButton
        function PlayConvolveSoundButtonPushed(app, event)
            

%Sample rate
Fs = 44100;

% Speed of sound in air m/s
Cair = 343;

%-------------------------------------------------------------------------%
                            %Room dimensions
%-------------------------------------------------------------------------%                            


% 
%DIMENSION 04 (UNCOMMENT TO USE THESE VALUES)
% Dimension in X direction (length) in meters
Lx = app.LengthEditField.Value;
% Dimension in Y direction (height) in meters
Ly = app.HeightEditField.Value
% Dimension in Z direction (width) in meters
Lz = app.WidthEditField.Value;

%Area of wall, sidewall and ceiling/floor 
A1 = Lx * Ly;
A2 = Lx * Lz;
A3 = Ly * Lz;

%Volume of the room
V = Lx * Ly * Lz;

%-------------------------------------------------------------------------%
            %Aborption coefficient and reflection coefficient
%-------------------------------------------------------------------------%


% Coefficients for wall 1
alpha1 = app.Wall1EditField.Value;
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R1 = sqrt(1-alpha1);

% coefficients for wall 2
alpha2 = app.Wall2EditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R2 = sqrt(1-alpha2);

% coefficients for floor
alpha3 = app.FloorEditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R3 = sqrt(1-alpha3);

% coefficients for ceiling
alpha4 = app.CeilingEditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R4 = sqrt(1-alpha4);

%coefficients for side wall 1
alpha5 = app.Sidewall1EditField.Value; 
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
end
R5 = sqrt(1-alpha5);

%coefficient for side wall 2
alpha6  = app.Sidewall2EditField.Value;  
if alpha1>1 || alpha1<-1
    error('Please enter the value for alpha 1 between -1 and 1');
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
  



[x,Fs] = audioread('castanetcut_44.wav');

% Checks if stereo and converts to mono
if size(x,2) == 2                            %detects if stereo or mono
    x = sum(x,2) / size(x,2);                %convert to mono
end


% Length of dry audio 
Lx = length(x);

%Length of impulse response sound
Li = length(impulse_resp);

%Summation of total lenght
N = Lx + Li;

%Initialising output vector
y = zeros(N,1);

%-------------------------------------------------------------------------%
          %Frequency domain convolution usinf fast convolution theorem
%-------------------------------------------------------------------------%
fftx = fft(x,N);
ffth = fft(impulse_resp,N);
y = ifft(fftx.*ffth);

soundsc(y,Fs);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create FloorEditField
            app.FloorEditField = uieditfield(app.UIFigure, 'numeric');
            app.FloorEditField.ValueChangedFcn = createCallbackFcn(app, @FloorEditFieldValueChanged, true);
            app.FloorEditField.Position = [125 322 100 22];
            app.FloorEditField.Value = 0.3;

            % Create CeilingEditFieldLabel
            app.CeilingEditFieldLabel = uilabel(app.UIFigure);
            app.CeilingEditFieldLabel.HorizontalAlignment = 'right';
            app.CeilingEditFieldLabel.Position = [68 280 42 22];
            app.CeilingEditFieldLabel.Text = 'Ceiling';

            % Create CeilingEditField
            app.CeilingEditField = uieditfield(app.UIFigure, 'numeric');
            app.CeilingEditField.ValueChangedFcn = createCallbackFcn(app, @CeilingEditFieldValueChanged, true);
            app.CeilingEditField.Position = [125 280 100 22];
            app.CeilingEditField.Value = 0.37;

            % Create Sidewall1EditFieldLabel
            app.Sidewall1EditFieldLabel = uilabel(app.UIFigure);
            app.Sidewall1EditFieldLabel.HorizontalAlignment = 'right';
            app.Sidewall1EditFieldLabel.Position = [46 230 64 22];
            app.Sidewall1EditFieldLabel.Text = 'Side wall 1';

            % Create Sidewall1EditField
            app.Sidewall1EditField = uieditfield(app.UIFigure, 'numeric');
            app.Sidewall1EditField.ValueChangedFcn = createCallbackFcn(app, @Sidewall1EditFieldValueChanged, true);
            app.Sidewall1EditField.Position = [125 230 100 22];
            app.Sidewall1EditField.Value = 0.25;

            % Create Sidewall2EditFieldLabel
            app.Sidewall2EditFieldLabel = uilabel(app.UIFigure);
            app.Sidewall2EditFieldLabel.HorizontalAlignment = 'right';
            app.Sidewall2EditFieldLabel.Position = [46 188 64 22];
            app.Sidewall2EditFieldLabel.Text = 'Side wall 2';

            % Create Sidewall2EditField
            app.Sidewall2EditField = uieditfield(app.UIFigure, 'numeric');
            app.Sidewall2EditField.ValueChangedFcn = createCallbackFcn(app, @Sidewall2EditFieldValueChanged, true);
            app.Sidewall2EditField.Position = [125 188 100 22];
            app.Sidewall2EditField.Value = 0.28;

            % Create Wall1EditFieldLabel
            app.Wall1EditFieldLabel = uilabel(app.UIFigure);
            app.Wall1EditFieldLabel.HorizontalAlignment = 'right';
            app.Wall1EditFieldLabel.Position = [70 142 38 22];
            app.Wall1EditFieldLabel.Text = 'Wall 1';

            % Create Wall1EditField
            app.Wall1EditField = uieditfield(app.UIFigure, 'numeric');
            app.Wall1EditField.ValueChangedFcn = createCallbackFcn(app, @Wall1EditFieldValueChanged, true);
            app.Wall1EditField.Position = [123 142 100 22];
            app.Wall1EditField.Value = 0.03;

            % Create Wall2EditFieldLabel
            app.Wall2EditFieldLabel = uilabel(app.UIFigure);
            app.Wall2EditFieldLabel.HorizontalAlignment = 'right';
            app.Wall2EditFieldLabel.Position = [70 100 38 22];
            app.Wall2EditFieldLabel.Text = 'Wall 2';

            % Create Wall2EditField
            app.Wall2EditField = uieditfield(app.UIFigure, 'numeric');
            app.Wall2EditField.ValueChangedFcn = createCallbackFcn(app, @Wall2EditFieldValueChanged, true);
            app.Wall2EditField.Position = [123 100 100 22];
            app.Wall2EditField.Value = 0.02;

            % Create LengthEditFieldLabel
            app.LengthEditFieldLabel = uilabel(app.UIFigure);
            app.LengthEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel.Position = [372 301 43 22];
            app.LengthEditFieldLabel.Text = 'Length';

            % Create LengthEditField
            app.LengthEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengthEditField.ValueChangedFcn = createCallbackFcn(app, @LengthEditFieldValueChanged, true);
            app.LengthEditField.Position = [430 301 100 22];
            app.LengthEditField.Value = 10;

            % Create WidthEditFieldLabel
            app.WidthEditFieldLabel = uilabel(app.UIFigure);
            app.WidthEditFieldLabel.HorizontalAlignment = 'right';
            app.WidthEditFieldLabel.Position = [378 230 37 22];
            app.WidthEditFieldLabel.Text = 'Width';

            % Create WidthEditField
            app.WidthEditField = uieditfield(app.UIFigure, 'numeric');
            app.WidthEditField.ValueChangedFcn = createCallbackFcn(app, @WidthEditFieldValueChanged, true);
            app.WidthEditField.Position = [430 230 100 22];
            app.WidthEditField.Value = 20;

            % Create HeightEditFieldLabel
            app.HeightEditFieldLabel = uilabel(app.UIFigure);
            app.HeightEditFieldLabel.HorizontalAlignment = 'right';
            app.HeightEditFieldLabel.Position = [374 153 41 22];
            app.HeightEditFieldLabel.Text = 'Height';

            % Create HeightEditField
            app.HeightEditField = uieditfield(app.UIFigure, 'numeric');
            app.HeightEditField.ValueChangedFcn = createCallbackFcn(app, @HeightEditFieldValueChanged, true);
            app.HeightEditField.Position = [430 153 100 22];
            app.HeightEditField.Value = 15;

            % Create PlayImpulseResponseButton
            app.PlayImpulseResponseButton = uibutton(app.UIFigure, 'push');
            app.PlayImpulseResponseButton.ButtonPushedFcn = createCallbackFcn(app, @PlayImpulseResponseButtonPushed, true);
            app.PlayImpulseResponseButton.Position = [70 33 195 28];
            app.PlayImpulseResponseButton.Text = 'Play Impulse Response';

            % Create PlayConvolveSoundButton
            app.PlayConvolveSoundButton = uibutton(app.UIFigure, 'push');
            app.PlayConvolveSoundButton.ButtonPushedFcn = createCallbackFcn(app, @PlayConvolveSoundButtonPushed, true);
            app.PlayConvolveSoundButton.Position = [350 33 180 28];
            app.PlayConvolveSoundButton.Text = 'Play Convolve Sound';

            % Create ImageSourceMethodReverberationLabel
            app.ImageSourceMethodReverberationLabel = uilabel(app.UIFigure);
            app.ImageSourceMethodReverberationLabel.HorizontalAlignment = 'center';
            app.ImageSourceMethodReverberationLabel.FontName = 'Avenir';
            app.ImageSourceMethodReverberationLabel.FontSize = 24;
            app.ImageSourceMethodReverberationLabel.FontWeight = 'bold';
            app.ImageSourceMethodReverberationLabel.Position = [78 397 428 35];
            app.ImageSourceMethodReverberationLabel.Text = 'Image Source Method Reverberation';

            % Create FloorLabel
            app.FloorLabel = uilabel(app.UIFigure);
            app.FloorLabel.HorizontalAlignment = 'right';
            app.FloorLabel.Position = [77 322 33 22];
            app.FloorLabel.Text = 'Floor';

            % Create AbsorptionCoefficientsLabel
            app.AbsorptionCoefficientsLabel = uilabel(app.UIFigure);
            app.AbsorptionCoefficientsLabel.HorizontalAlignment = 'right';
            app.AbsorptionCoefficientsLabel.FontWeight = 'bold';
            app.AbsorptionCoefficientsLabel.Position = [15 361 140 22];
            app.AbsorptionCoefficientsLabel.Text = 'Absorption Coefficients';
        end
    end

    methods (Access = public)

        % Construct app
        function app = ImageSource_ConvReverb_GUI_s1889125_sonawane

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
