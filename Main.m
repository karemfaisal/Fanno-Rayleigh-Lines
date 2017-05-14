%Draw a Fanno Line on the T - ?s diagram for :
%0Air flow, To= 288 K, Tx= 285K, Px= 98 kPa
%for chosen downstream pressure Py values: 40 38 36 34 32 30 28 26 24
%22 20 18 16 14 12 10 in KPa
%,,,,,
%Draw a Rayleigh Line on the T - ?s diagram for :
%Air flow, To= 288 K, Tx= 286K, Px= 97 kPa
%for chosen downstream pressure Py values:
%93 86 79 72 62 55 52 51 48 43 41 38 34 31 27 14 7 in KPa
%,,,,,
%Draw a Fanno line & Rayleigh Line on the same T - ?s diagram for :
%Air flow, To= 270 K, Tx= 148K, Px= 98 kPa
%for chosen downstream pressure Py values:
%For Fanno line Py=95:700 in KPa
%For Rayleigh line Py=95:500 in KPa
%...............................................

%Gas Dynamic Project ,, Mechanical Engineering 3rd year Alexandria
%University
%Coded By Karem Ali
k = 1.4;
r=.287;
c_p = 1.005; % Specific Heat Consumption for air
p_y=[40:-2:10]; % pressure for different Points
p_x = 98; % Pressure for a Specific Point x
t_x = 285; % Temp for Ther Specific Point x
t_0 = 288;  % stagnation temperature

% We Will Use The 'Solve' Function to get The Value of m_x
m_x = (2^(1/2)*(t_0 - t_x)^(1/2))/(t_x^(1/2)*(k - 1)^(1/2));  %MachNumber @ Point x

% Denisty * Velocity @ any point is Constant and its Value is 'c'
c = ( p_x / (r*t_x) ) * m_x * sqrt(k*r*t_x);

% We Will Use 'Solve' Function to get the value of t_y 
t_y = -(c_p.*p_y.*(p_y - ((2.*t_0.*c^2*r.^2 + c_p.*p_y.^2)/c_p).^(1/2)))/(c.^2*r.^2) ;  % Temp for points y

delta_s = (c_p .* log(t_y./t_x) - r.* log ( p_y./p_x)) ;  %Change in Entropy

%The Fanno Line   Plot 
subplot(3,1,1);
plot(delta_s,t_y) , xlabel('ds') , ylabel('Temperature') , title('Fanno Line  ');




%Begin Of Rayleigh Line  
p_y = [93 86 79 72 62 55 52 51 48 43 41 38 34 31 27 14 7];

% in Rayleigh Line   We have a Constant to fn(p,t)
constant = p_x+((c)^2*r*t_x)/p_x;

% We Will Use 'Solve' Function to get Ther Value of t_y
t_y =  (constant-p_y).*(p_y./(c.^2.*r)) ;
delta_s = (c_p .* log(t_y./t_x) - r.* log ( p_y./p_x)) ; %Change in Entropy

%Rayleigh Line   Plot
subplot(3,1,2);
plot(delta_s,t_y) , xlabel('ds') , ylabel('Temperature') , title('Rayleigh Line  ');



% Rayleigh Line   & Fanno Line   Together on The Same Plot
t_x = 148;
t_0 = 270;

%Rayleigh Line  
p_y=[98:10:700];
m_x=(2^(1/2)*(t_0 - t_x)^(1/2))/(t_x^(1/2)*(k - 1)^(1/2));
c = ( p_x / (r*t_x) ) * m_x * sqrt(k*r*t_x);
t_y =  -(c_p.*p_y.*(p_y - ((2.*t_0.*c.^2.*r.^2 + c_p.*p_y.^2)/c_p).^(1/2)))/(c.^2.*r.^2) ;
delta_s = (c_p .* log(t_y./t_x) - r.* log ( p_y./p_x)) ;

%Fanno Line  
p_y_2=[98:10:500];
c = ( p_x / (r*t_x) ) * m_x * sqrt(k*r*t_x);
constant=p_x+((c)^2*r*t_x)/p_x;
t_y_2 =  (constant-p_y_2).*(p_y_2./(c.^2.*r)) ;
delta_s_2 = (c_p .* log(t_y_2/t_x) - r.* log ( p_y_2./p_x)) ;

subplot(3,1,3);
plot(delta_s,t_y,delta_s_2,t_y_2) ,xlabel('ds') , ylabel('Temperature'),title('Fanno&Rayleigh Line'),legend('Rayleigh Line  ','Fanno Line  ');

%End of Code 
