% Настройка карты и передатчика
viewer = siteviewer(Buildings="map.osm", Basemap="topographic"); 

% Передатчик (tx)
tx = txsite(Name="Small cell transmitter", ...
    Latitude=55.01339898946948, ...
    Longitude=82.95073091983797, ...
    AntennaHeight=20, ...
    TransmitterPower=1, ...
    TransmitterFrequency=2.8e9);
show(tx)

rtpm = propagationModel("raytracing", ...
    Method="sbr", ...
    MaxNumReflections=1, ...
    BuildingsMaterial="concrete", ...
    TerrainMaterial="concrete");

% Добавляем явное моделирование отражения от земли
h_tx = tx.AntennaHeight;
h_rx = 0;
D = distance(tx.Latitude, tx.Longitude, 55.01314984273437, 82.94833302497865);
c = physconst("lightspeed");

tau_LOS = (sqrt(D^2 + (h_tx - h_rx)^2)) / c;
tau_GR = (sqrt(D^2 + (h_tx + h_rx)^2)) / c;

% Моделирование отраженного сигнала с учетом диэлектрической проницаемости
f = tx.TransmitterFrequency;
eps_r = 5.31; 
sigma = 0.001; 
eps_0 = 8.854e-12;

lambda = c / f;
k = 2 * pi / lambda;

% Комплексная диэлектрическая проницаемость
epsilon = eps_r - 1i * (sigma / (2 * pi * f * eps_0));

% Коэффициенты отражения для параллельной и перпендикулярной поляризации
R_parallel = (epsilon * cos(tau_GR) - sqrt(epsilon - sin(tau_GR)^2)) / ...
             (epsilon * cos(tau_GR) + sqrt(epsilon - sin(tau_GR)^2));
R_perpendicular = (cos(tau_GR) - sqrt(epsilon - sin(tau_GR)^2)) / ...
                   (cos(tau_GR) + sqrt(epsilon - sin(tau_GR)^2));

% Покрытие передатчика (показывает радиус действия)
coverage(tx, rtpm, ...
    SignalStrengths=-120:-5, ...
    MaxRange=250, ...
    Resolution=3, ...
    Transparency=0.6)

pause(10);

rx = rxsite(Name="Small cell receiver", ...
    Latitude=55.01314984273437, ...
    Longitude=82.94833302497865, ...
    AntennaHeight=0);

los(tx,rx)

rtpm.MaxNumReflections = 1;
clearMap(viewer)
raytrace(tx,rx,rtpm)

ss = sigstrength(rx,tx,rtpm);
disp("Received power using perfect reflection: " + ss + " dBm")

rtpm.BuildingsMaterial = "concrete";
rtpm.TerrainMaterial = "concrete";
raytrace(tx,rx,rtpm)

ss = sigstrength(rx,tx,rtpm);
disp("Received power using concrete materials: " + ss + " dBm")
pause(5);

rtpm.MaxNumReflections = 2;
rtpm.AngularSeparation = "low";

ss = sigstrength(rx,tx,rtpm);
disp("Received power with two-reflection paths: " + ss + " dBm")
pause(5);

clearMap(viewer)
raytrace(tx,rx,rtpm)

rtpm.MaxNumDiffractions = 1;

ss = sigstrength(rx,tx,rtpm);
disp("Received power with two-reflection and one-diffraction paths: " + ss + " dBm")
pause(5);

raytrace(tx,rx,rtpm)

rtpm.MaxNumReflections = 1;
rtpm.MaxNumDiffractions = 0;
clearMap(viewer)
show(tx)

coverage(tx, rtpm, ...
    SignalStrengths=-120:-5, ...
    MaxRange=250, ...
    Resolution=2, ...
    Transparency=0.4)

pause(10);


azvec = -180:180; % Azimuth angles (deg)
elvec = -90:90; % Elevation angles (deg)
SLA = 30; % Maximum side-lobe level attenuation (dB)
tilt = 0; % Tilt angle (deg)
az3dB = 65; % 3 dB beamwidth in azimuth (deg)
el3dB = 65; % 3 dB beamwidth in elevation (deg)
lambda = physconst("lightspeed")/tx.TransmitterFrequency; % Wavelength (m)

[az,el] = meshgrid(azvec,elvec);
azMagPattern = -min(12*(az/az3dB).^2,SLA);
elMagPattern = -min(12*((el-tilt)/el3dB).^2,SLA);
combinedMagPattern = -min(-(azMagPattern + elMagPattern), SLA); % Relative antenna gain (dB)

antennaElement = phased.CustomAntennaElement(MagnitudePattern=combinedMagPattern);
tx.Antenna = phased.URA(Size=[2 2], ...
    Element=antennaElement, ...
    ElementSpacing=[lambda/2 lambda/2]);

antennaDirectivity = pattern(tx.Antenna, tx.TransmitterFrequency);
antennaDirectivityMax = max(antennaDirectivity(:));
disp("Peak antenna directivity: " + antennaDirectivityMax + " dBi")
pause(5);

tx.AntennaAngle = -90;

clearMap(viewer)
show(rx)
pattern(tx, Transparency=0.6)
hide(tx)

rtpm.MaxNumReflections = 1;
rtpm.MaxNumDiffractions = 0;
ray = raytrace(tx,rx,rtpm);
disp(ray{1})

aod = ray{1}.AngleOfDeparture;
steeringaz = wrapTo180(aod(1)-tx.AntennaAngle(1));
steeringVector = phased.SteeringVector(SensorArray=tx.Antenna);
sv = steeringVector(tx.TransmitterFrequency,[steeringaz;aod(2)]);
tx.Antenna.Taper = conj(sv);

pattern(tx, Transparency=0.6)
raytrace(tx,rx,rtpm);
hide(tx)

ss = sigstrength(rx,tx,rtpm);
disp("Received power with beam steering: " + ss + " dBm")
pause(10);