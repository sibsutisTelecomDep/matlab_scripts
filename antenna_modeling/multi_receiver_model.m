% Настройка карты и передатчика
viewer = siteviewer(Buildings="map.osm", Basemap="topographic");

% Передатчик (tx)
tx = txsite(Name="Small cell transmitter", ...
    Latitude=55.01339898946948, ...
    Longitude=82.95073091983797, ...
    AntennaHeight=10, ...
    TransmitterPower=1, ...
    TransmitterFrequency=2.8e9);
show(tx)

rtpm = propagationModel("raytracing", ...
    Method="sbr", ...
    MaxNumReflections=3, ...
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

coverage(tx, rtpm, ...
    SignalStrengths=-120:-5, ...
    MaxRange=250, ...
    Resolution=3, ...
    Transparency=0.6)

pause(20);

rx1 = rxsite(Name="Small cell receiver 1", ...
    Latitude=55.01314984273437, ...
    Longitude=82.94833302497865, ...
    AntennaHeight=0);
los(tx, rx1)

rx2 = rxsite(Name="Small cell receiver 2", ...
    Latitude=55.0139680225036, ...
    Longitude=82.9491698741913, ...
    AntennaHeight=0);
los(tx, rx2)

rx3 = rxsite(Name="Small cell receiver 3", ...
    Latitude=55.012756126219216, ...
    Longitude=82.94897675514221, ...
    AntennaHeight=0);
los(tx, rx3)

rx4 = rxsite(Name="Small cell receiver 4", ...
    Latitude=55.01330698545722, ...
    Longitude=82.95013991951763, ...
    AntennaHeight=0);
los(tx, rx4)

% Моделирование для приемников
rtpm.MaxNumReflections = 1;
clearMap(viewer)
raytrace(tx, rx1, rtpm)
ss1 = sigstrength(rx1, tx, rtpm);
disp("Received power using perfect reflection for RX1: " + ss1 + " dBm")
pause(20);

raytrace(tx, rx2, rtpm)
ss2 = sigstrength(rx2, tx, rtpm);
disp("Received power using perfect reflection for RX2: " + ss2 + " dBm")
pause(20);

raytrace(tx, rx3, rtpm)
ss3 = sigstrength(rx3, tx, rtpm);
disp("Received power using perfect reflection for RX3: " + ss3 + " dBm")
pause(20);

raytrace(tx, rx4, rtpm)
ss4 = sigstrength(rx4, tx, rtpm);
disp("Received power using perfect reflection for RX4: " + ss4 + " dBm")
pause(20);

% Вычисление корреляции между приемниками
signal_strengths = [ss1, ss2, ss3, ss4];  % Мощности полученных сигналов от каждого приемника

% Расчет корреляции
corr_matrix = corrcoef(signal_strengths);
disp("Correlation matrix between received powers: ")
disp(corr_matrix)
pause(20);


% Моделирование с материалами бетона для всех приемников
rtpm.BuildingsMaterial = "concrete";
rtpm.TerrainMaterial = "concrete";

% Для RX1
raytrace(tx, rx1, rtpm)
ss1_concrete = sigstrength(rx1, tx, rtpm);
disp("Received power at RX1 with concrete materials: " + ss1_concrete + " dBm")
pause(20);

% Для RX2
raytrace(tx, rx2, rtpm)
ss2_concrete = sigstrength(rx2, tx, rtpm);
disp("Received power at RX2 with concrete materials: " + ss2_concrete + " dBm")
pause(20);

% Для RX3
raytrace(tx, rx3, rtpm)
ss3_concrete = sigstrength(rx3, tx, rtpm);
disp("Received power at RX3 with concrete materials: " + ss3_concrete + " dBm")
pause(20);

% Для RX4
raytrace(tx, rx4, rtpm)
ss4_concrete = sigstrength(rx4, tx, rtpm);
disp("Received power at RX4 with concrete materials: " + ss4_concrete + " dBm")
pause(20);

% Вычисление корреляции для материалов бетона
signal_strengths_concrete = [ss1_concrete, ss2_concrete, ss3_concrete, ss4_concrete];
corr_matrix_concrete = corrcoef(signal_strengths_concrete);
disp("Correlation matrix between received powers with concrete materials: ")
disp(corr_matrix_concrete)

pause(20);

% Добавление дифракции и отражений
rtpm.MaxNumReflections = 2;
rtpm.MaxNumDiffractions = 1;

% Моделирование с дифракцией и двумя отражениями
raytrace(tx, rx1, rtpm)
ss1_diffraction = sigstrength(rx1, tx, rtpm);
disp("Received power at RX1 with diffraction and two-reflections: " + ss1_diffraction + " dBm")
pause(20);

raytrace(tx, rx2, rtpm)
ss2_diffraction = sigstrength(rx2, tx, rtpm);
disp("Received power at RX2 with diffraction and two-reflections: " + ss2_diffraction + " dBm")
pause(20);

raytrace(tx, rx3, rtpm)
ss3_diffraction = sigstrength(rx3, tx, rtpm);
disp("Received power at RX3 with diffraction and two-reflections: " + ss3_diffraction + " dBm")
pause(20);

raytrace(tx, rx4, rtpm)
ss4_diffraction = sigstrength(rx4, tx, rtpm);
disp("Received power at RX4 with diffraction and two-reflections: " + ss4_diffraction + " dBm")
pause(20);

% Вычисление корреляции для модели с дифракцией
signal_strengths_diffraction = [ss1_diffraction, ss2_diffraction, ss3_diffraction, ss4_diffraction];
corr_matrix_diffraction = corrcoef(signal_strengths_diffraction);
disp("Correlation matrix between received powers with diffraction: ")
disp(corr_matrix_diffraction)
pause(20);

% === Расширенное моделирование погодных условий (влага, газы, туман) ===
disp("=== Weather Condition Modeling (Extended) ===");

% Расширенные параметры дождя
rain = propagationModel("rain", ...
    RainRate=20, ...                    % мм/ч — умеренный дождь
    Polarization="horizontal", ...
    RainHeight=3500);                  % Высота дождевого слоя (м)

% Газовые потери
gas = propagationModel("gas", ...
    Temperature=20, ...                % Температура воздуха (°C)
    Pressure=101.325, ...              % Давление (кПа)
    WaterVaporDensity=7.5);            % Влага в воздухе (г/м³)

% Туман/облачность
fog = propagationModel("fog", ...
    LiquidWaterDensity=0.2, ...        % Плотность воды в воздухе (г/м³)
    Visibility=500);                   % Видимость (м)

% Комбинируем модели
rtPlusWeather = rtpm + gas + rain + fog;

% Моделирование для всех приемников с погодными потерями
raytrace(tx, rx1, rtPlusWeather); ss1_weather = sigstrength(rx1, tx, rtPlusWeather);
disp("Received power at RX1 including weather loss: " + ss1_weather + " dBm")

raytrace(tx, rx2, rtPlusWeather); ss2_weather = sigstrength(rx2, tx, rtPlusWeather);
disp("Received power at RX2 including weather loss: " + ss2_weather + " dBm")

raytrace(tx, rx3, rtPlusWeather); ss3_weather = sigstrength(rx3, tx, rtPlusWeather);
disp("Received power at RX3 including weather loss: " + ss3_weather + " dBm")

raytrace(tx, rx4, rtPlusWeather); ss4_weather = sigstrength(rx4, tx, rtPlusWeather);
disp("Received power at RX4 including weather loss: " + ss4_weather + " dBm")

% Корреляционная матрица
signal_strengths_weather = [ss1_weather, ss2_weather, ss3_weather, ss4_weather];
corr_matrix_weather = corrcoef(signal_strengths_weather);
disp("Correlation matrix including weather impairments:")
disp(corr_matrix_weather)

pause(20);

% Отображение на карте
clearMap(viewer)
show(tx)
show(rx1)
show(rx2)
show(rx3)
show(rx4)

pause(20);

% Настройка антенны с лучевым управлением
rtpm.MaxNumReflections = 1;
rtpm.MaxNumDiffractions = 0;

show(tx)

coverage(tx, rtpm, ...
    SignalStrengths=-120:-5, ...
    MaxRange=250, ...
    Resolution=2, ...
    Transparency=0.4)

pause(20);

% Расчет теоретического RSRP для точек
Txpower = 40;  % Мощность передатчика (dBm)
antGain = 12;  % Коэффициент усиления антенны (dBi)
feedet = 3;    % Потери от фидера (dB)
BW = 10e6;     % Ширина полосы (Hz)
rsrp_freq = 15e3;  % Частота RSRP (Hz)

D1 = distance(tx.Latitude, tx.Longitude, rx1.Latitude, rx1.Longitude); % Расстояние от TX к RX1
D2 = distance(tx.Latitude, tx.Longitude, rx2.Latitude, rx2.Longitude); % Расстояние от TX к RX2
D3 = distance(tx.Latitude, tx.Longitude, rx3.Latitude, rx3.Longitude); % Расстояние от TX к RX3
D4 = distance(tx.Latitude, tx.Longitude, rx4.Latitude, rx4.Longitude); % Расстояние от TX к RX4

% Теоретический расчет RSRP
theoretical_rsrp1 = theoretical_RSRP(D1, Txpower, antGain, feedet, BW, rsrp_freq);
theoretical_rsrp2 = theoretical_RSRP(D2, Txpower, antGain, feedet, BW, rsrp_freq);
theoretical_rsrp3 = theoretical_RSRP(D3, Txpower, antGain, feedet, BW, rsrp_freq);
theoretical_rsrp4 = theoretical_RSRP(D4, Txpower, antGain, feedet, BW, rsrp_freq);

disp("Theoretical RSRP for RX1: " + theoretical_rsrp1 + " dBm");
disp("Theoretical RSRP for RX2: " + theoretical_rsrp2 + " dBm");
disp("Theoretical RSRP for RX3: " + theoretical_rsrp3 + " dBm");
disp("Theoretical RSRP for RX4: " + theoretical_rsrp4 + " dBm");

% Расчет отклонений между практическими и теоретическими значениями
delta1 = abs(ss1 - theoretical_rsrp1);
delta2 = abs(ss2 - theoretical_rsrp2);
delta3 = abs(ss3 - theoretical_rsrp3);
delta4 = abs(ss4 - theoretical_rsrp4);

% Нахождение минимального отклонения
[min_delta, min_idx] = min([delta1, delta2, delta3, delta4]);

disp(['Минимальное отклонение на точке: ', num2str(min_delta), ' dB']);

% Конфигурация мини-радиусов для каждой точки
radii1 = antennaConfiguration(antGain, Txpower, D1);
radii2 = antennaConfiguration(antGain, Txpower, D2);
radii3 = antennaConfiguration(antGain, Txpower, D3);
radii4 = antennaConfiguration(antGain, Txpower, D4);

% Проверка пересечений радиусов
threshold = 5;  % Пороговое значение для пересечения радиусов
if abs(radii1 - radii2) < threshold
    disp('Радиусы RX1 и RX2 пересекаются');
else
    disp('Радиусы RX1 и RX2 не пересекаются');
end

if abs(radii2 - radii3) < threshold
    disp('Радиусы RX2 и RX3 пересекаются');
else
    disp('Радиусы RX2 и RX3 не пересекаются');
end

if abs(radii3 - radii4) < threshold
    disp('Радиусы RX3 и RX4 пересекаются');
else
    disp('Радиусы RX3 и RX4 не пересекаются');
end

% Функция для расчета теоретического RSRP
function rsrp = theoretical_RSRP(D, Txpower, antGain, feedet, BW, rsrp_freq)
    % Считаем потери на расстоянии
    PL = pathLoss(D); % Потери на расстоянии
    rsrp = Txpower + antGain - PL - feedet + 10*log10(BW) - 10*log10(rsrp_freq);
end

% Функция для расчета потерь на расстоянии (FSPL)
function PL = pathLoss(D)
    freq = 2.8e9;  % Частота передатчика
    c = 3e8;       % Скорость света
    lambda = c / freq;
    FSPL = (4 * pi * D / lambda) ^ 2;
    PL = 10 * log10(FSPL); % Потери в dB
end

% Функция для конфигурации мини-радиусов антенны
function radii = antennaConfiguration(antGain, Txpower, D)
    theta = 30; % Угол лепестка антенны
    radius = (Txpower + antGain) - 20*log10(D);  % Простой расчет радиуса
    radii = radius;  % Минимальные радиусы для каждой точки
end

% Задание углов и характеристик антенны
azvec = -180:180; % Азимутальные углы (градусы)
elvec = -90:90; % Углы по высоте (градусы)
SLA = 30; % Максимальная угловая аттенюация боковых лепестков (dB)
tilt = 0; % Угол наклона (градусы)
az3dB = 65; % 3dB ширина в азимуте (градусы)
el3dB = 65; % 3dB ширина в вертикальной плоскости (градусы)
lambda = physconst("lightspeed") / tx.TransmitterFrequency; % Длина волны (м)

% Сеточная сетка углов
[az, el] = meshgrid(azvec, elvec);
azMagPattern = -min(12*(az/az3dB).^2, SLA);
elMagPattern = -min(12*((el-tilt)/el3dB).^2, SLA);
combinedMagPattern = -min(-(azMagPattern + elMagPattern), SLA); % Относительная усиление антенны (dB)

antennaElement = phased.CustomAntennaElement(MagnitudePattern=combinedMagPattern);
tx.Antenna = phased.URA(Size=[2 2], ...
    Element=antennaElement, ...
    ElementSpacing=[lambda/2 lambda/2]);

% Вычисление направленности антенны
antennaDirectivity = pattern(tx.Antenna, tx.TransmitterFrequency);
antennaDirectivityMax = max(antennaDirectivity(:));
disp("Peak antenna directivity: " + antennaDirectivityMax + " dBi")

pause(20);

tx.AntennaAngle = -90;

% Отображение на карте с учетом антенны
clearMap(viewer)
show(tx)
show(rx1)
show(rx2)
show(rx3)
show(rx4)
pattern(tx, Transparency=0.6)
hide(tx)

% Процесс лучевого трассирования с управлением направлением для каждого приемника
% Для RX1
ray = raytrace(tx, rx1, rtpm);
disp(ray{1})

aod = ray{1}.AngleOfDeparture;
steeringaz = wrapTo180(aod(1) - tx.AntennaAngle(1));
steeringVector = phased.SteeringVector(SensorArray=tx.Antenna);
sv = steeringVector(tx.TransmitterFrequency, [steeringaz; aod(2)]);
tx.Antenna.Taper = conj(sv);

% Отображение на карте с направленным лучом для RX1
pattern(tx, Transparency=0.6)
raytrace(tx, rx1, rtpm);
hide(tx)

% Вычисление мощности полученного сигнала с лучевым управлением для RX1
ss1 = sigstrength(rx1, tx, rtpm);
disp("Received power with beam steering at RX1: " + ss1 + " dBm")
pause(20);
% Для RX2
ray = raytrace(tx, rx2, rtpm);
disp(ray{1})

aod = ray{1}.AngleOfDeparture;
steeringaz = wrapTo180(aod(1) - tx.AntennaAngle(1));
steeringVector = phased.SteeringVector(SensorArray=tx.Antenna);
sv = steeringVector(tx.TransmitterFrequency, [steeringaz; aod(2)]);
tx.Antenna.Taper = conj(sv);

% Отображение на карте с направленным лучом для RX2
pattern(tx, Transparency=0.6)
raytrace(tx, rx2, rtpm);
hide(tx)

% Вычисление мощности полученного сигнала с лучевым управлением для RX2
ss2 = sigstrength(rx2, tx, rtpm);
disp("Received power with beam steering at RX2: " + ss2 + " dBm")

% Для RX3
ray = raytrace(tx, rx3, rtpm);
disp(ray{1})

aod = ray{1}.AngleOfDeparture;
steeringaz = wrapTo180(aod(1) - tx.AntennaAngle(1));
steeringVector = phased.SteeringVector(SensorArray=tx.Antenna);
sv = steeringVector(tx.TransmitterFrequency, [steeringaz; aod(2)]);
tx.Antenna.Taper = conj(sv);

% Отображение на карте с направленным лучом для RX3
pattern(tx, Transparency=0.6)
raytrace(tx, rx3, rtpm);
hide(tx)

% Вычисление мощности полученного сигнала с лучевым управлением для RX3
ss3 = sigstrength(rx3, tx, rtpm);
disp("Received power with beam steering at RX3: " + ss3 + " dBm")

% Для RX4
ray = raytrace(tx, rx4, rtpm);
disp(ray{1})

aod = ray{1}.AngleOfDeparture;
steeringaz = wrapTo180(aod(1) - tx.AntennaAngle(1));
steeringVector = phased.SteeringVector(SensorArray=tx.Antenna);
sv = steeringVector(tx.TransmitterFrequency, [steeringaz; aod(2)]);
tx.Antenna.Taper = conj(sv);

% Отображение на карте с направленным лучом для RX4
pattern(tx, Transparency=0.6)
raytrace(tx, rx4, rtpm);
hide(tx)

% Вычисление мощности полученного сигнала с лучевым управлением для RX4
ss4 = sigstrength(rx4, tx, rtpm);
disp("Received power with beam steering at RX4: " + ss4 + " dBm")
