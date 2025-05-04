clc; clear;

% Параметры
c = 3e8; % скорость света, м/с
f_c = 1.8725e9; % частота, Гц (Band 3 LTE)
f_c_ghz = f_c / 1e9; % частота, ГГц
h_BS = 25; % высота БС, м
h_UT = 1.5; % высота терминала, м
bsData = struct();
bsData.CI = [138256652];
BsLat = 55.013389;
BsLon = 82.950622;

data = readtable('138256652.csv');

if ~all(ismember({'lat', 'lon', 'ci', 'rsrp'}, data.Properties.VariableNames))
    error('CSV-файл должен содержать столбцы: lat, lon, ci, rsrp');
end

validCI = ismember(data.ci, bsData.CI);
if sum(validCI) == 0
    error('Нет данных с CI %d в файле 138256652.csv', bsData.CI(1));
end
data = data(validCI, :);

latArr = data.lat;
lonArr = data.lon;
ciArr = data.ci;
rsrpPhone = data.rsrp;

if ~isnumeric(latArr) || ~isnumeric(lonArr) || ~isnumeric(rsrpPhone)
    error('latArr, lonArr и rsrpPhone должны быть числовыми массивами');
end
validRows = ~isnan(latArr) & ~isnan(lonArr) & ~isnan(rsrpPhone);
if ~any(validRows)
    error('Все строки содержат NaN в lat, lon или rsrp');
end
latArr = double(latArr(validRows));
lonArr = double(lonArr(validRows));
ciArr = double(ciArr(validRows));
rsrpPhone = double(rsrpPhone(validRows));
numRows = length(latArr);

% Calc of angles and distances from the base station to the user
d_2D = arrayfun(@(lat,lon) haversine(BsLat, BsLon, lat, lon), latArr, lonArr);
d_3D = sqrt(d_2D.^2 + (h_BS - h_UT)^2);
az_angle = arrayfun(@(lat,lon) azimuth(BsLat, BsLon, lat, lon), latArr, lonArr);
relativeAzimuth = mod(az_angle - antennaAzimuth, 360); % Относительный азимут
relativeAzimuth(relativeAzimuth > 180) = relativeAzimuth(relativeAzimuth > 180) - 360; % Приведение к [-180, 180]
el_angle = rad2deg(atan2(h_BS - h_UT, d_2D));

%preliminary calculation for the Uma model, taking into account the distance
P_LOS = ones(numRows, 1);
idx = d_2D > 18;
P_LOS(idx) = (18 ./ d_2D(idx) + exp(-d_2D(idx)/63) .* (1 - 18 ./ d_2D(idx))) .* ...
    (1 + 5/4 * (h_UT - 1.5) * (d_2D(idx)/100).^3 .* exp(-d_2D(idx)/150));
    d_BP = (2 * pi * h_BS * h_UT * f_c) / c;

term = d_BP^2 + (h_BS - h_UT)^2;
PL_LOS_base = 28.0 + 20 * log10(f_c_ghz) + 22 * log10(d_3D);
idx_BP = d_2D > d_BP;
PL_LOS_base(idx_BP) = 28.0 + 20 * log10(f_c_ghz) + 40 * log10(d_3D(idx_BP)) - 9 * log10(term);

PL_NLOS_prime = 13.54 + 39.08 * log10(d_3D) + 20 * log10(f_c_ghz) - 0.6 * (h_UT - 1.5);
PL_NLOS_base = max(PL_LOS_base, PL_NLOS_prime);



PL_LOS_lib_base = zeros(numRows, 1);
PL_NLOS_lib_base = zeros(numRows, 1);

if exist('nrPathLoss', 'file')
    cfg = nrPathLossConfig('Scenario', 'UMa');
    for i = 1:numRows
        [x_bs, y_bs] = latlon2local(BsLat, BsLon, latArr(1), lonArr(1));
        [x_ue, y_ue] = latlon2local(latArr(i), lonArr(i), latArr(1), lonArr(1));
        if isnan(x_ue) || isnan(y_ue) || isnan(x_bs) || isnan(y_bs) || ~isnumeric(x_ue) || ~isnumeric(y_ue)
            warning('Некорректные координаты для строки %d: x_ue=%f, y_ue=%f, x_bs=%f, y_bs=%f', ...
                i, x_ue, y_ue, x_bs, y_bs);
            PL_LOS_lib_base(i) = NaN;
            PL_NLOS_lib_base(i) = NaN;
            continue;
        end
        posBS = double([x_bs; y_bs; h_BS]);
        posUE = double([x_ue; y_ue; h_UT]);
        if any(islogical(posUE)) || any(islogical(posBS))
            warning('Логический тип в posUE или posBS для строки %d', i);
            PL_LOS_lib_base(i) = NaN;
            PL_NLOS_lib_base(i) = NaN;
            continue;
        end
        [PL_LOS_lib_base(i), ~] = nrPathLoss(cfg, f_c, true, posBS, posUE); % LOS
        [PL_NLOS_lib_base(i), ~] = nrPathLoss(cfg, f_c, false, posBS, posUE); % NLOS
    end
else
    disp('5G Toolbox недоступен. Библиотечные расчеты RSRP будут пропущены.');
    PL_LOS_lib_base = NaN(numRows, 1);
    PL_NLOS_lib_base = NaN(numRows, 1);
end

P_Tx_dBm = 20:40;
L_cable = [3];
shadow_std = [6]; % from 3gpp TR38.901 - 7.4.1 page-28

[PTx, LC, SS] = ndgrid(P_Tx_dBm, L_cable, shadow_std);
combinations = [PTx(:), LC(:), SS(:)];
numComb = size(combinations, 1);

allFullTables = cell(numComb, 1);
meanRows = cell(numComb, 1);

if exist('parfor', 'file')
    parfor i = 1:numComb
        P_Tx_dBm_local = combinations(i, 1);
        L_cable_local = combinations(i, 2);
        shadow_std_local = combinations(i, 3);
        
        %main calc for Uma model

        rsrpLos = P_Tx_dBm_local - L_cable_local - PL_LOS_base;
        PL_NLOS = PL_NLOS_base + normrnd(0, shadow_std_local, numRows, 1);
        rsrpNlos = P_Tx_dBm_local - L_cable_local - PL_NLOS;
        
        if ~all(isnan(PL_LOS_lib_base))
            PL_NLOS_lib = PL_NLOS_lib_base + normrnd(0, shadow_std_local, numRows, 1);
            rsrpLosLib = P_Tx_dBm_local - L_cable_local - PL_LOS_lib_base;
            rsrpNlosLib = P_Tx_dBm_local - L_cable_local - PL_NLOS_lib;
        else
            rsrpLosLib = NaN(numRows, 1);
            rsrpNlosLib = NaN(numRows, 1);
        end
        
        diff_nlos = rsrpPhone - rsrpNlos; 
        diff_nlos_lib = rsrpPhone - rsrpNlosLib;
        
        tempFullTable = table(...
            repmat(P_Tx_dBm_local, numRows, 1), ...
            repmat(L_cable_local, numRows, 1), ...
            repmat(shadow_std_local, numRows, 1), ...
            latArr, ...
            lonArr, ...
            ciArr, ...
            d_2D, ...
            rsrpPhone, ...
            rsrpLos, ...
            rsrpNlos, ...
            rsrpLosLib, ...
            rsrpNlosLib, ...
            diff_nlos, ...
            diff_nlos_lib, ...
            'VariableNames', {'P_Tx_dBm', 'L_cable', 'shadow_std', 'lat', 'lon', 'ci', 'distance', 'rsrpPhone', 'rsrpLos', 'rsrpNlos', 'rsrpLosLib', 'rsrpNlosLib', 'diff_nlos', 'diff_nlos_lib'});
        
        allFullTables{i} = tempFullTable;
        
        mean_diff_nlos = mean(diff_nlos, 'omitnan');
        mean_diff_nlos_lib = mean(diff_nlos_lib, 'omitnan');
        
        meanRows{i} = [P_Tx_dBm_local, L_cable_local, shadow_std_local, mean_diff_nlos, mean_diff_nlos_lib];
    end
else
    for i = 1:numComb
        P_Tx_dBm_local = combinations(i, 1);
        L_cable_local = combinations(i, 2);
        shadow_std_local = combinations(i, 3);
        
        rsrpLos = P_Tx_dBm_local + G_tx + 0 - L_cable_local - PL_LOS_base;
        PL_NLOS = PL_NLOS_base + normrnd(0, shadow_std_local, numRows, 1);
        rsrpNlos = P_Tx_dBm_local + G_tx + 0 - L_cable_local - PL_NLOS;
        
        if ~all(isnan(PL_LOS_lib_base))
            PL_NLOS_lib = PL_NLOS_lib_base + normrnd(0, shadow_std_local, numRows, 1);
            rsrpLosLib = P_Tx_dBm_local + G_tx + 0 - L_cable_local - PL_LOS_lib_base;
            rsrpNlosLib = P_Tx_dBm_local + G_tx + 0 - L_cable_local - PL_NLOS_lib;
        else
            rsrpLosLib = NaN(numRows, 1);
            rsrpNlosLib = NaN(numRows, 1);
        end
        
        diff_nlos = rsrpPhone - rsrpNlos;
        diff_nlos_lib = rsrpPhone - rsrpNlosLib;
        
        tempFullTable = table(...
            repmat(P_Tx_dBm_local, numRows, 1), ...
            repmat(L_cable_local, numRows, 1), ...
            repmat(shadow_std_local, numRows, 1), ...
            latArr, ...
            lonArr, ...
            ciArr, ...
            d_2D, ...
            rsrpPhone, ...
            rsrpLos, ...
            rsrpNlos, ...
            rsrpLosLib, ...
            rsrpNlosLib, ...
            diff_nlos, ...
            diff_nlos_lib, ...
            'VariableNames', {'P_Tx_dBm', 'L_cable', 'shadow_std', 'lat', 'lon', 'ci', 'distance', 'rsrpPhone', 'rsrpLos', 'rsrpNlos', 'rsrpLosLib', 'rsrpNlosLib', 'diff_nlos', 'diff_nlos_lib'});
        
        allFullTables{i} = tempFullTable;
        
        mean_diff_nlos = mean(diff_nlos, 'omitnan');
        mean_diff_nlos_lib = mean(diff_nlos_lib, 'omitnan');
        
        meanRows{i} = [P_Tx_dBm_local, L_cable_local, shadow_std_local, mean_diff_nlos, mean_diff_nlos_lib];
    end
end

fullTable = vertcat(allFullTables{:});

meanMatrix = cell2mat(meanRows);
meanTable = array2table(meanMatrix, 'VariableNames', {'P_Tx_dBm', 'L_cable', 'shadow_std', 'mean_diff_nlos', 'mean_diff_nlos_lib'});

[min_mean_diff, idx] = min(abs(meanTable.mean_diff_nlos));

expected_Tx_dBm = meanTable.P_Tx_dBm(idx);
expected_L_min = meanTable.L_cable(idx);

fprintf("Smallest mean_diff_nlos - %.6f | Expected P_Tx_dBm = %.6f, Feeder_Loss - %.6f\n", min_mean_diff, expected_Tx_dBm, expected_L_min)

writetable(fullTable, '138256652_full_data.csv', 'Delimiter', ';');
writetable(meanTable, '138256652_mean_differences.csv', 'Delimiter', ';');

function dist = haversine(lat1, lon1, lat2, lon2)
    R = 6371e3;
    lat1 = deg2rad(double(lat1));
    lon1 = deg2rad(double(lon1));
    lat2 = deg2rad(double(lat2));
    lon2 = deg2rad(double(lon2));
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    dist = R * c;
end

function [x, y] = latlon2local(lat, lon, lat0, lon0)
    R = 6371e3;
    lat = deg2rad(double(lat));
    lon = deg2rad(double(lon));
    lat0 = deg2rad(double(lat0));
    lon0 = deg2rad(double(lon0));
    x = R * (lon - lon0) * cos(lat0);
    y = R * (lat - lat0);
    if ~isnumeric(x) || ~isnumeric(y)
        error('latlon2local вернул нечисловые значения: x=%s, y=%s', class(x), class(y));
    end
end

function az = azimuth(lat1, lon1, lat2, lon2)
    lat1 = deg2rad(double(lat1));
    lon1 = deg2rad(double(lon1));
    lat2 = deg2rad(double(lat2));
    lon2 = deg2rad(double(lon2));
    dlon = lon2 - lon1;
    y = sin(dlon) * cos(lat2);
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon);
    az = atan2(y, x);
    az = rad2deg(az);
    if az < 0
        az = az + 360;
    end
end