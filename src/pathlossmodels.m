% Параметры
c = 300000000; % скорость света, м/с
f_c = f_c_band; % частота, Гц (Band 3 LTE)
f_c_ghz = f_c / 1e9; % частота, ГГц
% h_BS = 25; % высота БС, м
% h_UT = 1.5; % высота терминала, м
bsData = struct();
bsData.CI = [138256652];
%bsData.bsLat = 55.013389;
%bsData.bsLon = 82.950622;
bsData.bsLat = bsLat; 
bsData.bsLon = bsLon; 

BsLat = bsData.bsLat;
BsLon = bsData.bsLon;
%antennaAzimuth = 180; % Азимут антенны (градусы)
%tilt = 5; % Наклон антенны (градусы)

P_Tx_dBm = p_tx_dBm_min:p_tx_dBm_max;
L_cable = [str2double(split(fider_loss, ',', 1))];

% Параметры для InF (Indoor Factory)
%h_BS_InF = 10; % Высота BS для InF, м
%h_UT_InF = 1.5; % Высота UE для InF, м
%clutterHeight_InF = 5; % Высота препятствий, м

data = readtable(data_path);

if ~all(ismember({'lat', 'lon', 'ci', 'rsrp'}, data.Properties.VariableNames))
    error('CSV-файл должен содержать столбцы: lat, lon, ci, rsrp');
end

validCI = ismember(data.ci, bsData.CI);
if sum(validCI) == 0
    error('Нет данных с CI %d в файле 138256652.csv', bsData.CI(1));
end
data = data(validCI, :);

validateData = validateDistance(data, BsLat, BsLon, 500);

latArr = validateData.lat;
lonArr = validateData.lon;

% Расчет углов и расстояний от базовой станции до пользователя
d_2D = arrayfun(@(lat,lon) haversine(BsLat, BsLon, lat, lon), latArr, lonArr);
d_3D = sqrt(d_2D.^2 + (h_BS - h_UT)^2);

ciArr = validateData.ci;
rsrpPhone = validateData.rsrp;

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

az_angle = arrayfun(@(lat,lon) azimuth(BsLat, BsLon, lat, lon), latArr, lonArr);
relativeAzimuth = mod(az_angle - antennaAzimuth, 360);
relativeAzimuth(relativeAzimuth > 180) = relativeAzimuth(relativeAzimuth > 180) - 360;
el_angle = rad2deg(atan2(h_UT - h_BS, d_2D));

angleToAntenna = zeros(numRows, 1);
for i = 1:numRows
    az_ant_rad = deg2rad(antennaAzimuth);
    tilt_rad = deg2rad(tilt);
    ant_vec = [cos(tilt_rad) * cos(az_ant_rad), cos(tilt_rad) * sin(az_ant_rad), sin(tilt_rad)];
    az_rad = deg2rad(az_angle(i));
    el_rad = deg2rad(el_angle(i));
    ue_vec = [cos(el_rad) * cos(az_rad), cos(el_rad) * sin(az_rad), sin(el_rad)];
    cos_theta = dot(ant_vec, ue_vec) / (norm(ant_vec) * norm(ue_vec));
    cos_theta = min(max(cos_theta, -1), 1);
    angleToAntenna(i) = rad2deg(acos(cos_theta));
end


azvec = -180:180; % Azimuth angles (deg)
elvec = -90:90; % Elevation angles (deg)
SLA = 30; % Maximum side-lobe level attenuation (dB)
tilt = 0; % Tilt angle (deg)
az3dB = 65; % 3 dB beamwidth in azimuth (deg)
el3dB = 65; % 3 dB beamwidth in elevation (deg)
lambda = physconst("lightspeed")/2.8e9; % Wavelength (m)

[az,el] = meshgrid(azvec,elvec);
azMagPattern = -min(12*(az/az3dB).^2,SLA);
elMagPattern = -min(12*((el-tilt)/el3dB).^2,SLA);
combinedMagPattern = -min(-(azMagPattern + elMagPattern), SLA); % Relative antenna gain (dB)


antennaElement = phased.CustomAntennaElement(...
    'MagnitudePattern', combinedMagPattern);

Antenna = phased.URA(Size=[2 2], ...
        Element=antennaElement, ...
        ElementSpacing=[lambda/2 lambda/2]);

directivePattern = pattern(Antenna, f_c);

maxdBi = max(directivePattern(:));
fprintf("Максимальное усиление: %.2f dBi\n", maxdBi);

% Расчет G_tx
G_tx = zeros(numRows, 1);
for i = 1:numRows
    G_tx(i) = pattern(antennaElement, f_c, relativeAzimuth(i), el_angle(i));
    % fprintf("G_tx(%d) = %f dBi\n", i, G_tx(i));
end

% Расчет path loss для всех сценариев (3GPP TR 38.901)

% 1. UMa (Urban Macrocell)
P_LOS_UMa = ones(numRows, 1);
idx = d_2D > 18;
P_LOS_UMa(idx) = (18 ./ d_2D(idx) + exp(-d_2D(idx)/63) .* (1 - 18 ./ d_2D(idx))) .* ...
    (1 + 5/4 * (h_UT - 1.5) * (d_2D(idx)/100).^3 .* exp(-d_2D(idx)/150));
d_BP_UMa = (4 * h_BS * h_UT * f_c) / c;
term_UMa = d_BP_UMa^2 + (h_BS - h_UT)^2;
PL_UMA_LOS = 28.0 + 20 * log10(f_c_ghz) + 22 * log10(d_3D);
idx_BP = d_2D > d_BP_UMa;
PL_UMA_LOS(idx_BP) = 28.0 + 20 * log10(f_c_ghz) + 40 * log10(d_3D(idx_BP)) - 9 * log10(term_UMa);
PL_UMA_NLOS_prime = 13.54 + 39.08 * log10(d_3D) + 20 * log10(f_c_ghz) - 0.6 * (h_UT - 1.5);
PL_UMA_NLOS = max(PL_UMA_LOS, PL_UMA_NLOS_prime);
% 2. RMa (Rural Macrocell)
h_BS_RMa = 35; % Высота BS для RMa, м
h_UT_RMa = 1.5; % Высота UE для RMa, м
d_3D_RMa = sqrt(d_2D.^2 + (h_BS_RMa - h_UT_RMa)^2);
d_BP_RMa = 2 * pi * h_BS_RMa * h_UT_RMa * f_c / c;
P_LOS_RMa = ones(numRows, 1);
idx_RMa = d_2D > 10;
P_LOS_RMa(idx_RMa) = exp(-(d_2D(idx_RMa) - 10) / 200);
PL_RMA_LOS = 20 * log10(40 * pi * d_3D_RMa * f_c_ghz / 3) + ...
    min(0.03 * 10^1.2, 10) * log10(d_3D_RMa) - ...
    min(0.044 * 10^1.2, 14.77) + 0.002 * log10(d_3D_RMa) .* d_3D_RMa;
idx_BP_RMa = d_2D > d_BP_RMa;
PL_RMA_LOS(idx_BP_RMa) = 20 * log10(40 * pi * d_BP_RMa * f_c_ghz / 3) + ...
    min(0.03 * 10^1.2, 10) * log10(d_BP_RMa) - ...
    min(0.044 * 10^1.2, 14.77) + 0.002 * log10(d_BP_RMa) .* d_BP_RMa + ...
    40 * log10(d_3D_RMa(idx_BP_RMa) / d_BP_RMa);
PL_RMA_NLOS = 161.04 - 7.1 * log10(10) + 7.5 * log10(10) - ...
    (24.37 - 3.7 * (10 / h_BS_RMa)^2) * log10(h_BS_RMa) + ...
    (43.42 - 3.1 * log10(h_BS_RMa)) * (log10(d_3D_RMa) - 3) + ...
    20 * log10(f_c_ghz) - (3.2 * (log10(11.75 * h_UT_RMa))^2 - 4.97);

% 3. UMi (Urban Microcell)
h_BS_UMi = 10; 
h_UT_UMi = 1.5;
d_3D_UMi = sqrt(d_2D.^2 + (h_BS_UMi - h_UT_UMi)^2);
P_LOS_UMi = ones(numRows, 1);
idx_UMi = d_2D > 18;
P_LOS_UMi(idx_UMi) = (18 ./ d_2D(idx_UMi) + exp(-d_2D(idx_UMi)/36) .* (1 - 18 ./ d_2D(idx_UMi)));
PL_UMI_LOS = 32.4 + 20 * log10(f_c_ghz) + 21 * log10(d_3D_UMi);
PL_UMI_NLOS = 22.4 + 35.3 * log10(d_3D_UMi) + 21.3 * log10(f_c_ghz) - 0.3 * (h_UT_UMi - 1.5);

% 4. InH (Indoor Hotspot)
d_3D_InH = d_2D;
P_LOS_InH = ones(numRows, 1);
idx_InH = d_2D > 1.2;
P_LOS_InH(idx_InH) = exp(-(d_2D(idx_InH) - 1.2) / 4.7);
PL_INH_LOS = 32.4 + 20 * log10(f_c_ghz) + 17 * log10(d_3D_InH);
PL_INH_NLOS = 147.4 - 43.3 * log10(10) + 20 * log10(f_c_ghz);

% 5. InF (Indoor Factory) - SL (Sparse Clutter, Low BS)
d_3D_InF = sqrt(d_2D.^2 + (h_BS_InF - h_UT_InF)^2);
P_LOS_InF = ones(numRows, 1);
idx_InF = d_2D > 1;
P_LOS_InF(idx_InF) = exp(-(d_2D(idx_InF) - 1) / 2.1);
PL_INF_SL_LOS = 31.84 + 20 * log10(f_c_ghz) + 21.5 * log10(d_3D_InF);
PL_INF_SL_NLOS = max(PL_INF_SL_LOS, 33.63 + 20 * log10(f_c_ghz) + 25.1 * log10(d_3D_InF) + ...
    20 * log10(1 + clutterHeight_InF / 10));

% % Библиотечные расчеты (если доступен 5G Toolbox)
% PL_UMA_LOS_lib = zeros(numRows, 1);
% PL_UMA_NLOS_lib = zeros(numRows, 1);
% PL_RMA_LOS_lib = zeros(numRows, 1);
% PL_RMA_NLOS_lib = zeros(numRows, 1);
% PL_UMI_LOS_lib = zeros(numRows, 1);
% PL_UMI_NLOS_lib = zeros(numRows, 1);
% PL_INH_LOS_lib = zeros(numRows, 1);
% PL_INH_NLOS_lib = zeros(numRows, 1);
% PL_INF_SL_LOS_lib = zeros(numRows, 1);
% PL_INF_SL_NLOS_lib = zeros(numRows, 1);

% if exist('nrPathLoss', 'file')
%     cfg_UMa = nrPathLossConfig('Scenario', 'UMa');
%     cfg_RMa = nrPathLossConfig('Scenario', 'RMa');
%     cfg_UMi = nrPathLossConfig('Scenario', 'UMi-StreetCanyon');
%     cfg_InH = nrPathLossConfig('Scenario', 'InH-Office');
%     cfg_InF = nrPathLossConfig('Scenario', 'InF-SL');
%     for i = 1:numRows
%         [x_bs, y_bs] = latlon2local(BsLat, BsLon, latArr(1), lonArr(1));
%         [x_ue, y_ue] = latlon2local(latArr(i), lonArr(i), latArr(1), lonArr(1));
%         if isnan(x_ue) || isnan(y_ue) || isnan(x_bs) || isnan(y_bs) || ~isnumeric(x_ue) || ~isnumeric(y_ue)
%             warning('Некорректные координаты для строки %d: x_ue=%f, y_ue=%f, x_bs=%f, y_bs=%f', ...
%                 i, x_ue, y_ue, x_bs, y_bs);
%             PL_UMA_LOS_lib(i) = NaN;
%             PL_UMA_NLOS_lib(i) = NaN;
%             PL_RMA_LOS_lib(i) = NaN;
%             PL_RMA_NLOS_lib(i) = NaN;
%             PL_UMI_LOS_lib(i) = NaN;
%             PL_UMI_NLOS_lib(i) = NaN;
%             PL_INH_LOS_lib(i) = NaN;
%             PL_INH_NLOS_lib(i) = NaN;
%             PL_INF_SL_LOS_lib(i) = NaN;
%             PL_INF_SL_NLOS_lib(i) = NaN;
%             continue;
%         end
%         posBS_UMa = double([x_bs; y_bs; h_BS]);
%         posUE_UMa = double([x_ue; y_ue; h_UT]);
%         posBS_RMa = double([x_bs; y_bs; h_BS_RMa]);
%         posUE_RMa = double([x_ue; y_ue; h_UT_RMa]);
%         posBS_UMi = double([x_bs; y_bs; h_BS_UMi]);
%         posUE_UMi = double([x_ue; y_ue; h_UT_UMi]);
%         posBS_InH = double([x_bs; y_bs; h_UT]); % BS и UE на одной высоте
%         posUE_InH = double([x_ue; y_ue; h_UT]);
%         posBS_InF = double([x_bs; y_bs; h_BS_InF]);
%         posUE_InF = double([x_ue; y_ue; h_UT_InF]);
        
%         if any(islogical(posUE_UMa))  any(islogical(posBS_UMa))
%             warning('Логический тип в posUE или posBS для строки %d', i);
%             PL_UMA_LOS_lib(i) = NaN;
%             PL_UMA_NLOS_lib(i) = NaN;
%             PL_RMA_LOS_lib(i) = NaN;
%             PL_RMA_NLOS_lib(i) = NaN;
%             PL_UMI_LOS_lib(i) = NaN;
%             PL_UMI_NLOS_lib(i) = NaN;
%             PL_INH_LOS_lib(i) = NaN;
%             PL_INH_NLOS_lib(i) = NaN;
%             PL_INF_SL_LOS_lib(i) = NaN;
%             PL_INF_SL_NLOS_lib(i) = NaN;
%             continue;
%         end
        
%         [PL_UMA_LOS_lib(i), ~] = nrPathLoss(cfg_UMa, f_c, true, posBS_UMa, posUE_UMa);
%         [PL_UMA_NLOS_lib(i), ~] = nrPathLoss(cfg_UMa, f_c, false, posBS_UMa, posUE_UMa);
%         [PL_RMA_LOS_lib(i), ~] = nrPathLoss(cfg_RMa, f_c, true, posBS_RMa, posUE_RMa);
%         [PL_RMA_NLOS_lib(i), ~] = nrPathLoss(cfg_RMa, f_c, false, posBS_RMa, posUE_RMa);
%         [PL_UMI_LOS_lib(i), ~] = nrPathLoss(cfg_UMi, f_c, true, posBS_UMi, posUE_UMi);
%         [PL_UMI_NLOS_lib(i), ~] = nrPathLoss(cfg_UMi, f_c, false, posBS_UMi, posUE_UMi);
%         [PL_INH_LOS_lib(i), ~] = nrPathLoss(cfg_InH, f_c, true, posBS_InH, posUE_InH);
%         [PL_INH_NLOS_lib(i), ~] = nrPathLoss(cfg_InH, f_c, false, posBS_InH, posUE_InH);
%         [PL_INF_SL_LOS_lib(i), ~] = nrPathLoss(cfg_InF, f_c, true, posBS_InF, posUE_InF);
%         [PL_INF_SL_NLOS_lib(i), ~] = nrPathLoss(cfg_InF, f_c, false, posBS_InF, posUE_InF);
%     end
% else
%     disp('5G Toolbox недоступен. Библиотечные расчеты RSRP будут пропущены.');
%     PL_UMA_LOS_lib = NaN(numRows, 1);
%     PL_UMA_NLOS_lib = NaN(numRows, 1);
%     PL_RMA_LOS_lib = NaN(numRows, 1);
%     PL_RMA_NLOS_lib = NaN(numRows, 1);
%     PL_UMI_LOS_lib = NaN(numRows, 1);
%     PL_UMI_NLOS_lib = NaN(numRows, 1);
%     PL_INH_LOS_lib = NaN(numRows, 1);
%     PL_INH_NLOS_lib = NaN(numRows, 1);
%     PL_INF_SL_LOS_lib = NaN(numRows, 1);
%     PL_INF_SL_NLOS_lib = NaN(numRows, 1);
% end

[PTx, LC] = ndgrid(P_Tx_dBm, L_cable);
combinations = [PTx(:), LC(:)];
numComb = size(combinations, 1);

allFullTables = cell(numComb, 1);
meanRows = cell(numComb, 1);

if exist('parfor', 'file')
    parfor i = 1:numComb
        P_Tx_dBm_local = combinations(i, 1);
        L_cable_local = combinations(i, 2);
        
        % Ширина полосы B (в Гц) для расчета RSRP
        B = 1.4e6; % Предположим, что B уже задано

        % Изменение формулы расчета RSRP
        rsrp_uma_los = P_Tx_dBm_local - L_cable_local - PL_UMA_LOS + G_tx - 10 * log10(B);
        rsrp_uma_nlos = P_Tx_dBm_local - L_cable_local - PL_UMA_NLOS + G_tx - 10 * log10(B);
        rsrp_rma_los = P_Tx_dBm_local - L_cable_local - PL_RMA_LOS + G_tx - 10 * log10(B);
        rsrp_rma_nlos = P_Tx_dBm_local - L_cable_local - PL_RMA_NLOS + G_tx - 10 * log10(B);
        rsrp_umi_los = P_Tx_dBm_local - L_cable_local - PL_UMI_LOS + G_tx - 10 * log10(B);
        rsrp_umi_nlos = P_Tx_dBm_local - L_cable_local - PL_UMI_NLOS + G_tx - 10 * log10(B);
        rsrp_inh_los = P_Tx_dBm_local - L_cable_local - PL_INH_LOS + G_tx - 10 * log10(B);
        rsrp_inh_nlos = P_Tx_dBm_local - L_cable_local - PL_INH_NLOS + G_tx - 10 * log10(B);
        rsrp_inf_sl_los = P_Tx_dBm_local - L_cable_local - PL_INF_SL_LOS + G_tx - 10 * log10(B);
        rsrp_inf_sl_nlos = P_Tx_dBm_local - L_cable_local - PL_INF_SL_NLOS + G_tx - 10 * log10(B);

                % if ~all(isnan(PL_UMA_LOS_lib))
        %     rsrp_uma_los_lib = P_Tx_dBm_local - L_cable_local - PL_UMA_LOS_lib + G_tx;
        %     rsrp_uma_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_UMA_NLOS_lib + G_tx;
        %     rsrp_rma_los_lib = P_Tx_dBm_local - L_cable_local - PL_RMA_LOS_lib + G_tx;
        %     rsrp_rma_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_RMA_NLOS_lib + G_tx;
        %     rsrp_umi_los_lib = P_Tx_dBm_local - L_cable_local - PL_UMI_LOS_lib + G_tx;
        %     rsrp_umi_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_UMI_NLOS_lib + G_tx;
        %     rsrp_inh_los_lib = P_Tx_dBm_local - L_cable_local - PL_INH_LOS_lib + G_tx;
        %     rsrp_inh_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_INH_NLOS_lib + G_tx;
        %     rsrp_inf_sl_los_lib = P_Tx_dBm_local - L_cable_local - PL_INF_SL_LOS_lib + G_tx;
        %     rsrp_inf_sl_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_INF_SL_NLOS_lib + G_tx;
        % else
        %     rsrp_uma_los_lib = NaN(numRows, 1);
        %     rsrp_uma_nlos_lib = NaN(numRows, 1);
        %     rsrp_rma_los_lib = NaN(numRows, 1);
        %     rsrp_rma_nlos_lib = NaN(numRows, 1);
        %     rsrp_umi_los_lib = NaN(numRows, 1);
        %     rsrp_umi_nlos_lib = NaN(numRows, 1);
        %     rsrp_inh_los_lib = NaN(numRows, 1);
        %     rsrp_inh_nlos_lib = NaN(numRows, 1);
        %     rsrp_inf_sl_los_lib = NaN(numRows, 1);
        %     rsrp_inf_sl_nlos_lib = NaN(numRows, 1);
        % end

        % Разницы между измеренным и расчетным RSRP
        diff_uma_los = rsrpPhone - rsrp_uma_los;
        diff_uma_nlos = rsrpPhone - rsrp_uma_nlos;
        diff_rma_los = rsrpPhone - rsrp_rma_los;
        diff_rma_nlos = rsrpPhone - rsrp_rma_nlos;
        diff_umi_los = rsrpPhone - rsrp_umi_los;
        diff_umi_nlos = rsrpPhone - rsrp_umi_nlos;
        diff_inh_los = rsrpPhone - rsrp_inh_los;
        diff_inh_nlos = rsrpPhone - rsrp_inh_nlos;
        diff_inf_sl_los = rsrpPhone - rsrp_inf_sl_los;
        diff_inf_sl_nlos = rsrpPhone - rsrp_inf_sl_nlos;
        
        tempFullTable = table(...
            repmat(P_Tx_dBm_local, numRows, 1), ...
            repmat(L_cable_local, numRows, 1), ...
            latArr, ...
            lonArr, ...
            ciArr, ...
            d_2D, ...
            relativeAzimuth, ...
            el_angle, ...
            G_tx, ...
            rsrpPhone, ...
            rsrp_uma_los, rsrp_uma_nlos, ...
            rsrp_rma_los, rsrp_rma_nlos, ...
            rsrp_umi_los, rsrp_umi_nlos, ...
            rsrp_inh_los, rsrp_inh_nlos, ...
            rsrp_inf_sl_los, rsrp_inf_sl_nlos, ...
            diff_uma_los, diff_uma_nlos, ...
            diff_rma_los, diff_rma_nlos, ...
            diff_umi_los, diff_umi_nlos, ...
            diff_inh_los, diff_inh_nlos, ...
                diff_inf_sl_los, diff_inf_sl_nlos, ...
                    'VariableNames', {...
                        'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                        'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone', ...
                        'rsrp_uma_los', 'rsrp_uma_nlos', ...
                        'rsrp_rma_los', 'rsrp_rma_nlos', ...
                        'rsrp_umi_los', 'rsrp_umi_nlos', ...
                        'rsrp_inh_los', 'rsrp_inh_nlos', ...
                        'rsrp_inf_sl_los', 'rsrp_inf_sl_nlos', ...
                        'diff_uma_los', 'diff_uma_nlos', ...
                        'diff_rma_los', 'diff_rma_nlos', ...
                        'diff_umi_los', 'diff_umi_nlos', ...
                        'diff_inh_los', 'diff_inh_nlos', ...
                        'diff_inf_sl_los', 'diff_inf_sl_nlos'});
                
        allFullTables{i} = tempFullTable;
        
        mean_diff_uma_nlos = mean(diff_uma_nlos, 'omitnan');
        mean_diff_rma_nlos = mean(diff_rma_nlos, 'omitnan');
        mean_diff_umi_nlos = mean(diff_umi_nlos, 'omitnan');
        mean_diff_inh_nlos = mean(diff_inh_nlos, 'omitnan');
        mean_diff_inf_sl_nlos = mean(diff_inf_sl_nlos, 'omitnan');
        mean_diff_uma_los = mean(diff_uma_los, 'omitnan');
        mean_diff_rma_los = mean(diff_rma_los, 'omitnan');
        mean_diff_umi_los = mean(diff_umi_los, 'omitnan');
        mean_diff_inh_los = mean(diff_inh_los, 'omitnan');
        mean_diff_inf_sl_los = mean(diff_inf_sl_los, 'omitnan');
        
        meanRows{i} = [P_Tx_dBm_local, L_cable_local, ...
            mean_diff_uma_nlos, mean_diff_rma_nlos, mean_diff_umi_nlos, ...
            mean_diff_inh_nlos, mean_diff_inf_sl_nlos, ...
            mean_diff_uma_los, mean_diff_rma_los, mean_diff_umi_los, ...
            mean_diff_inh_los, mean_diff_inf_sl_los];

            end
        else
            for i = 1:numComb
                P_Tx_dBm_local = combinations(i, 1);
                L_cable_local = combinations(i, 2);
                
                % Расчет RSRP для всех сценариев0
                rsrp_uma_los = P_Tx_dBm_local - L_cable_local - PL_UMA_LOS + G_tx; % для LOS
                rsrp_uma_nlos = P_Tx_dBm_local - L_cable_local - PL_UMA_NLOS + G_tx; % для NLOS
                rsrp_rma_los = P_Tx_dBm_local - L_cable_local - PL_RMA_LOS + G_tx;
                rsrp_rma_nlos = P_Tx_dBm_local - L_cable_local - PL_RMA_NLOS + G_tx;
                rsrp_umi_los = P_Tx_dBm_local - L_cable_local - PL_UMI_LOS + G_tx;
                rsrp_umi_nlos = P_Tx_dBm_local - L_cable_local - PL_UMI_NLOS + G_tx;
                rsrp_inh_los = P_Tx_dBm_local - L_cable_local - PL_INH_LOS + G_tx;
                rsrp_inh_nlos = P_Tx_dBm_local - L_cable_local - PL_INH_NLOS + G_tx;
                rsrp_inf_sl_los = P_Tx_dBm_local - L_cable_local - PL_INF_SL_LOS + G_tx;
                rsrp_inf_sl_nlos = P_Tx_dBm_local - L_cable_local - PL_INF_SL_NLOS + G_tx;

                % if ~all(isnan(PL_UMA_LOS_lib))
                %     rsrp_uma_los_lib = P_Tx_dBm_local - L_cable_local - PL_UMA_LOS_lib + G_tx;
                %     rsrp_uma_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_UMA_NLOS_lib + G_tx;
                %     rsrp_rma_los_lib = P_Tx_dBm_local - L_cable_local - PL_RMA_LOS_lib + G_tx;
                %     rsrp_rma_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_RMA_NLOS_lib + G_tx;
                %     rsrp_umi_los_lib = P_Tx_dBm_local - L_cable_local - PL_UMI_LOS_lib + G_tx;
                %     rsrp_umi_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_UMI_NLOS_lib + G_tx;
                %     rsrp_inh_los_lib = P_Tx_dBm_local - L_cable_local - PL_INH_LOS_lib + G_tx;
                %     rsrp_inh_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_INH_NLOS_lib + G_tx;
                %     rsrp_inf_sl_los_lib = P_Tx_dBm_local - L_cable_local - PL_INF_SL_LOS_lib + G_tx;
                %     rsrp_inf_sl_nlos_lib = P_Tx_dBm_local - L_cable_local - PL_INF_SL_NLOS_lib + G_tx;
                % else
                %     rsrp_uma_los_lib = NaN(numRows, 1);
                %     rsrp_uma_nlos_lib = NaN(numRows, 1);
                %     rsrp_rma_los_lib = NaN(numRows, 1);
                %     rsrp_rma_nlos_lib = NaN(numRows, 1);
                %     rsrp_umi_los_lib = NaN(numRows, 1);
                %     rsrp_umi_nlos_lib = NaN(numRows, 1);
                %     rsrp_inh_los_lib = NaN(numRows, 1);
                %     rsrp_inh_nlos_lib = NaN(numRows, 1);
                %     rsrp_inf_sl_los_lib = NaN(numRows, 1);
                %     rsrp_inf_sl_nlos_lib = NaN(numRows, 1);
                % end
                % Разницы между измеренным и расчетным RSRP
        diff_uma_los = rsrpPhone - rsrp_uma_los;
        diff_uma_nlos = rsrpPhone - rsrp_uma_nlos;
        diff_rma_los = rsrpPhone - rsrp_rma_los;
        diff_rma_nlos = rsrpPhone - rsrp_rma_nlos;
        diff_umi_los = rsrpPhone - rsrp_umi_los;
        diff_umi_nlos = rsrpPhone - rsrp_umi_nlos;
        diff_inh_los = rsrpPhone - rsrp_inh_los;
        diff_inh_nlos = rsrpPhone - rsrp_inh_nlos;
        diff_inf_sl_los = rsrpPhone - rsrp_inf_sl_los;
        diff_inf_sl_nlos = rsrpPhone - rsrp_inf_sl_nlos;
        
        tempFullTable = table(...
            repmat(P_Tx_dBm_local, numRows, 1), ...
            repmat(L_cable_local, numRows, 1), ...
            latArr, ...
            lonArr, ...
            ciArr, ...
            d_2D, ...
            relativeAzimuth, ...
            el_angle, ...
            G_tx, ...
            rsrpPhone, ...
            rsrp_uma_los, rsrp_uma_nlos, ...
            rsrp_rma_los, rsrp_rma_nlos, ...
            rsrp_umi_los, rsrp_umi_nlos, ...
            rsrp_inh_los, rsrp_inh_nlos, ...
            rsrp_inf_sl_los, rsrp_inf_sl_nlos, ...
            rsrp_inf_sl_los_lib, rsrp_inf_sl_nlos_lib, ...
            diff_uma_los, diff_uma_nlos, ...
            diff_rma_los, diff_rma_nlos, ...
            diff_umi_los, diff_umi_nlos, ...
            diff_inh_los, diff_inh_nlos, ...
            diff_inf_sl_los, diff_inf_sl_nlos, ...
            'VariableNames', {...
                'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone', ...
                'rsrp_uma_los', 'rsrp_uma_nlos', ...
                'rsrp_rma_los', 'rsrp_rma_nlos', ...
                'rsrp_umi_los', 'rsrp_umi_nlos', ...
                'rsrp_inh_los', 'rsrp_inh_nlos', ...
                'rsrp_inf_sl_los', 'rsrp_inf_sl_nlos', ...
                'diff_uma_los', 'diff_uma_nlos', ...
                'diff_rma_los', 'diff_rma_nlos', ...
                'diff_umi_los', 'diff_umi_nlos', ...
                'diff_inh_los', 'diff_inh_nlos', ...
                'diff_inf_sl_los', 'diff_inf_sl_nlos'});
        
        allFullTables{i} = tempFullTable;
        
        mean_diff_uma_nlos = mean(diff_uma_nlos, 'omitnan');
        mean_diff_rma_nlos = mean(diff_rma_nlos, 'omitnan');
        mean_diff_umi_nlos = mean(diff_umi_nlos, 'omitnan');
        mean_diff_inh_nlos = mean(diff_inh_nlos, 'omitnan');
        mean_diff_inf_sl_nlos = mean(diff_inf_sl_nlos, 'omitnan');
        mean_diff_uma_los = mean(diff_uma_los, 'omitnan');
        mean_diff_rma_los = mean(diff_rma_los, 'omitnan');
        mean_diff_umi_los = mean(diff_umi_los, 'omitnan');
        mean_diff_inh_los = mean(diff_inh_los, 'omitnan');
        mean_diff_inf_sl_los = mean(diff_inf_sl_los, 'omitnan');
        
        meanRows{i} = [P_Tx_dBm_local, L_cable_local, ...
            mean_diff_uma_nlos, mean_diff_rma_nlos, mean_diff_umi_nlos, ...
            mean_diff_inh_nlos, mean_diff_inf_sl_nlos, ...
            mean_diff_uma_los, mean_diff_rma_los, mean_diff_umi_los, ...
            mean_diff_inh_los, mean_diff_inf_sl_los];
    end
end

fullTable = vertcat(allFullTables{:});

meanMatrix = cell2mat(meanRows);
meanTable = array2table(meanMatrix, 'VariableNames', {...
    'P_Tx_dBm', 'L_cable', 'mean_diff_uma_nlos', 'mean_diff_rma_nlos', ...
    'mean_diff_umi_nlos', 'mean_diff_inh_nlos', 'mean_diff_inf_sl_nlos', ...
    'mean_diff_uma_los', 'mean_diff_rma_los', ...
    'mean_diff_umi_los', 'mean_diff_inh_los', 'mean_diff_inf_sl_los'});

fprintf('\nВсе значения mean_diff для каждой модели:\n');
fprintf('P_Tx_dBm | L_cable | UMa_NLOS | RMa_NLOS | UMi_NLOS | InH_NLOS | InF_SL_NLOS |  UMa_LOS |  RMa_LOS |  UMi_LOS |  InH_LOS |  InF_SL_LOS\n');
fprintf('---------|---------|----------|----------|----------|----------|-------------|----------|----------|----------|----------|-------------\n');
for i = 1:height(meanTable)
    fprintf('%.2f     | %.2f    | %.6f | %.6f | %.6f | %.6f | %.6f| %.6f | %.6f | %.6f | %.6f | %.6f\n', ...
        meanTable.P_Tx_dBm(i), meanTable.L_cable(i), ...
        meanTable.mean_diff_uma_nlos(i), meanTable.mean_diff_rma_nlos(i), ...
        meanTable.mean_diff_umi_nlos(i), meanTable.mean_diff_inh_nlos(i), ...
        meanTable.mean_diff_inf_sl_nlos(i), ...
        meanTable.mean_diff_uma_los(i), meanTable.mean_diff_rma_los(i), ...
        meanTable.mean_diff_umi_los(i), meanTable.mean_diff_inh_los(i), ...
        meanTable.mean_diff_inf_sl_los(i));
end


models = {'uma_nlos', 'rma_nlos', 'umi_nlos', 'inh_nlos', 'inf_sl_nlos'};
for i = 1:length(models)
    model_name = models{i};

    model_columns = {'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                     'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone'};
    model_columns = [model_columns, sprintf('rsrp_%s', model_name)];

    model_columns = [model_columns, sprintf('diff_%s', model_name)];
    
    model_full_table = fullTable(:, model_columns);
    
    filename = sprintf('full_data_%s.csv', model_name);
    writetable(model_full_table, filename, 'Delimiter', ';');
    fprintf('Сохранена полная таблица для модели %s в файл %s\n', upper(model_name), filename);
end


% Оптимальные параметры P_Tx_dBm и L_cable для каждой модели
optimal_params = struct();
for i = 1:length(models)
    model_name = models{i};
    col_name = sprintf('mean_diff_%s', model_name);
    if all(isnan(meanTable.(col_name))) || isempty(meanTable.(col_name))
        warning('Для модели %s нет валидных данных mean_diff. Используются значения по умолчанию.', upper(model_name));
        optimal_params.(model_name) = struct('P_Tx_dBm', NaN, 'L_cable', NaN);
        continue;
    end
    [~, idx] = min(abs(meanTable.(col_name)));
    if isempty(idx) || idx > height(meanTable)
        warning('Для модели %s не удалось найти минимальное mean_diff. Используются значения по умолчанию.', upper(model_name));
        optimal_params.(model_name) = struct('P_Tx_dBm', NaN, 'L_cable', NaN);
        continue;
    end
    optimal_params.(model_name) = struct('P_Tx_dBm', meanTable.P_Tx_dBm(idx), 'L_cable', meanTable.L_cable(idx));
end

% Создание и сохранение отфильтрованных таблиц для каждой модели
for i = 1:length(models)
    model_name = models{i};

    opt_P_Tx_dBm = optimal_params.(model_name).P_Tx_dBm;
    opt_L_cable = optimal_params.(model_name).L_cable;
    
    if isnan(opt_P_Tx_dBm) || isnan(opt_L_cable)
        warning('Пропуск сохранения таблицы для модели %s из-за отсутствия валидных параметров.', upper(model_name));
        continue;
    end
    
    idx = fullTable.P_Tx_dBm == opt_P_Tx_dBm & fullTable.L_cable == opt_L_cable;
    filtered_table = fullTable(idx, :);
    
    model_columns = {'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                     'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone'};
    model_columns = [model_columns, sprintf('rsrp_%s', model_name)];

    model_columns = [model_columns, sprintf('diff_%s', model_name)];
    
    model_filtered_table = filtered_table(:, model_columns);
    
    filename = sprintf('filtered_data_%s.csv', model_name);
    writetable(model_filtered_table, filename, 'Delimiter', ';');
    fprintf('Сохранена отфильтрованная таблица для модели %s в файл %s (P_Tx_dBm=%.2f, L_cable=%.2f)\n', ...
        upper(model_name), filename, opt_P_Tx_dBm, opt_L_cable);
end

% Поиск минимальных mean_diff для каждой модели
models = {'uma_nlos', 'rma_nlos', 'umi_nlos', 'inh_nlos', 'inf_sl_nlos', 'uma_los', 'rma_los', 'umi_los', 'inh_los', 'inf_sl_los'};
fprintf('\nМинимальные mean_diff для каждой модели:\n');
fprintf('Модель    | Min mean_diff | P_Tx_dBm | L_cable\n');
fprintf('----------|---------------|----------|---------\n');
for i = 1:length(models)
    model_name = models{i};
    col_name = sprintf('mean_diff_%s', model_name);
    if all(isnan(meanTable.(col_name))) || isempty(meanTable.(col_name))
        fprintf('%-9s | NaN           | NaN      | NaN\n', upper(model_name));
        continue;
    end
    [min_mean_diff, idx] = min(abs(meanTable.(col_name)));
    if isempty(idx) || idx > height(meanTable)
        fprintf('%-9s | NaN           | NaN      | NaN\n', upper(model_name));
        continue;
    end
    fprintf('%-9s | %.6f      | %.2f     | %.2f\n', ...
        upper(model_name), min_mean_diff, meanTable.P_Tx_dBm(idx), meanTable.L_cable(idx));

    optimal_params.model_name = struct('');

end

%%%%%

models = {'uma_los', 'rma_los', 'umi_los', 'inh_los', 'inf_sl_los'};
for i = 1:length(models)
    model_name = models{i};

    model_columns = {'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                     'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone'};
    model_columns = [model_columns, sprintf('rsrp_%s', model_name)];

    model_columns = [model_columns, sprintf('diff_%s', model_name)];
    
    model_full_table = fullTable(:, model_columns);
    
    % Сохраняем таблицу в CSV-файл
    filename = sprintf('full_data_%s.csv', model_name);
    writetable(model_full_table, filename, 'Delimiter', ';');
    fprintf('Сохранена полная таблица для модели %s в файл %s\n', upper(model_name), filename);
end


% Оптимальные параметры P_Tx_dBm и L_cable для каждой модели
optimal_params = struct();
for i = 1:length(models)
    model_name = models{i};
    col_name = sprintf('mean_diff_%s', model_name);
    if all(isnan(meanTable.(col_name))) || isempty(meanTable.(col_name))
        warning('Для модели %s нет валидных данных mean_diff. Используются значения по умолчанию.', upper(model_name));
        optimal_params.(model_name) = struct('P_Tx_dBm', NaN, 'L_cable', NaN);
        continue;
    end
    [~, idx] = min(abs(meanTable.(col_name)));
    if isempty(idx) || idx > height(meanTable)
        warning('Для модели %s не удалось найти минимальное mean_diff. Используются значения по умолчанию.', upper(model_name));
        optimal_params.(model_name) = struct('P_Tx_dBm', NaN, 'L_cable', NaN);
        continue;
    end
    optimal_params.(model_name) = struct('P_Tx_dBm', meanTable.P_Tx_dBm(idx), 'L_cable', meanTable.L_cable(idx));
end

for i = 1:length(models)
    model_name = models{i};

    opt_P_Tx_dBm = optimal_params.(model_name).P_Tx_dBm;
    opt_L_cable = optimal_params.(model_name).L_cable;
    
    if isnan(opt_P_Tx_dBm) || isnan(opt_L_cable)
        warning('Пропуск сохранения таблицы для модели %s из-за отсутствия валидных параметров.', upper(model_name));
        continue;
    end
    
    % Фильтруем fullTable по P_Tx_dBm и L_cable
    idx = fullTable.P_Tx_dBm == opt_P_Tx_dBm & fullTable.L_cable == opt_L_cable;
    filtered_table = fullTable(idx, :);
    
    % Выбираем релевантные столбцы
    model_columns = {'P_Tx_dBm', 'L_cable', 'lat', 'lon', 'ci', 'distance', ...
                     'relativeAzimuth', 'el_angle', 'gain', 'rsrpPhone'};
    model_columns = [model_columns, sprintf('rsrp_%s', model_name)];

    model_columns = [model_columns, sprintf('diff_%s', model_name)];
    
    model_filtered_table = filtered_table(:, model_columns);
    
    filename = sprintf('filtered_data_%s.csv', model_name);
    writetable(model_filtered_table, filename, 'Delimiter', ';');
    fprintf('Сохранена отфильтрованная таблица для модели %s в файл %s (P_Tx_dBm=%.2f, L_cable=%.2f)\n', ...
        upper(model_name), filename, opt_P_Tx_dBm, opt_L_cable);
end

%%%%%

los_files = {'filtered_data_inf_sl_los.csv', 'filtered_data_inh_los.csv', ...
             'filtered_data_rma_los.csv', 'filtered_data_uma_los.csv', ...
             'filtered_data_umi_los.csv'};
nlos_files = {'filtered_data_inf_sl_nlos.csv', 'filtered_data_inh_nlos.csv', ...
              'filtered_data_rma_nlos.csv', 'filtered_data_uma_nlos.csv', ...
              'filtered_data_umi_nlos.csv'};

function plot_rsrp(files, title_str, output_file, mode)
    
    model_names = containers.Map(...
    {'inf_sl_los', 'inf_sl_nlos', 'inh_los', 'inh_nlos', 'rma_los', 'rma_nlos', ...
     'uma_los', 'uma_nlos', 'umi_los', 'umi_nlos'}, ...
    {'InF-SL LOS', 'InF-SL NLOS', 'InH LOS', 'InH NLOS', 'RMa LOS', 'RMa NLOS', ...
     'UMa LOS', 'UMa NLOS', 'UMi LOS', 'UMi NLOS'});

    first_data = readtable(files{1}, 'Delimiter', ';');
    first_data = sortrows(first_data, 'distance');
    
    colors = {'r', 'g', 'm', 'c', 'k'};

    distance = first_data.distance;
    rsrp_phone = first_data.rsrpPhone;
    
    fig = figure('Position', [100, 100, 800, 480]);
    hold on;
    
    plot(distance, rsrp_phone, 'b-o', 'MarkerSize', 4, 'DisplayName', 'RSRP Phone');
    
    for i = 1:length(files)
        data = readtable(files{i}, 'Delimiter', ';');
        data = sortrows(data, 'distance');
        
        model_key = strrep(strrep(files{i}, 'filtered_data_', ''), '.csv', '');
        model_rsrp_col = sprintf('rsrp_%s', model_key);
        
        if ~ismember(model_rsrp_col, data.Properties.VariableNames)
            warning('Столбец %s не найден в файле %s', model_rsrp_col, files{i});
            continue;
        end
        
        rsrp_model = data.(model_rsrp_col);
        plot(distance, rsrp_model, sprintf('%s-x', colors{i}), 'MarkerSize', 4, ...
             'DisplayName', model_names(model_key));
    end
    
    xlim([min(distance) max(distance)]);
    xlabel('Distance (meters)');
    ylabel('RSRP (dBm)');
    title(sprintf('RSRP vs Distance (%s)', mode));
    grid on;
    
    legend('show');
    
    datacursormode on;
    
    hold off;
end

plot_rsrp(los_files, 'RSRP vs Distance for LOS Models', 'rsrp_plot_los.png', 'LOS');
plot_rsrp(nlos_files, 'RSRP vs Distance for NLOS Models', 'rsrp_plot_nlos.png', 'NLOS');


%%%%%

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

function validateTable = validateDistance(data, BsLat, BsLon, maxDistance)
    if ~all(ismember({'lat', 'lon'}, data.Properties.VariableNames))
        error('Таблица данных должна содержать столбцы ''lat'' и ''lon''.');
    end
    
    latArr = data.lat;
    lonArr = data.lon;
    
    d_2d = arrayfun(@(lat, lon) haversine(BsLat, BsLon, lat, lon), latArr, lonArr);
    
    validDistance = d_2d <= maxDistance;
    
    validateTable = data(validDistance, :);
end
