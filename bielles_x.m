%% Nettoyage de l’environnement
clear; clc;

%% Données matériaux : nom, E [GPa], sigma_D [MPa], G [GPa]
materiaux = {
    'Acier Böhler K190',    196, 800, 80;
    'Acier Maraging W720',  193, 735, 72;
    'Alu 7075',              72, 110, 27;
    'Alu Anticorodal',       69,  80, 26;
    'Alu Avional',           73, 100, 28;
    'Alu Contal',            72, 120, 27;
    'Titane 6Al-4V',        114, 500, 41;
    'Bronze Pfinodal',      127, 225, 44
};

%% Constantes globales
F = 5;                         % [N]
c_sec = 2;                     % Coefficient de sécurité
delta = 1.008e-3;                  % [m]
Ls = 200e-3;                   % [m]
alpha = delta / Ls;           % [rad]
R = 10e-3;                     % [m]
H1 = 45e-3;                    % [m]                  
coeff_geom = 12 * (2 * H1)^2;  % [m²] Facteur géométrique

k_target = 500;              % [N/m]
tol = 0.5;                   % ±10 %

k_min = 2e7;                  % Init raideur min [N/m]
k_max = 0;                    % Init raideur max [N/m]
min_name = ""; max_name = "";
min_d = 0; max_d = 0;
min_L = 0; max_L = 0;
min_rL = 0; max_rL = 0;
min_F = 0; max_F = 0;
min_delta = 0; max_delta = 0;
min_alpha = 0; max_alpha = 0;
min_sigma_d = 0; max_sigma_d = 0;
min_sigma_a = 0; max_sigma_a = 0;

%% Paramètres de balayage
d_vals = 100e-6 : 10e-6 : 400e-6;  % [m]
rL_vals = 20: 1 : 80;                     % [—] Ratio L/d

%% Préallocation des résultats
solutions = {};
idx = 1;

total_iter = size(materiaux,1) * numel(d_vals) * numel(rL_vals);
progress_step = round(total_iter / 20);
count = 0;

fprintf('📐 Dimensionnement tige 2 en cours...\n');

%% Balayage matériaux et géométries
for i = 1:size(materiaux, 1)
    nom = materiaux{i, 1};
    sigmaD = materiaux{i, 3} * 1e6;  % MPa → Pa
    E = materiaux{i, 2} * 1e9;       % GPa → Pa

    for d = d_vals
        I = (pi / 64) * d^4;
        EI = E * I;

        for rL = rL_vals
            count = count + 1;
            if mod(count, progress_step) == 0
                fprintf('Progression : %3.0f%%\n', 100 * count / total_iter);
            end

            L = rL * d;
            L = round(L, 3);  % [m]
            if L > 2 * R, continue, end
            if L ~= 12*10^-3, continue, end %Contrainte de modifications limité du 3D

            % Contrainte déplacement
            delta_adm = (2 * sigmaD * L^2) / (3 * E * d);
            if (2 * H1 * alpha) > delta_adm / c_sec, continue, end

            % Contrainte angulaire
            alpha_adm = (2 * sigmaD * L) / (pi^2 * E * d);
            if alpha > alpha_adm / c_sec, continue, end

            % Contrainte max rotation
            M_alpha = (EI * pi^2 * alpha) / (2 * L);
            sigma_alpha = (M_alpha * d) / (2 * I);
            if sigma_alpha > sigmaD / c_sec, continue, end

            % Contrainte max translation
            M_delta = (3 * EI / L^3) * (2 * H1 * alpha) * L;
            sigma_delta = (M_delta * d) / (2 * I);
            if sigma_delta > sigmaD / c_sec, continue, end

            % Flambage
            Fcr = (pi^2 * EI) / L^2;
            if F > Fcr / c_sec, continue, end

            % Calcul de la raideur
            k_geom = (coeff_geom * EI / L^3);
            k_rot = (EI * pi^2 / L^2);
            k = (k_geom + k_rot) / (Ls + R)^2;

            % Mise à jour des extrêmes
            if k < k_min
                k_min = k;
                min_name = nom;
                min_d = d * 1e3;  % [mm]
                min_L = L * 1e3;  % [mm]
                min_rL = rL;
                min_F = Fcr;
                min_delta = delta_adm * 1e3; % [mm]
                min_alpha = alpha_adm;
                min_sigma_d = sigma_delta / 1e6;   % [MPa]
                min_sigma_a = sigma_alpha / 1e6;   % [MPa]
            end
            if k > k_max
                k_max = k;
                max_name = nom;
                max_d = d * 1e3;
                max_L = L * 1e3;
                max_rL = rL;
                max_F = Fcr;
                max_delta = delta_adm * 1e3;   % [mm]
                max_alpha = alpha_adm;
                max_sigma_d = sigma_delta / 1e6; % [MPa]
                max_sigma_a = sigma_alpha / 1e6; % [MPa]
            end

            % Filtrage des solutions proches de la cible
            if abs(k - k_target) <= k_target * tol
                solutions{idx, 1} = nom;
                solutions{idx, 2} = d * 1e3;   % [mm]
                solutions{idx, 3} = L * 1e3;   % [mm]
                solutions{idx, 4} = rL;
                solutions{idx, 5} = k;
                solutions{idx, 6} = Fcr;   
                solutions{idx, 7} = delta_adm * 1e3;   % [mm]
                solutions{idx, 8} = alpha_adm;
                solutions{idx, 9} = sigma_delta / 1e6; % [MPa]
                solutions{idx, 10} = sigma_alpha / 1e6; % [MPa]
                idx = idx + 1;
            end
        end
    end
end

%% Affichage des résultats
if idx > 1
    T = cell2table(solutions(1:idx - 1, :), ...
        'VariableNames', {'Materiau', 'd_mm', 'L_mm', 'rL', 'k_N_per_m', ...
        'Fcr_N', 'delta_adm_mm', 'alpha_adm_rad', ...
        'sigma_delta_MPa', 'sigma_alpha_MPa'});

    disp('Configurations respectant la tolérance :');
    disp(T);
else
    disp('Aucune configuration ne respecte les contraintes.');
end

%% Résumé des extrêmes
fprintf('\n Statistiques extrêmes :\n');
fprintf('Raideur minimale : %.2f N/m [%s, d = %.2f mm, L = %.2f mm, Fcr = %.2f N, delta = %.2f mm, alpha = %.2f rad, sigma_d = %.2f MPa, sigma_a = %.2f MPa,]\n', ...
    k_min, min_name, min_d, min_L, min_F, min_delta, min_alpha, min_sigma_d, min_sigma_a);
fprintf('Raideur maximale : %.2f N/m [%s, d = %.2f mm, L = %.2f mm, Fcr = %.2f N, delta = %.2f mm, alpha = %.2f rad, sigma_d = %.2f MPa, sigma_a = %.2f MPa,]\n', ...
    k_max, max_name, max_d, max_L, max_F, max_delta, max_alpha, max_sigma_d, max_sigma_a);
fprintf('Raideur cible     : %.2f N/m\n', k_target);
