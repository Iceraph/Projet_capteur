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
delta = 1e-3;                  % [m]
Ls = 200e-3;                   % [m]
alpha = delta / Ls;           % [rad]
R = 10e-3;                     % [m]
C = 29e-3;                     % [m]
L10 = 23e-3;                   % [m]
c_sec = 2;                     % Coefficient de sécurité
k_tiges = 238.98;             % [N/m]
k_target = 5000 - k_tiges;    % [N/m]
tol = 0.10;                   % ±10 %

% Initialisation des extrêmes
k_min = 2e7; k_max = 0;
min_name = ""; max_name = "";
min_d = 0; max_d = 0;
min_L = 0; max_L = 0;
min_rL = 0; max_rL = 0;

%% Paramètres de balayage
d_vals = 100e-6 : 100e-6 : 10e-3; % [m]
rL_vals = 20 : 80;

%% Préallocation des résultats
solutions = {};
idx = 1;
total_iter = size(materiaux,1) * numel(d_vals) * numel(rL_vals);
count = 0;

solutions_par_materiau = cell(size(materiaux, 1), 1);
for i = 1:size(materiaux, 1)
    solutions_par_materiau{i} = {};
end

fprintf('🔧 Recherche en cours...\n');

%% Boucles principales
for i = 1:size(materiaux, 1)
    nom = materiaux{i, 1};
    E = materiaux{i, 2} * 1e9;     % [Pa]
    sigmaD = materiaux{i, 3} * 1e6;% [Pa]
    G = materiaux{i, 4} * 1e9;     % [Pa]

    solutions_mat_idx = 1;

    for d = d_vals
        I = (pi / 64) * d^4;
        EI = E * I;
        GI = G * I;

        for rL = rL_vals
            count = count + 1;
            if mod(count, round(total_iter / 20)) == 0
                fprintf('Progression : %3.0f%%\n', 100 * count / total_iter);
            end

            L = round(rL * d, 3); % [m]

            if L > sqrt(C^2 - R^2), continue, end

            % Contraintes de flexibilité
            if c_sec * R * alpha >= (sigmaD * L^2) / (3 * E * d), continue, end
            if c_sec * alpha >= (sigmaD * L) / (E * d), continue, end
            if c_sec * alpha >= (sigmaD * L) / (G * d), continue, end

            % Contraintes combinées
            sigma_1 = (E * alpha * d) / (2 * L);
            sigma_2 = (3 * E * R * alpha * d) / L^2;
            sigma_max = max(sigma_1, sigma_2);
            sigma_tau = (G * d * alpha) / (2 * L);
            sigma_VM = sqrt(sigma_1^2 + 3 * sigma_tau^2);

            if sigma_max > sigmaD / c_sec, continue, end

            % Flambage
            Fcr = (pi^2 * E * I) / (0.5 * L)^2;
            if F > Fcr / 2, continue, end

            % Raideur totale
            term1 = EI * (8 / L^2 + 9 * R^2 / (2 * L^3));
            term2 = GI * 16 / L;
            term3 = EI * ((3 * R^2 / L10^3) + 2 / L10^2);
            k = (term1 + term2 + term3) / (Ls + R)^2;

            % Mise à jour des extrêmes
            if k < k_min
                k_min = k; min_name = nom; min_d = d * 1e3; min_L = L * 1e3; min_rL = rL;
            end
            if k > k_max
                k_max = k; max_name = nom; max_d = d * 1e3; max_L = L * 1e3; max_rL = rL;
            end

            % Enregistrement si proche de la cible
            if abs(k - k_target) <= k_target * tol
                delta_k = abs(k - k_target);

                % Liste globale
                solutions{idx, 1} = nom;
                solutions{idx, 2} = d * 1e3;    % [mm]
                solutions{idx, 3} = L * 1e3;    % [mm]
                solutions{idx, 4} = rL;
                solutions{idx, 5} = k;
                solutions{idx, 6} = delta_k;
                idx = idx + 1;

                % Liste par matériau
                solutions_par_materiau{i}{solutions_mat_idx, 1} = nom;
                solutions_par_materiau{i}{solutions_mat_idx, 2} = d * 1e3;
                solutions_par_materiau{i}{solutions_mat_idx, 3} = L * 1e3;
                solutions_par_materiau{i}{solutions_mat_idx, 4} = rL;
                solutions_par_materiau{i}{solutions_mat_idx, 5} = k;
                solutions_par_materiau{i}{solutions_mat_idx, 6} = delta_k;
                solutions_mat_idx = solutions_mat_idx + 1;
            end
        end
    end
end

%% Résultats globaux
if idx > 1
    fprintf('\n Résultats globaux :\n');
    T = cell2table(solutions(1:idx - 1, :), ...
        'VariableNames', {'Materiau', 'd_mm', 'L_mm', 'rL', 'k_N_per_m', 'delta_k'});
    T = sortrows(T, 'delta_k', 'ascend');
    disp(T);

    %% Top 5 par matériau
    fprintf('\n Top 5 par matériau :\n');
    for i = 1:size(materiaux, 1)
        nom_materiau = materiaux{i, 1};
        sol_mat = solutions_par_materiau{i};
        if ~isempty(sol_mat)
            Tm = cell2table(sol_mat, ...
                'VariableNames', {'Materiau', 'd_mm', 'L_mm', 'rL', 'k_N_per_m', 'delta_k'});
            Tm = sortrows(Tm, 'delta_k', 'ascend');
            top_n = min(5, height(Tm));
            fprintf('\n🔹 %s — Meilleures %d configurations :\n', nom_materiau, top_n);
            disp(Tm(1:top_n, :));
        end
    end
else
    disp('Aucune configuration ne respecte les contraintes.');
end

%% Statistiques de raideur
fprintf('\n Statistiques des raideurs extrêmes :\n');
fprintf('Raideur minimale : %.2f N/m [%s — d = %.2f mm, L = %.2f mm]\n', k_min, min_name, min_d, min_L);
fprintf('Raideur maximale : %.2f N/m [%s — d = %.2f mm, L = %.2f mm]\n', k_max, max_name, max_d, max_L);
fprintf('Raideur cible     : %.2f N/m\n', k_target);
