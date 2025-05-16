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
delta = 1.008e-3;                  % [m]
Ls = 200e-3;                   % [m]
alpha = delta / Ls;           % [rad]
R = 8.5e-3;                     % [m]
C = 29e-3;                     % [m]
L12 = 23e-3;                   % [m]
c_sec = 2;                     % Coefficient de sécurité
k_tiges = 457.27;             % [N/m]
k_target = 5000 - k_tiges;    % [N/m]
tol = 0.10;                   % ±10 %

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

        % Contraintes combinées tige 12
        sigma_1_12 = (E * alpha * d) / (2 * L12);
        sigma_2_12 = (3 * E * R * alpha * d) / (L12)^2;
        sigma_max_12 = max(sigma_1_12, sigma_2_12);
        sigma_tau_12 = (G * d * alpha) / (2 * L12);
        sigma_VM_12 = sqrt(sigma_max_12^2 + 3 * sigma_tau_12^2);

        if sigma_VM_12 > sigmaD / c_sec, continue, end

        % Contraintes de flexibilité tige 12
        delta_adm_12 = (sigmaD * L12^2) / (3 * E * d);
        alpha_adm_12 = min((sigmaD * L12) / (E * d), (sigmaD * L12) / (G * d));
        if c_sec * R * alpha >= delta_adm_12, continue, end
        if c_sec * alpha >= alpha_adm_12, continue, end

        for rL = rL_vals
            count = count + 1;
            if mod(count, round(total_iter / 20)) == 0
                fprintf('Progression : %3.0f%%\n', 100 * count / total_iter);
            end

            L = round(rL * d, 3); % [m]

            if L > sqrt(C^2 - R^2), continue, end

            % Contraintes de flexibilité
            delta_adm = (sigmaD * L^2) / (3 * E * d);
            alpha_adm = min((sigmaD * L) / (E * d), (sigmaD * L) / (G * d));
            if c_sec * sqrt(2)/2 * R * alpha >= delta_adm, continue, end
            if c_sec * alpha >= alpha_adm, continue, end

            % Contraintes combinées
            sigma_1 = (E * alpha * d) / (2 * L);
            sigma_2 = (3 * E * R * alpha * d) / L^2;
            sigma_max = max(sigma_1, sigma_2);
            sigma_tau = (G * d * alpha) / (2 * L);
            sigma_VM = sqrt(sigma_max^2 + 3 * sigma_tau^2);

            if sigma_VM > sigmaD / c_sec, continue, end

            % Flambage
            Fcr = (pi^2 * E * I) / (0.5 * L)^2;
            if F > Fcr / 2, continue, end

            % Raideur totale
            term1 = EI * (8 / L^2 + 9 * R^2 / (2 * L^3));
            term2 = GI * 16 / L;
            term3 = EI * ((3 * R^2 / L12^3) + 2 / L12^2);
            k = (term1 + term2 + term3) / (Ls + R)^2;

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
                solutions{idx, 7} = sigma_VM / 1e6; % [MPa]
                solutions{idx, 8} = Fcr;
                solutions{idx, 9} = delta_adm * 1e6; % [µm]
                solutions{idx, 10} = alpha_adm;
                idx = idx + 1;

                % Liste par matériau
                solutions_par_materiau{i}{solutions_mat_idx, 1} = nom;
                solutions_par_materiau{i}{solutions_mat_idx, 2} = d * 1e3;
                solutions_par_materiau{i}{solutions_mat_idx, 3} = L * 1e3;
                solutions_par_materiau{i}{solutions_mat_idx, 4} = rL;
                solutions_par_materiau{i}{solutions_mat_idx, 5} = k;
                solutions_par_materiau{i}{solutions_mat_idx, 6} = delta_k;
                solutions_par_materiau{i}{solutions_mat_idx, 7} = sigma_VM / 1e6;
                solutions_par_materiau{i}{solutions_mat_idx, 8} = Fcr;
                solutions_par_materiau{i}{solutions_mat_idx, 9} = delta_adm * 1e6;
                solutions_par_materiau{i}{solutions_mat_idx, 10} = alpha_adm;
                solutions_mat_idx = solutions_mat_idx + 1;
            end
        end
    end
end

%% Résultats globaux
if idx > 1
    fprintf('\n Résultats globaux :\n');
    T = cell2table(solutions(1:idx - 1, :), ...
        'VariableNames', {'Materiau', 'd_mm', 'L_mm', 'rL', 'k_N_per_m', ...
        'delta_k', 'Sigma_VM_MPa', 'Fcr_N', 'delta_adm_µm', 'alpha_adm'});
    T = sortrows(T, 'delta_k', 'ascend');
    disp(T);

    %% Top 5 par matériau
    fprintf('\n Top 5 par matériau :\n');
    for i = 1:size(materiaux, 1)
        nom_materiau = materiaux{i, 1};
        sol_mat = solutions_par_materiau{i};
        if ~isempty(sol_mat)
            Tm = cell2table(sol_mat, ...
                'VariableNames', {'Materiau', 'd_mm', 'L_mm', 'rL', 'k_N_per_m', ...
                'delta_k', 'Sigma_VM_MPa', 'Fcr_N', 'delta_adm', 'alpha_adm'});
            Tm = sortrows(Tm, 'delta_k', 'ascend');
            top_n = min(5, height(Tm));
            fprintf('\n🔹 %s — Meilleures %d configurations :\n', nom_materiau, top_n);
            disp(Tm(1:top_n, :));
        end
    end
else
    disp('Aucune configuration ne respecte les contraintes.');
end

