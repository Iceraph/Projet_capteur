%% Nettoyage de l'environnement
clear; clc;

%% Données matériaux : nom, E [GPa], sigma_D [MPa]
materiaux = {
    'Acier Böhler',    196, 800;
    'Acier Maraging',  193, 735;
    'Alu 7075',         72, 110;
    'Alu Anticorodal',  69,  80;
    'Alu Avional',      73, 100;
    'Alu Contal',       72, 120;
    'Titane 6Al-4V',   114, 500;
    'Bronze',          127, 225
};

%% Constantes globales
k_tige = 1940;                 % [N/m]
k_table = 3026.5;              % [N/m]
F = 5;                         % [N]
delta = 1e-3;                  % [m]
r = 32e-3;                     % [m] Rayon de courbure
k_target = 5000 - k_tige - k_table;  % [N/m]
k_tol = 0.10;                  % ±10 %
k_min = k_target * (1 - k_tol);
k_max = k_target * (1 + k_tol);
h_err = 10e-6;                 % [m]
c_sec = 2;                     % Coefficient de sécurité
alpha = delta / r;            % [rad] Angle de rotation imposé

%% Pondération du score
w1 = 0.25;    % 1 / |L/h - 60|
w2 = 0.45;    % 1 / |k - k_target|
w3 = 0.15;    % 1 / |h - 200 µm|
w4 = 0.15;    % 1 / |b/h - 10|
epsilon = 1e-6;

%% Paramètres de balayage
h_vals = 130e-6 : 10e-6 : 1e-3;   % [m] Épaisseurs testées
rl_vals = 50 : 1 : 80;            % l = rl × h
rb_vals = 8 : 0.5 : 12;           % b = rb × h

%% Préallocation mémoire
max_solutions = size(materiaux, 1) * length(h_vals) * length(rl_vals) * length(rb_vals) / 10;
solutions = cell(ceil(max_solutions), 14);
sol_idx = 1;

%% Constante pour calcul raideur
coeff_k = 4 / r^2;

%% Balayage complet
for i = 1:size(materiaux, 1)
    nom = materiaux{i, 1};
    E = materiaux{i, 2} * 1e9;       % GPa → Pa
    sigmaD = materiaux{i, 3} * 1e6;  % MPa → Pa

    for h = h_vals
        h_cubed = h^3;

        for rl = rl_vals
            l = rl * h;
            if l > 0.040, continue, end  % Longueur maximale : 40 mm

            % Contraintes angulaires
            theta_max = (2 * sigmaD * l) / (E * h);
            if alpha > theta_max / c_sec, continue, end

            sigma_max = (alpha * h * E) / l;
            if sigma_max > sigmaD / c_sec, continue, end

            for rb = rb_vals
                b = rb * h;
                I = (b * h_cubed) / 12;

                F_crit = (pi^2 * E * I) / (0.7 *l)^2;
                if F > F_crit / c_sec, continue, end
            
                k = coeff_k * (E * I) / l;

                % Vérification de la raideur
                if k >= k_min && k <= k_max
                    delta_k = abs(k_target - k);

                    % Score pondéré
                    score = w1 / (abs(rl - 60) + epsilon) + ...
                            w2 / (delta_k + epsilon) + ...
                            w3 / (abs(h * 1e6 - 200) + epsilon) + ...
                            w4 / (abs(rb - 10) + epsilon);

                    % Stockage résultat
                    solutions(sol_idx, :) = {
                        nom, materiaux{i, 2}, sigmaD / 1e6, ...
                        h * 1e3, theta_max, b * 1e3, l * 1e3, ...
                        k, rl, rb, delta_k, score, sigma_max / 1e6, F_crit
                    };
                    sol_idx = sol_idx + 1;
                end
            end
        end
    end
end

%% Nettoyage des résultats vides
solutions = solutions(1:sol_idx - 1, :);

%% Affichage des résultats
if isempty(solutions)
    disp(' Aucun matériau ne satisfait les contraintes dans les plages définies.');
else
    % Création de la table
    T = cell2table(solutions, 'VariableNames', {
        'Alliage', 'E_GPa', 'SigmaD_MPa', ...
        'h_mm', 'theta_max', 'b_mm', 'l_mm', ...
        'k_N_per_m', 'rl', 'rb', ...
        'delta_k', 'Score', 'Sigma_max_MPa', 'Fcr_N'
    });

    % Tri par matériau puis score décroissant
    T_sorted = sortrows(T, {'Alliage', 'Score'}, {'ascend', 'descend'});

    % Sélection des 5 meilleures par matériau
    mat_list = unique(T_sorted.Alliage);
    meilleurs = table();

    for i = 1:length(mat_list)
        idx = strcmp(T_sorted.Alliage, mat_list{i});
        sous_table = T_sorted(idx, :);
        meilleurs = [meilleurs; sous_table(1:min(5, height(sous_table)), :)];
    end

    % Affichage
    disp('Meilleures solutions par matériau:');
    disp(meilleurs);

    fprintf('\nNombre total de solutions trouvées : %d\n', height(T));

    % Meilleur score global
    [best_score, best_idx] = max(T.Score);
    fprintf('\n Meilleure solution globale (score = %.4f):\n', best_score);
    disp(T(best_idx, :));
end
