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
k_tige = 1940;                  % [N/m] Raideur propre de la tige
F = 5;                          % [N] Force appliquée
nb_table = 4;                   % Nombre de tables (non utilisé)
delta = 1e-3;                   % [m] Déplacement visé
k_target = 5000 - k_tige;       % [N/m] Raideur à atteindre
k_tol = 0.10;                   % Tolérance ±10 %
k_min = k_target * (1 - k_tol); % Limite inférieure
k_max = k_target * (1 + k_tol); % Limite supérieure
h_err = 10e-6;                  % [m] Tolérance d'erreur sur h
c_sec = 2;                      % Coefficient de sécurité

%% Pondération des critères de score
w1 = 0.35;    % Poids : compacité (1 / (L/h))
w2 = 0.35;    % Poids : précision sur raideur (1 / |k - k_target|)
w3 = 0.20;    % Poids : encombrement (1 / L)
w4 = 0.10;    % Poids : proportion géométrique (1 / |b/h - 10|)
epsilon = 1e-6;  % Petite valeur pour éviter division par zéro

%% Paramètres de balayage
h_vals = 130e-6 : 10e-6 : 1e-3;  % [m] Épaisseurs testées
rl_vals = 50 : 1 : 80;           % l = rl × h
rb_vals = 8 : 0.5 : 12;          % b = rb × h

%% Préallocation mémoire
max_solutions = size(materiaux, 1) * length(h_vals) * length(rl_vals) * length(rb_vals) / 10;
solutions = cell(ceil(max_solutions), 13);
sol_idx = 1;

%% Constante pour le calcul de la raideur
coeff_k = 96;

%% Balayage complet des combinaisons
for i = 1:size(materiaux, 1)
    nom = materiaux{i, 1};
    E = materiaux{i, 2} * 1e9;       % Conversion GPa → Pa
    sigmaD = materiaux{i, 3} * 1e6;  % Conversion MPa → Pa
    
    h_min_coef = (3 * E * delta) / sigmaD;
    
    for h = h_vals
        h_cubed = h^3;
        
        for rl = rl_vals
            l = rl * h;
            if l > 0.020, continue, end  % Longueur maximale 20 mm

            h_min = h_min_coef / rl^2;
            if h <= h_min + h_err, continue, end  % Contraintes déplacement

            l_cubed = l^3;
            k_factor = coeff_k * E / l_cubed;

            for rb = rb_vals
                b = rb * h;
                I = (b * h_cubed) / 12;
                k = k_factor * I;

                % Contraintes mécaniques
                sigma_max = (3 * E * delta * h) / (2 * l^2);
                if sigma_max > sigmaD / c_sec, continue, end

                F_crit = (pi^2 * E * I) / l^2;
                if F/4 > F_crit / 2, continue, end

                % Vérification raideur
                if k >= k_min && k <= k_max
                    delta_k = abs(k - k_target);

                    % Score pondéré
                    score = w1 / rl + ...
                            w2 / (delta_k + epsilon) + ...
                            w3 / l + ...
                            w4 / (abs(rb - 10) + epsilon);
                    
                    % Stockage résultat
                    solutions(sol_idx, :) = {
                        nom, materiaux{i, 2}, sigmaD / 1e6, ...
                        h * 1e3, h_min * 1e3, b * 1e3, l * 1e3, ...
                        k, rl, rb, delta_k, score, sigma_max / 1e6
                    };
                    sol_idx = sol_idx + 1;

                    % Agrandissement dynamique si besoin
                    if sol_idx > size(solutions, 1)
                        solutions = [solutions; cell(max_solutions / 2, 13)];
                    end
                end
            end
        end
    end
end

%% Nettoyage des résultats vides
solutions = solutions(1:sol_idx - 1, :);

%% Affichage des résultats
if isempty(solutions)
    disp('Aucun matériau ne satisfait les contraintes dans les plages définies.');
else
    % Table finale avec noms explicites
    T = cell2table(solutions, 'VariableNames', {
        'Alliage', 'E_GPa', 'SigmaD_MPa', ...
        'h_mm', 'h_min_mm', 'b_mm', 'l_mm', ...
        'k_N_per_m', 'rl', 'rb', ...
        'delta_k', 'Score', 'Sigma_max_MPa'
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
    fprintf('\nMeilleure solution globale (score = %.4f):\n', best_score);
    disp(T(best_idx, :));
end
