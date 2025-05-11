%% Nettoyage de l'environnement
clear; clc;

%% Données matériaux : nom, E [GPa], sigma_D [MPa]
materiaux = {
    'Acier Böhler K190',    196, 800;
    'Acier Maraging W720',  193, 735;
    'Alu 7075',              72, 110;
    'Alu Anticorodal',       69,  80;
    'Alu Avional',           73, 100;
    'Alu Contal',            72, 120;
    'Titane 6Al-4V',        114, 500;
    'Bronze Pfinodal',      127, 225
};

%% Constantes globales
F = 5;                         % [N]
delta = 1e-3;                  % [m]
l = 10.8e-3;                   % [m] Longueur des tables
c_sec = 2;                     % Coefficient de sécurité
r_vals = 25e-3 : 1e-3 : 32e-3; % [m] Rayon de courbure (25 à 32 mm)
d_vals = 100e-6 : 10e-6 : 1e-3;% [m] Diamètre (100 µm à 1 mm)

%% Préallocation mémoire
max_solutions = length(materiaux) * length(d_vals) * length(r_vals) * 31; % pour Ld de 50 à 80
solutions = cell(max_solutions, 9);
idx = 1;

%% Constantes intermédiaires
k_const = (3 * delta / (5 * l))^2 * 3;
delta_l = (3 / (5 * l)) * delta^2;

%% Balayage
for i = 1:size(materiaux, 1)
    nom = materiaux{i, 1};
    E = materiaux{i, 2} * 1e9;       % GPa → Pa
    sigmaD = materiaux{i, 3} * 1e6;  % MPa → Pa

    for d = d_vals
        I = pi * d^4 / 64;   % Moment d'inertie
        EI = E * I;

        for r = r_vals
            for Ld = 50 : 80
                L = Ld * d;

                if L > 0.070, continue, end % Longueur max : 70 mm

                % Flambage
                Fcr = pi^2 * EI / (0.7 * L)^2;
                if F > Fcr / c_sec, continue, end

                % Déformations max
                if c_sec * delta_l >= (2 * sigmaD * L^2) / (3 * E * d), continue, end
                if c_sec * delta / r >= (sigmaD * L) / (E * d), continue, end

                % Contrainte maximale
                M = EI * (3 * delta_l / L^2 + 2 * delta / (L * r));
                sigma_max = (M * d / 2) / I;
                if sigma_max > sigmaD / c_sec, continue, end

                % Raideur
                k_z_in = 2 * (k_const * EI / L^3);
                if k_z_in > 50, continue, end

                k = k_z_in + (4 * EI) / ((r * L)^2);
                if k > 5000, continue, end

                % Stockage de la solution
                solutions(idx, :) = {
                    nom, d * 1e3, L * 1e3, r * 1e2, ...
                    Fcr, Ld, k, sigma_max / 1e6, k_z_in
                };
                idx = idx + 1;
            end
        end
    end
end

%% Nettoyage des solutions vides
solutions = solutions(1:idx - 1, :);

%% Création du tableau
T = cell2table(solutions, 'VariableNames', {
    'Matériau', 'd_mm', 'L_mm', 'r_cm', ...
    'Fcr_N', 'Ld', 'k_N_per_m', 'Sigma_max_MPa', 'k_(zin)_N_per_m'
});

%% Tri et statistiques
[min_k, min_idx] = min(T.k_N_per_m);
solution_min_k = T(min_idx, :);

% Tri par diamètre croissant
T = sortrows(T, 'd_mm');

%% Affichage
disp('Solutions triées par diamètre :');
disp(T);

fprintf('\nNombre total de solutions : %d\n', height(T));
fprintf('\n🏆 Solution avec raideur minimale (k = %.2f N/m):\n', min_k);
disp(solution_min_k);
