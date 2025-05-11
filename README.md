# Dimensionnement d’Éléments Flexibles pour un Capteur de Force

Projet académique de BA4 – *Conception de Mécanismes*  
**École Polytechnique Fédérale de Lausanne (EPFL)**

---

## 📌 Contexte

Ce projet s’inscrit dans le cadre du cours **Conception de Mécanismes (BA4)** à l’**EPFL**, dont l’objectif est de 
concevoir et dimensionner un **capteur de force à géométrie variable**. Les éléments flexibles sont calculés pour garantir la rigidité, 
et la sécurité mécanique dans ses différentes configurations.

Chaque fichier correspond à un sous-système ou un cas de charge particulier analysé et optimisé en autonomie.

---

## 🎯 Objectif

Ce projet contient une série d'outils MATLAB permettant de **dimensionner les éléments flexibles d’un capteur de force**. 
Ces éléments, soumis à différentes sollicitations mécaniques (translation, rotation, flexion, flambement...), 
sont optimisés en fonction de matériaux, géométries, et contraintes de sécurité pour répondre à un cahier des charges précis.

---

## 📁 Structure du projet

- `table_z.m`  
  → Dimensionnement d’une table à lame parallèle soumise à **flexion en translation** selon une raideur cible.

- `pivot_z.m`  
  → Dimensionnement d'un pivot à lame séparées soumis à une **rotation imposée** selon une raideur cible imposé par table_z.m et tige_z.m .

- `tige_z.m`  
  → Dimensionnement minimale d'une tige soumise à une **flexion et flambement**.

- `bielles_x.m`  
  → Dimensionnement minimale des bielles B et C, avec contraintes complexes de **moment, flambement et déformation**.

- `membranes_x.m`  
  → Dimensionnement complet d'une membrane  **avec flexion et torsion combinées**, dans un contexte géométrique imposé (rayon, hauteur de liaison, etc.).

---

## ⚙️ Technologies utilisées

- MATLAB R2022b ou ultérieur
- Scripts autonomes, ne nécessitant pas de dépendances externes
- Visualisation des résultats sous forme de `table`, triées par performance (raideur, contraintes, etc.)

---


## 📈 Résultats attendus

Les scripts génèrent :
- Une **liste triée** de configurations valides respectant toutes les contraintes.
- Une **analyse comparative** par matériau.
- Une **sélection automatique des meilleures solutions** en fonction d’un objectif de raideur.

---

## 🛠️ Utilisation

1. Ouvrir le script correspondant dans MATLAB.
2. Lancer le fichier (modifier les constantes selon les contraintes (géométriques et/ou sortie des autres programmes).
3. Lire les résultats dans la console (`disp`) ou adapter pour exporter en `.csv`.

---

## 📎 Licence

Projet académique – usage libre à des fins pédagogiques ou personnelles.  
Crédits requis en cas de réutilisation :  
> *“Projet EPFL – BA4 Conception de Mécanismes – 2025”*

---

