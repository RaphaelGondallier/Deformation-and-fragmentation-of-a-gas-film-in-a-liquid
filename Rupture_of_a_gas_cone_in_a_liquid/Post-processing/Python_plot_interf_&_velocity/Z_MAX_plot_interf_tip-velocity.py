import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
os.makedirs("Figures", exist_ok=True) #crée un dossier pour les images
os.makedirs("Archives", exist_ok=True) #crée un dossier pour les données de vitesse obtenues
# =============================================================================
# Simulation file parameters
# Carefully check these parameters before running!
# =============================================================================
Until_rupture = False # whether you want to stop at the rupture or not
L_dom =125
L0 = 100
R = 0.5 # = half of unit diameter
t_step = 0.01 
intervalle_sauvegarde = t_step
ti=0
tf=5
n = round((tf-ti)/intervalle_sauvegarde)+1 # nb de pas de temps enregistres avec outpufacets
max_refine = 12

parallel = False # whether the output_facets was run in parallel
nproc = 6 # if parallel how many cpus
# =============================================================================
# Post-processing
# =============================================================================
Z_list=[]
nb_compt_sup_max_iter=0
dx_min = L_dom/2**(max_refine)
iter_rupt = 0
rupture_time = 0

for i in range(n):
    t = i * intervalle_sauvegarde
    
    if not parallel:
        filename = f"pool_t_{t:06.3f}.dat"
        data = np.loadtxt(filename) #shape 2N, 2 [ [x1, y1]
                                    #              [x1', y1']
                                    #              [x2, y2]
                                    #              [x2', y2'] ]
    if parallel :
        valid_data = []
        for i in range(nproc):
            filename = f"pool_t_{t:06.3f}_pid{i}.dat"
        
        # Validation de l'existence de données (taille du fichier > 0 octet)
            if os.path.getsize(filename) > 0:
                d = np.loadtxt(filename, ndmin=2)
            # Sécurité supplémentaire au cas où le fichier ne contiendrait que des espaces
                if d.size > 0: 
                    valid_data.append(d)
    # Concaténation exclusive des tableaux de dimension (N, 2)
        data = np.vstack(valid_data)
        
    segments_haut = data.reshape(-1, 2, 2)  # assemble les segments
                                            #shape N, 2, 2 [ [[x1, y1]
                                             #               [x1', y1']]
                                             #              [[x2, y2]
                                             #               [x2', y2']] ]
    
    segments_bas = np.copy(segments_haut)
    segments_bas[:, :, 1] = -segments_bas[:, :, 1] # to have the bottom part y_bottom = -y_top
    
    segments = np.concatenate((segments_haut, segments_bas), axis=0)

# =============================================================================
#   Récupère le point de coord radiale nulle
# =============================================================================
    # Identification des index où l'ordonnée (rayon) est nulle avec une tolérance (dx_min/10 ici, inférieur à la limite de maille)
    # isclose() crée une liste de booléens indiquant si les rayons sont assez proches de 0
    # where() trouve l'index (ou les index si rupture) correspondant 
    liste_index_zero = np.where(np.isclose(data[:, 1], 0.0, atol=dx_min*10e-5))[0]
    # Évaluation du critère d'arrêt topologique
    if len(liste_index_zero) != 1:
        print(f"Arrêt du post-traitement à t={t:06.3f} : la liste contient {len(liste_index_zero)} index.")
        print(liste_index_zero + liste_index_zero//2 +1)
        if Until_rupture == True : # stop the simu iif this is true
            break
        if rupture_time == 0: # pour ne mettre à jour ce temps qu'une fois
            rupture_time = t
    elif len(liste_index_zero) == 1 and rupture_time == 0: 
        iter_rupt += 1 # Pour arrêter le calcul lors de la rupture
    
# =============================================================================
#   Récupère le point de coord axiale max 
# =============================================================================
    z_tip = np.max(data[:,0])
    id_tip = np.argmax(data[:,0])
    z_tip = data[id_tip, 0]
    r_tip = data[id_tip, 1]
    
    Z_list.append(z_tip)
# =============================================================================
#   Affichages interface
# =============================================================================
    fig, ax = plt.subplots()
    
    lc = LineCollection(segments, color='black', linewidth=0.2)
    ax.add_collection(lc)
    #ax.plot(z_tip, r_tip, marker='o', color='red', markersize=0.2, markeredgewidth=0.0)
    ax.scatter(z_tip, r_tip, color='red', s=0.2, edgecolors='none')
    
    ax.set_xlabel("z")
    ax.set_ylabel("r")

    ax.set_aspect('equal', adjustable='box')
    if Until_rupture == True :
        ax.set_xlim(L0*0.8, L0*1.1)
        ax.set_ylim(-3*R, 3*R)
    plt.title(f"t={t:06.3f}")
    plt.savefig(f"Figures/figure_t_{t:06.3f}.svg", bbox_inches="tight")  
    #plt.show()
    plt.close(fig)

# =============================================================================
# Vitesse de contraction
# =============================================================================
V=[0]
if Until_rupture == True : # ne plot qu'avant la rupture
    T=np.linspace(0, (iter_rupt-1)*intervalle_sauvegarde, iter_rupt)
if Until_rupture == False : # plot tout
    T=np.linspace(0, (n-1)*intervalle_sauvegarde, n)
Z = np.array(Z_list)
T = np.array(T)

# 2. Calcul de la vitesse (Schéma à 3 points d'ordre 2)
V = -np.gradient(Z, T)
#V_test = [-(Z[i-2]-Z[i+2]+8*Z[i+1]-8*Z[i-1])/(12*t_step) for i in range(2,n-1)]
# vitesse au moment de la rupture
if Until_rupture == False and rupture_time != 0: # la 2e condition fait qu'aucun point n'est tracé s'il n'y a pas de rupture
    rupture_velo = V[iter_rupt]
# Sauvegarde des variables cinématiques si a besoin de modifier les plots
    np.savez_compressed("Archives/kinematics_tip.npz", T=T, Z=Z, V=V, rupture_time=rupture_time, rupture_velo=rupture_velo)

np.savez_compressed("Archives/kinematics_tip.npz", T=T, Z=Z, V=V)

# Plots

# Graphique de la Position
plt.plot(T, Z, linestyle='-')
plt.xlabel("Time")
plt.ylabel("Tip axial coordinate")
if Until_rupture :
    plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)
plt.savefig("Figures/Tip_position.pdf")
plt.show()
plt.close()

# Graphique de la Vitesse
plt.plot(T, V, linestyle='-', color='darkred')
if Until_rupture == False and rupture_time != 0: # la 2e condition fait qu'aucun point n'est tracé s'il n'y a pas de rupture
    plt.plot(rupture_time, rupture_velo, marker='o', color='navy', markersize=4)
plt.xlabel("t*")
plt.ylabel("U*") 
if Until_rupture and rupture_time != 0:
    plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)
if Until_rupture and rupture_time != 0:
    print("No rupture! Increase t_f?")
plt.tight_layout()
plt.savefig("Figures/Tip_velocity.pdf")
plt.show()
plt.close()

plt.loglog(T, V, linestyle='-', color='darkred')
if Until_rupture == False and rupture_time != 0: # la 2e condition fait qu'aucun point n'est tracé s'il n'y a pas de rupture
    plt.plot(rupture_time, rupture_velo, marker='o', color='navy', markersize=4)
plt.xlabel("t*")
plt.ylabel("U*") 
if Until_rupture :
    plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde) # +10 pour voir la rupture sur le graphe
plt.savefig("Figures/loglog_tip_velocity.pdf")
plt.show()
plt.close()

#Vérification du déplacement minimal correctement capturé par le maillage
plt.plot(T, V*t_step/dx_min, label=r'$\xi(r=0) / dx_{min} $ (max refinement level = 12)', linestyle='-', color='darkred')
plt.plot(T, np.ones(len(T)), linestyle='--', label='Minimal mesh cell size limit')
plt.xlabel("Time")
plt.ylabel(r"Displacement of the tip normalized by" + "\n" + r"the minimal mesh size $\xi(r=0) / dx_{min}$")
if Until_rupture :
    plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)
plt.legend()
plt.tight_layout()
plt.savefig("Figures/Minimal_mesh_cell_size_limit.pdf")
plt.show()
plt.close()

plt.loglog(T, V*t_step/dx_min, label=r'$\xi(r=0) / dx_{min} $ (max refinement level = 12)', linestyle='-', color='darkred')
plt.loglog(T, np.ones(len(T)), linestyle='--', label='Minimal mesh cell size limit')
plt.xlabel("Time")
plt.ylabel(r"Displacement of the tip normalized by" + "\n" + r"the minimal mesh size $\xi(r=0) / dx_{min}$")
if Until_rupture :
    plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)
plt.legend()
plt.tight_layout()
plt.savefig("Figures/Minimal_mesh_cell_size_limit_loglog.pdf")
plt.show()
plt.close()
