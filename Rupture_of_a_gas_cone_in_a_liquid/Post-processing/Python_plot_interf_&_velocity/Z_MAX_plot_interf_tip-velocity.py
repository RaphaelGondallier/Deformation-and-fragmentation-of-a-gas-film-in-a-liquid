import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
os.makedirs("Figures", exist_ok=True) #crée un dossier pour les images

# =============================================================================
# Simulation file parameters
# =============================================================================
L_dom = 50
L0 = 40
R = 0.5 # = half of unit diameter
t_step = 0.01 
intervalle_sauvegarde = t_step
ti=0
tf=2
n= round((tf-ti)/intervalle_sauvegarde) # nb de pas de temps enregistres avec outpufacets
max_refine = 12

Z_list=[]
nb_compt_sup_max_iter=0

# =============================================================================
# Post-processing
# =============================================================================
dx_min = L_dom/2**(max_refine)
iter_rupt = 0

for i in range(n):
    t = i * intervalle_sauvegarde
    filename = f"pool_t_{t:06.3f}.dat"
    data = np.loadtxt(filename) #shape 2N, 2 [ [x1, y1]
                                #              [x1', y1']
                                #              [x2, y2]
                                #              [x2', y2'] ]

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
    liste_index_zero = np.where(np.isclose(data[:, 1], 0.0, atol=dx_min/10))[0].tolist()
    
    # Évaluation du critère d'arrêt topologique
    if len(liste_index_zero) != 1:
        print(f"Arrêt du post-traitement à t={t:06.3f} : la liste contient {len(liste_index_zero)} index.")
        break
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
    ax.plot(z_tip, r_tip, marker='o', color='red', markersize=0.5)
    
    ax.set_xlabel("z")
    ax.set_ylabel("r")

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(L0*0.8, L0*1.1)
    ax.set_ylim(-3*R, 3*R)
    plt.title(f"t={t:06.3f}")
    plt.savefig(f"Figures/figure_t_{t:06.3f}.svg", bbox_inches="tight")  
    plt.show()
    plt.close(fig)

# =============================================================================
# Vitesse de contraction
# =============================================================================
V=[0]
T=np.linspace(0, (iter_rupt-1)*intervalle_sauvegarde, iter_rupt)

Z = np.array(Z_list)
T = np.array(T)

# 2. Calcul de la vitesse (Schéma à 3 points d'ordre 2)
V = -np.gradient(Z, T)

plt.figure(figsize=(12, 5))

# Graphique de la Position
plt.subplot(1, 2, 1)
plt.plot(T, Z, linestyle='-')
plt.xlabel("Time")
plt.ylabel("Tip axial coordinate")
plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)

# Graphique de la Vitesse
plt.subplot(1, 2, 2)
plt.plot(T, V, linestyle='-', color='darkred')
plt.xlabel("Time")
plt.ylabel("Contraction velocity") 
plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)


plt.tight_layout()
plt.savefig("Figures/Tip_position_and_velocity.pdf")
plt.show()
plt.close()

plt.loglog(T, V, linestyle='-', color='darkred')
plt.xlabel("Time")
plt.ylabel("Contraction velocity") 
plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde) # +10 pour voir la rupture sur le graphe
plt.savefig("Figures/loglog_tip_velocity.pdf")
plt.show()
plt.close()

#Vérification du déplacement minimal correctement capturé par le maillage
plt.plot(T, V*t_step/dx_min, label=r'$\xi(r=0) / dx_{min} $ (max refinement level = 12)', linestyle='-', color='darkred')
plt.plot(T, np.ones(len(T)), linestyle='--', label='Minimal mesh cell size limit')
plt.xlabel("Time")
plt.ylabel(r"Displacement of the tip normalized by" + "\n" + r"the minimal mesh size $\xi(r=0) / dx_{min}$")
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
plt.xlim(t_step, (iter_rupt+10)*intervalle_sauvegarde)
plt.legend()
plt.tight_layout()
plt.savefig("Figures/Minimal_mesh_cell_size_limit_loglog.pdf")
plt.show()
plt.close()
