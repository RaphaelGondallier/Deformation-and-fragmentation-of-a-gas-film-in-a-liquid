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
R = L0/10
t_step = 0.01 
intervalle_sauvegarde = t_step
ti=0
tf=2
n= round((tf-ti)/intervalle_sauvegarde) # nb de pas de temps enregistres avec outpufacets
max_refine = 12

Z_list=[]
r_relat_error=[]
nb_compt_sup_max_iter=0
# =============================================================================
# Post-processing
# =============================================================================
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
#   Récupère le point de coord axiale max pour vérifier qu'il correspond au bout
# =============================================================================
    list_tip = data.copy()
    compt = 0
    max_iter = 100

    # Marges d'exploration initiales relatives (indépendantes de la rétraction)
    marge_z = 0.05 * L0  # Profondeur de recherche juste derrière le bout
    tol_r = 0.5 * R      # On commence large en rayon

    # On boucle tant que le point le plus en avant (max Z) 
    # n'est pas le point le plus proche de l'axe (min R)
    while (np.argmax(list_tip[:, 0]) != np.argmin(np.abs(list_tip[:, 1]))) and (compt < max_iter):
        compt += 1
        if compt==max_iter : nb_compt_sup_max_iter+=1
        # on trouve le Z maximum ACTUEL 
        z_max_actuel = np.max(list_tip[:, 0])
        
        # Mise à jour de la tolérance Z par rapport à ce nouveau bout
        tol_z = z_max_actuel - marge_z
        
        # Élimine la rupture au loin (ne garde que les points proches du bout)
        list_tip_new = data[data[:, 0] > tol_z]
        
        # Élimine les points de plus grand z que le bout, en cas de renfoncement du bout, en resserrant le rayon
        list_tip_new = list_tip_new[np.abs(list_tip_new[:, 1]) < tol_r]
        
        # si la tolérance devient trop fine et vide la liste, on s'arrête
        if len(list_tip_new) == 0:
            break
            
        list_tip = list_tip_new
        
        # Resserrement de la boîte pour forcer la convergence
        tol_r *= 0.95  # On referme doucement le cylindre autour de l'axe

    # Extraction finale 
    id_tip = np.argmax(list_tip[:, 0])
    z_tip = list_tip[id_tip, 0]
    Z_list.append(z_tip)
    r_tip = list_tip[id_tip, 1] 
    r_relat_error.append(abs(r_tip-0)/R)
    
# =============================================================================
#   Affichages interface
# =============================================================================
    fig, ax = plt.subplots()
    
    lc = LineCollection(segments, color='black', linewidth=0.2)
    ax.add_collection(lc)
    ax.plot(z_tip, r_tip, marker='o', color='red', markersize=2)
    
    ax.set_xlabel("z")
    ax.set_ylabel("r")

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, L_dom)
    ax.set_ylim(-R, R)
    plt.title(f"t={t}")
    plt.savefig(f"Figures/figure_t_{t:06.3f}.svg", bbox_inches="tight")
        
    plt.show()
    plt.close(fig)

plt.plot(np.linspace(0, (n-1)*intervalle_sauvegarde, n), r_relat_error)
plt.xlabel('t*')
plt.ylabel('Relative error')
print(f'{100*nb_compt_sup_max_iter/n}% des déterminations de la position du bout n ont pas convergé')
plt.savefig("Figures/Relative_error_on_tip_position.pdf")
plt.close()
# =============================================================================
# Vitesse de contraction
# =============================================================================
V=[0]
T=np.linspace(0, (n-1)*intervalle_sauvegarde, n)

Z = np.array(Z_list)
T = np.array(T)

# 2. Calcul de la vitesse (Schéma à 3 points d'ordre 2)
V = -np.gradient(Z, T)

plt.figure(figsize=(12, 5))

# Graphique de la Position
plt.subplot(1, 2, 1)
plt.plot(T, Z, marker='.', linestyle='-')
plt.xlabel("Time")
plt.ylabel("Tip axial coordinate")

# Graphique de la Vitesse
plt.subplot(1, 2, 2)
plt.plot(T, V, marker='.', linestyle='-', color='darkred')
plt.xlabel("Time")
plt.ylabel("Contraction velocity") 

plt.tight_layout()
plt.savefig("Figures/Tip_position_and_velocity.pdf")
plt.close()

plt.loglog(T, V, marker='.', linestyle='-', color='darkred')
plt.xlabel("Time")
plt.ylabel("Contraction velocity") 
plt.savefig("Figures/loglog_tip_velocity.pdf")
plt.close()

#Vérification du déplacement minimal correctement capturé par le maillage
Displ_list = V*T
dmin = L_dom/2**(max_refine)
t_min = np.argmax(Displ_list >= dmin)
print('Les résultats ne sont valables qu à partir de t = ', T[t_min])


plt.loglog(T, V*T, label='Displacement of the tip', marker='.', linestyle='-', color='darkred')
plt.loglog(T, dmin*np.ones(len(T)), linestyle='--', label='Minimal mesh cell size')
plt.axvline(x=T[t_min], color='orange', linestyle='--', zorder=0, label='Minimal valid time')
ax = plt.gca()
ax.text(T[t_min], -0.03, f"{T[t_min]:.4f}", 
        color='orange', 
        transform=ax.get_xaxis_transform(),
        ha='center', va='top',            
        fontsize=9)

plt.xlabel("Time")
plt.ylabel("Displacement") 
plt.legend()
plt.savefig("Figures/Minimal_valid_time.pdf")



# plt.loglog(T,V, linestyle='-', color='orange', label='Iterative method')
# plt.loglog(T,V1, linestyle='-', color='black', label='Max z method')
# plt.legend()
# plt.xlabel("Time")
# plt.ylabel("Contraction velocity") 
# plt.savefig("Figures/Compa_2methods_tip_tracking.pdf")
# plt.show()

# plt.loglog(T, abs(V1-V)/V, marker='.', linestyle='-', color='blue')
# plt.xlabel("Time")
# plt.ylabel("Relative error between the 2 methods") 
# plt.savefig("Figures/Relative_error_2methods_tip_tracking.pdf")
