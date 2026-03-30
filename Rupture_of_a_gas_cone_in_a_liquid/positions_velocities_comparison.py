import numpy as np
import matplotlib.pyplot as plt

# Configuration de la typographie
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "axes.formatter.use_mathtext": True
})

angles = [0, 15, 30, 34, 45, 60]

# Deux figures indépendantes
fig_Z, ax_Z = plt.subplots(figsize=(8, 6))
fig_V, ax_V = plt.subplots(figsize=(8, 6))

for angle in angles:
    filename = f"kinematics_tip{angle}.npz"
    
    try:
        data = np.load(filename)
    except FileNotFoundError:
        print(f"Fichier introuvable : {filename}")
        continue
        
    T = data['T']
    Z = data['Z']
    V = data['V']
    
    # Tracé sur la figure de position
    line_Z, = ax_Z.plot(T, Z, label=r"$\alpha = $" + f"{angle}")
    color = line_Z.get_color()
    
    # Tracé sur la figure de vitesse 
    ax_V.plot(T, V, color=color,label=r"$\alpha = $" + f"{angle}" + " mrad")
    
    # Traitement du point de rupture conditionnel
    if 'rupture_time' in data.files and 'rupture_velo' in data.files:
        t_rup = data['rupture_time']
        v_rup = data['rupture_velo']
        
        idx_rup = np.argmin(np.abs(T - t_rup))
        z_rup = Z[idx_rup]
        
        ax_Z.scatter(t_rup, z_rup, color=color, marker='o', s=30, zorder=3)
        ax_V.scatter(t_rup, v_rup, color=color, marker='o', s=30, zorder=3)

ax_Z.set_xlabel(r"t*")
ax_Z.set_ylabel(r"z*")
#ax_Z.grid(True, linestyle='--', alpha=0.5)
ax_Z.legend(loc='best')
ax_Z.set_xlim(0.0, 12.0) 

ax_V.set_xlabel(r"t*")
ax_V.set_ylabel(r"U*")
#ax_V.grid(True, linestyle='--', alpha=0.5)
ax_V.legend(loc='best')
ax_V.set_xlim(0.0, 12.0)  

fig_Z.tight_layout()
fig_V.tight_layout()

#fig_Z.savefig("positions.pdf", format='pdf', bbox_inches='tight')
#fig_V.savefig("velocities.pdf", format='pdf', bbox_inches='tight')

plt.show()