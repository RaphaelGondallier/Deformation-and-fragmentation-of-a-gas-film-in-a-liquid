// Saves the profile of the fields and phases as gifs and/or png sequences

#include <math.h> 
#include <stdlib.h> 
#include <stdio.h>
#include "axi.h" 
#include "navier-stokes/centered.h" 
#include "two-phase.h" 
#include "tension.h" 
#include "view.h" 

#define Oh 0.1118 
#define L_domain 62.5
#define L0 50

#define max_refine 12
#define min_refine 4
#define t_step 0.01
#define t_end 28
#define alpha 0.045

#define d_refine 10
#define time_factor 100

double t_post_proc = 8.0; 
double target_duration = 7.0;

// Generation or not of a GIF
int do_image = 1;
int do_gif = 1;

// Choose the desired fields
int do_ux = 0;
int do_uy = 0;
int do_p = 0;
int do_omega = 0;
int do_phases = 0;
int do_phi = 1; // Activation du champ de dissipation

// Déclaration globale pour éviter les conflits de libération mémoire (SIGABRT)
scalar omega[];
scalar inverted[];
scalar Phi[]; // Déclaration du nouveau champ

int main()
{
    run();
    return 0;
}

event init(i=0)
{
    char dump_filename[100];
    double t;
    
    p.nodump = false;

    omega.name = "omega";
    inverted.name = "inverted";
    Phi.name = "Phi";

    // Passage aux tableaux de dimension 6
    int active[6] = {do_ux, do_uy, do_p, do_omega, do_phases, do_phi};
    char * field_names[6] = {"u.x", "u.y", "p", "omega", "phases", "Phi"};
    char * dirs[6] = {"ux", "uy", "p", "omega", "phases", "Phi"}; 
    
    // Field scales (limite min de Phi fixée à 0.0 par définition positive)
    // La limite max (ici 50.0) devra potentiellement être ajustée selon vos valeurs extrêmes
    double limit_min[6] = {-1.0, -1.0, -2.0, -10.0, 0.0, -50.0}; // negative value for phi to have white at zero
    double limit_max[6] = { 1.0,  1.0,  2.0,  10.0, 0.0, 50.0};

    // Create the file of each field (and eventually remove what is inside)
    for (int j = 0; j < 6; j++) {
        if (active[j]) {
            char cmd[100];
            sprintf(cmd, "rm -rf %s && mkdir -p %s", dirs[j], dirs[j]);
            if(system(cmd)) {}
        }
    }
    
    // ---------------------------------------------------------
    // Images generation
    // ---------------------------------------------------------
    for(t = 0; t <= t_post_proc ; t += t_step * time_factor)
    { 
        sprintf(dump_filename, "../code/dumpfile_%06.3f", t);
        if (!restore(file = dump_filename)) continue;

        // Set reference pressure with the edge of the domain to remove fluctuations
        double p_ref = 0.;
        double n_ref = 0.;
        foreach(reduction(+:p_ref) reduction(+:n_ref)) {
            if (y > L_domain * 0.9) {
                p_ref += p[];
                n_ref += 1.0;
            }
        }
        
        if (n_ref > 0.) {
            p_ref /= n_ref;
            foreach() {
                p[] -= p_ref;
            }
        }

        if (do_omega) {
            foreach(){
                omega[] = (u.y[1,0]-u.y[-1,0]+u.x[0,-1]-u.x[0,1]) / (2.*Delta + SEPS);
            }
            boundary({omega}); 
        }

        // Dissipation rate calculation phi = D:D
        if (do_phi) {
            foreach() {
                double D_rr = (u.y[0,1] - u.y[0,-1]) / (2.*Delta);
                double D_zz = (u.x[1,0] - u.x[-1,0]) / (2.*Delta);
                double D_rz = 0.5 * ( (u.y[1,0] - u.y[-1,0]) / (2.*Delta) + (u.x[0,1] - u.x[0,-1]) / (2.*Delta) );
                
                // Tr(D) = 0 (incompressibility) implies D_tt = -(D_rr + D_zz)
                double D_tt = -(D_rr + D_zz); // avoid to calculate the 1/r term
                
                Phi[] = sq(D_rr) + sq(D_zz) + sq(D_tt) + 2.*sq(D_rz);
            }
            boundary({Phi});
        }

        if (do_p) {
            stats sp = statsf(p);
            fprintf(stderr, "t = %.2f | Pression min = %g, max = %g\n", t, sp.min, sp.max);
        }

        for (int j = 0; j < 6; j++) {
            if (!active[j]) continue;

            // Symmetry
            if (j == 1) {
                foreach() inverted[] = -u.y[];
                boundary({inverted});
            } else if (j == 3) {
                foreach() inverted[] = -omega[];
                boundary({inverted});
            }

            clear();
            view (tx=-0.5, fov=6.7, height=800, width=2400);
            
            // j=4 correspond au champ discret "phases". 
            // Tous les autres (j != 4) sont des champs continus.
            if (j != 4) {
                squares (field_names[j], min = limit_min[j], max = limit_max[j], map = blue_white_red);
                draw_vof ("f", lc = {0.0, 0.0, 0.0}, lw = 2.0);
                
                mirror(n={0,1,0}){
                    if (j == 1 || j == 3) {
                        squares ("inverted", min = limit_min[j], max = limit_max[j], map = blue_white_red);
                    } else {
                        squares (field_names[j], min = limit_min[j], max = limit_max[j], map = blue_white_red);
                    }
                    draw_vof ("f", lc = {0.0, 0.0, 0.0}, lw = 2.0);
                }
                
                colorbar (map = blue_white_red, min = limit_min[j], max = limit_max[j], size = 30, lscale = 1.5, lw = 5, fsize = 90);
            } 
            // RGB personnalised colors for the phases
            else {
                draw_vof ("f", filled = 1, fc = {0.0, 0.07, 0.38});
                draw_vof ("f", filled = -1, fc = {0.8, 0.8549, 0.941176});
                draw_vof ("f", lc = {0.0, 0.0, 0.0}, lw = 1.0);
                
                mirror(n={0,1,0}){
                    draw_vof ("f", filled = 1, fc = {0.0, 0.07, 0.38});
                    draw_vof ("f", filled = -1, fc = {0.8, 0.8549, 0.941176});
                    draw_vof ("f", lc = {0.0, 0.0, 0.0}, lw = 1.0);
                }
            }

            char file_path[100];
            if (do_image) {
                sprintf(file_path, "%s/%s-%06.3f.png", dirs[j], dirs[j], t);
                save(file = file_path);
            }
        }
    }

    // ---------------------------------------------------------
    // GIF generation if required
    // ---------------------------------------------------------
    if (do_gif && pid() == 0) {
        fprintf(stderr, "\n--- DÉBUT DE L'ASSEMBLAGE DES GIFS ---\n");
        
        int n_frames = 0;
        for(double t_count = 0; t_count <= t_post_proc ; t_count += t_step * time_factor) {
            n_frames++;
        }

        int delay = (int)((target_duration * 100.0) / n_frames);
        if (delay < 1) delay = 1; 

        for (int j = 0; j < 6; j++) {
            if (active[j]) {
                char cmd[500];
                
                sprintf(cmd, "cd %s && convert -density 150 -delay %d -loop 0 %s-*.png %s.gif", 
                        dirs[j], delay, dirs[j], dirs[j]);
              
                if(system(cmd)){}
            }
        }
    }
}
