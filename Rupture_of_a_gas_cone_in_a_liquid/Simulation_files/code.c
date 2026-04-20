/**
The two-phase Navier-Stokes solver is included.
*/
#include <math.h> 
#include "axi.h" //Axisymmetry around x axis
#include "navier-stokes/centered.h" 
#include "two-phase.h" //Two phased fluid with interface
#include "tension.h" //Surface tension
#include "view.h" //For exporting images

//Constants (unit distance D0 = 2e, the cone tip diameter)

#define Oh 0.1118 // /2 de manière a avoir le meme Oh que Prosper dans ses simu (pour voir si la diff d'angle crit vient de là)
#define L_domain 62.5
#define L0 50

#define max_refine 12
#define min_refine 4
#define t_step 0.01
#define t_end 28
#define alpha 0.0

#define d_refine 10
/** The geometrical shape of the filament is defined. This function is negative outside the filament and positive inside the filament.*/
double geometry(double x, double y)
{
    double Line_up = -y + L0*tan(alpha) - x*tan(alpha) + 0.5; // longueurs adimensionnées avec D
    double Line_right = -x + L0;
    double forme = min(Line_right,Line_up);
    double circle = -sq(y) - sq(L0 - 0.5*tan(alpha)- x) + sq(0.5/cos(alpha));
    double shape_final = max(forme,circle);
    return -shape_final;

}

int main()
{
    size(L_domain);
    origin(0.,0.);
    init_grid(1<<min_refine);

/** density ratio and viscosity ratio is fixed.

Basilisk solves Navier-Stokes equations with surface tension terms which read:

$$\rho(\partial_{t}\bm{u}+\bm{u}\cdot\nabla\bm{u})=-\nabla p+\nabla\cdot(2\mu\bm{D})+\sigma\kappa\delta_{s}\bm{n}$$

Diemnsionless equation is as below (details of non-dimensionalization can be found in the published paper):

$$\rho^{*}\frac{D\bm{u}}{Dt}= -\nabla p + Oh\mu^{*} \Delta \bm{u} + \delta^{*}_{s} \kappa^{*}\bm{n}$$
*/
    rho1 = 1.; //By definition rho1=1 since rho^* = rho / rho_w  
    //By definition rho2 is the ratio between the density of the gas and the density of the liquid
    //Both values for density are taken at 20°C and 1 atm from G. K. Batchelor - An introduction to fluid dynamics
    //Density of water is 0.9982 g/cm^3 while density of air is 1.205*10^-3 g/cm^3
    rho2 = 1./828.4;  

    /** By varying the Ohnesorge number $\mu_{1}$, the cases in the paper can be reproduced. */

    mu1 = Oh; // a retiré le D0 ici. mu1 is variable since mu^* = mu / mu_w, can be changed to change Ohnesorge number  
    //By definition mu2 is the ratio between the viscosity of the gas and the viscosity of the liquid times mu1
    //Both values for viscosity are taken at 20°C and 1 atm from G. K. Batchelor - An introduction to fluid dynamics
    //Viscosity of water is 1.002*10^-2 g/(cm s) while density of air is 1.81*10^-4 g/(cm s)
    mu2 = mu1/55.4;

    // surface tension is equal to 1 in dimensionless units
    f.sigma = 1.;

    //Boundary conditions for top, right, and left are set as symmetric (same as default)
    //Axisymmetric boundary conditions for the bottom boundary (set by axi.h) 
    u.n[top] = dirichlet(0.);
    u.t[top] = neumann(0.);
    u.n[right] = dirichlet(0.);
    u.t[right] = neumann(0.);
    u.n[left] = dirichlet(0.);
    u.t[left] = neumann(0.);

    run();
}

event init(t=0)
{
  /** The area around the interface is refined at the initial step.*/
    refine(abs(y-L0*tan(alpha) - 0.5 + x*tan(alpha)) < d_refine && level < max_refine); 
    fraction(f,geometry(x,y));
}

event adapt(i++)
{
   adapt_wavelet({f,u},(double []){1e-8,1e-3,1e-3},max_refine,min_refine); // a retiré sqrt(D0) ici. à chaque pas de temps, affine le maillage, en fonction de f et de u, le maillage est plus ou moins affiné s'il y a des variations importantes de u et de f ou non. Ceci dépend des trois nombres après qui s'appliquent respectivement à f, ux, uy
}

char file_name[100];

/** Simulation output and we do post-processing after the simulation*/
event outputfile(t+=t_step)
{
    scalar omega[];
    p.nodump = false;
    foreach(){
       omega[] = (u.y[1,0]-u.y[-1,0]+u.x[0,-1]-u.x[0,1]) / (2.*Delta + SEPS); //Vorticity definition for axisymmetry
    }
    sprintf(file_name,"dumpfile_%06.3f",t);
    dump(file=file_name);
}

//event image(t+=100*t_step)
//{
//    char png_filename[100];
//    clear();
//    view (tx=-0.5, fov=7, height=800, width=2400);
//    squares (color = "f", min = -1.5, max = 2.);
//    draw_vof (c = "f");   
//    mirror(n={0,1,0}){
//       squares (color = "f", min = -1.5, max = 2.);
//       draw_vof (c = "f");
//    }
//    sprintf(png_filename, "image-%.2f.png", t);
//    save(png_filename);
//}

event out_facets (t = 0; t +=t_step; t <= t_end) {
        char name[100];
        sprintf(name, "pool_t_%06.3f.dat", t);
        FILE * fp1 = fopen (name, "w");
        output_facets (f, fp1);

        fclose(fp1);
}

event end(t=t_end)
{
    
}
	
