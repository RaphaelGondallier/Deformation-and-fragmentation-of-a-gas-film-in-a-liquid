// Saves the profile of the interface as a DAT file 

#include <math.h> 
#include "axi.h" //Axisymmetry around x axis
#include "navier-stokes/centered.h" 
#include "two-phase.h" //Two phased fluid with interface
#include "tension.h" //Surface tension
#include "view.h" //For exporting images

//Constants	
//#define e 4.0 //épaisseur au bord (2/5 *D0)
//#define D0 2.5 // grandeurs deja adimensionnées "D0*" ici, de telle sorte à ce que e n'apparaisse pas dans les equations (là où avait D0/e, aura juste D0)
// Utilise desormais D0 = 2e pour adimensionner, donc pour le seuil de raffinement initial utilisera 1
#define Oh 0.001 // /2 de manière a avoir le meme Oh que Prosper dans ses simu (pour voir si la diff d'angle crit vient de là)
#define L_domain 50
#define L0 40
#define d_refine 2

//#define time_factor sqrt(pow(D0,3))
#define max_refine 10
#define min_refine 4
#define t_step 0.01
#define t_end 4.7
#define alpha 0.0

int main()
{
	run();
	return 1;
}

event init(i = 0)
{
	char filename[100];
	//double t_step 0.01;
	//double t_end 20.00;
	double t;
	char figurename[100];
	int i = 0;
	for(t = 0; t <= t_end; t += t_step)
	{ 
		sprintf(filename,"../code/dumpfile_%.2f", t);
		restore(file = filename);
		sprintf(figurename,"pool_t_%.2f.dat", t);
		FILE* fp = fopen (figurename, "w");
  		output_facets (f, fp); //f1 for water-oil; f2 for oil-air interface
		i++;
	}
}

event end(i=0)
{
}


