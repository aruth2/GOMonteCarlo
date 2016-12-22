#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gtk/gtk.h>
#include <pthread.h>

/*****************************************************************/
/* Kinetic Monte Carlo Simulation of Reduction in Graphene Oxide
 * Created By Anthony Ruth on July 16th, 2015
 * 
 * This program simulates the reduction of graphene oxide by the removal
 * of functional groups and sublimation of the carbon lattice. During the
 * reduction process the formation of individual graphenic clusters is monitored
 * to construct an emission and absorption profile.
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * */

//These 3 variable are used in subfile and need to be defined before the include statements
int debug = 0;
int verbosedebug = 0;  
struct sim *UIsim; //The main simulation object. During a nonensemble simulation this sim object stores the information associated with the simulation displayed in the UI. 
//During and ensemble simulation, this sim object holds the information associated with the ensembly averaged simulation.


//This file contains the Graphical User Interface
#include "aruthmath.h" //Additional math library
#include "aruthio.h" //Additional IO library

#include "GOSim.h" //The simulation library. Creates a simulation object and contains the operations for the simulation to proceed.
#include "GOplot.h"
#include "GODebug.h" //A debugging library. Prints debugging information to stdout. Used primarily in development.

GtkWidget *DirEntry; //The directory to store output files
GtkWidget *drawingarea; //Visualization of the graphene
GtkWidget *hydroxideremoval; //OH rate constant
GtkWidget *COremoval; //epoxide rate constant
GtkWidget *C3removal; //C3 rate constant
GtkWidget *Czigremoval; //Czig rate constant
GtkWidget *Carmremoval; //Carm rate constant
GtkWidget *C1removal; //C1 rate constant
GtkWidget *Ohop; //Oxygen Hopping rate constant
GtkWidget *Osource; //Oxygen absorption rate constant

GtkWidget *Ensemblespinbutton; //Number of simulations to run to create ensemble
GtkWidget *numbertoremove; //Number of Kinetic steps M, between checks of optical properties
GtkWidget *dimensionspin; //Height and Width of cell in units of graphene lattice constant
GtkWidget *Ospin; //Epoxide Coverage
GtkWidget *OHspin; //Hydroxide Coverage
GtkWidget *vacancyspin; //Vacancy Coverage
GtkWidget *TopLevel; //Toplevel Window Holds UI
GtkWidget *Periodiccheckbutton; //Debugging, checks if periodic boundary conditions are correct
GtkWidget *Covercheckbutton; //Debugging, Checks coverage of carbon atoms.
GtkWidget *savecheckbutton; //Debugging, 
GtkWidget *AtomNumberscheckbutton; //Debugging, prints adjacency list of carbon atoms, and functional group connections
GtkWidget *Debugcheckbutton; //Enables debugging
GtkWidget *VerboseDebugcheckbutton; //Enables verbose debugging


int continuous; //A single step or reduction till end
int maxthreads = 3; //Maximum number of threads to create at one time during ensemble calculation. Should be <= number of cores

char *workingdirectory = "working";//Directory which stores picture files during video creation


pthread_mutex_t mutex; //A mutex to lock the UIsim when an simulation in the ensemble finishes and needs to copy data over.
double *emissionvdegree;
double *absorptionvdegree;
int numbins = 1000; //The number of bins to use in distribution calculation
int maxpointsperbin = 100000;
int *numuseddegree;

void readbuttons(struct sim *s)
{
	//Loads the settings from the UI
	s->dimensions=50; 
	s->periodic=0;
	s->Oprobability = 0.3;

	if(numbertoremove != NULL)
	numremoves = gtk_spin_button_get_value(GTK_SPIN_BUTTON(numbertoremove));
	if(hydroxideremoval != NULL)
	hydroxideprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(hydroxideremoval));
	if(COremoval != NULL)
	COprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(COremoval));
	if(C1removal != NULL)
	C1probability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(C1removal));
	if(C3removal != NULL)
	C3probability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(C3removal));
	if(Czigremoval != NULL)
	Czigprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Czigremoval));
	if(Carmremoval != NULL)
	Carmprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Carmremoval));
	if(Ohop != NULL)
	Ohopprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Ohop));
	if(Osource != NULL)
	Osourceprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Osource));
	if(Debugcheckbutton != NULL)
	debug = gtk_toggle_button_get_active(Debugcheckbutton);
	if(VerboseDebugcheckbutton != NULL)
	verbosedebug = gtk_toggle_button_get_active(VerboseDebugcheckbutton);
	if(Periodiccheckbutton != NULL)
	s->periodic = gtk_toggle_button_get_active(Periodiccheckbutton);
	if(Covercheckbutton != NULL)
	cover = gtk_toggle_button_get_active(Covercheckbutton);
	
	if(dimensionspin != NULL)
	s->dimensions = gtk_spin_button_get_value(GTK_SPIN_BUTTON(dimensionspin));	
	if(Ospin != NULL)
	s->Oprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Ospin));
		
	if(OHspin != NULL)
	s->OHprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(OHspin));
	if(vacancyspin != NULL)
	s->vacancycoverage = gtk_spin_button_get_value(GTK_SPIN_BUTTON(vacancyspin));
	//if(reversibleratespin != NULL)
	//s->reversiblerate = gtk_spin_button_get_value(GTK_SPIN_BUTTON(reversibleratespin));
}

void savesettings()
{
	//Save the settings from the UI
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/settings",UIsim->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g\n",log(*(s->times+index)),*(s->QY+index));
	fclose(outfile);
	
	if(numbertoremove != NULL)
	numremoves = gtk_spin_button_get_value(GTK_SPIN_BUTTON(numbertoremove));
	if(hydroxideremoval != NULL)
	hydroxideprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(hydroxideremoval));
	if(COremoval != NULL)
	COprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(COremoval));
	if(C1removal != NULL)
	C1probability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(C1removal));
	if(C3removal != NULL)
	C3probability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(C3removal));
	if(Czigremoval != NULL)
	Czigprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Czigremoval));
	if(Carmremoval != NULL)
	Carmprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Carmremoval));
	if(Ohop != NULL)
	Ohopprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Ohop));
	if(Osource != NULL)
	Osourceprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Osource));
	if(Debugcheckbutton != NULL)
	debug = gtk_toggle_button_get_active(Debugcheckbutton);
	if(VerboseDebugcheckbutton != NULL)
	verbosedebug = gtk_toggle_button_get_active(VerboseDebugcheckbutton);
	if(Periodiccheckbutton != NULL)
	s->periodic = gtk_toggle_button_get_active(Periodiccheckbutton);
	if(Covercheckbutton != NULL)
	cover = gtk_toggle_button_get_active(Covercheckbutton);
	
	if(dimensionspin != NULL)
	s->dimensions = gtk_spin_button_get_value(GTK_SPIN_BUTTON(dimensionspin));	
	if(Ospin != NULL)
	s->Oprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(Ospin));
		
	if(OHspin != NULL)
	s->OHprobability = gtk_spin_button_get_value(GTK_SPIN_BUTTON(OHspin));
	if(vacancyspin != NULL)
	s->vacancycoverage = gtk_spin_button_get_value(GTK_SPIN_BUTTON(vacancyspin));
	//if(reversibleratespin != NULL)
	//s->reversiblerate = gtk_spin_button_get_value(GTK_SPIN_BUTTON(reversibleratespin));
}

void changedir()
{
	//Initializes the directory for output files
	if(DirEntry != NULL)
	UIsim->dir = gtk_entry_get_text(DirEntry);
	else
	UIsim->dir="default";
	char buffer[1000];
	sprintf(buffer,"mkdir %s\n",UIsim->dir);
	system(buffer);
}

void initializecallback(struct sim *s)
{
	//Reads the settings and initialized the GO sheet
	int olddimension = s->dimensions;
	readbuttons(UIsim);
	initialize(s);
}

gboolean continuousremoves(void *garbage)
{
	// This function is designed to be run by the g_timeout_add call inside the start cycle function. 
	// This function performs the simulation updating the depiction of the GO sheet along the way
	// This function also handles the creation of the video if that box is checked. 
	
	
	if((double)UIsim->numcarbonsold/UIsim->numcarbons > 1.02 && UIsim->numcarbonsold-UIsim->numcarbons > 100)//Provides periodic updates on the state of the system
	{
			system("date");
			printf("Carbons: %d, Oxygens: %d, Hydroxides: %d\n",UIsim->numcarbons,UIsim->numoxygens,UIsim->numhydroxides);
			UIsim->numcarbonsold = UIsim->numcarbons;
	}

	int save = gtk_toggle_button_get_active(savecheckbutton); //Whether a video should be created
	char buffer[1000];	
	removes(UIsim);
	
	clusterdetection(UIsim);
	gdk_threads_enter();
	draw();
	gdk_threads_leave();
	
	    double degree = *(UIsim->p1+UIsim->stepnumber) + 2 * *(UIsim->p2+UIsim->stepnumber) + 3 * *(UIsim->p3+UIsim->stepnumber);
		int degreeindex = degree/3 * numbins;
		if(numuseddegree[degreeindex] != maxpointsperbin-1)
		{
			emissionvdegree[maxpointsperbin*degreeindex+numuseddegree[degreeindex]] = *(UIsim->emission+UIsim->stepnumber);
			absorptionvdegree[maxpointsperbin*degreeindex+numuseddegree[degreeindex]] = *(UIsim->absorption+UIsim->stepnumber);
			numuseddegree[degreeindex]++;
		}
	
	if(save)//Create a single frame of the video
	{
		outputabsorption(UIsim);
		outputemission(UIsim);
		sprintf(buffer,"%s/absems%d.png",workingdirectory,UIsim->stepnumber-1);
		plotabsemsandcombine(buffer);
		sprintf(buffer,"montage -mode concatenate -tile x1 %s/image%d.png %s/absems%d.png %s/comboimage%d.png",workingdirectory,UIsim->stepnumber-1,workingdirectory,UIsim->stepnumber-1,workingdirectory,UIsim->stepnumber-1);
		system(buffer);
	}
	
	//Detects if the simulation is over
	if(UIsim->numcarbons >0 && UIsim->numsourcesites > 0 || UIsim->numhydroxides > 0 || UIsim->numhops > 0 && continuous)
	return TRUE;
	else
	{
	continuous = 0;
	if(save)
	{
		printf("Making video\n");
		char buffer[1000];
		system("rm test.mp4");
		sprintf(buffer,"avconv -y -framerate 30 -i %s/comboimage%%d.png -b:v 1000k %s/video.mp4\n",workingdirectory,UIsim->dir);
		printf(buffer);
		system(buffer);
		
		//sprintf(buffer,"rm %s/*\n",workingdirectory); //Removes the individual frames after the video is created
		//printf(buffer);
		//system(buffer);
	}
	outputemission(UIsim);
	outputabsorption(UIsim);
	outputQY(UIsim);
	outputprobabilities(UIsim);
	outputclusters(UIsim);
	outputrates(UIsim);
	//outputdistribution(UIsim); //Because the distribution data is enormous this function and the plot distribution function were disabled
	plotQY();
	plotabsorption();
	plotprobabilities();
	plotabsems();
	processdegreebins(1);
	plotemvsabs();
	plotclusters();
	plotrates();
	//plotdistribution();
	return FALSE;
	}
}

gboolean startcycle()
{
	//Either performs a single step of the simulation or starts running the simulation until the GO is completely gone
	continuous = !continuous;
	readbuttons(UIsim);
	if(emissionvdegree != NULL)
	free(emissionvdegree);
	if(absorptionvdegree != NULL)
	free(absorptionvdegree);
	if(numuseddegree != NULL)
	free(numuseddegree);
	emissionvdegree = malloc(numbins*maxpointsperbin*sizeof(double));
	absorptionvdegree = malloc(numbins*maxpointsperbin*sizeof(double));
	numuseddegree = malloc(numbins*sizeof(int));
	
	int index;
	for(index=0;index<numbins;index++)
	numuseddegree[index] = 0;
	if(continuous)
	{
	UIsim->numcarbonsold = UIsim->numcarbons;
	g_timeout_add(33, (GSourceFunc)continuousremoves,NULL);
	}
}



void draw()
{
	//Creates the depiction of the GO sheet.
	
	int shownumbers = gtk_toggle_button_get_active(AtomNumberscheckbutton);
	int save = gtk_toggle_button_get_active(savecheckbutton);
	GtkAllocation *allocation = (GtkAllocation *)malloc(sizeof(GtkAllocation));		//Create the drawable pixmap
	gtk_widget_get_allocation(drawingarea, allocation);
	guint width = allocation->width;
	guint height = allocation->height;
//	printf("Drawing Dimensions %d x %d\n",width,height);
	cairo_t *cr;
	cairo_surface_t *surface;
	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,width, height);
	cr = cairo_create (surface);
	
	
	cairo_set_source_rgb (cr, 1.0f, 1.0f, 1.0f); //draw white rectangle
	cairo_rectangle(cr,0,0,width,height);
	cairo_fill(cr);
	
	int index,carbonclusterindex,carbonindex;
	double cairox,cairoy,cairodiameter;
	cairodiameter = 0.5 * fmin(width,height)/UIsim->dimensions;
	
	char buffer[1000];
	for(index=0;index<UIsim->numcarbons;index++)
	{
	if(*(UIsim->carboncovered+index))
	{
	cairo_set_source_rgb (cr, 0.4, 0.4, 0.4); //grey if covered
	cairox = *(UIsim->carbonlocations+2*index) * width/UIsim->boundary[0]+cairodiameter/2;
	cairoy = *(UIsim->carbonlocations+2*index+1) * height/UIsim->boundary[1]+cairodiameter/2;
	cairo_arc(cr, cairox, cairoy, cairodiameter, 0, 2*M_PI);
	cairo_fill(cr);
	if(shownumbers)
	{
	sprintf(buffer,"%d",index);
	cairo_move_to(cr, cairox, cairoy);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_show_text(cr,buffer);
	}
	}
	}
	double r,g,b,h,s=1,v=1;
	for(index = 0;index<UIsim->numclusters;index++)
	{
		if((UIsim->clusters+index)->color != -1)
		{
		h = (double)((UIsim->clusters+index)->color)/360;
		gtk_hsv_to_rgb(h,s,v,&r,&g,&b);
		cairo_set_source_rgb (cr, r, g, b);
		}
		else
		cairo_set_source_rgb (cr, 0, 0, 0);
		for(carbonclusterindex=0;carbonclusterindex<(UIsim->clusters+index)->numcarbons;carbonclusterindex++)
		{
			carbonindex = *(((UIsim->clusters+index)->carbons)+carbonclusterindex);
			cairox = *(UIsim->carbonlocations+2*carbonindex) * width/UIsim->boundary[0]+cairodiameter/2;
			cairoy = *(UIsim->carbonlocations+2*carbonindex+1) * height/UIsim->boundary[1]+cairodiameter/2;
			cairo_arc(cr, cairox, cairoy, cairodiameter, 0, 2*M_PI);
			cairo_fill(cr);
			if(debug)
			if(carbonindex >=UIsim->numcarbons)
			printf("Drawing carbon number %d is greater than num carbons %d\n",carbonindex,UIsim->numcarbons);
			if(shownumbers)
			{
			sprintf(buffer,"%d",carbonindex);
			cairo_move_to(cr, cairox, cairoy);
			cairo_set_source_rgb (cr, 1, 1, 1);
			cairo_show_text(cr,buffer);
			cairo_set_source_rgb (cr, 0, 0, 0);
		}
			if(debug)
			if(*(UIsim->carboncovered+carbonindex))
			{
			
			printf("Carbon %d is covered but still part of cluster\n",carbonindex);
			int debugindex;
			for(debugindex = 0;debugindex<(UIsim->clusters+index)->numcarbons;debugindex++)
			printf("Carbon in this cluster: %d\n",*((UIsim->clusters+index)->carbons+debugindex));
			}
		}
		
	}
	
	
	for(index=0;index<UIsim->numoxygens;index++)
	{
	cairo_set_source_rgb (cr, 0.6, 0.6, 0.6); 
	cairox = *(UIsim->oxygenlocations+2*index) * width/UIsim->boundary[0]+cairodiameter/2;
	cairoy = *(UIsim->oxygenlocations+2*index+1) * height/UIsim->boundary[1]+cairodiameter/2;
	cairo_arc(cr, cairox, cairoy, cairodiameter/2, 0, 2*M_PI);
	cairo_fill(cr);
	if(shownumbers)
	{
	sprintf(buffer,"%d",index);
	cairo_move_to(cr, cairox, cairoy);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_show_text(cr,buffer);
	}
	}
	
	for(index=0;index<UIsim->numhydroxides;index++)
	{
	cairo_set_source_rgb (cr, 0.8, 0.8, 0.8); 
	cairox = *(UIsim->hydroxidelocations+2*index) * width/UIsim->boundary[0]+cairodiameter/2;
	cairoy = *(UIsim->hydroxidelocations+2*index+1) * height/UIsim->boundary[1]+cairodiameter/2;
	cairo_arc(cr, cairox, cairoy, cairodiameter/2, 0, 2*M_PI);
	cairo_fill(cr);
	if(shownumbers)
	{
	sprintf(buffer,"%d",index);
	cairo_move_to(cr, cairox, cairoy);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_show_text(cr,buffer);
	}
	}
	
	cairo_t *cr2;	
	cr2 = gdk_cairo_create (gtk_widget_get_window(drawingarea));
	cairo_set_source_surface(cr2, surface, 0, 0);
	cairo_paint(cr2);
	
	if(save)
	{
		char buffer[1000];
		sprintf(buffer,"%s/image%d.png",workingdirectory,UIsim->stepnumber-1);
		cairo_surface_write_to_png(surface,buffer);
	}
	
	cairo_destroy(cr2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void pausenow()
{
	//Enter a paused state until keyboard input relieves it
	scanf("%s\n");
}

void addsimdata(struct sim *s)
{
	//Adds the data from a finished simulation to the UIsim which will hold the ensemble simulation
	double *emission;
double *absorption;
double *p1, *p2, *p3;
double *times;
	pthread_mutex_lock(&mutex);
	
	int index;
	int degreeindex;
	double degree;
	int clustersizeindex;
	for(index=0;index < s->stepnumber;index++)
	{
		*(UIsim->times+index) += 1/ *(s->times+index);
		*(UIsim->emissions + index*(UIsim->numsims)+s->simnumber) = *(s->emission+index);
		*(UIsim->absorptions+index*(UIsim->numsims)+s->simnumber) = *(s->absorption+index);
		(*(UIsim->numsimsattime+index)) ++;
		
		*(UIsim->p1+index) += *(s->p1+index);
		*(UIsim->p2+index) += *(s->p2+index);
		*(UIsim->p3+index) += *(s->p3+index);
		*(UIsim->emittingclusters+index) += *(s->emittingclusters+index);
		
		*(UIsim->rates+numrates*index) += *(s->rates+numrates*index);
		*(UIsim->rates+numrates*index+1) += 	*(s->rates+numrates*index+1); 
		*(UIsim->rates+numrates*index+2) += 	*(s->rates+numrates*index+2);  
		*(UIsim->rates+numrates*index+3) += 	*(s->rates+numrates*index+3);  
		*(UIsim->rates+numrates*index+4) += 	*(s->rates+numrates*index+4); 	 	

		
		degree = *(s->p1+index) + 2 * *(s->p2+index) + 3 * *(s->p3+index);
		degreeindex = degree/3 * numbins;
		if(numuseddegree[degreeindex] != maxpointsperbin-1)
		{
			emissionvdegree[maxpointsperbin*degreeindex+numuseddegree[degreeindex]] = *(s->emission+index);
			absorptionvdegree[maxpointsperbin*degreeindex+numuseddegree[degreeindex]] = *(s->absorption+index);
			numuseddegree[degreeindex]++;
		}
		for(clustersizeindex=0;clustersizeindex<largestclustercounted;clustersizeindex++)
		*(UIsim->clusterdistribution+index*largestclustercounted + clustersizeindex) += *(s->clusterdistribution+index*largestclustercounted + clustersizeindex);
	}
	(UIsim->stepnumber) = fmax(UIsim->stepnumber,s->stepnumber);
	
	pthread_mutex_unlock(&mutex);
}
	

void * simthread(struct sim *s)
{
	//Controls a single simulation in the ensemble
	
	s->state = 1;
	printf("Started thread %s\n",s->dir);
	system("date");

	initialize(s);
	printf("Sim %s Carbons: %d, Oxygens: %d, Hydroxides: %d\n",s->dir,s->numcarbons,s->numoxygens,s->numhydroxides);

	while(s->numcarbons >0 || s->numsourcesites > 0 || s->numhydroxides > 0 || s->numhops > 0)
	{
		
	if((double)s->numcarbonsold/s->numcarbons > 1.02 && s->numcarbonsold-s->numcarbons > 100)
	{
			printf("Sim %s Carbons: %d, Oxygens: %d, Hydroxides: %d\n",s->dir,s->numcarbons,s->numoxygens,s->numhydroxides);
			s->numcarbonsold = s->numcarbons;
	}	
		
	removes(s);
	clusterdetection(s);	
	}
	addsimdata(s);
	freesim(s);
	printf("Finished sim %s\n",s->dir);
	s->state = 2;
}

void processdegreebins(double stochasticfactor)
{

	double *emaverages = malloc(numbins*sizeof(double));
	double *emstddevs = malloc(numbins*sizeof(double));
	double *abaverages = malloc(numbins*sizeof(double));
	double *abstddevs = malloc(numbins*sizeof(double));
	
	char buffer[1000];
	sprintf(buffer,"%s/VersusDegree.dat",UIsim->dir);
	FILE *outfile = fopen(buffer,"w");
	fprintf(outfile,"Degree Emission Dev(Emission) Absorption Dev(Absorption)");
	int index;
	double degree;
	for(index=0;index<numbins;index++)
	{
		degree = 3 * (double)index/numbins;
		emaverages[index] = listavg(emissionvdegree+maxpointsperbin*index,numuseddegree[index]);
		emstddevs[index] = liststddev(emissionvdegree+maxpointsperbin*index,numuseddegree[index])*stochasticfactor;
		abaverages[index] = listavg(absorptionvdegree+maxpointsperbin*index,numuseddegree[index]);
		abstddevs[index] = liststddev(absorptionvdegree+maxpointsperbin*index,numuseddegree[index])*stochasticfactor;
		fprintf(outfile,"%g %g %g %g %g\n",degree,emaverages[index],emstddevs[index],abaverages[index],abstddevs[index]);
	}
	fclose(outfile);
	plotemversusdegree();
	plotabsversusdegree();
}

void processtimes(double stochasticfactor)
{

	char buffer[1000];
	sprintf(buffer,"%s/TimeVariability.dat",UIsim->dir);
	FILE *outfile = fopen(buffer,"w");
	fprintf(outfile,"Time Emission Dev(Emission) Absorption Dev(Absorption)");
	int index;
	for(index=0;index<UIsim->stepnumber;index++)
	{
		*(UIsim->emission+index) = listavg(UIsim->emissions+index*UIsim->numsims,*(UIsim->numsimsattime+index));
		*(UIsim->emissionstddev+index) = liststddev(UIsim->emissions+index*UIsim->numsims,*(UIsim->numsimsattime+index))*stochasticfactor;
		*(UIsim->absorption+index) = listavg(UIsim->absorptions+index*UIsim->numsims,*(UIsim->numsimsattime+index));
		*(UIsim->absorptionstddev+index) = liststddev(UIsim->absorptions+index*UIsim->numsims,*(UIsim->numsimsattime+index))*stochasticfactor;
		*(UIsim->QY+index) = 274.5 * *(UIsim->emission+index) * 2.3 / *(UIsim->absorption+index)/UIsim->initialcarbons;		

		fprintf(outfile,"%g %g %g %g %g\n",*(UIsim->times+index),*(UIsim->emission+index),*(UIsim->emissionstddev+index),*(UIsim->absorption+index),*(UIsim->absorptionstddev+index));
	}
	fclose(outfile);
}


gboolean ensembleaverage()
{
	printf("Starting Ensemble Averaging in directory %s\n",UIsim->dir);
	
	int numsims = gtk_spin_button_get_value(Ensemblespinbutton);
	pthread_t *threads = malloc(numsims*sizeof(pthread_t));
	struct sim *sims = malloc(numsims*sizeof(struct sim));
	//states 0 not started, 1 started, 2 finished
	
	if(emissionvdegree != NULL)
	free(emissionvdegree);
	if(absorptionvdegree != NULL)
	free(absorptionvdegree);
	if(numuseddegree != NULL)
	free(numuseddegree);
	emissionvdegree = malloc(numbins*maxpointsperbin*sizeof(double));
	absorptionvdegree = malloc(numbins*maxpointsperbin*sizeof(double));
	numuseddegree = malloc(numbins*sizeof(int));
	
	UIsim->emissions = malloc(maxsteps*numsims*sizeof(double));
	UIsim->absorptions = malloc(maxsteps*numsims*sizeof(double));
	UIsim->numsimsattime = malloc(maxsteps*numsims*sizeof(double));
	UIsim->emissionstddev = malloc(maxsteps*sizeof(double));
	UIsim->absorptionstddev = malloc(maxsteps*sizeof(double));
	
	int numcarbonsflake = 160000000;
	
	int index;
	for(index=0;index<numbins;index++)
	numuseddegree[index] = 0;
	
	
	UIsim->numsims= numsims;
	char *dir;
	struct sim zz= {0};
	for(index=0;index<numsims;index++)
	{
		*(sims+index) = zz;
		dir = malloc(1000*sizeof(char));
		sprintf(dir,"%d",index);
		(sims+index)->dir = dir;
		readbuttons(sims+index);
		(sims+index)->state = 0;
		(sims+index)->simnumber = index;
	}
	
	int clustersizeindex;
	for(index=0;index<maxsteps;index++)
	{
	*(UIsim->times+index) = 0;
	*(UIsim->absorption+index) = 0;
	*(UIsim->emission+index) = 0;
	*(UIsim->p1+index) = 0;
	*(UIsim->p2+index) = 0;
	*(UIsim->p3+index) = 0;
	*(UIsim->emittingclusters+index) = 0;
	*(UIsim->rates+index*numrates) = *(UIsim->rates+index*numrates+1) = *(UIsim->rates+index*numrates+2) = *(UIsim->rates+index*numrates+3) = *(UIsim->rates+index*numrates+4) = 0;
	*(UIsim->numsimsattime+index) = 0;

	for(clustersizeindex=0;clustersizeindex<largestclustercounted;clustersizeindex++)
	*(UIsim->clusterdistribution+index*largestclustercounted+clustersizeindex) = 0;
	}
	
	
	
	int nextsim = 0;
	int numactivesims = 0;
	void *retval;
	int returnvalue;
	
	
	while(numactivesims > 0 || nextsim <numsims)
	{
		numactivesims = 0;
		for(index=0;index<nextsim;index++)
		{
			if((sims+index)->state == 1)
			numactivesims++;
		}
	//	printf("Num active threads %d\n",numactivesims);
		while(numactivesims < maxthreads && nextsim<numsims)
		{
			pthread_create(threads + nextsim,NULL,simthread,sims+nextsim);
			nextsim++;
			numactivesims++;
		}
		sleep(1);
	}
	for(index=0;index<UIsim->stepnumber;index++)
	{
	*(UIsim->times+index) = 1/ *(UIsim->times+index);
	*(UIsim->p1+index) /= numsims;
	*(UIsim->p2+index) /= numsims;
	*(UIsim->p3+index) /= numsims;
	}
	double stochasticfactor = (double)numcarbonsflake/(sims->initialcarbons);

	stochasticfactor = 1/sqrt(stochasticfactor);
	
	processtimes(stochasticfactor);
	outputemission(UIsim);
	outputabsorption(UIsim);
	outputQY(UIsim);
	outputprobabilities(UIsim);
	outputclusters(UIsim);
	outputrates(UIsim);
	//outputdistribution(UIsim);
	//double stochasticfactor = (double)numcarbonsflake/(sims->initialcarbons)/numsims;

	//stochasticfactor = 1;
	printf("The stochastic factor is %g\n",stochasticfactor);
	
	
	plotQY();
	plotabsorption();
	plotprobabilities();
	plotabsems();
	plotemvsabs();
	plotclusters();
	plotrates();
	processdegreebins(stochasticfactor);
	//plotdistribution();
	printf("Simulations done\n");
	return FALSE;
}


int main ( int argc, char **argv )
{
	//Main function initializes connects and packs the UI
	gdk_threads_init();
	gtk_init(NULL,NULL);
	pthread_mutex_init(&mutex,NULL);
	
	srand(0);
	
	UIsim = malloc(sizeof(struct sim));
	struct sim zz = {0};
	*UIsim = zz;
	
	changedir();
	
	initializecallback(UIsim);
	
	
	DirEntry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(DirEntry),"default");
	g_signal_connect_swapped(DirEntry, "activate", G_CALLBACK (changedir),NULL);

	
	drawingarea = gtk_drawing_area_new ();
  gtk_widget_set_size_request (drawingarea, 480,480);

	GtkWidget *drawingbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(drawingbox),drawingarea,FALSE,FALSE,10);//Used to control dimensions of drawing area


	Periodiccheckbutton = gtk_check_button_new_with_label ("Periodic Boundary Conditions");
	g_signal_connect_swapped(Periodiccheckbutton, "toggled", G_CALLBACK (initializecallback),UIsim);
		
	Covercheckbutton = gtk_check_button_new_with_label ("Cover With Impurities");
	g_signal_connect_swapped(Periodiccheckbutton, "toggled", G_CALLBACK (initializecallback),UIsim);

	
	GtkWidget *dimensionstext = gtk_label_new("Dimensions");
	dimensionspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(50,0,1000, 1,1,0)),.01,0);
	GtkWidget *dimensionsbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(dimensionsbox),dimensionstext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(dimensionsbox),dimensionspin, FALSE, TRUE, 20);
	g_signal_connect_swapped(dimensionspin, "value-changed", G_CALLBACK (initializecallback),UIsim);
	g_signal_connect_swapped(dimensionspin, "value-changed", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(dimensionspin, "value-changed", G_CALLBACK (draw),NULL);
	
	
	GtkWidget *Otext = gtk_label_new("Oxygen Coverage");
	Ospin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.175,0,1, 0.001,0.001,0)),.001,3);
	GtkWidget *Obox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(Obox),Otext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(Obox),Ospin, FALSE, TRUE, 20);
	g_signal_connect_swapped(Ospin, "value-changed", G_CALLBACK (initializecallback),UIsim);
	g_signal_connect_swapped(Ospin, "value-changed", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(Ospin, "value-changed", G_CALLBACK (draw),NULL);
	
	
	GtkWidget *OHtext = gtk_label_new("Hydroxide Coverage");
	OHspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.3,0,1, 0.001,0.001,0)),.001,3);
	GtkWidget *OHbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(OHbox),OHtext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(OHbox),OHspin, FALSE, TRUE, 20);
	g_signal_connect_swapped(OHspin, "value-changed", G_CALLBACK (initializecallback),UIsim);
	g_signal_connect_swapped(OHspin, "value-changed", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(OHspin, "value-changed", G_CALLBACK (draw),NULL);
	
	
	GtkWidget *vacancytext = gtk_label_new("Vacancy Coverage");
	vacancyspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.1,0,1, 0.001,0.001,0)),.001,3);
	GtkWidget *vacancybox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(vacancybox),vacancytext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(vacancybox),vacancyspin, FALSE, TRUE, 20);
	g_signal_connect_swapped(vacancyspin, "value-changed", G_CALLBACK (initializecallback),UIsim);
	g_signal_connect_swapped(vacancyspin, "value-changed", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(vacancyspin, "value-changed", G_CALLBACK (draw),NULL);
	
	//GtkWidget *reversibleratetext = gtk_label_new("Reversible Rate");
	//reversibleratespin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0,0,100, 0.001,0.001,0)),.001,3);
	//GtkWidget *reversibleratebox = gtk_hbox_new(FALSE,0);
	//gtk_box_pack_start (GTK_BOX(reversibleratebox),reversibleratetext, FALSE, TRUE, 20);
	//gtk_box_pack_start (GTK_BOX(reversibleratebox),reversibleratespin, FALSE, TRUE, 20);
	
	
	GtkWidget *hydroxideremovaltext = gtk_label_new("Hydroxide Removal Coef.");
	hydroxideremoval = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(1000,0,1000000, 0.11,1,0)),.01,3);
	GtkWidget *hydroxideremovalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(hydroxideremovalbox),hydroxideremovaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(hydroxideremovalbox),hydroxideremoval, FALSE, TRUE, 20);
	
	
	GtkWidget *COremovaltext = gtk_label_new("CO Removal Coef.");
	COremoval = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(10,0,1000000, 0.11,1,0)),.01,3);
	GtkWidget *COremovalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(COremovalbox),COremovaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(COremovalbox),COremoval, FALSE, TRUE, 20);
	
	
	GtkWidget *C1removaltext = gtk_label_new("Singly-bonded Carbon Removal Coef.");
	C1removal = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(1,0,1000000, 0.11,1,0)),.001,3);
	GtkWidget *C1removalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(C1removalbox),C1removaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(C1removalbox),C1removal, FALSE, TRUE, 20);
	
	GtkWidget *Czigremovaltext = gtk_label_new("Zig-Zag Edge Removal Coef.");
	Czigremoval = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.01,0,1000000, 0.11,1,0)),.001,3);
	GtkWidget *Czigremovalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(Czigremovalbox),Czigremovaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(Czigremovalbox),Czigremoval, FALSE, TRUE, 20);
	
	GtkWidget *Carmremovaltext = gtk_label_new("Armchair Edge Removal Coef.");
	Carmremoval = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.1,0,1000000, 0.11,1,0)),.001,3);
	GtkWidget *Carmremovalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(Carmremovalbox),Carmremovaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(Carmremovalbox),Carmremoval, FALSE, TRUE, 20);
	
	GtkWidget *C3removaltext = gtk_label_new("Triply-bonded Carbon Removal Coef.");
	C3removal = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0.001,0,1000000, 0.11,1,0)),.001,3);
	GtkWidget *C3removalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(C3removalbox),C3removaltext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(C3removalbox),C3removal, FALSE, TRUE, 20);
	
	
	GtkWidget *Ohoptext = gtk_label_new("Oxygen Hopping Coef.");
	Ohop = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0,0,1000000, 0.1,1,0)),.01,3);
	GtkWidget *Ohopbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(Ohopbox),Ohoptext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(Ohopbox),Ohop, FALSE, TRUE, 20);
	
	GtkWidget *OSourcetext = gtk_label_new("Oxygen Source Coef.");
	Osource= gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(0,0,1000000, 0.11,1,0)),.01,3);
	GtkWidget *OSourcebox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(OSourcebox),OSourcetext, FALSE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(OSourcebox),Osource, FALSE, TRUE, 20);
	
	
	savecheckbutton = gtk_check_button_new_with_label ("Make Video");
	
	
	GtkWidget *playbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(playbutton),gtk_image_new_from_stock(GTK_STOCK_REFRESH,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(playbutton, "clicked", G_CALLBACK (initializecallback),UIsim);
	g_signal_connect_swapped(playbutton, "clicked", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(playbutton, "clicked", G_CALLBACK (draw),NULL);

	
	GtkWidget *removebutton = gtk_button_new_with_label("TimeStep");
	GtkWidget *continuousremovebutton = gtk_button_new_with_label("ContinuousSteps");
	numbertoremove = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(1,0,100000, 1,1,0)),.01,0);
	GtkWidget *removalbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(removalbox),removebutton, TRUE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(removalbox),continuousremovebutton, TRUE, TRUE, 20);
	gtk_box_pack_start (GTK_BOX(removalbox),numbertoremove, TRUE, TRUE, 20);
	g_signal_connect_swapped(removebutton, "clicked", G_CALLBACK (readbuttons),UIsim);
	g_signal_connect_swapped(removebutton, "clicked", G_CALLBACK (removes),UIsim);
	g_signal_connect_swapped(removebutton, "clicked", G_CALLBACK (clusterdetection),UIsim);
	g_signal_connect_swapped(removebutton, "clicked", G_CALLBACK (draw),NULL);
	g_signal_connect_swapped(continuousremovebutton, "clicked", G_CALLBACK (startcycle),NULL);

		
	GtkWidget *ensemblebox = gtk_hbox_new(FALSE,0);
	Ensemblespinbutton = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new(100,1,100000, 1,10,0)),1,0);
	GtkWidget *ensemblebutton = gtk_button_new_with_label("Begin Ensemble Calculation");
	g_signal_connect_swapped(ensemblebutton, "clicked", G_CALLBACK (ensembleaverage),NULL);
	gtk_box_pack_start (GTK_BOX(ensemblebox),Ensemblespinbutton, TRUE, TRUE, 20);	
	gtk_box_pack_start (GTK_BOX(ensemblebox),ensemblebutton, TRUE, TRUE, 20);	
	
	
	GtkWidget *adjacencybutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(adjacencybutton),gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (printadjacency),UIsim);
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (printhops),UIsim);
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (printcarbons),UIsim);
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (printoxygens),UIsim);
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (checknumcovered),UIsim);
	g_signal_connect_swapped(adjacencybutton, "clicked", G_CALLBACK (hopcheck),UIsim);

	AtomNumberscheckbutton = gtk_check_button_new_with_label ("Display Numbers on Atoms");
	Debugcheckbutton = gtk_check_button_new_with_label ("Debugging Mode");
	VerboseDebugcheckbutton = gtk_check_button_new_with_label ("Very Verbose Debugging");

	
	g_signal_connect_swapped(drawingarea,"draw",G_CALLBACK (draw),NULL);
	
	GtkWidget *vbox = gtk_vbox_new(FALSE,0);
	GtkWidget *hbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(vbox),DirEntry, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Periodiccheckbutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Covercheckbutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),dimensionsbox, FALSE, TRUE, 2);	
	gtk_box_pack_start (GTK_BOX(vbox),Obox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),OHbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),vacancybox, FALSE, TRUE, 2);
	//gtk_box_pack_start (GTK_BOX(vbox),reversibleratebox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),playbutton, FALSE, TRUE, 2);

	gtk_box_pack_start (GTK_BOX(vbox),hydroxideremovalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),COremovalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),C1removalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Carmremovalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Czigremovalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),C3removalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Ohopbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),OSourcebox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),savecheckbutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),removalbox, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),ensemblebox, FALSE, TRUE, 2);
	

	gtk_box_pack_start (GTK_BOX(vbox),adjacencybutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),AtomNumberscheckbutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),Debugcheckbutton, FALSE, TRUE, 2);
	gtk_box_pack_start (GTK_BOX(vbox),VerboseDebugcheckbutton, FALSE, TRUE, 2);

	
	gtk_box_pack_start (GTK_BOX(hbox),drawingbox, FALSE, FALSE, 10);
	gtk_box_pack_start (GTK_BOX(hbox),vbox, FALSE, TRUE, 10);

	
	TopLevel = gtk_window_new (GTK_WINDOW_TOPLEVEL); 
	gtk_container_add (GTK_CONTAINER (TopLevel), hbox);
	
	gtk_widget_show_all(TopLevel);
	g_signal_connect_swapped (TopLevel, "delete_event", G_CALLBACK (exit), 0);
		
	gtk_main();
}
