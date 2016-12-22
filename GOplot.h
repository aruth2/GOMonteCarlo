void printpm3dfile(char *filename,int *x, double *y,double *z, double *matrix1,int *matrix2, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	printf("Printing pm3d file\n");
	FILE *outfile = fopen(filename,"w");
	int i;
	int step, time;
	double epsilon = 1e-9;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
			step = i%numberOfColumns;
			time = i/numberOfColumns;
				if((time*clusterlinesamples/numberOfRows - (time-1)*clusterlinesamples/numberOfRows == 1)
			 || ((step)*clusterlinesamples/numberOfColumns - (step-1)*clusterlinesamples/numberOfColumns == 1))
			{
			 fprintf(outfile,"%d\t%g\t%g\t%g\t%d\n",*(x+i),*(y+i),*(z+i),emissionofcluster(i%numberOfColumns) * *(matrix1+i),150);			
			//printf("step %d time %d\n",step,time);
		 }
			 else
			fprintf(outfile,"%d\t%g\t%g\t%g\t%d\n",*(x+i),*(y+i),*(z+i),emissionofcluster(i%numberOfColumns) * *(matrix1+i),*(matrix2+i));
//			if(*(matrix1+i) != 0)
//			printf("Stepnumber %d cluster size %d\n",i/numberOfColumns,i%numberOfColumns);
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
	fclose(outfile);
}

void printsurfacefile(char *filename,int *matrix1, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	FILE *outfile = fopen(filename,"w");
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			fprintf(outfile,"%d %d %d\n",i%numberOfColumns,i/numberOfColumns,*(matrix1+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
	fclose(outfile);
}



void outputdistribution(struct sim *s)
{
	printf("Printing distribution to file for %d steps and %d largest cluster\n",s->stepnumber,largestclustercounted);
	char buffer[1000];
	sprintf(buffer,"%s/distribution.dat",s->dir);
	int *sizematrix = malloc((s->stepnumber) * largestclustercounted*sizeof(int));
	double *timematrix = malloc((s->stepnumber) * largestclustercounted*sizeof(double));
	double *degreematrix = malloc((s->stepnumber) * largestclustercounted*sizeof(double));
	int i,j;
	
	if(sizematrix == NULL || timematrix == NULL || s->times == NULL)
	printf("NULL MATRIX\n");
	for(i = 0;i<s->stepnumber;i++)
	for(j=0;j<largestclustercounted;j++)
	{
	*(sizematrix + i*largestclustercounted+j) = j;
	*(timematrix+i*largestclustercounted+j) = *(s->times+i);
	*(degreematrix+i*largestclustercounted+j) = *(s->p1 +i) + *(s->p2+i)*2 + *(s->p3+i)*3;
	}
	printf("Printing file\n");
	printpm3dfile(buffer,sizematrix,timematrix,degreematrix,s->clusterdistribution,sizematrix,s->stepnumber,largestclustercounted);
	//sprintf(buffer,"%s/surface.dat",s->dir);
	//printsurfacefile(buffer,s->clusterdistribution,s->stepnumber,largestclustercounted);
	//free(sizematrix);
	free(sizematrix);
	free(timematrix);
	free(degreematrix);
}

void outputemission(struct sim *s)
{
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/emission.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g\n",log(*(s->times+index)),*(s->emission+index));
	fclose(outfile);

}

void outputQY(struct sim *s)
{
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/QY.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g\n",log(*(s->times+index)),*(s->QY+index));
	fclose(outfile);

}

void outputabsorption(struct sim *s)
{
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/absorption.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g\n",log(*(s->times+index)),*(s->absorption+index));
	fclose(outfile);

}

void outputprobabilities(struct sim *s)
{
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/probs.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g %g %g %g\n",log(*(s->times+index)),*(s->p1+index),*(s->p2+index), *(s->p3+index),*(s->p1+index) + 2 * *(s->p2+index) + 3 * *(s->p3+index));
	fclose(outfile);
}

void outputclusters(struct sim *s)
{
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/clusters.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g\n",log(*(s->times+index)),*(s->emittingclusters+index));
	fclose(outfile);

}

void normalizerates(struct sim *s)
{
	int rateindex,stepindex;
	double sum;
	for(rateindex=0;rateindex<numrates;rateindex++)
	{
		sum = 0;
	for(stepindex=0;stepindex<s->stepnumber;stepindex++)
	sum+= *(s->rates+stepindex*numrates+rateindex);
	printf("Sum %g\n",sum);
	for(stepindex=0;stepindex<s->stepnumber;stepindex++)
	*(s->rates+stepindex*numrates+rateindex) /= sum;
	}
}

void outputrates(struct sim *s)
{
	printf("Printing Rates\n");
	//normalizerates(s);
	int index;
	char buffer[1000];
	sprintf(buffer,"%s/rates.dat",s->dir);
	FILE *outfile = fopen(buffer,"w");
	for(index = 0;index<s->stepnumber;index++)
	fprintf(outfile,"%g %g %g %g %g %g %g %g %g\n",log(*(s->times+index)),*(s->rates+index*numrates), *(s->rates+index*numrates+1), *(s->rates+index*numrates+2), *(s->rates+index*numrates+3), *(s->rates+index*numrates+4),*(s->rates+index*numrates+5),*(s->rates+index*numrates+6),*(s->rates+index*numrates+7));
	fclose(outfile);
}

void plotQY()
{
	printf("Plotting QY\n");
	char *gnuplotscript = 
	"set term postscript eps color enhanced \"Helvetica\" 24\n"
	"set output \"QY.eps\"\n"
	"set ylabel \"QY(%)\"\n"
	"set xlabel \"Time(log scale)\"\n"
	"set title \"Quantum Yield\"\n"
	"set key noautotitles\n"
	"plot 'QY.dat' u 1:2 w lines lc 'black'\n";
	
	
char buffer[1000];
sprintf(buffer,"%s/QYplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript);	
	
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot QYplotscript",UIsim->dir);
system(buffer);
}

void plotabsorption()
{
		printf("Plotting absorption with %d steps\n",UIsim->stepnumber);
	char *gnuplotscript = 
	"set term postscript eps color enhanced \"Helvetica\" 24\n"
	"set output \"absorption.eps\"\n"
	"set ylabel \"Absorption(arb. units)\"\n"
	"set xlabel \"Time(log scale)\"\n"
	"set title \"Absorption\"\n"
	"unset ytics\n"
	"set key noautotitles\n";
	
char buffer[1000];
sprintf(buffer,"%s/absplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript);
	
	fprintf(outfile,"plot \"absorption.dat\" u 1:2 w lines lc 'black'\n");
	
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot absplotscript",UIsim->dir);
system(buffer);
}

void plotdistribution()
{
printf("Plotting distribution\n");
	
char *gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 24\n"
"set output 'distribution.eps'\n"
"set ylabel 'Time(Arb. Units)'\n"
"set xlabel 'Cluster Size'\n"
"set xlabel  rotate parallel\n"
"set xlabel offset -3\n"
"set zlabel 'Cluster Distribution'\n"
"set zlabel rotate parallel\n"
"set zlabel offset +3\n"
"set format z ''\n"
"set view 30,80\n"
"set cbrange [54:275]\n"
"set logscale y\n"
"unset key\n"
"set contour surface\n"
"set cntrparam levels auto 10\n"
"set hidden3d\n"
"set palette defined (-0.001 'grey', 0 'violet', 1 'blue', 2 'cyan', 3 'green',4 'yellow', 5 'orange', 6 'red', 6.001 'black') \n"
"splot 'distribution.dat' u ($1):($2*exp(13)):($4):($5) w pm3d";

char buffer[1000];
sprintf(buffer,"%s/distributionplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");

fprintf(outfile,gnuplotscript);
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot distributionplotscript",UIsim->dir);
system(buffer);
	
	gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 24\n"
"set output 'distributionvsdegree.eps'\n"
"set ylabel 'Coordination Number <d>'\n"
"set xlabel 'Cluster Size'\n"
"set xlabel  rotate parallel\n"
"set xlabel offset -3\n"
"set zlabel 'Cluster Distribution'\n"
"set zlabel rotate parallel\n"
"set zlabel offset +3\n"
"set format z ''\n"
"set view 30,80\n"
"set cbrange [54:275]\n"
"unset key\n"
"set contour surface\n"
"set cntrparam levels auto 10\n"
"set hidden3d\n"
"set palette defined (-0.001 'grey', 0 'violet', 1 'blue', 2 'cyan', 3 'green',4 'yellow', 5 'orange', 6 'red', 6.001 'black') \n"
"splot 'distribution.dat' u ($1):($3):($4):($5) w pm3d";

sprintf(buffer,"%s/distributionvsdegreeplotscript",UIsim->dir);
outfile = fopen(buffer,"w");

fprintf(outfile,gnuplotscript);
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot distributionvsdegreeplotscript",UIsim->dir);
system(buffer);
}

void plotemvsabs()
{
	char *gnuplotscript = 
	"set term postscript eps color enhanced 'Helvetica' 36\n"
	"set output 'emvsabs.eps'\n"
	"set ylabel 'Emission (Arb. Units)'\n"
	"set xlabel 'Absorption(%)'\n"
	"set format y ''\n"
	"set key noautotitles\n"
	"plot 'VersusDegree.dat' u 4:2 w points pt 7 ps 1.5 lc 'black' \n";
	
	char buffer[1000];
	sprintf(buffer,"%s/emvsabsplotscript",UIsim->dir);
	FILE *outfile = fopen(buffer,"w");

	fprintf(outfile,gnuplotscript);
	fclose(outfile);
	sprintf(buffer,"cd %s && gnuplot emvsabsplotscript",UIsim->dir);
	system(buffer);
	
}

void plotabsems()
{
printf("Plotting absems\n");
	
char *gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 36\n"
"set output 'absems.eps'\n"
"set logscale x\n"
"set xlabel 'Time(Arb. Units)'\n"
"set ytics nomirror\n"
"set ylabel 'Absorption(%)'\n"
"set y2label 'Emission'\n"
"set y2tics\n"
"set key noautotitles\n"
"plot 'absorption.dat' u (exp($1+13)):($2) axes x1y1 w lines lt 3, 'emission.dat' u (exp($1+13)):($2) axes x1y2 w lines lt 1";

char buffer[1000];
sprintf(buffer,"%s/absemsplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");

fprintf(outfile,gnuplotscript);
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot absemsplotscript",UIsim->dir);
system(buffer);
	
}

void plotabsemsandcombine(char *outfilename)
{
//This function creates an image for a single frame for the video. Currently it requires hardcoding the xrange and yrange.
//It should be updated to allow for a passed variable to set the range and to run the simulation twice.

char *gnuplotscript = 
"set term png font Helvetica 24\n"
"set output '%s'\n"
"set logscale x\n"
"set xlabel 'Time(a. u.)'\n"
"set xrange [1000000:3000000]\n"
"set yrange [0:1.5]\n"
"set y2range [0:7.5]\n"
"set ytics nomirror\n"
"set label 1 'Absorption(%)' at screen 0.05,0.33 rotate by 90 front nopoint tc rgb 'blue'\n"
"set label 2 'Emission' at screen 0.95,0.4 rotate by 90 front nopoint tc rgb 'red'\n"
"unset xtics\n"
"set ytics nomirror\n"
"set y2tics\n"
"set key noautotitles\n"
"plot 'absorption.dat' u (exp($1+9)):($2) axes x1y1 w lines lt 3 lw 5 lc rgb 'blue', 'emission.dat' u (exp($1+9)):($2) axes x1y2 w lines lt 1 lw 5 lc rgb 'red'";

char buffer[1000];
sprintf(buffer,"%s/absemscombineplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript,outfilename);
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot absemscombineplotscript",UIsim->dir);
system(buffer);
}

void plotclusters()
{
printf("Plotting clusters\n");
	
char *gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 36\n"
"set output 'clusters.eps'\n"
"set logscale x\n"
"set xlabel 'Time(Arb. Units)'\n"
"set ylabel 'NumClusters'\n"
"set key noautotitles\n"
"plot 'clusters.dat' u (exp($1+13)):($2) w lines\n";

char buffer[1000];
sprintf(buffer,"%s/clustersplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");

fprintf(outfile,gnuplotscript);
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot clustersplotscript",UIsim->dir);
system(buffer);	
}

void plotprobabilities()
{
	printf("Plotting probabilities\n");	
	char *gnuplotscript = 
	"set term postscript eps color enhanced \"Helvetica\" 36\n"
	"set output \"probs.eps\"\n"
	"set ylabel \"Number\"\n"
	"set xlabel \"Time(log scale)\"\n"
	"set title \"Area/Edge Comparison\"\n"
	"set key noautotitles\n";
	
char buffer[1000];
sprintf(buffer,"%s/probsplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
	fprintf(outfile,gnuplotscript);
	
	fprintf(outfile,"plot \"probs.dat\" u 1:2 w lines title \"P1\",\"probs.dat\" u 1:3 w lines title \"P2\",\"probs.dat\" u 1:4 w lines title \"P3\",\"probs.dat\" u 1:5 w lines title \"<d>\" \n");
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot probsplotscript",UIsim->dir);
system(buffer);
}

void plotemversusdegree()
{
printf("Plotting Em v degree\n");
		
char *gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 36\n"
"set output 'emversusdegree.eps'\n"
"set xlabel 'Coordination Number <d>'\n"
"set ylabel 'Emission(Arb. Units)'\n"
"set key noautotitles\n"
"plot 'VersusDegree.dat' u 1:2:3 every 10 w errorbars ps 1 lc 'black'\n";

char buffer[1000];
sprintf(buffer,"%s/emvsdegreeplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript);
	
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot emvsdegreeplotscript",UIsim->dir);
system(buffer);	
}
void plotabsversusdegree()
{
	printf("Plotting Abs v degree\n");

	char *gnuplotscript = 
"set term postscript eps color enhanced \"Helvetica\" 36\n"
"set output 'absversusdegree.eps'\n"
"set xlabel 'Coordination Number <d>'\n"
"set ylabel 'Absorption(%)'\n"
"set key noautotitles\n"
"plot 'VersusDegree.dat' u 1:4:5 every 10 w errorbars ps 1 lc 'black'\n";
	
char buffer[1000];
sprintf(buffer,"%s/absvsdegreeplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript);
	
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot absvsdegreeplotscript",UIsim->dir);
system(buffer);
		
}

void plotrates()
{
	printf("Plotting Rates\n");
	char *gnuplotscript = 

"set term postscript eps color enhanced 'Helvetica' 24\n"
"set output 'rates.eps'\n"
"set xlabel 'Time (Arb. Units)'\n"
"set ylabel 'Rate'\n"
"set logscale y\n"
"set logscale x\n"
"set format y '%%.0e'\n"
"set xrange [0.025:3000]\n"
"set key noautotitles\n"
"set key bottom left\n"
"t = 15.5\n"
"plot 'rates.dat' u (exp($1+t)):($2) w lines title 'OH removal' dashtype 1 lw 5"
", 'rates.dat' u (exp($1+t)):($3) w lines title 'CO removal' dashtype 2 lw 7"
", 'rates.dat' u (exp($1+t)):($4) w lines title 'C1 removal' dashtype 3 lw 7"
", 'rates.dat' u (exp($1+t)):($5) w lines title 'C3 removal' dashtype 4 lw 7"
", 'rates.dat' u (exp($1+t)):($6) w lines title 'Czig removal' dashtype 5 lw 7"
", 'rates.dat' u (exp($1+t)):($7) w lines title 'Carm removal' dashtype 6 lw 7"
", 'rates.dat' u (exp($1+t)):($8) w lines title 'O hopping' dashtype 7 lw 7"
", 'rates.dat' u (exp($1+t)):($9) w lines title 'C hopping' dashtype 8 lw 7";

char buffer[1000];
sprintf(buffer,"%s/ratesplotscript",UIsim->dir);
FILE *outfile = fopen(buffer,"w");
fprintf(outfile,gnuplotscript);
	
fclose(outfile);
sprintf(buffer,"cd %s && gnuplot ratesplotscript",UIsim->dir);
system(buffer);
	
	
}
