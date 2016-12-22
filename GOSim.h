#define maxconnections 3
#define largestclustercounted 400
#define clusterminsize 20
#define clusterlinesamples 10
double lv1[2] = {1.5,sqrt(3)/2};
double lv2[2] = {0,sqrt(3)};
double basis[2] = {1,0};
double oxygendistance = sqrt(1.76);//lv1 - basis along different direction/2 
double epsilon = 0.01;

double hydroxideprobability;
double COprobability;
double C1probability;
double Czigprobability;
double Carmprobability;
double C3probability;
double Ohopprobability;
double Osourceprobability;
int numremoves;
int numrates = 8;
int cover;

struct cluster
{
	int *carbons;
	int numcarbons;
	int color;
	
};

struct hop
{
	double Olocnew[2];
	double Olocold[2];
	int Ccovered;
	int Cuncovered;
	int Cunchanged;
	int oxygen;
};
#define maxsteps 1000000
struct sim
{
int state; // 0 not started, 1 running, 2 finished
double *carbonlocations;
double *oxygenlocations;
double *hydroxidelocations;
double *emission;
double *absorption;

double *emissions;
double *emissionstddev;
double *absorptions;
double *absorptionstddev;
int *numsimsattime;
int numsims;
int simnumber;

double *p1, *p2, *p3;
double *times;
double *emittingclusters;
int numcarbons;
int numcovered;
int numoxygens;
int numhydroxides;
int dimensions;
int *carboncovered;
int initialcarbons;
int periodic;
int *adjacencylist;
int *numadjacent;
int *numoxidizedadjacent;
int stepnumber;
int numcarbonsold;
char *dir;
struct cluster *clusters;
int numclusters;
struct hop *hops;
int numhops;
int numsourcesites;
double boundary[2];

double Oprobability;
double OHprobability;
double vacancycoverage;
double reversiblerate;

double *rates;
double *clusterdistribution;
double *QY;

//double *originalcarbonlocations;
//int *originaladjacencylist;
//int *availableunoccupiedsites;
//int *adjacencylistmap;
};//Split up this way so that multiple simulations can be ran in parallel

double emissionofcluster(int clustersize)
{
	return clustersize * exp(-0.44 * sqrt(clustersize));
}

void freeclusters(struct sim *s,int dimension)
{
	int index1;
	if(s->clusters != NULL)
	{
	for(index1=0;index1<s->initialcarbons;index1++)//Must use old dimensions
	free((s->clusters + index1)->carbons);
	free(s->clusters);
	}
}

void freesim(struct sim *s)
{
	
	
	
	if(s->adjacencylist != NULL)
	free(s->adjacencylist);
//	if(s->originaladjacencylist != NULL)
//	free(s->originaladjacencylist);
	if(s->numadjacent != NULL)
	free(s->numadjacent);
	if(s->hops != NULL)
	{
		free(s->hops);
	}
	if(s->emission !=NULL)
	free(s->emission);
	if(s->absorption !=NULL)
	free(s->absorption);
	if(s->p1 !=NULL)
	free(s->p1);
	if(s->p2 !=NULL)
	free(s->p2);
	if(s->p3 !=NULL)
	free(s->p3);
	if(s->times != NULL)
	free(s->times);
	if(s->rates != NULL)
	free(s->rates);
	if(s->emittingclusters != NULL)
	free(s->emittingclusters);
	if(s->clusters != NULL)
	freeclusters(s,s->dimensions);
	if(s->clusterdistribution != NULL)
	free(s->clusterdistribution);
}

void initialize(struct sim *s)
{
	freesim(s);
	printf("Making new board with dimensions %d and O prob. %g\n",s->dimensions,s->Oprobability);

	int index1,index2,index3;
	
	int n = 2*s->dimensions/3;//truncated
	int m = s->dimensions/sqrt(3);//truncated
	s->boundary[0] = 3 * n / 2;
	s->boundary[1] = m * sqrt(3);
	
	
	double x,y;
	s->numcarbons = 0;
	for(index1 = -s->dimensions;index1<=s->dimensions;index1++)
		for(index2 = -s->dimensions;index2<=s->dimensions;index2++)
			for(index3 = 0;index3<=1;index3++)
			{
				x = index1 * *lv1 + index2* *lv2 + index3 * *basis;
				y = index1 * *(lv1+1) + index2 * *(lv2+1) + index3* *(basis+1);
				//printf("x %g y %g\n",x,y);
				if(x>= 0 && x <s->boundary[0]-epsilon && y >= 0 && y <s->boundary[1]-epsilon)
				s->numcarbons++;
			}
	
	s->carbonlocations = malloc(2 * s->numcarbons * sizeof(double));
	s->oxygenlocations = malloc(2 * s->numcarbons * sizeof(double));
	s->hydroxidelocations = malloc(2 * s->numcarbons * sizeof(double)); //give them same amount of memory
	s->carboncovered = malloc(s->numcarbons*sizeof(int));


	s->adjacencylist = malloc(maxconnections * s->numcarbons*sizeof(int));
//	s->originaladjacencylist = malloc(maxconnections * s->numcarbons*sizeof(int));
	s->numadjacent = malloc(s->numcarbons*sizeof(int));
	s->numoxidizedadjacent = malloc(s->numcarbons*sizeof(int));
	for(index1 = 0;index1<s->numcarbons;index1++)
	{
	*(s->numadjacent+index1)=0;
	*(s->numoxidizedadjacent+index1)=0;
	}
	printf("Board will have %d carbons\n",s->numcarbons);
	
	s->numcarbons = 0;
	double progresscheck;
	for(index1 = -s->dimensions;index1<=s->dimensions;index1++)
		for(index2 = -s->dimensions;index2<=s->dimensions;index2++)
			for(index3 = 0;index3<=1;index3++)
			{
				x = index1 * *lv1 + index2* *lv2 + index3 * *basis;
				y = index1 * *(lv1+1) + index2 * *(lv2+1) + index3* *(basis+1);
				if(x>= 0 && x <s->boundary[0]-epsilon && y >= 0 && y <s->boundary[1]-epsilon)
				{
					*(s->carbonlocations+2*s->numcarbons) = x;
					*(s->carbonlocations+2*s->numcarbons+1) = y;
					*(s->carboncovered + s->numcarbons) = 0;
					s->numcarbons++;
					//addadjacency(s,s->numcarbons-1);
					addcarbon(s,s->numcarbons-1);
					if(debug)
					{
					progresscheck = log2(s->numcarbons); //Just here to keep user informed of progress
					if((int)progresscheck == progresscheck)
					printf("Num carbons %d\n",s->numcarbons);
					}
				}
			}
			if(verbosedebug)
			printadjacency(s);
	//printf("readying adjacency\n");
	//for(index1=0;index1<numcarbons;index1++)
		if(debug)
		checkadjacency(s);
	s->initialcarbons = s->numcarbons;
	
	int selectedcarbon;
	int selectedneighbor;
	
	double distance;
	s->numcovered = 0;
	s->numoxygens = 0;
	s->numhydroxides = 0;
	s->stepnumber=0;
	int searchstart;
	double rand;
	//checknumcovered();
	
	//for(index=0;index<s->numcarbons;index++)
	//{
	//	*(s->originaladjacencylist+index) = *(s->adjacencylist+index);
	//}
	
	
	if(s->vacancycoverage+s->Oprobability+s->OHprobability >1)
	{
	printf("Too many carbons are covered or vacant. Reducing OH probability to match\n");
	s->OHprobability = 1- s->vacancycoverage-s->Oprobability;
	}
	
	int numvacancies = s->numcarbons * s->vacancycoverage;
	int index;
	for(index=0;index<numvacancies;index++)//Create Vacancies
	{
		selectedcarbon = generaterandom() * s->numcarbons;
		removeC(s,selectedcarbon);
	}
	
	//s->Oprobability = s->Oprobability/(1-s->Oprobability);
	double Oselectionchance = s->Oprobability/(s->Oprobability+2*s->OHprobability);
	
	if(cover)
	while(s->numcovered <s->initialcarbons * (s->Oprobability + s->OHprobability))
	{
		//printf("Carbons: %d, Oxygens: %d, Hydroxides: %d\n",numcarbons,numoxygens,numhydroxides);
		//checknumcovered();
		searchstart = generaterandom()*s->numcarbons;
		rand = generaterandom();
		selectedcarbon = -1;
		for(;selectedcarbon != searchstart;selectedcarbon++,selectedcarbon %= s->numcarbons)
			{
				if(selectedcarbon == -1)
				selectedcarbon = searchstart;
				//printf("Checking carbon %d of %d using %d\n",selectedcarbon,numcarbons,searchstart);
				
				//if(selectedcarbon >numcarbons)
				//printf("error\n");
		//if(rand<s->Oprobability)//add oxygen
		if(rand<Oselectionchance)
		{

				if(!*(s->carboncovered+selectedcarbon) && *(s->numadjacent+selectedcarbon) -*(s->numoxidizedadjacent+selectedcarbon) > 0)
				{
					selectedneighbor = generaterandom() * (*(s->numadjacent+selectedcarbon)-*(s->numoxidizedadjacent+selectedcarbon));
					selectedneighbor = *(s->adjacencylist+ maxconnections*selectedcarbon + selectedneighbor);//Assumes the lowest entries in the list are the unoxidized ones
					
					placeoxygen(s,selectedcarbon,selectedneighbor,s->oxygenlocations+2*s->numoxygens);
					//*(oxygenlocations + 2*numoxygens) = (*(carbonlocations + 2*selectedcarbon) + *(carbonlocations+2*selectedneighbor))/2;
					//*(oxygenlocations + 2*numoxygens+1) = (*(carbonlocations + 2*selectedcarbon+1) + *(carbonlocations+2*selectedneighbor+1))/2;
					s->numoxygens++;
					oxidizecarbon(s,selectedneighbor);
					oxidizecarbon(s,selectedcarbon);
					//removeadjacency(s,selectedneighbor);
					//removeadjacency(s,selectedcarbon);
					s->numcovered += 2;
					break;
				}
			
		}
		else//add hydrogen
		{
					if(!*(s->carboncovered+selectedcarbon))
					{
					*(s->hydroxidelocations + 2*s->numhydroxides) = *(s->carbonlocations + 2*selectedcarbon);
					*(s->hydroxidelocations + 2*s->numhydroxides+1) = *(s->carbonlocations + 2*selectedcarbon+1);
					s->numhydroxides++;
					//removeadjacency(s,selectedcarbon);
					oxidizecarbon(s,selectedcarbon);
					s->numcovered++;
					break;
					}
			
			}
		}
	}
	
	printf("Carbons: %d, Oxygens: %d, Hydroxides: %d\n",s->numcarbons,s->numoxygens,s->numhydroxides);
	
	s->numhops = 0;
	s->clusters = malloc(s->initialcarbons*sizeof(struct cluster));
	for(index1=0;index1<s->initialcarbons;index1++)
	(s->clusters + index1)->carbons = malloc(s->numcarbons*sizeof(int));
	s->hops = malloc(4*(s->numoxygens+10)*sizeof(struct hop));
	s->emission = malloc(maxsteps*sizeof(double));
	s->absorption = malloc(maxsteps*sizeof(double));
	s->p1 = malloc(maxsteps*sizeof(double));
	s->p2 = malloc(maxsteps*sizeof(double));
	s->p3 = malloc(maxsteps*sizeof(double));
	s->times = malloc(maxsteps*sizeof(double));
	s->emittingclusters = malloc(maxsteps*sizeof(double));
	s->rates = malloc(maxsteps*numrates*sizeof(double));
	s->clusterdistribution = malloc(maxsteps * largestclustercounted * sizeof(double));
	s->QY = malloc(maxsteps*sizeof(double));
	}
	
	void placeoxygen(struct sim *s,int C1,int C2, double *Oarray)
{
	if(s->periodic)
	{
		double xdis = fabs(*(s->carbonlocations+2*C1)-*(s->carbonlocations+2*C2));
		double xsum = *(s->carbonlocations+2*C1)+*(s->carbonlocations+2*C2);
		if(xdis > 1.5)
			if(xsum > s->boundary[0])
				*(Oarray) = (xsum - s->boundary[0])/2;
			else
				*(Oarray) = (xsum + s->boundary[0])/2;
		else
		{
			*(Oarray) = xsum/2;
		}
		double ydis = fabs(*(s->carbonlocations+2*C1+1)-*(s->carbonlocations+2*C2+1));
		double ysum = *(s->carbonlocations+2*C1+1)+*(s->carbonlocations+2*C2+1);
		if(ydis > 1.5)
			if(ysum > s->boundary[1])
				*(Oarray+1) = (ysum - s->boundary[1])/2;
			else
				*(Oarray+1) = (ysum + s->boundary[1])/2;
		else
		{
			*(Oarray+1) = ysum/2;
		}
		
	}
	else
	{
		*(Oarray) = (*(s->carbonlocations+2*C1)+*(s->carbonlocations+2*C2))/2;
		*(Oarray+1) = (*(s->carbonlocations+2*C1+1)+*(s->carbonlocations+2*C2+1))/2;
	}
}

void addcarbon(struct sim *s,int carbonindex)
{
	if(debug)
	{	
		if(verbosedebug)
		printf("adding carbon %d to adjacency lists\n",carbonindex);
		if(carbonindex >= s->numcarbons)
		printf("Carbon %d was added to adjacency even though its greater than number of carbons %d\n",carbonindex,s->numcarbons);
		//if(!*(s->carboncovered+carbonindex))
		//printf("Carbon %d was not covered prior to this event\n",carbonindex);
	}
	
	int searchindex;
	*(s->carboncovered+carbonindex) = 0;
	double distance;
	for(searchindex = 0;searchindex<s->numcarbons;searchindex++)
		{
			if(!*(s->carboncovered+searchindex) && searchindex != carbonindex && !intcontains(s->adjacencylist+maxconnections*carbonindex,searchindex,*(s->numadjacent+carbonindex)))
			{
				distance = distancebetweenpoints(s->carbonlocations+2*carbonindex, s->carbonlocations+2*searchindex,2,s->periodic,s->boundary);
				if(distance < 1.01 )//Add each to others adjacency list
				{
				*(s->adjacencylist+maxconnections*carbonindex + *(s->numadjacent+carbonindex)) = searchindex;
				(*(s->numadjacent+carbonindex))++;
				*(s->adjacencylist+maxconnections*searchindex + *(s->numadjacent+searchindex)) = carbonindex;
				(*(s->numadjacent+searchindex))++;
				}
			}
		}
}

void reducecarbon(struct sim *s,int carbonindex)
{
	if(debug)
	{	
		if(verbosedebug)
		printf("adding carbon %d to adjacency lists\n",carbonindex);
		if(carbonindex >= s->numcarbons)
		printf("Carbon %d was added to adjacency even though its greater than number of carbons %d\n",carbonindex,s->numcarbons);
		if(!*(s->carboncovered+carbonindex))
		printf("Carbon %d was not covered prior to this event\n",carbonindex);
	}
	
	int searchindex;
	*(s->carboncovered+carbonindex) = 0;
	int index,index2,index3,adjacencyindex;
	int swap;
	
	if(*(s->numadjacent+carbonindex) > 0) 
	for(adjacencyindex = 0;adjacencyindex< *(s->numadjacent+carbonindex);adjacencyindex++)
	{
		index = *(s->adjacencylist+maxconnections*carbonindex+adjacencyindex);//pull removal index out of this ones list
		if(verbosedebug)
		printf("connecting carbon %d to carbon %d\n",carbonindex,index);
		for(index2 = 0;index2<*(s->numadjacent+index);index2++)
		{
			if(*(s->adjacencylist+maxconnections*index+index2) == carbonindex)
			{
			for(index3 = index2;index3>=*(s->numadjacent+index)-*(s->numoxidizedadjacent+index)+1; index3--)
			{
				swap = *(s->adjacencylist+maxconnections*index+index3);
				*(s->adjacencylist+maxconnections*index+index3) = *(s->adjacencylist+maxconnections*index+index3-1);
				*(s->adjacencylist+maxconnections*index+index3-1) = swap;
			}	
			(*(s->numoxidizedadjacent+index))--;

			break;
			}	
		}
	}
	
}
/*
void removecarbon(struct sim *s,int carbonindex)
{
	if(verbosedebug)
	printf("removing carbon %d from adjacency lists\n",carbonindex);
	if(debug)
	if(*(s->carboncovered+carbonindex) == 1)
	printf("Trying to remove carbon %d from adjacency lists, but it is already covered\n",carbonindex);
	*(s->carboncovered+carbonindex) = 1;
	int index,index2,index3,adjacencyindex;
	if(*(s->numadjacent+carbonindex) > 0) 
	for(adjacencyindex = 0;adjacencyindex< *(s->numadjacent+carbonindex);adjacencyindex++)
	{
		index = *(s->adjacencylist+maxconnections*carbonindex+adjacencyindex);//pull removal index out of this ones list
		if(verbosedebug)
		printf("disconnecting carbon %d from carbon %d\n",carbonindex,index);
		for(index2 = 0;index2<*(s->numadjacent+index);index2++)
		{
			if(*(s->adjacencylist+maxconnections*index+index2) == carbonindex)
			{
				if(index2 <*(s->numadjacent+index)-1)
			for(index3 = index2;index3<*(s->numadjacent+index)-1;index3++)
			*(s->adjacencylist+maxconnections*index+index3) = *(s->adjacencylist+maxconnections*index+index3+1);
			(*(s->numadjacent+index))--;
			break;
			}	
		}
	}
	*(s->numadjacent+carbonindex) = 0;
	
}*/

void oxidizecarbon(struct sim *s,int carbonindex)
{
	if(verbosedebug)
	printf("removing carbon %d from adjacency lists\n",carbonindex);
	if(debug)
	if(*(s->carboncovered+carbonindex) == 1)
	printf("Trying to remove carbon %d from adjacency lists, but it is already covered\n",carbonindex);
	*(s->carboncovered+carbonindex) = 1;
	int index,index2,index3,adjacencyindex;
	int swap;
	if(*(s->numadjacent+carbonindex) > 0) 
	for(adjacencyindex = 0;adjacencyindex< *(s->numadjacent+carbonindex);adjacencyindex++)
	{
		index = *(s->adjacencylist+maxconnections*carbonindex+adjacencyindex);//pull removal index out of this ones list
		if(verbosedebug)
		printf("disconnecting carbon %d from carbon %d\n",carbonindex,index);
		for(index2 = 0;index2<*(s->numadjacent+index);index2++)
		{
			if(*(s->adjacencylist+maxconnections*index+index2) == carbonindex)
			{
			for(index3 = index2;index3<*(s->numadjacent+index)-*(s->numoxidizedadjacent+index)-1; index3++)
			{
				swap = *(s->adjacencylist+maxconnections*index+index3);
				*(s->adjacencylist+maxconnections*index+index3) = *(s->adjacencylist+maxconnections*index+index3+1);
				*(s->adjacencylist+maxconnections*index+index3+1) = swap;
			}	
			(*(s->numoxidizedadjacent+index))++;

			break;
			}	
		}
	}
	
}

void removeadjacency(struct sim *s,int carbonindex)
{
	if(verbosedebug)
	printf("removing carbon %d from adjacency lists\n",carbonindex);
	//if(debug)
	//if(*(s->carboncovered+carbonindex) == 1)
	//printf("Trying to remove carbon %d from adjacency lists, but it is already covered\n",carbonindex);
	//*(s->carboncovered+carbonindex) = 1;
	int covered = *(s->carboncovered+carbonindex);
	int index,index2,index3,adjacencyindex;
	if(*(s->numadjacent+carbonindex) > 0) 
	for(adjacencyindex = 0;adjacencyindex< *(s->numadjacent+carbonindex);adjacencyindex++)
	{
		index = *(s->adjacencylist+maxconnections*carbonindex+adjacencyindex);//pull removal index out of this ones list
		if(verbosedebug)
		printf("disconnecting carbon %d from carbon %d\n",carbonindex,index);
		for(index2 = 0;index2<*(s->numadjacent+index);index2++)
		{
			if(*(s->adjacencylist+maxconnections*index+index2) == carbonindex)
			{
				if(index2 <*(s->numadjacent+index)-1)
			for(index3 = index2;index3<*(s->numadjacent+index)-1;index3++)
			*(s->adjacencylist+maxconnections*index+index3) = *(s->adjacencylist+maxconnections*index+index3+1);
			(*(s->numadjacent+index))--;
			if(covered)
			(*(s->numoxidizedadjacent+index))--;
			break;
			}	
		}
	}
	*(s->numadjacent+carbonindex) = 0;
	
}

void Chopcheck(struct sim *s,int carbonindex)
{
	int oxygenindex;
	double distance,distance2;
	struct hop *h;
	int carbonindex2;
	for(oxygenindex = 0;oxygenindex<s->numoxygens;oxygenindex++)
	{
		distance = distancebetweenpoints(s->carbonlocations+2*carbonindex,s->oxygenlocations+2*oxygenindex,2,s->periodic,s->boundary);
		if(distance<oxygendistance)
		{
			h = (s->hops+s->numhops);
			h->Cunchanged = -1;
			h->Cuncovered = -1;
			h->oxygen = oxygenindex;
			h->Ccovered = carbonindex;
			h->Olocold[0] = *(s->oxygenlocations+2*oxygenindex);
			h->Olocold[1] = *(s->oxygenlocations+2*oxygenindex+1);
			for(carbonindex2 = 0;carbonindex2<s->numcarbons;carbonindex2++)
			if(carbonindex2 != carbonindex)
			{
				distance = distancebetweenpoints(s->carbonlocations+2*carbonindex,s->carbonlocations+2*carbonindex2,2,s->periodic,s->boundary);
				distance2 = distancebetweenpoints(s->carbonlocations+2*carbonindex2,s->oxygenlocations+2*oxygenindex,2,s->periodic,s->boundary);
				if(distance <1.01 && distance2 < 0.51)
				h->Cunchanged = carbonindex2;
				else
				if(distance < 2.01  && distance2 < 0.51)
				h->Cuncovered = carbonindex2;
			}
			
			placeoxygen(s,h->Cunchanged,h->Ccovered,h->Olocnew);
			//h->Olocnew[0] = (*(carbonlocations + 2* h->Cunchanged) + *(carbonlocations + 2* h->Ccovered))/2;
			//h->Olocnew[1] = (*(carbonlocations + 2* h->Cunchanged+1) + *(carbonlocations + 2* h->Ccovered+1))/2;
//			if(h->Cuncovered != -1 && h->Cunchanged != -1)//This is a really crude "fix" 
			s->numhops++;
			if(verbosedebug)
			printhop(s,h);
		}
	}
	
}

void Ohopcheck(struct sim *s,int oxygenindex, int carbonexclusion)
{
	int carbonindex;
	double distance,distance2;
	struct hop *h;
	int carbonindex2;
	for(carbonindex = 0;carbonindex<s->numcarbons;carbonindex++)
	{
		distance = distancebetweenpoints(s->carbonlocations+2*carbonindex,s->oxygenlocations+2*oxygenindex,2,s->periodic,s->boundary);
		if(distance<oxygendistance && !*(s->carboncovered+carbonindex) && carbonindex != carbonexclusion)
		{
			h = (s->hops+s->numhops);
			h->Cunchanged = -1;
			h->Cuncovered = -1;
			h->oxygen = oxygenindex;
			h->Ccovered = carbonindex;
			h->Olocold[0] = *(s->oxygenlocations+2*oxygenindex);
			h->Olocold[1] = *(s->oxygenlocations+2*oxygenindex+1);
			for(carbonindex2 = 0;carbonindex2<s->numcarbons;carbonindex2++)
			if(carbonindex2 != carbonindex)
			{
				distance = distancebetweenpoints(s->carbonlocations+2*carbonindex,s->carbonlocations+2*carbonindex2,2,s->periodic,s->boundary);
				distance2 = distancebetweenpoints(s->carbonlocations+2*carbonindex2, s->oxygenlocations+2*oxygenindex,2,s->periodic,s->boundary);
				if(distance <1.01 && distance2 < 0.51)
				h->Cunchanged = carbonindex2;
				else
				if(distance < 2.01  && distance2 < 0.51)
				h->Cuncovered = carbonindex2;
			}
			placeoxygen(s,h->Cunchanged,h->Ccovered,h->Olocnew);

			//h->Olocnew[0] = (*(carbonlocations + 2* h->Cunchanged) + *(carbonlocations + 2* h->Ccovered))/2;
			//h->Olocnew[1] = (*(carbonlocations + 2* h->Cunchanged+1) + *(carbonlocations + 2* h->Ccovered+1))/2;
//			if(h->Cuncovered != -1 && h->Cunchanged != -1)//This is a really crude "fix" 
			s->numhops++;
			if(verbosedebug)
			printhop(s,h);
		}
	}
}

void removehop(struct sim *s,int removalindex)
{
	int swapindex;
	if(removalindex <s->numhops-1)
	for(swapindex = removalindex;swapindex<s->numhops-1;swapindex++)
	*(s->hops+swapindex) = *(s->hops+swapindex+1);
	s->numhops--;
}

void performhop(struct sim *s,int hopindex)
{
	struct hop *h =malloc(sizeof(struct hop));
	*h = *(s->hops+hopindex);//deep copy
	int index;
	
	*(s->oxygenlocations + 2 * (h->oxygen)) = h->Olocnew[0];
	*(s->oxygenlocations + 2 * (h->oxygen)+1) = h->Olocnew[1];
	
	if(debug)
	{
	double distance;
		distance = distancebetweenpoints(s->oxygenlocations+2*(h->oxygen),h->Olocnew,2,s->periodic,s->boundary);		
	if(distance > 1.5)
	{
	printf("Hop Distance %g too large, hOold %g x %g, hOnew %g x %g, Oold %g x %g\n",distance,h->Olocold[0], h->Olocold[1], h->Olocnew[0], h->Olocnew[1], *(s->oxygenlocations+2*h->oxygen),*(s->oxygenlocations+2*h->oxygen+1)); 
	printhop(s,h);
	}
	
	if(h->Ccovered >= s->numcarbons)
	printf("Ccovered is too high\n");
	
	if(h->Cuncovered >= s->numcarbons)
	printf("Cuncovered is too high\n");
	
	if(h->Cunchanged >= s->numcarbons)
	printf("Cunchanged is too high\n");
	
	if(*(s->carboncovered+h->Ccovered))
	printf("Carbon to cover is already covered\n");
	
	if(!*(s->carboncovered+h->Cuncovered))
	printf("Carbon to uncovere is already uncovered\n");
	
	if(!*(s->carboncovered+h->Cunchanged))
	printf("Carbon unchanged is uncovered\n");
	}
	//removeadjacency(s,h->Ccovered);
	oxidizecarbon(s,h->Ccovered);
	reducecarbon(s,h->Cuncovered);
	//addadjacency(s,h->Cuncovered);
	for(index = 0;index<s->numhops;index++)
	if((s->hops+index)->oxygen == h->oxygen && index != hopindex)
	{
	if(verbosedebug)
	printf("removing hop because oxygen moved\n");
	removehop(s,index);
	index--;
	hopindex--;
	}
	
	Chopcheck(s,h->Cuncovered);
	Ohopcheck(s,h->oxygen,h->Cuncovered);
	
	for(index = 0;index<s->numhops;index++)
	if(*(s->carboncovered+(s->hops+index)->Ccovered))
	{
		if(verbosedebug)
		printf("removing hop because carbon covered\n");//This will remove the performed hop as well
	removehop(s,index);
	index--;
	}
	free(h);
	h = NULL;
}
void countsourcesites(struct sim *s)
{
	int index;
	int sum=0;
	
	for(index=0;index<s->numcarbons;index++)
	{
		if(!*(s->carboncovered+index))
		sum += *(s->numadjacent+index)-*(s->numoxidizedadjacent+index);
	}
	//sum /= 2;//sites get double counted but its okay for placement
	s->numsourcesites = sum;
	//return sum;
	if(verbosedebug)
	printf("Number of source sites %d\n",sum);
}

void oxygensource(struct sim *s,int sourceindex)
{
	int index;
	int sum=0;
	
	int selectedneighbor,selectedcarbon;
	for(index=0;index<s->numcarbons;index++)
	{
		if(!*(s->carboncovered+index))
		sum += *(s->numadjacent+index)-*(s->numoxidizedadjacent+index);
		if(sum >= sourceindex) 
		{
			
			selectedcarbon = index;
			selectedneighbor = *(s->adjacencylist+maxconnections*selectedcarbon + sum-sourceindex);
			//if(*(s->carboncovered+selectedcarbon))
			//printf("Selected Carbon already covered\n");
			//if(*(s->carboncovered+selectedneighbor))
			//printf("Selected neighbor already covered\n");
			placeoxygen(s,selectedcarbon,selectedneighbor,s->oxygenlocations+2*s->numoxygens);
			s->numoxygens++;
			oxidizecarbon(s,selectedneighbor);
			oxidizecarbon(s,selectedcarbon);
			//removeadjacency(s,selectedneighbor);
			//removeadjacency(s,selectedcarbon);
			s->numcovered += 2;
			
			break;
		}
	}
	//sum /= 2;//sites get double counted but its okay for placement
	//s->numsourcesites = sum;
	for(index = 0;index<s->numhops;index++)
	if((s->hops+index)->Ccovered == selectedcarbon || (s->hops+index)->Ccovered == selectedneighbor)
	{
		if(verbosedebug)
		printf("removing hop because carbon covered\n");//This will remove the performed hop as well
	removehop(s,index);
	index--;
	}
	Ohopcheck(s,s->numoxygens-1,-1);
	
	//return sum;
}

void clusterdetection(struct sim *s)
{
	s->numclusters = 0;
	double distance;
	int index1,index2;
		
	int clusterindex = 0, carbonclusterindex;	
	int usedcarbons[s->numcarbons];
	int numusedcarbons = 0;
	int adjacencyindex;
	int carbonindex,searchindex;
	for(index1=0;index1<s->numcarbons;index1++)
		if(!*(s->carboncovered+index1) && !intcontains(usedcarbons,index1,numusedcarbons))//must have at least 1 nieghbor and not have been used before
		{
			
			*((s->clusters + clusterindex)->carbons) = index1;
			(s->clusters + clusterindex)->numcarbons = 1;
			*(usedcarbons+numusedcarbons) = index1;
			numusedcarbons++; 
			s->numclusters++;	
			for(carbonclusterindex=0;carbonclusterindex<(s->clusters+clusterindex)->numcarbons;carbonclusterindex++)
			{
				carbonindex = *(((s->clusters+clusterindex)->carbons)+carbonclusterindex);
				if(debug)
				if(*(s->numadjacent + carbonindex)-*(s->numoxidizedadjacent+carbonindex) < 0 || *(s->numadjacent + carbonindex) > maxconnections)
				printf("error in numadjacent carbon %d has %d adjacent %d oxidized\n",carbonindex,*(s->numadjacent + carbonindex),*(s->numadjacent + carbonindex));
				if(*(s->numadjacent+carbonindex)-*(s->numoxidizedadjacent+carbonindex))
				for(adjacencyindex=0;adjacencyindex< *(s->numadjacent + carbonindex)-*(s->numoxidizedadjacent+carbonindex);adjacencyindex++)
				{
						searchindex = *(s->adjacencylist + maxconnections*carbonindex + adjacencyindex);
						/*The carboncovered check and usedcarbon check should not be necessary*/
						if(/* !*(carboncovered+searchindex) && !contains(usedcarbons,searchindex,numusedcarbons) && */!intcontains((s->clusters+clusterindex)->carbons,searchindex,(s->clusters+clusterindex)->numcarbons))
						{
						*(((s->clusters + clusterindex)->carbons)+((s->clusters + clusterindex)->numcarbons)) = searchindex;
						((s->clusters + clusterindex)->numcarbons)++;
						*(usedcarbons+numusedcarbons) = searchindex;
						numusedcarbons++;
						}
				}
			}
			clusterindex++;
		}
				double energy;
	
	*(s->emission+s->stepnumber) = 0;
	*(s->absorption+s->stepnumber) = 0;
	*(s->p1+s->stepnumber) = 0;
	*(s->p2+s->stepnumber) = 0;
	*(s->p3+s->stepnumber) = 0;
	*(s->emittingclusters+s->stepnumber) = 0;
	int num;
	for(index1 = 0;index1<s->numcarbons;index1++)
	{
		if(*(s->numadjacent+index1)-*(s->numoxidizedadjacent+index1) == 1)
		(*(s->p1+s->stepnumber))+= (double)1/s->initialcarbons;
		if(*(s->numadjacent+index1)-*(s->numoxidizedadjacent+index1) == 2)
		(*(s->p2+s->stepnumber))+= (double)1/s->initialcarbons;
		if(*(s->numadjacent+index1)-*(s->numoxidizedadjacent+index1) == 3)
		(*(s->p3+s->stepnumber))+= (double)1/s->initialcarbons;
	}
	for(index1 = 0;index1<largestclustercounted;index1++)
	*(s->clusterdistribution+(s->stepnumber)*largestclustercounted+index1) = 0;
	
		for(clusterindex=0;clusterindex<s->numclusters;clusterindex++)
		{
			//
			if((s->clusters+clusterindex)->numcarbons < largestclustercounted && (s->clusters+clusterindex)->numcarbons >= clusterminsize )
			{
			(*(s->clusterdistribution+(s->stepnumber)*largestclustercounted+(s->clusters+clusterindex)->numcarbons))++;
			//printf("Stepnumber %d clustersize %d\n",s->stepnumber,(s->clusters+clusterindex)->numcarbons);
			}
			energy = 12.9 * sqrt(M_PI/(s->clusters+clusterindex)->numcarbons);
			if(energy > 1.38 && energy <3.1)
			//if(energy <3.1)
			{
			(s->clusters+clusterindex)->color = (int)((energy - 1.38) * 360/(3.1-1.38) + 320) % 360;
			//*(s->emission+s->stepnumber) += (double)(s->clusters+clusterindex)->numcarbons/s->initialcarbons;
			*(s->emission+s->stepnumber) += emissionofcluster((s->clusters+clusterindex)->numcarbons);
			(*(s->emittingclusters+s->stepnumber))++;
			}
			else
			 (s->clusters+clusterindex)->color = -1;
			
			if(energy < 3.1)
			*(s->absorption+s->stepnumber) += 2.3 *  (s->clusters+clusterindex)->numcarbons/s->initialcarbons;
			// printf("Cluster of size %d energy %g(eV) color %d\n",(clusters+clusterindex)->numcarbons,energy,(clusters+clusterindex)->color);
		}
		
		*(s->QY+s->stepnumber) = 274.5 * *(s->emission+s->stepnumber) * 2.3 / *(s->absorption+s->stepnumber)/s->initialcarbons;		
}

void removehydroxide(struct sim *s,int removalindex)
{
	int index;
	double distance;
	for(index = 0;index<s->numcarbons;index++)
	{
		distance = distancebetweenpoints(s->carbonlocations+2*index,s->hydroxidelocations+2*removalindex,2,s->periodic,s->boundary);
		if(distance < 0.5)
		{
			//addadjacency(s,index);
			reducecarbon(s,index);
			Chopcheck(s,index);
			break;
		}
	}
	
	for(index = removalindex;index<s->numhydroxides-1;index++)
	{
		*(s->hydroxidelocations +2*index) = *(s->hydroxidelocations+2*(index+1));
		*(s->hydroxidelocations + 2*index+1) = *(s->hydroxidelocations +2*(index+1)+1);
	}
	s->numhydroxides--;
	s->numcovered--;
}

void removeCO(struct sim *s,int removalindex)
{
	if(verbosedebug)
	printf("Removing CO %d\n",removalindex);
	int index,index2,uncoverindex = -1;
	int Ctoremove = 2*generaterandom();//whether to remove first or second C we come accross
	int skipped = 0;
	double distance;
	for(index = 0;index<s->numcarbons;index++)
	{
		distance = distancebetweenpoints(s->carbonlocations+2*index,s->oxygenlocations+2*removalindex,2,s->periodic,s->boundary);
		if(distance < 1)
		{
			if(Ctoremove ==0 && skipped == 0 || skipped == 1 && Ctoremove == 1)
			{
			removeC(s,index);
			index--;
			//if(uncoverindex > index)
			//uncoverindex--;
			}
			else
			{
			uncoverindex = index;
			}
			skipped++;
			
		}
	}
//	if(debug)
	if(uncoverindex == -1)
	printf("Never found the carbon to uncover\n");
	
	for(index = removalindex;index<s->numoxygens-1;index++)
	{
		*(s->oxygenlocations +2*index) = *(s->oxygenlocations+2*(index+1));
		*(s->oxygenlocations + 2*index+1) = *(s->oxygenlocations +2*(index+1)+1);
	}
	for(index = 0;index<s->numhops;index++)
	{
		if((s->hops+index)->oxygen == removalindex)
		{
		removehop(s,index);
		index--;
		if(verbosedebug)
		printf("Removed hop because oxygen gone\n");	
		}
		else
		if((s->hops+index)->oxygen > removalindex)
		{
		((s->hops+index)->oxygen)--;
		if(verbosedebug)
		printf("Reduced hop number because oxygen removed\n");
		}
	}
	s->numoxygens--;
	s->numcovered-=2;
//	if(uncoverindex != -1)//Bad correction
//	{
	//addadjacency(s,uncoverindex);
	reducecarbon(s,uncoverindex);
	Chopcheck(s,uncoverindex);
//}
}

void removeO(struct sim *s,int removalindex)
{
	if(verbosedebug)
	printf("Removing O %d\n",removalindex);
	int index,index2,uncoverindex1 = -1,uncoverindex2 = -1;
	//int Ctoremove = 2*generaterandom();//whether to remove first or second C we come accross
	int found = 0;
	double distance;
	for(index = 0;index<s->numcarbons;index++)
	{
		distance = distancebetweenpoints(s->carbonlocations+2*index,s->oxygenlocations+2*removalindex,2,s->periodic,s->boundary);
		if(distance < 1)
		{
			if(found)
			{
			uncoverindex2 = index;
			break;
			}
			else
			uncoverindex1 = index;
			
			found++;
		}
	}
	if(uncoverindex1 == -1 || uncoverindex2 == -1)
	printf("Never found the carbon to uncover\n");
	
	for(index = removalindex;index<s->numoxygens-1;index++)
	{
		*(s->oxygenlocations +2*index) = *(s->oxygenlocations+2*(index+1));
		*(s->oxygenlocations + 2*index+1) = *(s->oxygenlocations +2*(index+1)+1);
	}
	for(index = 0;index<s->numhops;index++)
	{
		if((s->hops+index)->oxygen == removalindex)
		{
		removehop(s,index);
		index--;
		if(verbosedebug)
		printf("Removed hop because oxygen gone\n");	
		}
		else
		if((s->hops+index)->oxygen > removalindex)
		{
		((s->hops+index)->oxygen)--;
		if(verbosedebug)
		printf("Reduced hop number because oxygen removed\n");
		}
	}
	s->numoxygens--;
	s->numcovered-=2;

	reducecarbon(s,uncoverindex1);
	reducecarbon(s,uncoverindex2);
	Chopcheck(s,uncoverindex1);
	Chopcheck(s,uncoverindex2);

}

	
void removeC(struct sim *s,int removalindex)
{
	int index,adjacencyindex,index2,index3;
		removeadjacency(s,removalindex);
	
	for(index = 0;index<s->numcarbons;index++)
	{
		for(adjacencyindex = 0;adjacencyindex <maxconnections;adjacencyindex++)
		{
			if(*(s->adjacencylist+maxconnections*index+adjacencyindex) > removalindex)
			(*(s->adjacencylist+maxconnections*index+adjacencyindex))--;
		}
	}
//	printf("here I am\n");
	for(index = removalindex;index<s->numcarbons-1;index++)
	{
		*(s->carbonlocations +2*index) = *(s->carbonlocations+2*(index+1));
		*(s->carbonlocations + 2*index+1) = *(s->carbonlocations +2*(index+1)+1);
		*(s->carboncovered+index) = *(s->carboncovered+index+1);
		*(s->numadjacent+index) = *(s->numadjacent+index+1);
		*(s->numoxidizedadjacent+index) = *(s->numoxidizedadjacent+index+1);
		
		for(adjacencyindex = 0;adjacencyindex <maxconnections;adjacencyindex++)
		{
			*(s->adjacencylist+maxconnections*index+adjacencyindex) = *(s->adjacencylist+maxconnections*(index+1)+adjacencyindex);
		}
	}
//	printf("listen\n");
	for(index = 0;index<s->numhops;index++)
	{
		if((s->hops+index)->Ccovered == removalindex)
		{
		removehop(s,index);
		index--;
	}
	else
	{
		if((s->hops+index)->Ccovered > removalindex)
		((s->hops+index)->Ccovered)--;
		if((s->hops+index)->Cuncovered > removalindex)
		((s->hops+index)->Cuncovered)--;
		if((s->hops+index)->Cunchanged > removalindex)
		((s->hops+index)->Cunchanged)--;
	}
	}
	s->numcarbons--;
}

	
gboolean removes(struct sim *s)
{
	
	int index;
	double accruedtime = 0;
	for(index=0;index<numremoves;index++)
	{
	removal(s);
	if(s->stepnumber != 0)
	accruedtime += *(s->times+s->stepnumber) - *(s->times+s->stepnumber-1);
	}
	if(s->stepnumber != 0)
	*(s->times+s->stepnumber) = *(s->times+s->stepnumber-1) + accruedtime;
	
	s->stepnumber++;
	return FALSE;
}
	
	
gboolean removal(struct sim *s)
{
	//debug = gtk_toggle_button_get_active(Debugcheckbutton);
	//verbosedebug = gtk_toggle_button_get_active(VerboseDebugcheckbutton);
	countsourcesites(s);
	int index;
	if(debug)
	{
	checknumcovered(s);
	hopcheck(s);
	}
	double OH,CO,C1, Czig, Carm, C3,Ohop,Osource;
	
		
	OH = s->numhydroxides*hydroxideprobability;
	CO = s->numoxygens*COprobability;
	Ohop = Ohopprobability* s->numhops;
	Osource = Osourceprobability * s->numsourcesites;
	
	C1=Czig=Carm=C3=0;
	int C1C[s->numcarbons],CzigC[s->numcarbons],CarmC[s->numcarbons],C3C[s->numcarbons];
	int numC1=0,numCzig=0,numCarm=0,numC3=0;
	for(index =0;index<s->numcarbons;index++)
	{
		if(!*(s->carboncovered+index))
		switch(*(s->numadjacent+index))
		{
			case 0:
			C1+=C1probability;
			C1C[numC1] = index;
			numC1++;
			break;
			
			case 1:
			C1+=C1probability;
			C1C[numC1] = index;
			numC1++;
			break;
			
			case 3:
			C3+=C3probability;
			C3C[numC3] = index;
			numC3++;
			break;
			
			case 2:
			{
				/*printf("neighbor 1 is %d and has %d neighbors, neighbor 2 is %d and has %d neighbors\n",
				*(s->adjacencylist + maxconnections*index+0),*(s->numadjacent+ *(s->adjacencylist + maxconnections*index+0)),
				*(s->adjacencylist + maxconnections*index+1),*(s->numadjacent+ *(s->adjacencylist + maxconnections*index+1)));
				*/if(*(s->numadjacent+ *(s->adjacencylist + maxconnections*index+0)) <= 2 ||
				 *(s->numadjacent+ *(s->adjacencylist + maxconnections*index+1)) <= 2)
				 //Armchair
				 {
				 Carm+=Carmprobability;
				 CarmC[numCarm] = index;
				 numCarm++;
				}
				 else
				{
				 Czig+=Czigprobability;
				 CzigC[numCzig] = index;
				 numCzig++;
				}
			}
			
			break;
		}
			
		
	}
	/*
	*(s->rates+numrates*s->stepnumber) = OH;
	*(s->rates+numrates*s->stepnumber+1) = CO;
	*(s->rates+numrates*s->stepnumber+2) = C1;
	*(s->rates+numrates*s->stepnumber+3) = C3;
	*(s->rates+numrates*s->stepnumber+4) = Czig;
	*(s->rates+numrates*s->stepnumber+5) = Carm;	
	*(s->rates+numrates*s->stepnumber+6) = Ohop;	
	*(s->rates+numrates*s->stepnumber+7) = Osource;	
	*/
	*(s->rates+numrates*s->stepnumber) = s->numhydroxides;
	*(s->rates+numrates*s->stepnumber+1) = s->numoxygens;
	*(s->rates+numrates*s->stepnumber+2) = s->numcarbons;
	*(s->rates+numrates*s->stepnumber+3) = s->numcarbons;
	*(s->rates+numrates*s->stepnumber+4) = s->numcarbons;
	*(s->rates+numrates*s->stepnumber+5) = s->numcarbons;	
	*(s->rates+numrates*s->stepnumber+6) = Ohop;	
	*(s->rates+numrates*s->stepnumber+7) = Osource;	
	
	
	
	
	double sum = OH+CO+C1 + C3 + Czig + Carm+Ohop + Osource;
	OH /= sum;
	CO = CO/sum+OH;
	C1 = C1/sum+CO;
	C3 = C3/sum+C1;
	Czig = Czig/sum+C3;
	Carm = Carm/sum+Czig;
	Ohop = Ohop/sum+Carm;
	Osource = Osource/sum+Ohop;
	if(verbosedebug)
	{
	printf("OH %g CO %g C1 %g C3 %g Czig %g Carm %g Hop %g SOurce %g\n",OH,CO,C1,C3,Czig,Carm,Ohop,Osource);
	printf("num C1 %d C3 %d Czig %d Carm %d O %d OH %d\n",numC1,numC3,numCzig,numCarm,s->numoxygens,s->numhydroxides);
	}
	
	
	
	//printf("A\n");
	if(s->stepnumber != 0)
	*(s->times+s->stepnumber) = *(s->times+s->stepnumber-1) + 1/sum;
	else
	*(s->times) = 0;
	//printf("B\n");
	
	
	double rand = generaterandom(); 
	int removalindex;
	if(rand <= OH)
	{
		removalindex = generaterandom() * s->numhydroxides;
		if(verbosedebug)
		printf("removing hydroxide %d of %d\n",removalindex,s->numhydroxides);
		removehydroxide(s,removalindex);
	}
	if(rand >OH && rand <= CO)//Changing this to O not CO
	{
		removalindex = generaterandom() * s->numoxygens;
		if(verbosedebug)
		printf("removing CO %d of %d\n",removalindex,s->numoxygens);
		//removeCO(s,removalindex);
		removeO(s,removalindex);
	}
	if(rand > CO && rand <=  C1)
	{
		removalindex = C1C[(int)(generaterandom() * numC1)];
		
		if(verbosedebug)
		printf("removing carbon %d of %d\n",removalindex,s->numcarbons);
		removeC(s,removalindex);
	}
	if(rand > C1 && rand <=  C3)
	{
		removalindex = C3C[(int)(generaterandom() * numC3)];
		
		
		if(verbosedebug)
		printf("removing carbon %d of %d\n",removalindex,s->numcarbons);
		removeC(s,removalindex);
	}
	if(rand > C3 && rand <=  Czig)
	{
		removalindex = CzigC[(int)(generaterandom() * numCzig)];
		
		if(verbosedebug)
		printf("removing carbon %d of %d\n",removalindex,s->numcarbons);
		removeC(s,removalindex);
	}
	if(rand > Czig && rand <=  Carm)
	{
		removalindex = CarmC[(int)(generaterandom() * numCarm)];
		
		if(verbosedebug)
		printf("removing carbon %d of %d\n",removalindex,s->numcarbons);
		removeC(s,removalindex);
	}
	if(rand > Carm && rand <=  Ohop)
	{
		removalindex = generaterandom() * s->numhops;
		if(verbosedebug)
		printhop(s->hops+removalindex);
		performhop(s,removalindex);
	}
	if(rand >  Ohop)
	{
		removalindex = generaterandom() * s->numsourcesites;
		if(verbosedebug)
		printf("Sourcing oxygen %d\n",removalindex);
		oxygensource(s,removalindex);
	}
	if(debug)
	{
	checknumcovered(s);
	hopcheck(s);
	}
	return FALSE;
}

