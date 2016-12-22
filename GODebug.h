//This area is for debugging functions, These functions add significantly to run time and are not necessary for the operation of the
//program. They are however useful for the detection of bugs.


void printadjacency(struct sim *s)
{
	int index1;
	
	for(index1=0;index1<s->numcarbons;index1++)
	{
		printf("\nLooking at carbon %d its coverage is %d numadjacent %d numoxidizedadjacent %d\n",index1,*(s->carboncovered+index1),*(s->numadjacent+index1),*(s->numoxidizedadjacent+index1));
		printneighbors(s,index1);
	}
	//printf("\n");
}

void printneighbors(struct sim *s,int index1)
{
	int index2;
	for(index2=0;index2<*(s->numadjacent+index1);index2++)
		printf("%d\t",*(s->adjacencylist+maxconnections*index1 +index2));
		printf("\n");
}

void printhops(struct sim *s)
{
	int index;
	struct hop *h;
	for(index = 0;index<s->numhops;index++)
	{
		h = (s->hops+index);
		printhop(s,h);
	}
}

void printhop(struct sim *s,struct hop *h)
{
		printf("numhops %d Hop is C1 C2 C3 O %d %d %d %d %g x %g -> %g x %g\n",s->numhops,h->Ccovered,h->Cunchanged,h->Cuncovered,h->oxygen,h->Olocold[0],h->Olocold[1],h->Olocnew[0],h->Olocnew[1]);
}

void printoxygens(struct sim *s)
{
	int index;
	for(index = 0;index<s->numoxygens;index++)
	printf("Oxygen %d is at %g %g\n",index,*(s->oxygenlocations+2*index),*(s->oxygenlocations+2*index+1));
}

void printcarbons(struct sim *s)
{
	int index;
	for(index = 0;index<s->numcarbons;index++)
	printf("Carbon %d is at %g %g\n",index,*(s->carbonlocations+2*index),*(s->carbonlocations+2*index+1));
}

void checknumcovered(struct sim *s)
{
	int index,count;
	for(index = 0,count = 0;index<s->numcarbons;index++)
	count += *(s->carboncovered+index);
	
	if(s->numcovered != s->numhydroxides+2*s->numoxygens || count != s->numcovered)
	printf("%d carbons, %d covered, %d O, %d OH Num covered %d num counted %d\n",s->numcarbons,s->numcovered,s->numoxygens,s->numhydroxides,s->numcovered,count);
	
}
void checkadjacency(struct sim *s)
{
	//Checks if numadjacent is <= 3 and if covered carbons have any adjacent
	int index;
	for(index=0;index<s->numcarbons;index++)
	{
		if(*(s->carboncovered+index) && *(s->numadjacent+index) > 0)
		{
		printf("Carbon %d is covered but has %d neighbors\n",index,*(s->numadjacent+index));
		printneighbors(s,index);
		}
		if(*(s->numadjacent+index) > 3)
		{
		printf("Carbon %d has %d neighbors\n",index,*(s->numadjacent+index));
		printneighbors(s,index);
		}
	}
}


void hopcheck(struct sim *s)
{
	//printf("Checking hops\n");
	int index;
	struct hop *h;
	double distance1,distance2;
	for(index = 0;index<s->numhops;index++)
	{
		h = s->hops+index;
		if(h->Ccovered >= s->numcarbons || h->Cunchanged >= s->numcarbons || h->Cuncovered >= s->numcarbons)
		{
		printf("Carbon in hop too high numcarbons %d covered %d unchanged %d uncovered %d\n",s->numcarbons,h->Ccovered,h->Cunchanged,h->Cuncovered);
		pausenow();
		}
		if(h->oxygen >= s->numoxygens)
		{
		printf("Oxygen in hop too high numoxygens %d in hop %d",s->numoxygens,h->oxygen);
		pausenow();
	}
		distance1 = distancebetweenpoints(s->oxygenlocations+2*(h->oxygen),s->carbonlocations+2*(h->Cuncovered),2,s->periodic,s->boundary);
		distance2 = distancebetweenpoints(s->oxygenlocations+2*(h->oxygen),s->carbonlocations+2*(h->Cunchanged),2,s->periodic,s->boundary);
		if(distance1 > 1)
		{
		printf("Oxygen and uncovered carbons too far\n");
		printhop(s,h);
		pausenow();
		}
		if(distance2 > 1)
		{
		printf("Oxygen and unchanged carbons too far\n");
		printhop(s,h);
		pausenow();
		}
	}
	//printf("Done checking hops\n");
}
