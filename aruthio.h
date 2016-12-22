

int intcontains(int *list, int search, int numitems)
{
	if(numitems == 0)
	return 0;
	
	int index;
	for(index = 0;index<numitems;index++)
	if(*(list+index) == search)
	return 1;
	
	return 0;
}

int doublecontains(double *list, double search,double epsilon, int numitems)
{
	
	
	if(numitems == 0)
	return 0;
	
	int index;
	for(index = 0;index<numitems;index++)
	if(fabs(*(list+index) - search) <= epsilon)
	return index+1;//Returns index +1 to distinguish 0 and 1
	
	return 0;
}


int containstwice(double *list1, double search1, double *list2, double search2, int numitems)
{
	//Double search function
	double epsilon = 1e-10;
	if(numitems == 0)
	return 0;
	
	int index;
	for(index = 0;index<numitems;index++)
	if(fabs(*(list1+index)- search1)<epsilon && fabs(*(list2+index)- search2)<epsilon)
	return 1;
	
	return 0;
}


int countcolumns(char *filename)
{
	char *garbage1 = malloc(20000*sizeof(char));
	char *garbage2;	
	FILE *infile = fopen(filename,"r");
	
	int counter = 0;
	fscanf(infile,"%[^\n]",garbage1);
	garbage2 = strtok (garbage1," \t");
	while(garbage2 != NULL)
	{
		garbage2 = strtok(NULL," \t");		
		counter++;
	}
	fclose(infile);	
	free(garbage1);
	free(garbage2);
	return counter;
}
int countlines(char *filename)
{
	char *garbage1 = malloc(20000*sizeof(char));
	char *garbage2 = malloc(20000*sizeof(char));
	char *discard = malloc(10*sizeof(char));	
	FILE *infile = fopen(filename,"r");
	printf("Counting lines in %s\n",filename);
	sprintf(garbage1,".*?");
	sprintf(garbage2,"|@#");

	int counter = -1;
//printf("Result of strcmp %d !strcmp %d\n",strcmp(garbage1,garbage2),!strcmp);
	while(strcmp(garbage1,garbage2))
	{
	strcpy(garbage2,garbage1);
	fscanf(infile,"%[^\n]%[\n]",garbage1,discard);
	//printf("got %s %s\n",garbage1,garbage2);
	counter++;
	}
	fclose(infile);	
	free(garbage1);
	free(garbage2);
	free(discard);
	return counter;
	
}

void printmatrix(double *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			printf("%.03g \t",*(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				printf("\n");
			}
			
	}
}

void makematrix(char *output,double *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	char *buffer = malloc(1000*sizeof(char));
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
			sprintf(buffer,"%s",output);
			if ((i+1) % numberOfColumns == 0)
			{ 
				sprintf(output,"%s%.16g",buffer,*(matrix+i));  
				sprintf(buffer,"%s",output);
				sprintf(output,"%s\n",buffer);
			}
			else
			sprintf(output,"%s%.16g \t",buffer,*(matrix+i));
			
	}
}

void printintmatrix(int *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			printf("%d \t",*(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				printf("\n");
			}
			
	}
}

void printintmatrixfile(FILE *outfile,int *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			fprintf(outfile,"%d \t",*(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
}


void printdoublematrixfile(FILE *outfile,double *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			fprintf(outfile,"%g \t",*(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
}

void readdoublematrixfile(FILE *infile,double *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			fscanf(infile,"%lf \t",(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fscanf(infile,"\n");
			}
			
	}
}

void readintmatrixfile(FILE *infile,int *matrix, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	//printf("Reading int matrix from file\n");
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		//	printf("Matrix is %d\n",*(matrix+i));
			fscanf(infile,"%d \t",(matrix+i));
			//printf("Scanned index %d and found %d\n",i,(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fscanf(infile,"\n");
			}
			
	}
}


void printlabeledintmatrixfile(FILE *outfile,int *matrix, int numberOfRows, int numberOfColumns, char *labels, int stride)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		if(i%numberOfColumns == 0)
		fprintf(outfile,"%s\t",labels + (int)(i/numberOfColumns)*stride);
			fprintf(outfile,"%d \t",*(matrix+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
}

void strip(FILE *infile, char *search, char *output,char *delimeter, int fill)
{
	/**************************************************************/
	/* finds an option that fits the mold                         */
	/* search(white space)delimeter(white space)ouput(white space)
	 * if fill = 1, fill in output. Else return after search is found
	 * */
	 
	 //int index;
	 //printf("Stripping\n");
	 char *line = malloc(1000*sizeof(char));
	 char *garbage = malloc(1000*sizeof(char));
	 char buffer[1000];
	 size_t size = 1000;
	 while ( 1 == fscanf( infile, " %999[^\n]%*[\n]", line ) )
	 //while(!feof(infile) )
	 {
		 //getline(&line,&size,infile);
		
		 //fscanf(infile,"%[^\n]%[\n]",line,garbage);
		 //printf("Line %s\n",line);
		 if(strstr(line,search))
		 {
			 if(fill)
			 {
				 sprintf(buffer,"\%[^%s]\%[%s %s]\%[^%s %s]",delimeter,delimeter,"\\t",delimeter,"\\t");
				// printf("scan string %s\n",buffer);
			 sscanf(strstr (line, search),buffer,garbage,garbage,output);
			 //printf("%s\n",output);
			 //sscanf(strstr (line, search),"%[^:]%[: \t]%s",garbage,garbage,output);
			 //printf("%s\n",output);
			}
			 free(line);
			 free(garbage);
			 return;
		 }
	 }
	 if(fill)
	 sprintf(buffer,"");
}

void scanmatrix(FILE *infile, int *matrix, int numrows, int numcolumns)
{
	/****************************************************/
	/* This program reads the contents of a matrix.     */
	/* 								                    */
	/****************************************************/
	int i;
	for (i = 0;i<numrows*numcolumns ;i++ )
	{
		
			fscanf(infile,"%d ",(matrix+i));
			if ((i+1) % numcolumns == 0)
			{
				fscanf(infile,"\n");
			}
		
	
	}
} 

void stripchar(char *string,char removal)
{
	int index;
	char *copy = malloc(strlen(string));
	for(index=0;*(string+index)!= '\0';index++)
	{
		while(*(string+index) == removal)
		{
			strncpy(copy,string+index+1,strlen(string+index+1));
			strncpy(string+index,copy,strlen(string+index+1));
			*(string+strlen(string)-1) ='\0';
		}
	}
}
