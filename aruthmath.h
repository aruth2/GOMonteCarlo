int truthmathverbose = 0;
#include <stdlib.h>

//#define pi 3.14159265359

/*    From: "A SURVEY OF COMPUTATIONAL PHYSICS" 
   by RH Landau, MJ Paez, and CC BORDEIANU 
   Copyright Princeton University Press, Princeton, 2008.
   Electronic Materials copyright: R Landau, Oregon State Univ, 2008;
   MJ Paez, Univ Antioquia, 2008; & CC Bordeianu, Univ Bucharest, 2008
   Support by National Science Foundation
*/
// gauss.c: Points and weights for Gaussian quadrature  
// comment: this file has to reside in the same directory as integ.c  
void gauss(int npts, int job, double a, double b, 
                                           double x[], double w[]) {
	
//     npts     number of points                                 
//     job= 0  rescaling uniformly between (a,b)                  
//           1  for integral (0,b) with 50% pts inside (0, ab/(a+b)
//           2  for integral (a,inf) with 50% inside (a,b+2a)        
//     x, w     output grid points and weights.                     

  int     m, i, j; 
  double  t, t1, pp, p1, p2, p3;
	double eps= 1e-15;      // eps= accuracy to adjust //changed accuracy from 1e-10 to 1e-15
  m= (npts+1)/2;
      for (i = 1; i<= m; i++) {  
        t = cos(M_PI*(i-0.25)/(npts+0.5));
        t1= 1;
        while((fabs(t-t1))>= eps) { 
          p1= 1.0;
          p2= 0.0;
          for (j = 1; j<= npts; j++) {
            p3= p2;
            p2= p1;
            p1= ((2*j-1)*t*p2-(j-1)*p3)/j;
         }
         pp= npts*(t*p1-p2)/(t*t-1);
         t1= t;
         t = t1 - p1/pp;
       }   
       x[i-1]= -t;
       x[npts-i]= t;
       w[i-1]   = 2.0/((1-t*t)*pp*pp);
       w[npts-i]= w[i-1];
    }
    if (job == 0)   {
      for (i = 0; i<npts ; i++)  {
            x[i]= x[i]*(b-a)/2.0+(b+a)/2.0;
            w[i]= w[i]*(b-a)/2.0;
      }
    }
      if (job == 1)  {
        for (i = 0; i<npts; i++)  {
          x[i]= a*b*(1+x[i]) / (b+a-(b-a)*x[i]);
          w[i]= w[i]*2*a*b*b /((b+a-(b-a)*x[i])*(b+a-(b-a)*x[i]));
       }
    }
      if (job == 2)   {
        for (i = 0; i<npts; i++)  {
          x[i]= (b*x[i]+b+a+a) / (1-x[i]);
          w[i]=  w[i]*2*(a+b)  /((1-x[i])*(1-x[i]));
         }
      }
}
double generaterandom()
{	
/******************************************/
/* Generates a random double from 0       */
/* to 1. The number of different numbers  */
/* allowed is equal to RAND_MAX           */
/******************************************/
	double random;
	random = ((double)rand())/((double)RAND_MAX);	
	return(random);
}

double mindistance(double *x, double *y, int numpoints)
{
	//Given a list of points
	//Returns the smallest distance between 2 points
	
	int index1,index2;
	double distance,min;
	if(numpoints < 2)
	return 0;
	min = pow(pow(*x - *(x+1),2) + pow(*y-*(y+1),2),0.5);
	for(index1=0;index1<numpoints-1;index1++)
	for(index2=index1+1;index2<numpoints;index2++)
	{
		distance = pow(pow(*(x+index1) - *(x+index2),2) + pow(*(y+index1)-*(y+index2),2),0.5);
		if(distance < min)
		min = distance;
	}
	return min;
	
}

double dotproduct(double *vec1,double *vec2,int numpoints)
{
	//numpoints should be length of arrays vec1 and vec2
	int index;
	double sum;
	for(index=0,sum=0;index<numpoints;index++)
	{
		sum += *(vec1+index) * *(vec2+index);
	}
	return sum;
}

void getParabola(double *xcoords, double *ycoords, double *A, double *B, double *C)
{
	//This functions solves for a parabola of the form y = Ax^2 + Bx +C
	//Using the first 3 x,y pairs in xcoords and ycoords

	double denom = (*(xcoords)-*(xcoords+1)) * (*(xcoords)-*(xcoords+2)) * (*(xcoords+1)-*(xcoords+2));

	*A = (*(xcoords+2) * (*(ycoords+1) - *(ycoords))+ *(xcoords+1) * (*(ycoords) - *(ycoords+2)) + *(xcoords) * (*(ycoords+2) - *(ycoords+1) ))/denom;
	*B = (pow(*(xcoords+2),2.0) * (*(ycoords) - *(ycoords+1)) + pow(*(xcoords+1),2.0) * (*(ycoords+2) - *(ycoords)) + pow(*(xcoords),2.0) * (*(ycoords+1) - *(ycoords+2) ))/denom;
	*C = (*(xcoords+1) * *(xcoords+2) * (*(xcoords+1)- *(xcoords+2)) * *(ycoords) + *(xcoords+2) * *(xcoords) * (*(xcoords+2)-*(xcoords)) * *(ycoords+1) + *(xcoords) * *(xcoords+1) * (*(xcoords)- *(xcoords +1)) * *(ycoords+2) )/denom;
if(truthmathverbose)	
printf("Coords (%g,%g) (%g,%g) (%g,%g) got parabola y = %gx^2 + %gx + %g\n",*(xcoords),*(ycoords),*(xcoords+1),*(ycoords+1),*(xcoords+2),*(ycoords+2),*A,*B,*C); 
}

double trapezoidalintegration(double *xcoords, double *ycoords,int numpoints)
{
	double increment = (*(xcoords + numpoints-1)- *xcoords)/(numpoints-1);
	double *weights = malloc(numpoints*sizeof(double));
	int index;
	for (index = 0; index < numpoints; index++)
	{
		if (index == 0 || index == (numpoints-1))
			*(weights+index)= .5 * increment; //Endpoints are weighted 1/2
		else
			*(weights+index) = increment;
	}
	return dotproduct(ycoords,weights,numpoints);
}

double simpsonsruleintegration(double *xcoords, double *ycoords,int numpoints)
{
	double increment = (*(xcoords + numpoints-1)- *xcoords)/(numpoints-1);
	double *weights = malloc(numpoints*sizeof(double));
	int index;
	for (index = 0; index < numpoints; index++)
	{
		if (index % 2 == 0)
		{
			if (index == 0 ||index == numpoints-1)
			{
				*(weights+index)= increment/3;//left and right endpoints are at even indices and weighted 1/3
			}
			else
			{
				*(weights+index)= 2*increment/3;//Even indices that are not endpoints are weighted 2/3
			}
		}
		else
		{
			*(weights+index)= 4 * increment/3; //Odd indices are weighted 4/3
		}
	}
		return dotproduct(ycoords,weights,numpoints);
}

double projection(double x1, double y1, double x2, double y2)
{
	//proj(a->b) = |a|cos(q)/|b| = |a.b|/|b|^2
if(truthmathverbose)
printf("Vector 1 %g %g Vector 2 %g %g\n",x1,y1,x2,y2);
	double proj = (x1*x2+y1*y2)/(pow(x2,2.0)+pow(y2,2.0));
if(truthmathverbose)
printf("Projections %g\n",proj);
	if(fabs(proj)>10e10)
		return 1;
	else
	return proj;
}

void crossproduct(double *vector, double *v1, double *v2)
{
	int i1,i2,i3;
	//double *vector=malloc(3*sizeof(double));
	//vector[0]=vector[1]=vector[2]=0;
	
	vector[0] = v1[1] * v2[2] - v1[2] * v2[1];
	vector[1] = -v1[0] * v2[2] + v1[2] * v2[0];
	vector[2] = v1[0] * v2[1] - v1[1] * v2[0];
	
	//return vector;
}

void unitvector(double *vector, double phi, double theta)
{
	//Given 2 polar angles, returns a 3 dimensional vector with unit length
	//double *vector = malloc(3*sizeof(double));
	vector[0] = sin(theta) * cos(phi);
	vector[1] = sin(theta) * sin(phi);
	vector[2] = cos(theta);
	//return vector;
}

void attenuate(double *v1, double *v2, double attenuation)
{
	//Attenuates vector v1 in the direction of vector v2
	//so that |v1| = constant
	//but the direction is closer to that of v2
	//v1 = (1-a) v1 + a * v2 * |v1|/|v2|
	
	double mag1 = pow(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2),0.5);
	double mag2 = pow(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2),0.5);
	//printf("Original %g %g %g vector %g %g %g mag1 %g mag2 %g\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],mag1,mag2);
	
	v1[0] = (1-attenuation) * v1[0] + attenuation * v2[0] * mag1/mag2;
	v1[1] = (1-attenuation) * v1[1] + attenuation * v2[1] * mag1/mag2;
	v1[2] = (1-attenuation) * v1[2] + attenuation * v2[2] * mag1/mag2;
	//printf("Attenuated Vector %g %g %g\n",v1[0],v1[1],v1[2]);
	
}


void toPolar(double x, double y, double z, double *r, double *phi, double *theta)
{
	*(r) = sqrt(x*x+y*y+z*z);
	*(phi) = atan2(y,x);
	if(*(phi) < 0)
	*(phi) += 2 *M_PI;
    double s = sqrt(x*x+y*y);
    *(theta) = acos(z/ *r);
    if(*theta < 0)
    *theta = M_PI/2 - *theta; 
}

void bubblesort(double *array, void *data, int numberOfRows, int datasize)
{
	/****************************************************************************************/
	/* This function uses the bubble sorting algorithm to sort the array and swaps the data */
	/* in the same way the array is manueveured. This was created             */
	/* for an array of eigenvalues and a matrix of eigenvectors, but may extend to other    */
	/* problems as well. The array is sorted from lowest to highest.						*/
	/****************************************************************************************/
	//printf("Sorting Eigenvalues\n");
	
	double swapSpace;
	void *swapdata = malloc(datasize);
	int outerIndex, innerIndex;
	
	for (outerIndex = 0;outerIndex<(numberOfRows-1) ;outerIndex++ )
	{
		for (innerIndex = 0;innerIndex<(numberOfRows-outerIndex-1) ;innerIndex++ )
		{
			if (*(array+innerIndex) > *(array+innerIndex+1))
			{
				swapSpace = *(array+innerIndex);
				*(array+innerIndex) = *(array+innerIndex+1);
				*(array+innerIndex+1) = swapSpace;

				memcpy(swapdata,data+datasize*innerIndex,datasize);
				memcpy(data+datasize*innerIndex,data+datasize*(innerIndex+1),datasize);
				memcpy(data+datasize*(innerIndex+1),swapdata,datasize);
			}
		}
	}
}

void bubblesorttwodatas(double *array, void *data1, void *data2, int numberOfRows, int datasize1,int datasize2)
{
	/****************************************************************************************/
	/* This function uses the bubble sorting algorithm to sort the array and swaps the data */
	/* in the same way the array is manueveured. This was created             */
	/* for an array of eigenvalues and a matrix of eigenvectors, but may extend to other    */
	/* problems as well. The array is sorted from lowest to highest.						*/
	/****************************************************************************************/
	//printf("Sorting Eigenvalues\n");
	
	double swapSpace;
	void *swapdata1 = malloc(datasize1);
	void *swapdata2 = malloc(datasize2);

	double epsilon = 0.0001;

	int outerIndex, innerIndex;
	
	for (outerIndex = 0;outerIndex<(numberOfRows-1) ;outerIndex++ )
	{
		for (innerIndex = 0;innerIndex<(numberOfRows-outerIndex-1) ;innerIndex++ )
		{
			if (*(array+innerIndex) - *(array+innerIndex+1) > epsilon)
			{
				swapSpace = *(array+innerIndex);
				*(array+innerIndex) = *(array+innerIndex+1);
				*(array+innerIndex+1) = swapSpace;

				memcpy(swapdata1,data1+datasize1*innerIndex,datasize1);
				memcpy(data1+datasize1*innerIndex,data1+datasize1*(innerIndex+1),datasize1);
				memcpy(data1+datasize1*(innerIndex+1),swapdata1,datasize1);
				
				memcpy(swapdata2,data2+datasize2*innerIndex,datasize2);
				memcpy(data2+datasize2*innerIndex,data2+datasize2*(innerIndex+1),datasize2);
				memcpy(data2+datasize2*(innerIndex+1),swapdata2,datasize2);
			}
		}
	}
}

void bubblesortthreedatas(double *array, void *data1, void *data2, void *data3, int numberOfRows, int datasize1,int datasize2, int datasize3)
{
	/****************************************************************************************/
	/* This function uses the bubble sorting algorithm to sort the array and swaps the data */
	/* in the same way the array is manueveured. This was created             */
	/* for an array of eigenvalues and a matrix of eigenvectors, but may extend to other    */
	/* problems as well. The array is sorted from lowest to highest.						*/
	/****************************************************************************************/
	//printf("Sorting Eigenvalues\n");
	
	double swapSpace;
	void *swapdata1 = malloc(datasize1);
	void *swapdata2 = malloc(datasize2);
	void *swapdata3 = malloc(datasize3);
	double epsilon = 0.0001;

	int outerIndex, innerIndex;
	
	for (outerIndex = 0;outerIndex<(numberOfRows-1) ;outerIndex++ )
	{
		for (innerIndex = 0;innerIndex<(numberOfRows-outerIndex-1) ;innerIndex++ )
		{
			if (*(array+innerIndex) - *(array+innerIndex+1) < epsilon)
			{
				swapSpace = *(array+innerIndex);
				*(array+innerIndex) = *(array+innerIndex+1);
				*(array+innerIndex+1) = swapSpace;

				memcpy(swapdata1,data1+datasize1*innerIndex,datasize1);
				memcpy(data1+datasize1*innerIndex,data1+datasize1*(innerIndex+1),datasize1);
				memcpy(data1+datasize1*(innerIndex+1),swapdata1,datasize1);
				
				memcpy(swapdata2,data2+datasize2*innerIndex,datasize2);
				memcpy(data2+datasize2*innerIndex,data2+datasize2*(innerIndex+1),datasize2);
				memcpy(data2+datasize2*(innerIndex+1),swapdata2,datasize2);
				
				memcpy(swapdata3,data3+datasize3*innerIndex,datasize3);
				memcpy(data3+datasize3*innerIndex,data3+datasize3*(innerIndex+1),datasize3);
				memcpy(data3+datasize3*(innerIndex+1),swapdata3,datasize3);
			}
		}
	}
}

void bubbleintsort(int *array, void *data, int numberOfRows, int datasize)
{
	/****************************************************************************************/
	/* This function uses the bubble sorting algorithm to sort the array and swaps the data */
	/* in the same way the array is manueveured. This was created             */
	/* for an array of eigenvalues and a matrix of eigenvectors, but may extend to other    */
	/* problems as well. The array is sorted from lowest to highest.						*/
	/****************************************************************************************/
	printf("Sorting Eigenvalues\n");
	
	int swapSpace;
	void *swapdata = malloc(datasize);
	int outerIndex, innerIndex;
	
	for (outerIndex = 0;outerIndex<(numberOfRows-1) ;outerIndex++ )
	{
		for (innerIndex = 0;innerIndex<(numberOfRows-outerIndex-1) ;innerIndex++ )
		{
			if (*(array+innerIndex) > *(array+innerIndex+1))
			{
				swapSpace = *(array+innerIndex);
				*(array+innerIndex) = *(array+innerIndex+1);
				*(array+innerIndex+1) = swapSpace;

				memcpy(swapdata,data+datasize*innerIndex,datasize);
				memcpy(data+datasize*innerIndex,data+datasize*(innerIndex+1),datasize);
				memcpy(data+datasize*(innerIndex+1),swapdata,datasize);
			}
		}
	}
}


void circlecenter(double x1, double y1, double x2, double y2, double x3, double y3, double *x0, double *y0)
{
	double r1sq = pow(x1,2) + pow(y1,2);
	double r2sq = pow(x2,2) + pow(y2,2);
	double r3sq = pow(x3,2) + pow(y3,2);

	double s = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
	
	*x0 = (r1sq*(y2-y3)+r2sq*(y3-y1)+r3sq*(y1-y2))/s/2;
	*y0 = (r1sq*(x3-x2)+r2sq*(x1-x3)+r3sq*(x2-x1))/s/2;
	
}

int findfactor(int number,int lower, int upper)
   {
	   //Tries to find a factor between lower and upper inclusive
	   //returns the lowest number found that is a factor
	   //if no factor is found returns 0
	   int index;
	   for(index=lower;index<=upper;index++)
		   if(number % index == 0)
			   return index;
	   
	   return 0;
   }

void ctofortran(double *matrix,double *fortranMatrix,int numberOfRows,int numberOfColumns)
{
	/******************************************************************/
	/* This function converts a row-major C matrix into a column-major*/
	/* fortran "array".												  */
	/******************************************************************/
	int columnNumber,rowNumber;
	for (columnNumber = 0;columnNumber<numberOfColumns ;columnNumber++ )
	{
		for (rowNumber = 0;rowNumber<numberOfRows ;rowNumber++ )
		{
			fortranMatrix[rowNumber+numberOfRows*columnNumber] = *(matrix + columnNumber+numberOfColumns*rowNumber);
		}
	}
}



void copydoubles(double *target, double *input, int numcopy)
{
	int index;
	for(index = 0;index<numcopy;index++)
	{
	 *(target+index) = *(input+index);
	}
}

void copyints(int *target, int *input, int numcopy)
{
	int index;
	for(index = 0;index<numcopy;index++)
	{
	 *(target+index) = *(input+index);
	}
}

double listmin(double *list, int num)
{
	double min = *list;
	int index;
	for(index = 1;index<num;index++)
	if(*(list+index) < min)
	min = *(list+index);
	return min;
}

int listintmin(int *list, int num)
{
	int min = *list;
	int index;
	for(index = 1;index<num;index++)
	if(*(list+index) < min)
	min = *(list+index);
	return min;
}

double listmax(double *list, int num)
{
	double max = *list;
	int index;
	for(index = 1;index<num;index++)
	if(*(list+index) > max)
	max = *(list+index);
	return max;
}

int listintmax(int *list, int num)
{
	int max = *list;
	int index;
	for(index = 1;index<num;index++)
	if(*(list+index) > max)
	max = *(list+index);
	return max;
}

double listavg(double *list, int num)
{
	int index;
	double sum = 0;
	for(index=0;index<num;index++)
	sum += *(list+index);
	return sum/num;
}

double listintavg(int *list, int num)
{
	int index;
	double sum = 0;
	for(index=0;index<num;index++)
	sum += *(list+index);
	return sum/num;
}




double liststddev(double *list, int num)
{
	double avg = listavg(list,num);
	double variance=0;
	int index;
	for(index=0;index<num;index++)
	variance += pow(*(list+index)-avg,2)/num;
	
	return sqrt(variance);
}

/*
void listcat(double *list1, int num1, double *list2, int num2)
{
	/*************************************************************/
	/* Concatenates list2 to the end of list1, assumes list1 is large enough to accept the data
	 * 
	 * */
	/* 
	 int index;
	 for(index = 0;index<num2;index++)
	 {
		 *(list1+index+num1) = *(list2+index);
	 }
}*/

double * matrixpoint(double *matrix, double *point, int dimensions, int num)
{
	/*
	 * Repeatedly applies matrix to a point and generates a list of points created.
	 * */
	 double *list = malloc(dimensions*num*sizeof(double));
	 
	 
	 int index;
	 for(index = 0;index<dimensions;index++)
	 *(list+index) = *(point+index);
	 for(index = 1;index<num;index++)
	 matrixdoublemultiplication(list+dimensions*index,matrix,list + dimensions*(index-1),dimensions,dimensions,1);
	return list;
}

void linearinterpolation(double *start, double *end, double *points, int numdimensions, int numpoints)
{
	double *startcopy = malloc(numdimensions*sizeof(double));
	double *endcopy = malloc(numdimensions*sizeof(double));
	copydoubles(startcopy,start,numdimensions);
	copydoubles(endcopy,end,numdimensions);
	int pointindex,dimensionindex;
if(truthmathverbose)
printf("Start %g %g end %g %g\n",*start,*(start+1),*end,*(end+1));
	for(pointindex=0;pointindex<numpoints;pointindex++)
	{
		for(dimensionindex=0;dimensionindex<numdimensions;dimensionindex++)
		{
		*(points+pointindex*numdimensions+dimensionindex) = (*(endcopy+dimensionindex)-*(startcopy+dimensionindex))	* (double)pointindex/(double)(numpoints) + *(startcopy+dimensionindex);
if(truthmathverbose)
printf("%g\t",*(points+pointindex*numdimensions+dimensionindex));
		}
if(truthmathverbose)
printf("\n");
	}
if(truthmathverbose)
printf("\n\n\n");
free(startcopy);
free(endcopy);
}
void parabolicinterpolation(double *start,double *middle, double *end, double *points, int numpoints)
{
	int pointindex;
if(truthmathverbose)
printf("Start %g %g end %g %g\n",*start,*(start+1),*end,*(end+1));
	double *xvalues = malloc(3*sizeof(double));
	double *yvalues = malloc(3*sizeof(double));
	*xvalues = *start;
	*(xvalues + 1) = *middle;
	*(xvalues + 2) = *end;
	*yvalues = *(start+1);
	*(yvalues+1) = *(middle+1);
	*(yvalues+2) = *(end+1);
	double A,B,C, x;
	getParabola(xvalues,yvalues,&A,&B,&C);
	for(pointindex=0;pointindex<numpoints;pointindex++)
	{
		
		x = *(points+pointindex*2) = (*(end)-*(start))	* (double)pointindex/(double)(numpoints-1) + *(start);
		*(points+pointindex*2+1) = A*pow(x,2)+B*x+C;
if(truthmathverbose)
printf("%g %g\n",*(points+pointindex*2),*(points+pointindex*2+1));
	}
if(truthmathverbose)
printf("\n\n\n");
	free(xvalues);
	free(yvalues);
}

void circle(double *x, double *y, int numpoints, double radius, double xcenter, double ycenter)
{

		int index;
		for (index = 0;index<numpoints;index++)
		{
		*(x+index)= xcenter + radius * cos(2*M_PI/index);
		*(y+index) = ycenter + radius * sin(2*M_PI/index);
	}
}
	
void getsquare(int *x, int *y, int radius)
    {
    	//creates the perimeter of a square in no particular order;
    	//the size of x and y should be 4 * radius -2
    	int index;
    	for(index=0;index<radius;index++)
    	{
    		x[index]=index;
    		y[index]=0;
    		if(index != 0)
    		{
    		x[index+radius-1] = radius-1;
    		y[index+radius-1] = index;
    		x[index+3*radius-2] = 0;
    		y[index+3*radius-2] = index;
    		}
    		x[index+2*radius-1] = index;
    		y[index+2*radius-1] = radius-1;
    		
    	}
  
    	
    }
    
    double towardzero(double number1, double number2)
{
	double sign = number1/abs(number1);
	double value = number1 - sign * abs(number2);
	if(value * number1 >0 )
		return value;
	else 
		return 0;
}
	long gcd(long number1, long number2)
	{
		if(number1 == 0)
		return number2;
		if(number2 == 0)
		return number1;
		
		return gcd((long)fmax(number1,number2) % (long)fmin(number1,number2),(long)fmin(number1,number2));
	}
	
	long lcm(long number1, long number2)
	{
		int multiple1 = number1,multiple2 = number2;
		
		while(multiple1 != multiple2)
		{
			if(multiple1<multiple2)
			multiple1+=number1;
			else
			multiple2+=number2;
		}
		return multiple1;
	}
	long arithmaticmean(long num1, long num2)
	{
		if((num1 + num2) % 2 == 0)
		return (num1+num2)/2;
		else
		return 0;
	}
	
	long geometricmean(long num1, long num2)
	{
		long product = num1 * num2;
		if((long)pow(product,0.5) == pow(product,0.5))
		return pow(product,0.5);
		else
		return 0;
	}
	
	void insertsort(int *list, int numitems,int newnumber, int ascending)
	{
		int lower = 0, upper = numitems;
		int middle = (lower + upper)/2;
		for(;lower != upper;middle = (lower + upper)/2)
		{
			if(ascending)
			{
			if(newnumber > *(list+middle) && middle != lower)
			lower = middle;
			else
			upper = middle;
			}
			else
			{
			if(newnumber <= *(list+middle) && middle != lower)
			lower = middle;
			else
			upper = middle;
			}
		}
		int index;
		for(index = numitems;index>middle;index--)
		*(list+index) = *(list+index-1);
		*(list + middle) = newnumber;

	}
	
	void printintlist(int *list, int numprint)
	{
		int index;
		for(index = 0;index<numprint;index++)
		printf("%d\n",*(list+index));
	}
	
	void printdoublelist(double *list, int numprint)
	{
		int index;
		for(index = 0;index<numprint;index++)
		printf("%g\n",*(list+index));
	}
	
	void matrixintmultiplication(int *output, int *a, int *b, int dim1, int dim2, int dim3)
	{
		/*****************************************************/
		/* Perform O(dim1*dim3) = A(dim1*dim2) X B(dim2*dim3)*/
		/*****************************************************/
		
		int index1,index2,index3;
		
		int *values = malloc(dim1*dim3*sizeof(double));
	
	
		for(index1=0;index1<dim1;index1++)
		for(index3=0;index3<dim3;index3++)
		for(index2=0,*(values+index1*dim3+index3)=0;index2<dim2;index2++)
		*(values+index1*dim3+index3) += *(a + index1*dim2+index2) * *(b+index2*dim3+index3);
	
		for(index1=0;index1<dim1;index1++)
		for(index3=0;index3<dim3;index3++)
		*(output+index1*dim3+index3) = *(values+index1*dim3+index3);
		
		free(values);
	}
	double pointplanedistance(double *point, double *plane)
	{
		//Distance between point (x,y,z) and plane Ax + By + Cz = D
		return fabs(dotproduct(point,plane,3) + plane[3])/sqrt(pow(plane[0],2) + pow(plane[1],2) + pow(plane[2],2));
	}
	
	int withinbox(double *position, double *latticevectors)
	{
		//Determines whether the position is within a box produced by the 3 lattice vectors
		
		double n1[4],n2[4],n3[4];
		crossproduct(n1,latticevectors,latticevectors+3);
		crossproduct(n2,latticevectors+6,latticevectors);
		crossproduct(n3,latticevectors+3,latticevectors+6);
		double outter[3];
		int i,j;
		for(i=0;i<3;i++)
		{
			outter[i] = 0;
		for(j=0;j<3;j++)
		outter[i] += latticevectors[3*j+i];
		}
		n1[3] = dotproduct(n1,outter,3);
		n2[3] = dotproduct(n2,outter,3);
		n3[3] = dotproduct(n3,outter,3);
		if(pointplanedistance(position,n1) > n1[3])
		return 0;
		if(pointplanedistance(position,n2) > n2[3])
		return 0;
		if(pointplanedistance(position,n3) > n3[3])
		return 0;
		double d1 = n1[3],d2 = n2[3], d3 = n3[3];
		n1[3] = n2[3] = n3[3] = 0;
		if(pointplanedistance(position,n1) >= d1)
		return 0;
		if(pointplanedistance(position,n2) >= d2)
		return 0;
		if(pointplanedistance(position,n3) >= d3)
		return 0;
		
		return 1;
	}
	
	void matrixdoublemultiplication(double *output, double *a, double *b, int dim1, int dim2, int dim3)
	{
		/*****************************************************/
		/* Perform O(dim1*dim3) = A(dim1*dim2) X B(dim2*dim3)*/
		/*****************************************************/
		
		int index1,index2,index3;
		
		double *values = malloc(dim1*dim3*sizeof(double));
	
	
		for(index1=0;index1<dim1;index1++)
		for(index3=0;index3<dim3;index3++)
		for(index2=0,*(values+index1*dim3+index3)=0;index2<dim2;index2++)
		*(values+index1*dim3+index3) += *(a + index1*dim2+index2) * *(b+index2*dim3+index3);
	
		for(index1=0;index1<dim1;index1++)
		for(index3=0;index3<dim3;index3++)
		*(output+index1*dim3+index3) = *(values+index1*dim3+index3);
		
		free(values);
	}
	double determinant(double *matrix, int size);

double *transpose(double *matrix, int a, int b)
{
	double *transposedmatrix = malloc(a * b * sizeof(double));
	int i,j;
	for(i = 0;i<a;i++)
	for(j = 0;j<b;j++)
	*(transposedmatrix + j*a + i) = *(matrix +i * b + j);
	return transposedmatrix;
}

double singlecofactor(double *matrix, int size, int a, int b)
{
	//printf("Finding single cofactor %d x %d of matrix\n",a,b);
	//printmatrix(matrix,size,size);
			int k,l, kprime, lprime; //The indices of the submatrix
			double *submatrix = malloc((size-1)*(size-1)*sizeof(double));
			double value;
			
			for(k = 0;k<size-1;k++)
			{
				if(k>=a)
				kprime = k +1;
				else
				kprime = k;
			for(l = 0;l<size-1;l++)
			{
				if(l>=b)
				lprime = l+1;
				else
				lprime = l;
				*(submatrix+k*(size-1) +l) = *(matrix+kprime*(size) + lprime);
				}
			}
			value =determinant(submatrix,size-1);
			free(submatrix); 
			return value;
}

double *cofactor (double *matrix, int size)
{
	//Constructs a cofactor matrix from a matrix
	
	
		double *outputmatrix = malloc(size*size*sizeof(double));
		int i,j; //The indices of the matrix to use
		for(i = 0;i<size;i++)
		for(j = 0;j<size;j++)
		{
			*(outputmatrix + i*size+j) = pow(-1,i+j) * singlecofactor(matrix,size,i,j);
		}
	//printf("Cofactor Matrix\n");
	//printmatrix(outputmatrix,size,size);
	return outputmatrix;
}

double determinant(double *matrix, int size)
{
	//printf("Finding determinant of matrix\n");
	//printmatrix(matrix,size,size);
	if(size == 1)
	return *matrix;
	
	else
	{
		int i;
		double value = 0;
		for(i = 0;i<size;i++)
		value += pow(-1,i) * *(matrix+i) * singlecofactor(matrix,size,0,i);
		return  value;

	}

}

double *inversematrix(double *matrix, int size)
{
	//printf("Inverting matrix \n");
	//printmatrix(matrix,size,size);
	double *cofact = cofactor(matrix,size);
	double *outputmatrix = transpose(cofact,size,size);
	free(cofact);
	int i,j;
	double factor = determinant(matrix,size);
	for(i = 0;i<size;i++)
	for(j = 0;j<size;j++)
	*(outputmatrix +i*size+j) /= factor;
	return outputmatrix;
}

void maskcompress(double *input, double *output, int numdim, int N, int M, int *mask)
{
	/* Compresses a NxNxNxNxN... array/matrix/hypermatrix into a
	 * N-MxN-MxN-MxN-M.... array/matrix/hypermatrix
	 * by leaving out rows/columns/hypercolumns what have a mask value of 1
	 * 
	 * */
	 
	 int inputindices[numdim];//Ordered with least significant first
	 int inputindex;
	 int skipped[numdim];
	 int dimindex;
	 int outputindex;
	 
	 int counter;
	 int skipable;
	 for(inputindex = 0;inputindex<pow(N,numdim);inputindex++)
	 {
		 skipable = 0;
		 for(counter = inputindex,dimindex = 0;dimindex<numdim;dimindex++)
		 {
			 *(inputindices+dimindex) = counter % N;
			 counter /= N;
			 if(*(inputindices+dimindex) == 0)
			 *(skipped+dimindex) = 0;
			 
			if(*(mask+*(inputindices+dimindex)) == 1)
			 {
				 inputindex += pow(N,dimindex) - 1;
				 skipable = 1;
				 (*(skipped+dimindex))++;
				 break;
			 }
			 
		 }
		 	 if(skipable)
			 continue;
			 else
			 {
				 //printf("\n");		
				 for(dimindex = numdim-1,outputindex = 0;dimindex>=0;dimindex--)
				 {
				//	 printf("Dimension %d in %d out %d\t",dimindex,*(inputindices+dimindex),*(inputindices+dimindex) - *(skipped+dimindex));
					 
					 outputindex = outputindex * (N-M) + *(inputindices+dimindex) - *(skipped+dimindex);
				 }
				// printf("Input index %d output index %d\n",inputindex,outputindex);
				 *(output + outputindex) = *(input + inputindex);
			 }
	 }
	 
}

void maskdecompress(double *output, double *input,  int numdim, int N, int M, int *mask)
{
	/* Compresses a N=MxN-MxN-MxN-MxN-M... array/matrix/hypermatrix into a
	 * NxNxNxN.... array/matrix/hypermatrix
	 * by skipping over rows/columns/hypercolumns what have a mask value of 1 (and setting them to 0
	 * 
	 * */
	 
	 int outputindices[numdim];//Ordered with least significant first
	 int outputindex;
	 int skipped[numdim];
	 int dimindex;
	 int inputindex;
	 
	 int counter;
	 int skipable;
	 
	 int next;
	 for(outputindex = 0;outputindex<pow(N,numdim);outputindex++)
	 {
		 skipable = 0;
		 for(counter = outputindex,dimindex = 0;dimindex<numdim;dimindex++)
		 {
			 *(outputindices+dimindex) = counter % N;
			 counter /= N;
			 if(*(outputindices+dimindex) == 0)
			 *(skipped+dimindex) = 0;
			 
			if(*(mask+*(outputindices+dimindex)) == 1)
			 {
				 next = outputindex + pow(N,dimindex) - 1;
				 for(;outputindex<=next && outputindex<pow(N,numdim);outputindex++)
				 {
				 (*(output+outputindex)) = 0;
				}
				 outputindex--;
				 //outputindex + outputindex += pow(N,dimindex) - 1; // fill zeros
				 skipable = 1;
				 (*(skipped+dimindex))++;
				 break;
			 }
			 
		 }
		 	 if(skipable)
			 continue;
			 else
			 {
				 //printf("\n");		
				 for(dimindex = numdim-1,inputindex = 0;dimindex>=0;dimindex--)
				 {
				//	 printf("Dimension %d in %d out %d\t",dimindex,*(inputindices+dimindex),*(inputindices+dimindex) - *(skipped+dimindex));
					 
					 inputindex = inputindex * (N-M) + *(outputindices+dimindex) - *(skipped+dimindex);
				 }
				// printf("Input index %d output index %d\n",inputindex,outputindex);
				 *(output + outputindex) = *(input + inputindex);
			 }
	 }
	 
}

		
	int inttrace(int *matrix, int dimension)
	{
		/******************************************/
		/* Returns the trace of a square matrix   */
		/******************************************/
		
		int sum =0;
		int index;
		for(index=0;index<dimension;index++)
		sum+= *(matrix+index*dimension+index);
		printf("Trace is %d\n",sum);
		return sum;
	}
	
	/***************************************/
	/* GRAPHS                              */
	/***************************************/
	
	void graphdistance(int *adjmatrix, int numnodes, int *output)
	{
		/*Computes the distance between each node in a graph. Stops when all offdiagonal entries have been filled */
		int i,j;
		for(i=0;i<numnodes;i++)
		for(j=0;j<numnodes;j++)
		output[i*numnodes+j] = 0;
		
		int changed = 1;
		int distance;
		int matrixpower[numnodes*numnodes];
		int maxdistance = 100;
		copyints(matrixpower,adjmatrix,numnodes*numnodes);
		
		for(distance = 1;changed && distance < maxdistance;distance++)
		{
			changed = 0;
			if(distance != 1)
			matrixintmultiplication(matrixpower,adjmatrix,matrixpower,numnodes,numnodes,numnodes);
			for(i=0;i<numnodes;i++)
			for(j=0;j<numnodes;j++)
			{
				if(i != j && !output[i*numnodes+j] && matrixpower[i*numnodes+j])
				{
					changed = 1;
					output[i*numnodes+j] = distance;
				}
			}
		}
		
	}
	double globalclustering(int *adjmatrix, int numnodes)
	{
		int *asquared = malloc(numnodes*numnodes*sizeof(int));
		int *acubed = malloc(numnodes*numnodes*sizeof(int));
		
		matrixintmultiplication(asquared,adjmatrix,adjmatrix,numnodes,numnodes,numnodes);
		matrixintmultiplication(acubed,asquared,adjmatrix,numnodes,numnodes,numnodes);
		
		int t1 = inttrace(acubed,numnodes);
		int t2 = inttrace(asquared,numnodes);
		
		int t3;
		int index;
		for(index=0,t3=0;index<numnodes*numnodes;index++)
		t3+=*(asquared+index);
		
		return (double)t1/((double)t3-(double)t2);
	}
	
	void clustering(int *adjmatrix, int numnodes, double *clustering)
	{
		/**********************************************************/
		/* Compute the clustering coefficient for each node       */
		/* using the formula c(i) = 2(A^3)ii/ki(ki-1)             */
		/* where ki = (A^2)ii                                     */
		/**********************************************************/
		
		int *asquared = malloc(numnodes*numnodes*sizeof(int));
		int *acubed = malloc(numnodes*numnodes*sizeof(int));
		
		matrixintmultiplication(asquared,adjmatrix,adjmatrix,numnodes,numnodes,numnodes);
		matrixintmultiplication(acubed,asquared,adjmatrix,numnodes,numnodes,numnodes);
		
		int index,element;
		for(index=0;index<numnodes;index++)
		{
			element = index*numnodes+index;
		*(clustering+index) = 2* *(acubed + element) / ( *(asquared+ element) * (*(asquared+ element)-1));
		}
	}


double distancebetweenpoints(double *vector1, double *vector2, int numdimensions, int periodic, double *boundary)
	{
		int index;
		double distance = 0;
		double displacement;
		for(index = 0;index<numdimensions;index++)
		{
		displacement = *(vector1+index) - *(vector2+index);
		if(periodic && fabs(displacement) > (*(boundary+index))/2)
		displacement = *(boundary+index) - fabs(displacement);
		distance += pow(displacement,2);
	}
		distance = pow(distance,0.5);
		return distance;
	}

int factorial(int i)
{
if(i == 0)
return 1;
else
return i * factorial(i-1);
}

double norm(double *vector, int numdimensions)
{
	double value = 0;
	int i;
	for(i = 0;i<numdimensions;i++)
	value += *(vector+i) * *(vector+i);
	
	return sqrt(value);
}

double directionalcosine(double *a, double *b, int numdimensions)
{
	return dotproduct(a,b,numdimensions)/norm(a,numdimensions)/norm(b,numdimensions);
}

double doubletrace(double *matrix, int size)
{
	double sum=0;
	int index;
	for(index=0;index<size;index++)
	sum+= matrix[index*size+index];
	
	return sum;
}

int antitrace(int *matrix, int size)
{
	int sum = 0;
	int i,j;
	for(i=0;i<size;i++)
	for(j=0;j<size;j++)
	if(i != j)
	sum += matrix[i*size+j];
	
	return sum;
	
}

void graphmeasures(int *adjacencymatrix, int size)
{
	printf("Taking metrics of graph for adjacency matrix:\n");
	printintmatrix(adjacencymatrix,size,size);

	
	int degrees[size];
	int i,j;
	for(i = 0;i<size;i++)
	{
		degrees[i] = 0;
	for(j = 0;j<size;j++)
	{
		degrees[i] += adjacencymatrix[i*size+j];
		}
	}
	bubbleintsort(degrees,NULL,size,0);
	printf("Degree sequence is :\n");
	printintlist(degrees,size);
	
	//int paths, polygons;
	int maxsize = 10; //Max path/side length
	int matrixpower[size*size];
	copyints(matrixpower,adjacencymatrix,size*size);
	
	for(i=2;i<maxsize;i++)
	{
		
		matrixintmultiplication(matrixpower,adjacencymatrix,matrixpower,size,size,size);
		//printf("Adjacency matrix to power %d\n",i);
		//printintmatrix(matrixpower,size,size);
		printf("Paths of length %d %d Polygons of size %d (including repeated edges)  %d\n",i,antitrace(matrixpower,size),i,inttrace(matrixpower,size));
		//printf("Paths of length %d %d\n",i,antitrace(matrixpower,size));
	}
	
	int distancematrix[size*size];
	graphdistance(adjacencymatrix,size,distancematrix);
	printf("Distance matrix is:\n");
	printintmatrix(distancematrix,size,size);
	
	int eccentricity[size];
	for(i=0;i<size;i++)
	eccentricity[i] = listintmax(distancematrix+i*size,size);
	
	bubbleintsort(eccentricity,NULL,size,0);
	printf("Eccentricity list\n");
	printintlist(eccentricity,size);
	
	printf("Radius %d diameter %d\n",listintmin(eccentricity,size),listintmax(eccentricity,size));
	printf("Average distance %g\n",listintavg(distancematrix,size*size));
}
