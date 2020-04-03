//Author: Ajay Singh Pawar
//This file contains a program which implements floyd's algorithm parallely
//Partioning:-We do domain decomposition by defining our primitive task as the single element of the matrix. So in other words work done by the primitive task is to find the correct value of the element it is associated with.
//Communication-In the outer most loop for a given value of k we need to communicate a[k,m] to all tasks associated with column m and a[m,k] to all tasks associated with row m
//Agglomeration- We agglomerate in terms of rows because in C matrix is generally stored in row major form so it makes computation faster 
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
int minimum(int a,int b)
{
	if(a>=b)
		return b;
	else
		return a;
}
int main(int argc, char *argv[])
{
int id,p,low,high,root,n=atoi(argv[2])*atoi(argv[3]),l,h,dim=atoi(argv[2]),value,value1;
int array[n];
char *s;
FILE *file;
MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&p);		//Initializing MPI to run
	low= (n*id)/p;
	high= (((id+1)*n)/p)-1;
	if(id==(p-1))
	{
		file=fopen(argv[1],"r");
		for(int i=0 ;i<n ;i++)					//Reading argv[2]*2 elements from the file named argv[1]
		{
			fscanf(file,"%d",&array[i]);
		}
		for(int i=0 ;i<p-1; i++)
		{
			l=i*n/p;
			h=(i+1)*n/p-1;						 	//According to the partition based in low and high sending the elements to the resp.proc
			MPI_Send(&array[l],h-l+1,MPI_INT,i,0,MPI_COMM_WORLD);
		}
		fclose(file);					
	}
	else
	{
		MPI_Recv(&array[low],high-low+1,MPI_INT,p-1,0,MPI_COMM_WORLD,&status);	//Receiving the elements in the different processes
	}


	for(int k=0;k<dim;k++)
	{	

		for(int t=0;t<dim;t++)			
		{
		
			for(int z=0;z<dim;z++)
			{
		
						value=t*dim+k;						//with respect to row
						value1=k*dim+z;					//with respect to column
						for(int i=0;i<p;i++)
							if((((i*n)/p)<=value)&&(((((i+1)*n)/p)-1)>=value))
								root=i;

						MPI_Bcast(&array[value],1,MPI_INT,root,MPI_COMM_WORLD); //broadcasting elements in row
							for(int i=0;i<p;i++)
								if((((i*n)/p)<=value1)&&(((((i+1)*n)/p)-1)>=value1))
									root=i;
						MPI_Bcast(&array[value1],1,MPI_INT,root,MPI_COMM_WORLD);		//broadcasting elements of column
				if((low<=t*dim+z)&&(high>=t*dim+z))
				{
					array[t*dim+z]=minimum(array[t*dim+z], array[t*dim+k]+array[k*dim+z]);
				}
			}
		}
	}
	if(id!=0)
	{
		for(int t=0; t<n; t++)
		{
				if((low<=t)&&(high>=t))
					MPI_Send(&array[low],high-low+1,MPI_INT,0,0,MPI_COMM_WORLD);
					t+=high-low;			
		}
	}
	if(id==0)
	{
		for(int t=1; t<p; t++)
			MPI_Recv(&array[(t*n)/p],((t+1)*n)/p-(t*n)/p+1,MPI_INT,t,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		for(int i=0;i<dim;i++)
		{
			for(int j=0;j<dim;j++)
			{
				printf("%d\t",array[i*dim+j]);
			}
			printf("\n");
		}
	}
MPI_Finalize();
return 0;
}
