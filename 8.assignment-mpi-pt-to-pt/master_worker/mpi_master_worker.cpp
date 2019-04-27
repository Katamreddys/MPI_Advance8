#include <mpi.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

float f1(float x, int intensity);
float f2(float x, int intensity);
float f3(float x, int intensity);
float f4(float x, int intensity);

#ifdef __cplusplus
}
#endif


int main (int argc, char* argv[]) {

  if (argc < 6) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <functionid> <a> <b> <n> <intensity>"<<std::endl;
    return -1;
  }
  
    MPI_Init( &argc, &argv );
	
	int rankMPI,sizeMpi;
        int fid = atoi(argv[1]);
	float a  = atof(argv[2]);
	float b = atof(argv[3]);
	int n = atoi(argv[4]);
	int intensity = atoi(argv[5]);
	float (*func)(float,int);
	switch(fid)
  {
     case 1 : func = &f1;
              break;

     case 2 : func = &f2;
                break;

     case 3 : func = &f3;
                break;

     case 4 : func = &f4;
                break;

     default :  cerr<<"Functionid does not exists, Enter 1, 2, 3 or 4" <<endl; 
                return -1;                

   };
    MPI_Comm_rank( MPI_COMM_WORLD, &rankMPI );
	MPI_Comm_size( MPI_COMM_WORLD, &sizeMpi );
	float IntSumm;
	int p = n/sizeMpi;
	
    //MPI_Comm_split( MPI_COMM_WORLD, rank == 0, 0, &new_comm );
    if (rankMPI == 0) 
	{
		double Time_Start = MPI_Wtime(); 
		float ValInt = 0.0;
		for(int i=1;i<sizeMpi;i++)
		{
			
			int start = (i-1)*p,end = i*p;
			if(i == (sizeMpi-1))
				end = n;
			MPI_Send(&start,1,MPI_INT,i,0,MPI_COMM_WORLD);
			MPI_Send(&end,1,MPI_INT,i,1,MPI_COMM_WORLD);
			
		}
		int nr =1;
		while(nr != sizeMpi)
		{
			float par_inte;
			int cnt;
			MPI_Status status;
			MPI_Recv(&par_inte,1,MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
			MPI_Get_count(&status, MPI_INT, &cnt);
			int value = -1;
			MPI_Send(&value,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
			nr++;
			ValInt += par_inte;
		}
           	double tend = MPI_Wtime(); 
		cout<<((b-a)/n)*ValInt<<endl;
		cerr<<tend-Time_Start<<endl;
	}
    else
	{ 	
		int start,end;
		while(true){
			MPI_Recv(&start,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if(start == -1)
				break; 
			MPI_Recv(&end,1,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			IntSumm = 0.0;
			for(int i=start;i<end;i++)
			{
				float x = (i+0.5)*((b-a)/n);
				IntSumm += func(a+x,intensity);
			}
			MPI_Send(&IntSumm,1,MPI_FLOAT,0,2,MPI_COMM_WORLD);
		
		}
	}

    MPI_Finalize( );

  return 0;
}