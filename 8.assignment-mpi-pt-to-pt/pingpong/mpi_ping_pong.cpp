#include <mpi.h>
#include <iostream>
using namespace std;

int main (int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <element>"<<std::endl;
    return -1;
  }
  int element;
  MPI_Init(&argc,&argv);
  element = stoi(argv[1]);
	int rankMPI, sizeMPI;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankMPI);
  MPI_Comm_size(MPI_COMM_WORLD, &sizeMPI);

  int Val;
  if (rankMPI == 0) {
      //element
    Val = element;
      //send and receive
    MPI_Send(&Val, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    cout<<"Send from 0 to 1"<<endl;
    MPI_Recv(&Val, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"After adding Recive"<<endl;
    cout<<"After adding: " <<Val<<endl;
  } else if (rankMPI == 1) {
      
      //receive add and send
    MPI_Recv(&Val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"Received to 1 from 0"<<endl;

    Val = Val + 2;
    MPI_Send(&Val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    cout<<"Sending to 0 from 1 after adding"<<endl;

  }


  MPI_Finalize();
  return 0;
}