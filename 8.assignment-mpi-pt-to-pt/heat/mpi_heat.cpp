#include <mpi.h>
#include <math.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
 #endif

  double generate2DHeat(long N, long global_i, long global_j);

  int check2DHeat(double** H, long N, long rank, long P, long k); //this assumes array of array and grid block decomposition

 #ifdef __cplusplus
}
#endif

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1;
  }

  MPI_Init(&argc,&argv);
  long N, K;
  N = atol(argv[1]);
  K = atol(argv[2]);

  int wd,NoP;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&wd);
   MPI_Comm_size(MPI_COMM_WORLD,&NoP);

   int p = sqrt(NoP);
   long ED = N/p; 
   
   int RD = wd/p,CD = wd%p;
  // use double for heat 2d 

   double**  Harr = new double*[ED];
   double**  GArr = new double*[ED];
   for(long i=0;i<ED;i++)
     {
     Harr[i] = new double[ED];
     GArr[i] = new double[ED];
     
   }
   
   long RS = (RD*ED),CS = (CD*ED);
   long RE = RS+ED,CE = CS+ED;

  for (long global_i = RS,RST=0; global_i<RE; global_i++,RST++) {
    for (long global_j= CS,CST=0; global_j<CE; global_j++,CST++) {
       Harr[RST][CST] = generate2DHeat(N,global_i, global_j);
    }
  } 
  // write code here
  
  double *TArr = new double[ED];
  double *BArr = new double[ED];
  double *RArr = new double[ED];
  double *LArr = new double[ED];
  double *TLArr = new double[ED];
  double *BLArr = new double[ED];
  double *RLArr = new double[ED];
  double *LLArr = new double[ED];
  
  int top,bottom,left,right;
  left =CD?wd-1:-1;
  right = (CD == (p-1))?-1:wd+1;
  top = wd-p;
  bottom = wd+p;
  MPI_Status status[4];
  MPI_Request requests[8];
  long count = 0;
  double start = MPI_Wtime();
  for (long it = 0; it<K; it++) 
  {
     for(long ind =0,ite = 0;ind < ED;ind++)
       {
   LArr[ite] = Harr[ind][0];
   RArr[ite] = Harr[ind][ED-1];
   TArr[ite]  = Harr[0][ind];
   BArr[ite] = Harr[ED-1][ind];
   LLArr[ite] = Harr[ind][0];
   RLArr[ite] = Harr[ind][ED-1];
   TLArr[ite]  = Harr[0][ind];
   BLArr[ite] = Harr[ED-1][ind];
   ite++;
       }
     count = 0;
     //  cout<<"I am rank "<<wd<<endl; 
     if(top>=0){
       MPI_Isend(TLArr,ED,MPI_DOUBLE,top,0,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(right != -1){
       MPI_Isend(RLArr,ED,MPI_DOUBLE,right,1,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(bottom < NoP){
       MPI_Isend(BLArr,ED,MPI_DOUBLE,bottom,2,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(left != -1){
       MPI_Isend(LLArr,ED,MPI_DOUBLE,left,3,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     for(long i=1;(i+1)<ED;i++)
       {
   for(long j=1;(j+1)<ED;j++)
     {
       GArr[i][j] = (Harr[i][j] + Harr[i-1][j] + Harr[i+1][j] + Harr[i][j-1] + Harr[i][j+1])/5;
     }
       }
     // cout<<"I am rank "<<wd<<" and sent "<<count<<" messages"<<endl; 
     // MPI_Waitall(count,requests,status);
     if(top >= 0){
       MPI_Recv(TArr,ED,MPI_DOUBLE,top,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(right != -1){
       MPI_Recv(RArr,ED,MPI_DOUBLE,right,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(bottom < NoP){
       MPI_Recv(BArr,ED,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     if(left != -1){
       MPI_Recv(LArr,ED,MPI_DOUBLE,left,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     // cout<<"I am rank "<<wd<<" and received "<<count<<" messages"<<endl;
     GArr[0][0] = (Harr[0][0] + TArr[0] + LArr[0] + Harr[0][1] + Harr[1][0])/5;
     GArr[0][ED-1] = (Harr[0][ED-1] + TArr[ED-1] + RArr[0]+Harr[0][ED-2]+Harr[1][ED-1])/5;
     GArr[ED-1][0] = (Harr[ED-1][0] + BArr[0] + LArr[ED-1]+Harr[ED-1][1]+Harr[ED-2][0])/5;
     GArr[ED-1][ED-1] = (Harr[ED-1][ED-1] + BArr[ED-1] +
          RArr[ED-1]+Harr[ED-1][ED-2]+Harr[ED-2][ED-1])/5;
     for(long i=1,j=1;i<ED-1;i++,j++)
   {
      GArr[0][i] = (Harr[0][i]+Harr[1][i]+TArr[i]+Harr[0][i-1]+Harr[0][i+1])/5;
      GArr[ED-1][i] = (Harr[ED-1][i]+Harr[ED-2][i]+BArr[i]+Harr[ED-1][i-1]+Harr[ED-1][i+1])/5;
      GArr[j][0] = (Harr[j][0]+LArr[j]+Harr[j][1]+Harr[j-1][0]+Harr[j+1][0])/5;
      GArr[j][ED-1] = (Harr[j][ED-1]+Harr[j-1][ED-1]+Harr[j+1][ED]+RArr[j]+Harr[j][ED-2])/5;
   }
      Harr = GArr;
      // check2DHeat(double** H, long n, long rank, long P, long k)
      check2DHeat(Harr,N,wd,NoP,it);
      MPI_Waitall(count,requests,status);
             
   }
  if(wd == 0)
    {
       double end = MPI_Wtime();
       cerr<<end-start<<endl; 
    }
 for(long i=0;i<ED;i++)
   delete[] Harr[i];
  delete[] Harr;
  delete[] TArr;
  delete[] BArr;
  delete[] RArr;
  delete[] LArr;

  MPI_Finalize();

  return 0;
}