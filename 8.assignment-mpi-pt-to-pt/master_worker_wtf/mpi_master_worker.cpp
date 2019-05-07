#include <mpi.h>
  #include <iostream>
  #include <stdio.h>
  #include <stdlib.h>
  #include <cstdlib>
  #include <string.h>
  #include <chrono>
  #include <cmath>
using namespace std;
using seconds = chrono::seconds;
using check_time = std::chrono::high_resolution_clock;

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

float numerical_integration(int function_id, float a, int limitArr[], float intensity,  float pre_caculation)
{
  float sum = 0;
  
  for (int i=limitArr[0]; i<limitArr[1]; i++){

    float x = (a + (i + 0.5) * pre_caculation);
    
    switch(function_id){      
      case 1:
      sum +=f1(x, intensity);
      break;

      case 2:
      sum +=f2(x, intensity);
      break;

      case 3:
      sum +=f3(x, intensity);
      break;

      case 4:
      sum +=f4(x, intensity);
      break;
    }
    
  }
  return sum;
}

void GetArrValues(int *arr, int value){
  for(int i = 0; i<3; i++){
    arr[i] = value;
    
  }
}

int main (int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  if (argc < 6) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <function_id> <a> <b> <n> <intensity>"<<std::endl;
    return -1;
  }
  int function_id = stoi(argv[1]);
  float a = stof(argv[2]);
  float b = stof(argv[3]);
  int n = stoi(argv[4]);
  float intensity = stoi(argv[5]);
  float pre_caculation = ((b - a) / n );
  float parSum = 0;
  
  int mpiRank, Mpisize;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Mpisize);
  MPI_Status status;
  float totalSum = 0;
  int par = (n/(Mpisize - 1))/8;
  auto timeStart = check_time::now();  
  

  if(mpiRank == 0){

    int sentChunk = 0,workGiven = 0, limitArr[2], Poin;
    
    for(int i=1; i< Mpisize; i++){
     for(int chunk = 0; chunk < 3; chunk++){
      limitArr[0] = workGiven * par;
      limitArr[1] = (workGiven + 1)* par;
      Poin = limitArr[1];
      sentChunk++;
      workGiven++;

      if(limitArr[0] >= n){
        limitArr[0] = -1;
        sentChunk--;
      }else if(limitArr[1]>n){
        limitArr[1] = n;
      }

      MPI_Send(limitArr, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
    }      
  }
  while(sentChunk != 0){

    MPI_Recv(&parSum, 1, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    sentChunk--;
    totalSum+=parSum;
    limitArr[0] = Poin;
    limitArr[1] = Poin + par;
    Poin = limitArr[1];
    sentChunk++;

    if(limitArr[0] >= n){
      limitArr[0] = -1;
      sentChunk--;
    }else if(limitArr[1]>n){
      limitArr[1] = n;
    }

    MPI_Send(limitArr, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
  }
}else{
  int WRArr[2], INArr[3], BArr1[2], BArr2[2], Barr3[2], RC = 3, SArr[3], send = 0, REArr[3], outcome;
  float sum = 0.0;
  MPI_Status stat[3];
  MPI_Request request[3];

  GetArrValues(SArr, 1);
  GetArrValues(REArr, 1);
  

  while(true){
    if(SArr[0] == 1 && REArr[0] == 1){
      MPI_Irecv(BArr1, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &request[0]);
    }
    if(SArr[1] == 1 && REArr[1] == 1){
      MPI_Irecv(BArr2, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &request[1]);
    }
    if(SArr[2] == 1 && REArr[2] == 1){
      MPI_Irecv(Barr3, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &request[2]);
    }
    GetArrValues(REArr, 0);
    
    while(true){
      MPI_Waitsome(3, request, &outcome, INArr, stat);
      if(outcome > 0){
        break;
      }
    }

    for(int i=0; i<outcome; i++){
      send = 1;
      if(INArr[i] == 0){
        
        WRArr[0] = BArr1[0];
        WRArr[1] = BArr1[1];

        if(WRArr[0]==-1){
          SArr[0] = 0;
          RC--;
          send = 0;
        }else{
          REArr[0] = 1;
        }
      }if(INArr[i] == 1){
        
        WRArr[0] = BArr2[0];
        WRArr[1] = BArr2[1];

        if(WRArr[0]==-1){
          SArr[1] = 0;
          RC--;
          send = 0;
        }else{
          REArr[1] = 1;
        }
      }if(INArr[i] == 2){
        
        WRArr[0] = Barr3[0];
        WRArr[1] = Barr3[1];

        if(WRArr[0]==-1){
          SArr[2] = 0;
          RC--;
          send = 0;
        }else{
          REArr[2] = 1;
        }
      }

      if(send == 1){
        sum = numerical_integration(function_id, a, WRArr, intensity, pre_caculation);
        MPI_Send(&sum, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      }

    }
    if(RC == 0){
      break;
    }
  }
}

if(mpiRank == 0){

  cout<<totalSum*pre_caculation;
  auto timeEnd = check_time::now();
  auto TotalTime = timeEnd - timeStart;
  auto secs = std::chrono::duration_cast<std::chrono::duration<float>>(TotalTime);
  cerr<<secs.count();

}
