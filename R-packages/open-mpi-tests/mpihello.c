# include <stdio.h>
# include "mpi.h"
int main (int argc , char** argv )
{
      int rank, size , nameLen ;
      char processorName [ MPI_MAX_PROCESSOR_NAME ] ;
      MPI_Init (&argc , &argv ) ;
      MPI_Comm_rank ( MPI_COMM_WORLD, &rank ) ;
      MPI_Comm_size ( MPI_COMM_WORLD, &size ) ;
      MPI_Get_processor_name ( processorName , &nameLen ) ;
      printf( "Hello , rank %d , size %d on processor %s \n" ,
                  rank , size , processorName ) ;
      MPI_Finalize( ) ;
      return 0;
}
