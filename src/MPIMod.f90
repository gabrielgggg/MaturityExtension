MODULE MPIMod
  USE iso_Fortran_env, ONLY: wp => real64
  USE mpi
  IMPLICIT NONE

  INTEGER :: mpiErr
  INTEGER :: workerId
  INTEGER :: workerNo
  INTEGER :: chunk, chunkCh, ssSzCh

  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffV, bufflV, buffVd, bufflVd
  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffQS, bufflQS, buffQL, bufflQL
  
  REAL(wp), DIMENSION(:), ALLOCATABLE :: buffdpr, buffldpr, buffchpr, bufflchpr
  INTEGER, DIMENSION(:), ALLOCATABLE :: buffchS, bufflchS, buffchL, bufflchL
  
CONTAINS

  SUBROUTINE mpiSetup(ssSz, chPr, okFlag)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ssSz, chPr
    INTEGER, INTENT(OUT) :: okFlag

    CALL MPI_Init(mpiErr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, workerNo, mpiErr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, workerId, mpiErr)

    chunk = ssSz / workerNo
    chunkCh = chunk * chPr
    ssSzCh = ssSz * chPr

    IF (MOD(ssSz, workerNo) == 0) THEN
      okFlag = 1

      !
      ! Buffers for MPI syncronization
      !
      ALLOCATE(buffV(chunk))
      ALLOCATE(buffLV(ssSz))
      
      ALLOCATE(buffVd(chunk))
      ALLOCATE(buffLVd(ssSz))
      
      ALLOCATE(buffQS(chunk))
      ALLOCATE(buffLQS(ssSz))
      
      ALLOCATE(buffQL(chunk))
      ALLOCATE(buffLQL(ssSz))
      
      ALLOCATE(buffdpr(chunk))
      ALLOCATE(buffLdpr(ssSz))
      
      ALLOCATE(buffchpr(chunkCh))
      ALLOCATE(buffLchpr(ssSzCh))
      ALLOCATE(buffchS(chunkCh))
      ALLOCATE(bufflchS(ssSzCh))
      ALLOCATE(buffchL(chunkCh))
      ALLOCATE(bufflchL(ssSzCh))
    ELSE
      okFlag = 0
    END IF
  END SUBROUTINE mpiSetup

  SUBROUTINE mpiFinalize()
    IMPLICIT NONE

    CALL MPI_Finalize(mpiErr)

    DEALLOCATE(buffV, bufflV, buffQS, bufflQS, buffQL, bufflQL)
    DEALLOCATE(buffVd, bufflVd)
    DEALLOCATE(buffdpr, buffldpr, buffchpr, bufflchpr)
    DEALLOCATE(buffchS, bufflchS, buffchL, bufflchL)
  END SUBROUTINE mpiFinalize

END MODULE MPIMod
