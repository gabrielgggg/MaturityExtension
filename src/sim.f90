MODULE sim
  USE iso_Fortran_env, ONLY: wp => real64
  IMPLICIT NONE

  REAL(wp) :: piwp = 4.0 * ATAN(1.0)
CONTAINS

  SUBROUTINE fixSeed()
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
    INTEGER :: n
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = 1989
    CALL RANDOM_SEED(put = seed)
  END SUBROUTINE fixSeed

  SUBROUTINE simMarkov(sz, piMat, newIx, oldIx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sz, oldIx
    REAL(wp), DIMENSION(sz, sz), INTENT(IN) :: piMat
    INTEGER, INTENT(OUT) :: newIx
!    REAL(wp) :: rno, ssum
!    INTEGER :: ix

!    CALL simUniform(rno)
    CALL simDiscrete(sz, piMat(oldIx, :), newIx)

!    ssum = 0
!    DO ix = 1,sz
!      ssum = ssum + piMat(oldIx, ix)
!      IF ( rno < ssum .AND. ssum - piMat(oldIx, ix) <= rno ) THEN
!        newIx = ix
!        EXIT
!      END IF
!    END DO
  END SUBROUTINE simMarkov

  SUBROUTINE simDiscrete(sz, pmf, drawIx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sz
    REAL(wp), DIMENSION(sz), INTENT(IN) :: pmf
    INTEGER, INTENT(OUT) :: drawIx
    REAL(wp) :: rno, ssum
    INTEGER :: ix

    CALL simUniform(rno)

    ssum = 0
    DO ix = 1,sz
      ssum = ssum + pmf(ix)
      IF ( rno < ssum .AND. rno >= ssum - pmf(ix) ) THEN
        drawIx = ix
        EXIT
      END IF
    END DO
  END SUBROUTINE simDiscrete

  SUBROUTINE simStdNormal(val)
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: val
    REAL(wp) :: v1, v2
 
    CALL RANDOM_NUMBER(v1)
    CALL RANDOM_NUMBER(v2)
    val = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * piwp * v2 )
  END SUBROUTINE simStdNormal

  SUBROUTINE simUniform(val)
    IMPLICIT NONE
    REAL(wp), INTENT(OUT) :: val

    CALL RANDOM_NUMBER(val)
  END SUBROUTINE simUniform

  SUBROUTINE simDiscreteUniform(sz, retIx)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sz
    INTEGER, INTENT(OUT) :: retIx
    REAL(wp) :: rno

    CALL simUniform(rno)
    retIx = FLOOR(sz * rno) + 1
  END SUBROUTINE simDiscreteUniform

END MODULE sim
