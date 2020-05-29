MODULE jmpMod
  USE iso_Fortran_env, ONLY: wp => real64
  USE mpi
  USE MPIMod
  USE Param
  USE NL
  USE sim
  IMPLICIT NONE

  ! IO:
  CHARACTER (LEN=*), PARAMETER :: ioDir = "./"
  CHARACTER (LEN=*), PARAMETER :: outDir = ioDir // "results/"
  CHARACTER (LEN=*), PARAMETER :: shockDir = ioDir // "shocks/"
  CHARACTER, PARAMETER:: tab = CHAR(9)

  ! Grids
  REAL(wp), DIMENSION(ySz) :: y, EdCre
  REAL(wp), DIMENSION(ySz, ySz) :: exoTran
  REAL(wp), DIMENSION(bSSz) :: bS
  REAL(wp), DIMENSION(bLSz) :: bL

  ! ArrayIxs: values, prices, policies
  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: V0, V1
  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: Va0, Va1

  REAL(wp), DIMENSION(ySz) :: Vd0, Vd1
  REAL(wp), DIMENSION(ySz) :: hy, uhy, Vaut

  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: qS0, qS1, qL0, qL1
  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: qSa0, qSa1, qLa0, qLa1
  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: xiS, xiL, xiSkeep, xiLkeep

  REAL(wp), DIMENSION(ySz, bSSz, bLSz) :: defPr
  REAL(wp), DIMENSION(:, :, :, :), ALLOCATABLE :: chPr
  INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: chS, chL

  REAL(wp), DIMENSION(:, :, :, :), ALLOCATABLE :: chPrA
  INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: chSa, chLa

  REAL(wp), DIMENSION(ySz, chPrSz) :: gammaPr
  INTEGER, DIMENSION(ySz, chPrSz) :: gammaS, gammaL

  TYPE stateSp
    INTEGER :: yIx, bSIx, bLIx
  END TYPE stateSp
  
  TYPE(stateSp), DIMENSION(:), ALLOCATABLE :: toSS

CONTAINS
  !
  !
  !

  PURE FUNCTION uu(cc) RESULT(uuval)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: cc
    REAL(wp) :: uuval

    IF (cc <= 0.0_wp) THEN
      uuval = veryNeg
    ELSEIF (crra == 1.0) THEN
      uuval = LOG(cc)
    ELSEIF (crra == 2.0) THEN
      uuval = 1.0_wp - 1.0_wp / cc
    ELSE
      uuval = (cc**(1.0_wp - crra) - 1.0_wp) / (1.0_wp - crra)
    END IF
  END FUNCTION uu
  
  SUBROUTINE loadShocks()
    IMPLICIT NONE

    OPEN(unit=100, file=shockDir // "y.tab", status="old", action="read", position="rewind")
    READ (100, *) y
    CLOSE(100)

    OPEN(unit=100, file=shockDir // "PI.tab", status="old", action="read", position="rewind")
    READ (100, *) exoTran
    CLOSE(100)
    exoTran = TRANSPOSE(exoTran)
    exoTran(:, ySz) = MAX( 0.0_wp , 1.0_wp - SUM(exoTran(:, 1:(ySz-1)), 2) )
  END SUBROUTINE loadShocks
  
  SUBROUTINE generateGrids()
    IMPLICIT NONE
    INTEGER :: ix, yIx, bSIx, bLIx

    ALLOCATE(toSS(spaceSz))
    ix = 1
    DO bLIx = 1,bLSz
    DO bSIx = 1,bSSz
    DO yIx = 1,ySz
      toSS(ix)%yIx = yIx
      toSS(ix)%bSIx = bSIx
      toSS(ix)%bLIx = bLIx
      ix = ix + 1
    END DO
    END DO
    END DO

    CALL linspace(bS(1:bSsplit), bSmin, bSmid, bSsplit)
    CALL linspace(bS(bSsplit:bSSz), bSmid, bSmax, bSSz - bSsplit + 1)
    
    CALL linspace(bL(1:bLsplit), bLmin, bLmid, bLsplit)
    CALL linspace(bL(bLsplit:bLSz), bLmid, bLmax, bLSz - bLsplit + 1)

    DO yIx = 1,ySz
      hy(yIx) = y(yIx) - MAX( 0.0_wp, lbd0 * y(yIx) + lbd1 * y(yIx) * y(yIx) )
      ! hy(yIx) = MIN( y(yIx), 1.0_wp - lbd )
      uhy(yIx) = uu( hy(yIx) )
    END DO
  END SUBROUTINE generateGrids

  SUBROUTINE findVaut()
    IMPLICIT NONE
    INTEGER :: yIx
    REAL(wp), DIMENSION(ySz) :: Vaut0
    REAL(wp) :: errV

    Vaut = 0.0_wp

    errV = 1.0_wp
    DO WHILE (errV > epsV)
      Vaut0 = Vaut
      Vaut = (1.0_wp - beta) * uhy + beta * MATMUL( exoTran, Vaut0 )
      errV = MAXVAL(ABS(Vaut - Vaut0))
    END DO
  END SUBROUTINE

  SUBROUTINE saveEquilibrium()
    IMPLICIT NONE
    OPEN(unit=100, file=outDir // "parameters.txt")
    WRITE (100, "(I4)") ySz
    WRITE (100, "(I4)") bSSz
    WRITE (100, "(I4)") bLSz
    WRITE (100, "(I4)") chPrSz
    WRITE (100, "(F9.6)") crra
    WRITE (100, "(F9.6)") beta
    WRITE (100, "(F9.6)") rf
    WRITE (100, "(F9.6)") eta
    WRITE (100, "(F9.6)") etaA
    WRITE (100, "(F9.6)") deltaS
    WRITE (100, "(F9.6)") deltaL
    WRITE (100, "(F9.6)") alpha
    WRITE (100, "(F9.6)") muS
    WRITE (100, "(F9.6)") lbd0
    WRITE (100, "(F9.6)") lbd1
    WRITE (100, "(F15.10)") rho
    WRITE (100, "(F15.10)") rhoDef
    WRITE (100, "(F15.10)") rhoNash
    WRITE (100, "(F9.6)") adj
    WRITE (100, "(F9.6)") bSshare
    CLOSE(100)

    OPEN(unit=100, file=outDir // "bS.txt")
    WRITE (100, "(F9.6)") bS
    CLOSE(100)

    OPEN(unit=100, file=outDir // "bL.txt")
    WRITE (100, "(F9.6)") bL
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "V.txt")
    WRITE (100, "(F16.6)") V1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "Va.txt")
    WRITE (100, "(F16.6)") Va1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qL.txt")
    WRITE (100, "(F9.6)") qL1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qS.txt")
    WRITE (100, "(F9.6)") qS1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qLa.txt")
    WRITE (100, "(F9.6)") qLa1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qSa.txt")
    WRITE (100, "(F9.6)") qSa1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "xiL.txt")
    WRITE (100, "(F9.6)") xiL
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "xiS.txt")
    WRITE (100, "(F9.6)") xiS
    CLOSE(100)

    OPEN(unit=100, file=outDir // "defPr.txt")
    WRITE (100, "(F9.6)") defPr
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "Vd.txt")
    WRITE (100, "(F16.6)") Vd1
    CLOSE(100)

    OPEN(unit=100, file=outDir // "Vaut.txt")
    WRITE (100, "(F16.6)") Vaut
    CLOSE(100)

    OPEN(unit=100, file=outDir // "gammaPr.txt")
    WRITE (100, "(F9.6)") gammaPr
    CLOSE(100)

    OPEN(unit=100, file=outDir // "gammaS.txt")
    WRITE (100, "(I4)") gammaS
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "gammaL.txt")
    WRITE (100, "(I4)") gammaL
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "chPr.txt")
    WRITE (100, "(F9.6)") chPr
    CLOSE(100)

    OPEN(unit=100, file=outDir // "chS.txt")
    WRITE (100, "(I4)") chS
    CLOSE(100)    
    
    OPEN(unit=100, file=outDir // "chL.txt")
    WRITE (100, "(I4)") chL
    CLOSE(100)    

    OPEN(unit=100, file=outDir // "chPrA.txt")
    WRITE (100, "(F9.6)") chPrA
    CLOSE(100)

    OPEN(unit=100, file=outDir // "chSa.txt")
    WRITE (100, "(I4)") chSa
    CLOSE(100)    
    
    OPEN(unit=100, file=outDir // "chLa.txt")
    WRITE (100, "(I4)") chLa
    CLOSE(100)    
    
    OPEN(unit=100, file=outDir // "EdCre.txt")
    WRITE (100, "(F9.6)") EdCre
    CLOSE(100)
  END SUBROUTINE saveEquilibrium

  ! Load arrays from .txt files
  SUBROUTINE loadEquilibrium()
    IMPLICIT NONE

    IF (workerId == 0) WRITE (*, *) "Loading policies and values from files."
    IF (workerId == 0) WRITE (*, *) "Make sure the grids and parameters are the same."

    OPEN(unit=100, file=outDir // "V.txt", status="old", position="rewind")
    READ (100, "(F16.6)") V1
    V0 = V1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "Va.txt", status="old", position="rewind")
    READ (100, "(F16.6)") Va1
    Va0 = Va1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qL.txt", status="old", position="rewind")
    READ (100, "(F9.6)") qL1
    qL0 = qL1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qS.txt", status="old", position="rewind")
    READ (100, "(F9.6)") qS1
    qS0 = qS1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qLa.txt", status="old", position="rewind")
    READ (100, "(F9.6)") qLa1
    qLa0 = qLa1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "qSa.txt", status="old", position="rewind")
    READ (100, "(F9.6)") qSa1
    qSa0 = qSa1
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "xiL.txt", status="old", position="rewind")
    READ (100, "(F9.6)") xiL
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "xiS.txt", status="old", position="rewind")
    READ (100, "(F9.6)") xiS
    CLOSE(100)

    OPEN(unit=100, file=outDir // "defPr.txt", status="old", position="rewind")
    READ (100, "(F9.6)") defPr
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "Vd.txt", status="old", position="rewind")
    READ (100, "(F16.6)") Vd1
    Vd0 = Vd1
    CLOSE(100)

    OPEN(unit=100, file=outDir // "gammaPr.txt", status="old", position="rewind")
    READ (100, "(F9.6)") gammaPr
    CLOSE(100)

    OPEN(unit=100, file=outDir // "gammaS.txt", status="old", position="rewind")
    READ (100, "(I4)") gammaS
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "gammaL.txt", status="old", position="rewind")
    READ (100, "(I4)") gammaL
    CLOSE(100)
    
    OPEN(unit=100, file=outDir // "chPr.txt", status="old", position="rewind")
    READ (100, "(F9.6)") chPr
    CLOSE(100)

    OPEN(unit=100, file=outDir // "chS.txt", status="old", position="rewind")
    READ (100, "(I4)") chS
    CLOSE(100)    
    
    OPEN(unit=100, file=outDir // "chL.txt", status="old", position="rewind")
    READ (100, "(I4)") chL
    CLOSE(100)    

    OPEN(unit=100, file=outDir // "chPrA.txt", status="old", position="rewind")
    READ (100, "(F9.6)") chPrA
    CLOSE(100)

    OPEN(unit=100, file=outDir // "chSa.txt", status="old", position="rewind")
    READ (100, "(I4)") chSa
    CLOSE(100)    
    
    OPEN(unit=100, file=outDir // "chLa.txt", status="old", position="rewind")
    READ (100, "(I4)") chLa
    CLOSE(100)    

    OPEN(unit=100, file=outDir // "EdCre.txt", status="old", position="rewind")
    READ (100, "(F9.6)") EdCre
    CLOSE(100)    
  END SUBROUTINE loadEquilibrium

  ! Simulate
  SUBROUTINE simulate()
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: yPath, bSpath, bLpath, dPath
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: cPath, gdpPath, qSpath, qLpath
    INTEGER, ALLOCATABLE, DIMENSION(:) :: yEv, yDefEv, bSev, bLev, gammaSev, gammaLev
    INTEGER :: tt, evCnt, ix, drawCh, yDefaulted, returnS, returnL
    REAL(wp) :: rno, ddPr
    LOGICAL :: inDef, inArr, justReturn
    
    ALLOCATE(yPath(simSz), bSpath(simSz), bLpath(simSz), dPath(simSz))
    ALLOCATE(cPath(simSz), gdpPath(simSz), qSpath(simSz), qLpath(simSz))
    ALLOCATE(yEv(simSz), yDefEv(simSz),bSev(simSz), bLev(simSz), gammaSev(simSz), gammaLev(simSz))

    CALL fixSeed()

    yPath = 1
    bSpath = 1
    bLpath = 1
    dPath = 0
    cPath = 0.0_wp
    gdpPath = 0.0_wp
    qSpath = 0.0_wp
    qLpath = 0.0_wp

    yDefaulted = 1
    evCnt = 1

    DO tt = 2,simSz-1

      CALL simMarkov(ySz, exoTran, yPath(tt), yPath(tt-1))
      justReturn = .FALSE.
      returnS = 1
      returnL = 1

      IF ( dPath(tt-1) == 2 ) THEN ! Prev in arr
        CALL simUniform(rno)
        IF (rno < etaA) THEN ! To market
          dPath(tt) = 0
          ! Check Default
          ! Use chPr later
        ELSE ! Stay
          dPath(tt) = 2
          ! Use chPrA
          CALL simDiscrete(chPrSz, chPrA(yPath(tt), bSpath(tt-1), bLpath(tt-1), :), drawCh)
          bSpath(tt) = chSa(yPath(tt), bSpath(tt-1), bLpath(tt-1), drawCh)
          bLpath(tt) = chLa(yPath(tt), bSpath(tt-1), bLpath(tt-1), drawCh)
        END IF
      ELSE IF ( dPath(tt-1) == 1 ) THEN ! Prev in default
        CALL simUniform(rno)
        IF (rno < eta) THEN ! To arr
          dPath(tt) = 2
          ! Use gamma-s
          CALL simDiscrete(chPrSz, gammaPr(yPath(tt), :), drawCh)
          ! bSpath(tt) = gammaS(yPath(tt), drawCh)
          ! bLpath(tt) = gammaL(yPath(tt), drawCh)
          returnS = gammaS(yPath(tt), drawCh)
          returnL = gammaL(yPath(tt), drawCh)

          yEv(evCnt) = yPath(tt)
          yDefEv(evCnt) = yDefaulted
          bSev(evCnt) = bSpath(tt-1)
          bLev(evCnt) = bLpath(tt-1)
          gammaSev(evCnt) = returnS
          gammaLev(evCnt) = returnL
          evCnt = evCnt + 1
          
          justReturn = .TRUE.
          
          CALL simDiscrete(chPrSz, chPrA(yPath(tt), returnS, returnL, :), drawCh)
          bSpath(tt) = chSa(yPath(tt), returnS, returnL, drawCh)
          bLpath(tt) = chLa(yPath(tt), returnS, returnL, drawCh)
        ELSE ! Stay
          dPath(tt) = 1
          bSpath(tt) = bSpath(tt-1)
          bLpath(tt) = bLpath(tt-1)
        END IF
      ELSE
        ! CALL simUniform(rno)
        dPath(tt) = 0
        ! Check Default
        ! Use chPr later
      END IF

      IF (dPath(tt) == 0) THEN
          CALL simUniform(rno)
          ddPr = defPr(yPath(tt), bSpath(tt-1), bLpath(tt-1))
          IF (rho < ddPr) THEN
            ! Choose default
            dPath(tt) = 1
            yDefaulted = yPath(tt)
            bSpath(tt) = bSpath(tt-1)
            bLpath(tt) = bLpath(tt-1)
          ELSE
            CALL simDiscrete(chPrSz, chPr(yPath(tt), bSpath(tt-1), bLpath(tt-1), :) / (1.0_wp - ddPr), drawCh)
            bSpath(tt) = chS(yPath(tt), bSpath(tt-1), bLpath(tt-1), drawCh)
            bLpath(tt) = chL(yPath(tt), bSpath(tt-1), bLpath(tt-1), drawCh)
          END IF
      END IF

      IF (dPath(tt) == 1) THEN ! Default
        gdpPath(tt) = hy(yPath(tt))
        cPath(tt) = gdpPath(tt)
        qSpath(tt) = xiS(yPath(tt), bSpath(tt), bLpath(tt))
        qLpath(tt) = xiL(yPath(tt), bSpath(tt), bLpath(tt))
      ELSE IF (dPath(tt) == 2) THEN ! Arrears
        gdpPath(tt) = y(yPath(tt))
        qSpath(tt) = qSa1(yPath(tt), bSpath(tt), bLpath(tt)) 
        qLpath(tt) = qLa1(yPath(tt), bSpath(tt), bLpath(tt)) 
        IF (justReturn) THEN
          cPath(tt) = gdpPath(tt) - kS * bS(returnS) - kL * bL(returnL) &
            + qSpath(tt) * ( bS(bSpath(tt)) - (1.0_wp - deltaS) * bS(returnS) ) &
            + qLpath(tt) * ( bL(bLpath(tt)) - (1.0_wp - deltaL) * bL(returnL) )
        ELSE 
          cPath(tt) = gdpPath(tt) - kS * bS(bSpath(tt-1)) - kL * bL(bLpath(tt-1)) &
            + qSpath(tt) * ( bS(bSpath(tt)) - (1.0_wp - deltaS) * bS(bSpath(tt-1)) ) &
            + qLpath(tt) * ( bL(bLpath(tt)) - (1.0_wp - deltaL) * bL(bLpath(tt-1)) )
          !IF ( bS(bSpath(tt)) + bL(bLpath(tt)) > 0.0_wp ) THEN
          !    cPath(tt) = cPath(tt) - adj * ( &
          !      bS(bSpath(tt)) / (bS(bSpath(tt)) + bL(bLpath(tt))) - bSshare )**2
          !END IF
        END IF
      ELSE
        gdpPath(tt) = y(yPath(tt))
        qSpath(tt) = qS1(yPath(tt), bSpath(tt), bLpath(tt)) 
        qLpath(tt) = qL1(yPath(tt), bSpath(tt), bLpath(tt)) 
        cPath(tt) = gdpPath(tt) - kS * bS(bSpath(tt-1)) - kL * bL(bLpath(tt-1)) &
          + qSpath(tt) * ( bS(bSpath(tt)) - (1.0_wp - deltaS) * bS(bSpath(tt-1)) ) &
          + qLpath(tt) * ( bL(bLpath(tt)) - (1.0_wp - deltaL) * bL(bLpath(tt-1)) )
        IF ( bS(bSpath(tt)) + bL(bLpath(tt)) > 0.0_wp ) THEN
            cPath(tt) = cPath(tt) - adj * ( &
              bS(bSpath(tt)) / (bS(bSpath(tt)) + bL(bLpath(tt))) - bSshare )**2
        END IF
      END IF
    END DO

    OPEN(unit=100, file=outDir // "events.tab")
    DO ix = 1,evCnt-1
      WRITE (100, "(7(I10,A))") ix, TAB, yEv(ix), TAB, bSev(ix), TAB, &
        bLev(ix), TAB, gammaSev(ix), TAB, gammaLev(ix), TAB, yDefEv(ix), TAB
    END DO
    CLOSE(100)

    OPEN(unit=100, file=outDir // "simulation.tab")
    WRITE (100, "(A9,A,3(A3,A,A10,A),A1,A,4(A10,A))") "T", TAB, "yIx", TAB, "y", TAB, "bSix", TAB, &
      "bS", TAB, "bLix", TAB, "bL", TAB, "d", TAB, "qS", TAB, "qL", TAB, "gdp", TAB, "c", TAB
    DO ix = 5000,simSz-500
      WRITE (100, "(I9,A,3(I3,A,F10.6,A),I1,A,4(F10.6,A))") ix, TAB, yPath(ix), TAB, y(yPath(ix)), TAB, &
        bSpath(ix), TAB, bS(bSpath(ix)), TAB, bLpath(ix), TAB, bL(bLpath(ix)), TAB, &
        dPath(ix), TAB, qSpath(ix), TAB, qLpath(ix), TAB, gdpPath(ix), TAB, cPath(ix), TAB
    END DO
    CLOSE(100)

!    OPEN(unit=100, file=outDir // "qSpath.txt")
!    WRITE (100, "(F9.6)") qSpath(10:simSz-1)
!    CLOSE(100)
!
!    OPEN(unit=100, file=outDir // "qLpath.txt")
!    WRITE (100, "(F9.6)") qLpath(10:simSz-1)
!    CLOSE(100)
!
!    OPEN(unit=100, file=outDir // "cPath.txt")
!    WRITE (100, "(F9.6)") cPath(10:simSz-1)
!    CLOSE(100)
!
!    OPEN(unit=100, file=outDir // "gdpPath.txt")
!    WRITE (100, "(F9.6)") gdpPath(10:simSz-1)
!    CLOSE(100)
!
!    OPEN(unit=100, file=outDir // "yPath.txt")
!    WRITE (100, "(I4)") yPath(10:simSz-1)
!    CLOSE(100)    
!
!    OPEN(unit=100, file=outDir // "bSpath.txt")
!    WRITE (100, "(I4)") bSpath(10:simSz-1)
!    CLOSE(100)    
!
!    OPEN(unit=100, file=outDir // "bLpath.txt")
!    WRITE (100, "(I4)") bLpath(10:simSz-1)
!    CLOSE(100)    
!
!    OPEN(unit=100, file=outDir // "dPath.txt")
!    WRITE (100, "(I4)") dPath(10:simSz-1)
!    CLOSE(100)    
  END SUBROUTINE simulate

END MODULE jmpMod
