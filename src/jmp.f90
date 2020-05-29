!
! Mihalache 2018
! Sovereign Default Resolution Through Maturity Extension
!
PROGRAM jmp
  USE iso_Fortran_env, ONLY: wp => real64
  USE mpi
  USE Param
  USE MPIMod
  USE NL
  USE jmpMod
  IMPLICIT NONE

  INTEGER, DIMENSION(8) :: timeVals
  INTEGER :: ix, yIx, bSIx, bLIx, ggSIx, ggLIx, yPrIx, bSprIx, bLprIx
  INTEGER :: iter, okFlag, ii, jj, ffS, ffL
  REAL(wp) :: errQ, errV, errVd, errQS, errQL, errValues, ggPr, ffPr
  TYPE(stateSp) :: st

  REAL(wp) :: netRes, cons, vHere, maxTmp, vdTmp, Ib, divd, maxTmp2
  REAL(wp) :: qStmp, qLtmp, dd, ellS, ellL, bSum
  REAL(wp), DIMENSION(bSSz, bLSz) :: contVal, dSov, dCre, nash
  REAL(wp), DIMENSION(chPrSz) :: chTmp, dCreTmp
  INTEGER, DIMENSION(chPrSz) :: chStmp, chLtmp
  INTEGER, DIMENSION(3) :: locQS, locQL, locV
  INTEGER, DIMENSION(1) :: locVd
  
  !
  ! Start MPI
  !
  CALL mpiSetup(spaceSz, chPrSz, okFlag)
  IF (okFlag == 0) THEN
    IF (workerId == 0) WRITE (*, "(A,I10,A,I3)") "MPI okFlag == 0. ", spaceSz, " ", workerNo
    STOP
  ELSE
    IF (workerId == 0) WRITE (*, *) "Start..."
  END IF

  ALLOCATE(chPr(ySz, bSSz, bLSz, chPrSz))
  ALLOCATE(chS(ySz, bSSz, bLSz, chPrSz))
  ALLOCATE(chL(ySz, bSSz, bLSz, chPrSz))
  ALLOCATE(chPrA(ySz, bSSz, bLSz, chPrSz))
  ALLOCATE(chSa(ySz, bSSz, bLSz, chPrSz))
  ALLOCATE(chLa(ySz, bSSz, bLSz, chPrSz))

  ! Shocks and grids
  CALL loadShocks()
  CALL generateGrids()
  CALL findVaut()

  IF (loadGuess) THEN
    CALL loadEquilibrium()
  ELSE
    ! Initialisation
    qS0 = 1.0_wp
    qS1 = 1.0_wp
    qL0 = 1.0_wp
    qL1 = 1.0_wp
    qSa0 = 1.0_wp
    qSa1 = 1.0_wp
    qLa0 = 1.0_wp
    qLa1 = 1.0_wp
    xiS = 1.0_wp
    xiL = 1.0_wp
    gammaS = 1
    gammaL = 1
    gammaPr = 1.0_wp / (0.0_wp + chPrSz)
    defPr = 0.0_wp
    chPr = 1.0_wp / (0.0_wp + chPrSz)
    chS = 1
    chL = 1
    chPrA = 1.0_wp / (0.0_wp + chPrSz)
    chSa = 1
    chLa = 1
    Vd0 = Vaut
    Vd1 = Vaut
    FORALL (bSIx = 1:bSSz, bLIx = 1:bLSz)
      V0(:, bSIx, bLIx) = 0.0_wp ! Vaut
      V1(:, bSIx, bLIx)  = 0.0_wp ! Vaut
      Va0(:, bSIx, bLIx) = 0.0_wp ! Vaut
      Va1(:, bSIx, bLIx)  = 0.0_wp ! Vaut
    END FORALL
  END IF

  IF (computeVFI) THEN
    errV = 1.0_wp
    errQ = 1.0_wp
    iter = 1
    DO WHILE ((errValues > epsV .OR. errQ > epsQ) .AND. iter <= maxIter)

      IF (workerId == 0 .AND. verbose) WRITE (*, *) "   V..."
      !
      ! Update V
      !
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx,contVal,bSprIx,bLprIx, &
      !$OMP& netRes,chTmp,chStmp,chLtmp,ii,cons,vHere,maxTmp,ellS,ellL,bSum, &
      !$OMP& Ib,maxTmp2,divd)
      DO ix = 1,chunk
        st = toSS(workerId * chunk + ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        netRes = y(yIx) - kS * bS(bSIx) - kL * bL(bLIx)

        FORALL (bSprIx = 1:bSSz, bLprIx = 1:bLSz)
          contVal(bSprIx, bLprIx) = DOT_PRODUCT( exoTran(yIx, :), V0(:, bSprIx, bLprIx) )
        END FORALL
        contVal = beta * contVal

        chTmp = veryNeg
        chStmp = 1
        chLtmp = 1

        DO bSprIx = 1,bSSz
          DO bLprIx = 1,bLSz
            ellS = bS(bSprIx) - (1.0_wp - deltaS) * bS(bSIx)
            ellL = bL(bLprIx) - (1.0_wp - deltaL) * bL(bLIx)
            cons = netRes + qS0(yIx, bSprIx, bLprIx) * ellS + qL0(yIx, bSprIx, bLprIx) * ellL
            IF ( bS(bSprIx) + bL(bLprIx) > 0.0_wp ) THEN
              cons = cons - adj * ( bS(bSprIx) / (bS(bSprIx) + bL(bLprIx)) - bSshare )**2
            END IF
            IF (cons <= 0.0_wp .OR. &
              ! DOT_PRODUCT( exoTran(yIx, :), defPr(:, bSprIx, bLprIx) ) > maxDefPr) CYCLE
              ! qL0(yIx, bSprIx, bLprIx) < minq .OR. qS0(yIx, bSprIx, bLprIx) < minq) CYCLE
              qL0(yIx, bSprIx, bLprIx) < minq) CYCLE
            vHere = (1.0_wp - beta) * uu(cons) + contVal(bSprIx, bLprIx)

            ! No better than the current worst, discard
            IF (vHere < chTmp(chPrSz)) CYCLE

            DO ii = chPrSz-1,1,-1
              IF (vHere >= chTmp(ii+1) .AND. vHere < chTmp(ii)) THEN
                ! Insert at position ii+1 and
                ! Move down positions ii+1 and higher
                IF(ii+2 <= chPrSz) THEN
                  chTmp(ii+2:chPrSz) = chTmp(ii+1:chPrSz-1)
                  chStmp(ii+2:chPrSz) = chStmp(ii+1:chPrSz-1)
                  chLtmp(ii+2:chPrSz) = chLtmp(ii+1:chPrSz-1)
                END IF
                chTmp(ii+1) = vHere
                chStmp(ii+1) = bSprIx
                chLtmp(ii+1) = bLprIx
                EXIT
              END IF
            END DO

            IF (vHere > chTmp(1)) THEN
              ! Best one so far
              ! Move everything else down
              chTmp(2:chPrSz) = chTmp(1:chPrSz-1)
              chStmp(2:chPrSz) = chStmp(1:chPrSz-1)
              chLtmp(2:chPrSz) = chLtmp(1:chPrSz-1)
              chTmp(1) = vHere
              chStmp(1) = bSprIx
              chLtmp(1) = bLprIx
            END IF

          END DO ! bLprIx
        END DO ! bSprIx

        maxTmp = MAXVAL(chTmp)
        bSum = SUM(EXP((chTmp - maxTmp) / rho))
        Ib = maxTmp + rho * LOG( bSum / (0.0_wp + chPrSz) )

        maxTmp2 = MAX(Vd0(yIx), Ib)
        divd = EXP((Ib - maxTmp2) / rhoDef) + EXP((Vd0(yIx) - maxTmp2) / rhoDef)

        buffdpr(ix) = EXP((Vd0(yIx) - maxTmp2) / rhoDef) / divd
        buffchS( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chStmp
        buffchL( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chLtmp
        buffchpr( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = (1.0_wp - buffdpr(ix)) & 
          * EXP((chTmp - maxTmp) / rho) / bSum

        buffV(ix) = maxTmp2 + rhoDef * LOG(divd / 2.0_wp)

!        maxTmp = MAX( Vd0(yIx), MAXVAL(chTmp) )
!        bSum = SUM(EXP((chTmp - maxTmp) / rho)) + EXP((Vd0(yIx) - maxTmp) / rho)

!        buffdpr(ix) = EXP((Vd0(yIx) - maxTmp) / rho) / bSum
!        buffchS( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chStmp
!        buffchL( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chLtmp
!        buffchpr( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = EXP((chTmp - maxTmp) / rho) / bSum

!        buffV(ix) = maxTmp + rho * LOG(bSum / (chPrSz + 1.0_wp))
      END DO ! ix
      !$OMP END PARALLEL DO

      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      CALL MPI_ALLGATHER(buffV, chunk, MPI_DOUBLE_PRECISION, bufflV, chunk, &
        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      CALL MPI_ALLGATHER(buffdpr, chunk, MPI_DOUBLE_PRECISION, buffldpr, chunk, &
        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)

      CALL MPI_ALLGATHER(buffchpr, chunkCh, MPI_DOUBLE_PRECISION, bufflchpr, chunkCh, &
        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)

      CALL MPI_ALLGATHER(buffchS, chunkCh, MPI_INTEGER, buffLchS, chunkCh, &
        MPI_INTEGER, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      CALL MPI_ALLGATHER(buffchL, chunkCh, MPI_INTEGER, bufflchL, chunkCh, &
        MPI_INTEGER, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      ! From buffers to arrays
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx)
      DO ix = 1,spaceSz
        st = toSS(ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        V1(yIx, bSIx, bLIx) = bufflV(ix)
        defPr(yIx, bSIx, bLIx) = buffldpr(ix)
        chPr(yIx, bSIx, bLIx, :) = bufflchpr( ((ix-1)*chPrSz+1):(ix*chPrSz) )
        chS(yIx, bSIx, bLIx, :) = bufflchS( ((ix-1)*chPrSz+1):(ix*chPrSz) )
        chL(yIx, bSIx, bLIx, :) = bufflchL( ((ix-1)*chPrSz+1):(ix*chPrSz) )
      END DO
      !$OMP END PARALLEL DO


      IF (workerId == 0 .AND. verbose) WRITE (*, *) "   Va..."
      !
      ! Update V arrears
      !
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx,contVal,bSprIx,bLprIx, &
      !$OMP& netRes,chTmp,chStmp,chLtmp,ii,cons,vHere,maxTmp,ellS,ellL,bSum)
      DO ix = 1,chunk
        st = toSS(workerId * chunk + ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        netRes = y(yIx) - kS * bS(bSIx) - kL * bL(bLIx)

        FORALL (bSprIx = 1:bSSz, bLprIx = 1:bLSz)
          contVal(bSprIx, bLprIx) = DOT_PRODUCT( exoTran(yIx, :), &
            etaA * V0(:, bSprIx, bLprIx) + (1.0_wp - etaA) * Va0(:, bSprIx, bLprIx) )
        END FORALL
        contVal = beta * contVal

        chTmp = veryNeg
        chStmp = 1
        chLtmp = 1

        DO bSprIx = 1,bSSz
          DO bLprIx = 1,bLSz
            ellS = bS(bSprIx) - (1.0_wp - deltaS) * bS(bSIx)
            ellL = bL(bLprIx) - (1.0_wp - deltaL) * bL(bLIx)
            cons = netRes + qSa0(yIx, bSprIx, bLprIx) * ellS + qLa0(yIx, bSprIx, bLprIx) * ellL
            !IF ( bS(bSprIx) + bL(bLprIx) > 0.0_wp ) THEN
            !  cons = cons - adj * ( bS(bSprIx) / (bS(bSprIx) + bL(bLprIx)) - bSshare )**2
            !END IF
            IF (ellS > 0.0_wp .OR. ellL > 0.0_wp .OR. cons <= 0.0_wp) CYCLE
            vHere = (1.0_wp - beta) * uu(cons) + contVal(bSprIx, bLprIx)

            ! No better than the current worst, discard
            IF (vHere < chTmp(chPrSz)) CYCLE

            DO ii = chPrSz-1,1,-1
              IF (vHere >= chTmp(ii+1) .AND. vHere < chTmp(ii)) THEN
                ! Insert at position ii+1 and
                ! Move down positions ii+1 and higher
                IF(ii+2 <= chPrSz) THEN
                  chTmp(ii+2:chPrSz) = chTmp(ii+1:chPrSz-1)
                  chStmp(ii+2:chPrSz) = chStmp(ii+1:chPrSz-1)
                  chLtmp(ii+2:chPrSz) = chLtmp(ii+1:chPrSz-1)
                END IF
                chTmp(ii+1) = vHere
                chStmp(ii+1) = bSprIx
                chLtmp(ii+1) = bLprIx
                EXIT
              END IF
            END DO

            IF (vHere > chTmp(1)) THEN
              ! Best one so far
              ! Move everything else down
              chTmp(2:chPrSz) = chTmp(1:chPrSz-1)
              chStmp(2:chPrSz) = chStmp(1:chPrSz-1)
              chLtmp(2:chPrSz) = chLtmp(1:chPrSz-1)
              chTmp(1) = vHere
              chStmp(1) = bSprIx
              chLtmp(1) = bLprIx
            END IF

          END DO ! bLprIx
        END DO ! bSprIx

        maxTmp = MAXVAL(chTmp)
        bSum = SUM(EXP((chTmp - maxTmp) / rho))

        buffchS( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chStmp
        buffchL( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = chLtmp
        buffchpr( ((ix-1)*chPrSz+1):(ix*chPrSz) ) = EXP((chTmp - maxTmp) / rho) / bSum

        buffV(ix) = maxTmp + rho * LOG(bSum / (chPrSz + 0.0_wp))
      END DO ! ix
      !$OMP END PARALLEL DO

      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      CALL MPI_ALLGATHER(buffV, chunk, MPI_DOUBLE_PRECISION, bufflV, chunk, &
        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)

      CALL MPI_ALLGATHER(buffchpr, chunkCh, MPI_DOUBLE_PRECISION, bufflchpr, chunkCh, &
        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)

      CALL MPI_ALLGATHER(buffchS, chunkCh, MPI_INTEGER, buffLchS, chunkCh, &
        MPI_INTEGER, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      CALL MPI_ALLGATHER(buffchL, chunkCh, MPI_INTEGER, bufflchL, chunkCh, &
        MPI_INTEGER, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      ! From buffers to arrays
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx)
      DO ix = 1,spaceSz
        st = toSS(ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        Va1(yIx, bSIx, bLIx) = bufflV(ix)
        chPrA(yIx, bSIx, bLIx, :) = bufflchpr( ((ix-1)*chPrSz+1):(ix*chPrSz) )
        chSa(yIx, bSIx, bLIx, :) = bufflchS( ((ix-1)*chPrSz+1):(ix*chPrSz) )
        chLa(yIx, bSIx, bLIx, :) = bufflchL( ((ix-1)*chPrSz+1):(ix*chPrSz) )
      END DO
      !$OMP END PARALLEL DO

      IF (workerId == 0 .AND. verbose) WRITE (*, *) "   q..."
      !
      ! Update q
      !
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,yPrIx,bSprIx,bLprIx,dd,ii,qStmp,qLtmp,ffS,ffL,ffPr)
      DO ix = 1,chunk
        st = toSS(workerId * chunk + ix)
        yIx = st%yIx
        bSprIx = st%bSIx
        bLprIx = st%bLIx
        
        qStmp = 0.0_wp
        qLtmp = 0.0_wp
        DO yPrIx = 1,ySz
          dd = defPr(yPrIx, bSprIx, bLprIx)
          qStmp = qStmp + exoTran(yIx, yPrIx) * (dd * xiS(yPrIx, bSprIx, bLprIx) + (1.0_wp - dd) * kS)
          qLtmp = qLtmp + exoTran(yIx, yPrIx) * (dd * xiL(yPrIx, bSprIx, bLprIx) + (1.0_wp - dd) * kL)
          DO ii = 1,chPrSz
            ffS = chS(yPrIx, bSprIx, bLprIx, ii)
            ffL = chL(yPrIx, bSprIx, bLprIx, ii)
            ffPr = chPr(yPrIx, bSprIx, bLprIx, ii)
            qStmp = qStmp + exoTran(yIx, yPrIx) * (1.0_wp - deltaS) * qS0(yPrIx, ffS, ffL) * ffPr
            qLtmp = qLtmp + exoTran(yIx, yPrIx) * (1.0_wp - deltaL) * qL0(yPrIx, ffS, ffL) * ffPr
          END DO
        END DO
        
        buffQS(ix) = qStmp / (1.0_wp + rf)
        buffQL(ix) = qLtmp / (1.0_wp + rf)
      END DO
      !$OMP END PARALLEL DO

      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      CALL MPI_ALLGATHER(buffQS, chunk, MPI_DOUBLE_PRECISION, bufflQS, &
        chunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      CALL MPI_ALLGATHER(buffQL, chunk, MPI_DOUBLE_PRECISION, bufflQL, &
        chunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      ! From buffers to arrays
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx)
      DO ix = 1,spaceSz
        st = toSS(ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        qS1(yIx, bSIx, bLIx) = bufflQS(ix)
        qL1(yIx, bSIx, bLIx) = bufflQL(ix)
      END DO
      !$OMP END PARALLEL DO

      IF (workerId == 0 .AND. verbose) WRITE (*, *) "   qa..."
      !
      ! Update q arrears
      !
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,yPrIx,bSprIx,bLprIx,dd,ii,qStmp,qLtmp,ffS,ffL,ffPr)
      DO ix = 1,chunk
        st = toSS(workerId * chunk + ix)
        yIx = st%yIx
        bSprIx = st%bSIx
        bLprIx = st%bLIx
        
        qStmp = 0.0_wp
        qLtmp = 0.0_wp
        DO yPrIx = 1,ySz
          DO ii = 1,chPrSz
            ffS = chSa(yPrIx, bSprIx, bLprIx, ii)
            ffL = chLa(yPrIx, bSprIx, bLprIx, ii)
            ffPr = chPrA(yPrIx, bSprIx, bLprIx, ii)
            qStmp = qStmp + exoTran(yIx, yPrIx) * qSa0(yPrIx, ffS, ffL) * ffPr
            qLtmp = qLtmp + exoTran(yIx, yPrIx) * qLa0(yPrIx, ffS, ffL) * ffPr
          END DO
        END DO
        
        buffQS(ix) = etaA * qS0(yIx, bSprIx, bLprIx) &
          + (1.0_wp - etaA) * (kS + (1.0_wp - deltaS) * qStmp) / (1.0_wp + rf)
        buffQL(ix) = etaA * qL0(yIx, bSprIx, bLprIx) &
          + (1.0_wp - etaA) * (kL + (1.0_wp - deltaL) * qLtmp) / (1.0_wp + rf)
      END DO
      !$OMP END PARALLEL DO

      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      CALL MPI_ALLGATHER(buffQS, chunk, MPI_DOUBLE_PRECISION, bufflQS, &
        chunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      CALL MPI_ALLGATHER(buffQL, chunk, MPI_DOUBLE_PRECISION, bufflQL, &
        chunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpiErr)
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      
      ! From buffers to arrays
      !$OMP PARALLEL DO PRIVATE(ix,st,yIx,bSIx,bLIx)
      DO ix = 1,spaceSz
        st = toSS(ix)
        yIx = st%yIx
        bSIx = st%bSIx
        bLIx = st%bLIx

        qSa1(yIx, bSIx, bLIx) = bufflQS(ix)
        qLa1(yIx, bSIx, bLIx) = bufflQL(ix)
      END DO
      !$OMP END PARALLEL DO


      IF (recovery) THEN
        IF (workerId == 0 .AND. verbose) WRITE (*, *) "   Nash, Xi..."
      !
      ! Nash and Xi-s
      !
      !$OMP PARALLEL DO PRIVATE(yIx,bSIx,bLIx,ggSIx,ggLIx,dSov,dCre,dCreTmp,&
      !$OMP& nash,chTmp,chStmp,chLtmp,ii,vHere,maxTmp,bSum,qStmp,qLtmp,jj,ffS,ffL,ffPr)
      DO yIx = 1,ySz
        DO ggSIx = 1,bSSz
        DO ggLIx = 1,bLSz
          dSov(ggSIx, ggLIx) = Va1(yIx, ggSIx, ggLIx) - Vaut(yIx)

          qStmp = 0.0_wp
          qLtmp = 0.0_wp
          DO jj = 1,chPrSz
            ffS = chSa(yIx, ggSIx, ggLIx, jj)
            ffL = chLa(yIx, ggSIx, ggLIx, jj)
            ffPr = chPrA(yIx, ggSIx, ggLIx, jj)
            qStmp = qStmp + qSa1(yIx, ffS, ffL) * ffPr
            qLtmp = qLtmp + qLa1(yIx, ffS, ffL) * ffPr
          END DO

          dCre(ggSIx, ggLIx) = (kS + (1.0_wp - deltaS) * qStmp) * bS(ggSIx) &
            + (kL + (1.0_wp - deltaL) * qLtmp) * bL(ggLIx)
        END DO
        END DO
        WHERE (dSov > 0.0_wp .AND. dCre > 0.0_wp)
          nash = alpha * LOG(dSov) + (1.0_wp - alpha) * LOG(dCre)
        ELSEWHERE
          nash = veryNeg
        END WHERE

        chTmp = veryNeg
        dCreTmp = veryNeg
        chStmp = 1
        chLtmp = 1 
        DO ggSIx = 1,bSSz
        DO ggLIx = 1,bLSz
          vHere = nash(ggSIx, ggLIx)
          IF (vHere < veryNeg + 0.01_wp) CYCLE
          ! No better than the current worst, discard
          IF (vHere < chTmp(chPrSz)) CYCLE

          DO ii = chPrSz-1,1,-1
            IF (vHere >= chTmp(ii+1) .AND. vHere < chTmp(ii)) THEN
              ! Insert at position ii+1 and
              ! Move down positions ii+1 and higher
              IF(ii+2 <= chPrSz) THEN
                chTmp(ii+2:chPrSz) = chTmp(ii+1:chPrSz-1)
                dCreTmp(ii+2:chPrSz) = dCreTmp(ii+1:chPrSz-1)
                chStmp(ii+2:chPrSz) = chStmp(ii+1:chPrSz-1)
                chLtmp(ii+2:chPrSz) = chLtmp(ii+1:chPrSz-1)
              END IF
              chTmp(ii+1) = vHere
              dCreTmp(ii+1) = dCre(ggSIx, ggLIx)
              chStmp(ii+1) = ggSIx
              chLtmp(ii+1) = ggLIx
              EXIT
            END IF
          END DO

          IF (vHere > chTmp(1)) THEN
            ! Best one so far
            ! Move everything else down
            chTmp(2:chPrSz) = chTmp(1:chPrSz-1)
            dCreTmp(2:chPrSz) = dCreTmp(1:chPrSz-1)
            chStmp(2:chPrSz) = chStmp(1:chPrSz-1)
            chLtmp(2:chPrSz) = chLtmp(1:chPrSz-1)
            chTmp(1) = vHere
            dCreTmp(1) = dCre(ggSIx, ggLIx)
            chStmp(1) = ggSIx
            chLtmp(1) = ggLIx
          END IF
        END DO
        END DO

        maxTmp = MAXVAL(chTmp)
        bSum = SUM(EXP((chTmp - maxTmp) / rhoNash))

        gammaS(yIx, :) = chStmp
        gammaL(yIx, :) = chLtmp
        gammaPr(yIx, :) = EXP((chTmp - maxTmp) / rhoNash) / bSum
        EdCre(yIx) = DOT_PRODUCT(gammaPr(yIx, :), dCreTmp)
      END DO
      !$OMP END PARALLEL DO 

      !
      ! deltaCreditors
      !
      !DO yIx = 1,ySz
      !  dd = 0.0_wp
      !  !$OMP PARALLEL DO PRIVATE(ii,jj,ggSIx,ggLIx,ggPr,ffS,ffL,ffPr) REDUCTION(+:dd)
      !  DO ii = 1,chPrSz
      !    ggPr = gammaPr(yIx, ii)
      !    IF (ggPr > 1.0D-8) THEN
      !      ggSIx = gammaS(yIx, ii)
      !      ggLIx = gammaL(yIx, ii)
      !      dd = dd + ( kS * bS(ggSIx) + kL * bL(ggLIx) ) * ggPr
      !      DO jj = 1,chPrSz
      !        ffS = chSa(yIx, ggSIx, ggLIx, jj)
      !        ffL = chLa(yIx, ggSIx, ggLIx, jj)
      !        ffPr = chPrA(yIx, ggSIx, ggLIx, jj)
      !        dd = dd + (1.0_wp - deltaS) * qSa1(yIx, ffS, ffL) * ffPr * bS(ggSIx) * ggPr
      !        dd = dd + (1.0_wp - deltaL) * qLa1(yIx, ffS, ffL) * ffPr * bL(ggLIx) * ggPr
      !      END DO
      !    ELSE
      !      dd = dd
      !    END IF
      !  END DO
      !  !$OMP END PARALLEL DO 
      !  EdCre(yIx) = dd
      !END DO
      
      xiSkeep = xiS
      xiLkeep = xiL
      
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(bLIx,bSIx,yIx,divd)
      DO bLIx = 1,bLSz
        DO bSIx = 1,bSSz
          IF (bLIx == 1 .AND. bSIx == 1) THEN
            xiS(:, bSIx, bLIx) = 1.0_wp
            xiL(:, bSIx, bLIx) = 1.0_wp
            CYCLE
          END IF
          
         divd = muS * bS(bSIx) + bL(bLIx)
          
          DO yIx = 1,ySz
            IF (bSIx > 1) THEN
              xiS(yIx, bSIx, bLIx) = ((1.0_wp - eta) * DOT_PRODUCT(exoTran(yIx, :), xiSkeep(:, bSIx, bLIx)) &
                + eta * muS * DOT_PRODUCT(exoTran(yIx, :), EdCre) / divd) / (1.0_wp + rf)
            ELSE
              xiS(yIx, bSIx, bLIx) = 0.0_wp
            END IF
            
            IF (bLIx > 1) THEN
              xiL(yIx, bSIx, bLIx) = ((1.0_wp - eta) * DOT_PRODUCT(exoTran(yIx, :), xiLkeep(:, bSIx, bLIx)) &
                + eta * DOT_PRODUCT(exoTran(yIx, :), EdCre) / divd ) / (1.0_wp + rf)
            ELSE
              xiL(yIx, bSIx, bLIx) = 0.0_wp
            END IF
          END DO

        END DO
      END DO
      !$OMP END PARALLEL DO 

      ELSE
        xiS = 0.0_wp
        xiL = 0.0_wp
        gammaS = 1
        gammaL = 1
        gammaPr = 0.0_wp
        gammaPr(:, 1) = 1.0_wp
        EdCre = 0.0_wp
      END IF ! Recovery
      
      DO yIx = 1,ySz
        vdTmp = (1.0_wp - beta) * uhy(yIx) + beta * (1.0_wp - eta) * DOT_PRODUCT(exoTran(yIx, :), Vd0)
        dd = 0.0_wp
        !$OMP PARALLEL DO PRIVATE(yPrIx,ggSIx,ggLIx,ii) REDUCTION(+:dd) COLLAPSE(2)
        DO yPrIx = 1,ySz
          DO ii = 1,chPrSz
            ggSIx = gammaS(yPrIx, ii)
            ggLIx = gammaL(yPrIx, ii)
            dd = dd + exoTran(yIx, yPrIx) * gammaPr(yPrIx, ii) * Va1(yPrIx, ggSIx, ggLIx)
          END DO
        END DO 
        !$OMP END PARALLEL DO
        Vd1(yIx) = vdTmp + beta * eta * dd
      END DO


      IF (.FALSE. .AND. workerId == 0) THEN
        WRITE (*, *) MAXVAL(xiS), MAXVAL(xiL), MAXVAL(qS0), MAXVAL(qL0), MAXVAL(qSa0), MAXVAL(qLa0)
        WRITE (*, *) MINVAL(xiS), MINVAL(xiL), MINVAL(qS0), MINVAL(qL0), MINVAL(qSa0), MINVAL(qLa0)
        WRITE (*, *) MINVAL(EdCre), MAXVAL(EdCre)
      END IF

      !
      ! Compute errors; Weighted average q
      !
      !$OMP PARALLEL SECTIONS
      !$OMP SECTION
      errQS = MAX( MAXVAL(ABS(qS1 - qS0)), MAXVAL(ABS(qSa1 - qSa0)) )
      locQS = MAXLOC(ABS(qS1 - qS0))
      qS0 = (1.0_wp - updateQ) * qS0 + updateQ * qS1
      qSa0 = (1.0_wp - updateQ) * qSa0 + updateQ * qSa1

      !$OMP SECTION
      errQL = MAX( MAXVAL(ABS(qL1 - qL0)), MAXVAL(ABS(qLa1 - qLa0)) )
      locQL = MAXLOC(ABS(qL1 - qL0))
      qL0 = (1.0_wp - updateQ) * qL0 + updateQ * qL1
      qLa0 = (1.0_wp - updateQ) * qLa0 + updateQ * qLa1

      !$OMP SECTION
      errV = MAX( MAXVAL(ABS(V1 - V0)), MAXVAL(ABS(Va1 - Va0)) )
      locV = MAXLOC(ABS(V1 - V0))
      V0 = V1
      Va0 = Va1
      
      !$OMP SECTION
      errVd = MAXVAL(ABS(Vd1 - Vd0))
      locVd = MAXLOC(ABS(Vd1 - Vd0))
      Vd0 = Vd1
      !$OMP END PARALLEL SECTIONS
      
      errQ = MAX( errQS, errQL )
      errValues = MAX(errV, errVd)

      ! Iterate
      IF (workerId == 0) THEN 
        WRITE (*, "(I6,A,ES10.3,A,ES10.3)") iter, " errValues = ", errValues, ", errQ = ", errQ
        IF (.TRUE.) WRITE (*, "(A,ES10.3,A,3I5)") "   qS: ", errQS, " at ", locQS
        IF (.TRUE.) WRITE (*, "(A,ES10.3,A,3I5)") "   qL: ", errQL, " at ", locQL
        IF (.TRUE.) WRITE (*, "(A,ES10.3,A,3I5)") "    V: ", errV, " at ", locV
        IF (.TRUE.) WRITE (*, "(A,ES10.3,A,3I5)") "   Vd: ", errVd, " at ", locVd
        CALL timestamp()
      END IF
      iter = iter + 1
    END DO
  ELSE IF (workerId == 0 .AND. .NOT. loadGuess) THEN
    WRITE (*, *) "Not computing, loading results from .txt-s"
    CALL loadEquilibrium()
  END IF

  CALL mpiFinalize()

  IF (workerId == 0) THEN
    IF (computeVFI) THEN
      WRITE (*, *) "Done. Saving equilibrium..."
      CALL saveEquilibrium()
    END IF

    WRITE (*, *) "Simulate..."
    CALL simulate()
    WRITE (*, *) "End."
  END IF

CONTAINS

END PROGRAM jmp
