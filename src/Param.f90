MODULE Param
  USE iso_Fortran_env, ONLY: wp => real64
  IMPLICIT NONE

  ! Calibration
  REAL(wp), PARAMETER :: crra = 2.0_wp ! CRRA
  REAL(wp), PARAMETER :: beta = 0.94_wp ! Discounting
  REAL(wp), PARAMETER :: rf = 0.032_wp ! Risk-free rate
  REAL(wp), PARAMETER :: deltaS = (1.0_wp + rf) / 1.0_wp - rf ! Maturity short
  REAL(wp), PARAMETER :: deltaL = (1.0_wp + rf) / 10.0_wp - rf ! Maturity long
  REAL(wp), PARAMETER :: alpha = 0.945_wp ! Sov bargain power
  REAL(wp), PARAMETER :: eta = 0.33_wp ! Pr negotiations
  REAL(wp), PARAMETER :: etaA = 0.25_wp ! Pr return to market
  REAL(wp), PARAMETER :: muS = 0.5_wp ! Assignment
  REAL(wp), PARAMETER :: lbd0 = -0.85_wp, lbd1 = 1.0_wp
  REAL(wp), PARAMETER :: adj = 0.02_wp, bSshare = 0.33_wp

  REAL(wp), PARAMETER :: kS = deltaS + rf, kL = deltaL + rf
  
  ! State space
  INTEGER, PARAMETER :: ySz = 21
  INTEGER, PARAMETER :: bSSz = 74, bLSz = 74, bSsplit = 37, bLsplit = 37
  REAL(wp), PARAMETER :: bSmin = 0.0, bSmid = 0.3, bSmax = 0.6
  REAL(wp), PARAMETER :: bLmin = 0.0, bLmid = 0.3, bLmax = 0.6

  INTEGER, PARAMETER :: spaceSz = ySz * bSSz * bLSz

  ! Choice Prs Limit
  INTEGER, PARAMETER :: chPrSz = 300

  ! Numerical
  REAL(wp), PARAMETER :: veryNeg = -1.0D+4, updateQ = 1.0_wp
  REAL(wp), PARAMETER :: epsV = 1.0D-6
  REAL(wp), PARAMETER :: epsQ = 1.0D-4
  ! REAL(wp), PARAMETER :: maxDefPr = 0.15_wp
  REAL(wp), PARAMETER :: minq = 0.7_wp

  REAL(wp), PARAMETER :: rho = 1.0D-5
  REAL(wp), PARAMETER :: rhoDef = rho
  REAL(wp), PARAMETER :: rhoNash = rho

  INTEGER, PARAMETER :: maxIter = 500

  LOGICAL, PARAMETER :: loadGuess = .FALSE.
  LOGICAL, PARAMETER :: computeVFI = .TRUE.
  LOGICAL, PARAMETER :: recovery = .TRUE.
  LOGICAL, PARAMETER :: verbose = .FALSE.
  
  ! Simulation
  INTEGER, PARAMETER :: simSz = 150002 ! > 5000 only
END MODULE Param
