C ---------------------------------------------------------------------- 
C  det: Calculate determinant of 3x3 matrix
C  Input: A (3x3 matrix)
C ---------------------------------------------------------------------- 
      function eq(STRESS)
          implicit none
          REAL eq,STRESS
          dimension STRESS(6)
          eq=sqrt(STRESS(1)*STRESS(1)
     .      +STRESS(2)*STRESS(2)
     .      +STRESS(3)*STRESS(3)
     .      -STRESS(1)*STRESS(2)
     .      -STRESS(2)*STRESS(3)
     .      -STRESS(3)*STRESS(1)
     .      +3.0*(STRESS(4)*STRESS(4)
     .      +STRESS(5)*STRESS(5)
     .      +STRESS(6)*STRESS(6)))
      end


C ----------------------------------------------------------------------    
C  matinv: Invert matrix
C  Input: n1    : Lower dimension of square matrix
C         n2    : Upper dimension of square matrix
C          A    : Square matrix of dimension n1:n2 x n1:n2 to be inverted
C  Output: Ainv : Inverted matrix
C
C Taken from http://wp.me/p61TQ-zb
C Based on a lecture by Prof. McFarland
C http://math.uww.edu/~mcfarlat/inverse.htm
C ----------------------------------------------------------------------
      subroutine matinv(n1,n2,A,AINV)
      implicit none
      integer i,j,k,n1,n2
      real*8 A(n1:n2,n1:n2),AINV(n1:n2,n1:n2),B(n1:n2,n1:2*n2-n1+1)
     1,pivot,xnum

! Initialization
      do i=n1,n2
      do j=n1,n2
      AINV(i,j)=0.0D0
      end do
      end do

! Make augmented matrix
      do i=n1,n2
          do j=n1,n2
              B(i,j)=0.0D0
              B(i,j+n2-n1+1)=0.0D0

              B(i,j)=A(i,j)
                  if(i.EQ.j) then
                      B(i,j+n2-n1+1)=1.0D0
                  end if
        end do
      end do

      do i=n1,n2

! Choose the leftmost non-zero element as pivot
        do j=n1,n2
            if(dabs(B(i,j)).gt.0)then
                pivot=B(i,j)
                exit
            end if
        end do

! Step 1: Change the chosen pivot into "1" by dividing
! the pivot's row by the pivot number
          do j=n1,2*n2-n1+1
              B(i,j)=B(i,j)/pivot
          end do
          pivot=B(i,i) !Update pivot value

! Step 2: Change the remainder of the pivot's column into 0's
! by adding to each row a suitable multiple of the pivot row
        do k=n1,n2 !row
              if(k.ne.i) then
                  xnum=B(k,i)/pivot !Same column with current pivot
                  do j=n1,2*n2-n1+1 !column
                      B(k,j)=B(k,j)-xnum*B(i,j)
                  end do
              end if
        end do

      end do

! Prepare the final inverted matrix
      do i=n1,n2
          do j=n1,n2
              AINV(i,j)=B(i,j+n2-n1+1)
          end do
      end do

      return
      end







       SUBROUTINE VUMAT(
! READ ONLY - DO NOT MODIFY
     . NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME,
     . TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY,
     . STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     . STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW,
     . STRETCHNEW, DEFGRADNEW, FIELDNEW,
! WRITE ONLY - DO NOT READ
     . STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
!-----------------------------------------------------------------------
!     ABAQUS implicit variable declaration included in VABA_PARAM.INC
!     states the following:
!     a to h are real variables
!     o to z are real variables
!     i to n are integer variables
!-----------------------------------------------------------------------
      INCLUDE 'VABA_PARAM.INC'
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     ABAQUS variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     . CHARLENGTH(*), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(*),
     . TEMPOLD(*), STRETCHOLD(*), DEFGRADOLD(*), FIELDOLD(*),
     . STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV),
     . ENERINTERNOLD(NBLOCK),  ENERINELASOLD(NBLOCK), TEMPNEW(*),
     . STRETCHNEW(*), DEFGRADNEW(*), FIELDNEW(*),
     . STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     . ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C
      CHARACTER*80 CMNAME
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION STRESS(NDIR+NSHR) ! Stress tensor inside UMAT
      DIMENSION DFDS(NDIR+NSHR)   ! Derivative of the yield function
      DIMENSION STRESSK(NDIR+NSHR)

      DIMENSION EP(NDIR+NSHR)
      DIMENSION DEP(NDIR+NSHR)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables could be declared as follow
!     but it is not required due to the VABA_PARAM.INC file
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!   Material parameters
      REAL YOUNG       ! Young's modulus 
      REAL POISS       ! Poisson's ratio
      REAL C11,C12,C44 ! Elasticity matrix C components
      REAL SIGMA0      ! Initial yield stress
      REAL ET          ! Tangent modulus
!   Internal variables
      REAL P           ! Equivalent plastic strain
!   Plasticity variables
      REAL DLAMBDA     ! Plastic multiplier
      REAL DDLAMBDA    ! Increment in plastic multiplier
!   Yield function variables
      REAL F           ! Yield function
      REAL PHI         ! Equivalent stress
      REAL SIGMAY      ! Yield stress
!   Computational variables
      REAL RESNOR      ! Convergence criterion
      REAL TOL         ! Tolerance for the RMAP algorithm
      INTEGER MXITER   ! Maximum number of iteration for the RMAP
      INTEGER ITER     ! Number of iteration for the RMAP
!   Added by newbie
      REAL DFDP        ! Gradient of "Yield" criterion
      REAL P0
      REAL STRESSK
      REAL DPDT

      REAL DT

      REAL P0DOT
      REAL S
      REAL SV
      REAL DSVDP
      
      REAL EP
      REAL DEP

      PARAMETER(one=1.d0,zero=0.d0)

!-----------------------------------------------------------------------
!     Coefficients for viscous stress
!-----------------------------------------------------------------------
      S = 5e6
      P0DOT = 0.0001


      EP = 0.0
      DEP = 0.0

!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      YOUNG  = PROPS(1)
      POISS  = PROPS(2)
      SIGMA0 = PROPS(3)
      ET     = PROPS(4)
      TOL    = PROPS(5)
      MXITER = PROPS(6)
!-----------------------------------------------------------------------
!     Compute elasticity matrix
!-----------------------------------------------------------------------
      C11    = YOUNG*(1.0-POISS)/((1.0+POISS)*(1.0-2.0*POISS)) !la+2nu
      C12    = POISS*C11/(1.0-POISS) !Lame lambda
      C44    = 0.5*YOUNG/(1.0+POISS) !Lame nu
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Loop over integration points
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO i=1,NBLOCK
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        If time = 0 then pure elastic computation
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         IF(TOTALTIME.eq.0.0)THEN
            STRESSNEW(i,1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                       +C11*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C11*STRAININC(i,3)
            STRESSNEW(i,4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESSNEW(i,5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESSNEW(i,6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0
         ELSE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        Elastic-predictor-corrector scheme
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Elastic prediction
!-----------------------------------------------------------------------
            STRESS(1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                                +C11*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C11*STRAININC(i,3)
            STRESS(4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESS(5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESS(6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0   
!-----------------------------------------------------------------------
!           Equivalent stress
!-----------------------------------------------------------------------
            PHI = eq(STRESS)
!-----------------------------------------------------------------------
!           Equivalent plastic strain from previous time step
!-----------------------------------------------------------------------
            P         = STATEOLD(i,1) 
            P0        = STATEOLD(i,1)

!-----------------------------------------------------------------------
!           Initialize the plastic multiplier
!-----------------------------------------------------------------------
            DLAMBDA  = 0.00001

            SV = S*log(one+DLAMBDA/(DT*P0DOT))
            DSVDP = S/(DT*P0DOT+DLAMBDA)
!-----------------------------------------------------------------------
!           Update yield stress
!-----------------------------------------------------------------------
            SIGMAY   = SIGMA0+ET*P

!-----------------------------------------------------------------------
!           Compute yield function
!-----------------------------------------------------------------------
            F        = PHI-SIGMAY - 3*C44*DLAMBDA - SV
            DFDP     = -ET -3.0*C44 -DSVDP

!-----------------------------------------------------------------------
!           Compute the derivative of the yield function
!-----------------------------------------------------------------------
            IF(PHI.eq.0)THEN
               DENOM = 1.0
            ELSE
               DENOM = PHI
            ENDIF        
c
            DFDS(1) = (STRESS(1)-0.5*(STRESS(2)+STRESS(3)))/DENOM
            DFDS(2) = (STRESS(2)-0.5*(STRESS(3)+STRESS(1)))/DENOM
            DFDS(3) = (STRESS(3)-0.5*(STRESS(1)+STRESS(2)))/DENOM
            DFDS(4) = 3.0*STRESS(4)/DENOM
            DFDS(5) = 3.0*STRESS(5)/DENOM
            DFDS(6) = 3.0*STRESS(6)/DENOM
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!           Check for plasticity
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Plastic flow
!-----------------------------------------------------------------------
            IF(F.GT.0.0)THEN ! Plastic flow
               DO ITER=1,MXITER
!-----------------------------------------------------------------------
!                 Compute increment in plastic multiplier
!-----------------------------------------------------------------------
                  DDLAMBDA  = -F/DFDP
!-----------------------------------------------------------------------
!                 Update equivalent plastic strain
!-----------------------------------------------------------------------
                  DLAMBDA = DLAMBDA+DDLAMBDA
                  P      = P0 +DLAMBDA

!                  DEP(1) = -DDLAMBDA*DFDS(1)
!                  DEP(2) = -DDLAMBDA*DFDS(2)
!                  DEP(3) = -DDLAMBDA*DFDS(3)
!                  DEP(4) = -DDLAMBDA*DFDS(4)
!                  DEP(5) = -DDLAMBDA*DFDS(5)
!                  DEP(6) = -DDLAMBDA*DFDS(6)
!
!                  print *,DFDS(1)
!
!                  EP(1) = EP(1) + DEP(1)
!                  EP(2) = EP(2) + DEP(2)
!                  EP(3) = EP(3) + DEP(3)
!                  EP(4) = EP(4) + DEP(4)
!                  EP(5) = EP(5) + DEP(5)
!                  EP(6) = EP(6) + DEP(6)
!
!
!            STRESS(1) = STRESSOLD(i,1)+C11*(STRAININC(i,1)-EP(1))
!     .                                +C12*(STRAININC(i,2)-EP(2))
!     .                                +C12*(STRAININC(i,3)-EP(3))
!            STRESS(2) = STRESSOLD(i,2)+C12*(STRAININC(i,1)-EP(1))
!     .                                +C11*(STRAININC(i,2)-EP(2))
!     .                                +C12*(STRAININC(i,3)-EP(3))
!            STRESS(3) = STRESSOLD(i,3)+C12*(STRAININC(i,1)-EP(1))
!     .                                +C12*(STRAININC(i,2)-EP(2))
!     .                                +C11*(STRAININC(i,3)-EP(3))
!            STRESS(4) = STRESSOLD(i,4)+C44*(STRAININC(i,4)-EP(4))*2.0
!            STRESS(5) = STRESSOLD(i,5)+C44*(STRAININC(i,5)-EP(5))*2.0
!            STRESS(6) = STRESSOLD(i,6)+C44*(STRAININC(i,6)-EP(6))*2.0



                                    STRESSK(1)= STRESS(1)-DLAMBDA*(C11*DFDS(1)
     .                             +C12*DFDS(2)
     .                        +C12*DFDS(3))
                  STRESSK(2)= STRESS(2)-DLAMBDA*(C12*DFDS(1)
     .                        +C11*DFDS(2)
     .                        +C12*DFDS(3))
                  STRESSK(3)= STRESS(3)-DLAMBDA*(C12*DFDS(1)
     .                        +C12*DFDS(2)
     .                        +C11*DFDS(3))
                  STRESSK(4)= STRESS(4)-DLAMBDA*C44*DFDS(4)
                  STRESSK(5)= STRESS(5)-DLAMBDA*C44*DFDS(5)
                  STRESSK(6)= STRESS(6)-DLAMBDA*C44*DFDS(6)
                  
                  SIGMAY   = SIGMA0+ET*P
                  PHI = eq(STRESSK)

!                  F        = PHI-SIGMAY - 3*C44*DLAMBDA - SV
                  F        = PHI-SIGMAY - SV

                  DFDP     = -ET -3.0*C44 - DSVDP

!-----------------------------------------------------------------------
!                 Update viscous stress
!-----------------------------------------------------------------------
                  IF(DLAMBDA.GT.zero)THEN
                        SV = S*log(one+DLAMBDA/(DT*P0DOT))
                        DSVDP = S/(DT*P0DOT+DLAMBDA)
                  ELSE
                        SV = zero
                        DSVDP = zero
                  ENDIF

!-----------------------------------------------------------------------
!                 Update yield stress
!-----------------------------------------------------------------------
                  
!-----------------------------------------------------------------------
!           Compute yield function
!-----------------------------------------------------------------------




!                  PHI = eq(STRESSK)


                  IF(PHI.eq.0)THEN
                     DENOM = 1.0
                  ELSE
                     DENOM = PHI
                  ENDIF        
      
!                  DFDS(1) = (STRESS(1)-0.5*(STRESS(2)+STRESS(3)))/DENOM
!                  DFDS(2) = (STRESS(2)-0.5*(STRESS(3)+STRESS(1)))/DENOM
!                  DFDS(3) = (STRESS(3)-0.5*(STRESS(1)+STRESS(2)))/DENOM
!                  DFDS(4) = 3.0*STRESS(4)/DENOM
!                  DFDS(5) = 3.0*STRESS(5)/DENOM
!                  DFDS(6) = 3.0*STRESS(6)/DENOM
                 DFDS(1) = (STRESSK(1)-0.5*(STRESSK(2)+STRESSK(3)))/DENOM
                 DFDS(2) = (STRESSK(2)-0.5*(STRESSK(3)+STRESSK(1)))/DENOM
                 DFDS(3) = (STRESSK(3)-0.5*(STRESSK(1)+STRESSK(2)))/DENOM
                 DFDS(4) = 3.0*STRESSK(4)/DENOM
                 DFDS(5) = 3.0*STRESSK(5)/DENOM
                 DFDS(6) = 3.0*STRESSK(6)/DENOM



!-----------------------------------------------------------------------
!                 Compute convergence criterion
!-----------------------------------------------------------------------
                  RESNOR = ABS(F/SIGMAY)
                  

!-----------------------------------------------------------------------
!                 Check for convergence
!-----------------------------------------------------------------------
                  IF(RESNOR.LE.TOL)THEN ! RMAP has converged
                  PRINT *,ITER
!-----------------------------------------------------------------------
!                    Update the stress tensor
!-----------------------------------------------------------------------
                        STRESSNEW(i,1)= STRESS(1)-DLAMBDA*(C11*DFDS(1)
     .                             +C12*DFDS(2)
     .                        +C12*DFDS(3))
                        STRESSNEW(i,2)= STRESS(2)-DLAMBDA*(C12*DFDS(1)
     .                        +C11*DFDS(2)
     .                        +C12*DFDS(3))
                        STRESSNEW(i,3)= STRESS(3)-DLAMBDA*(C12*DFDS(1)
     .                        +C12*DFDS(2)
     .                        +C11*DFDS(3))
                        STRESSNEW(i,4)= STRESS(4)-DLAMBDA*C44*DFDS(4)
                        STRESSNEW(i,5)= STRESS(5)-DLAMBDA*C44*DFDS(5)
                        STRESSNEW(i,6)= STRESS(6)-DLAMBDA*C44*DFDS(6)
!-----------------------------------------------------------------------
!                    Update the history variables
!-----------------------------------------------------------------------
                     STATENEW(i,1) = P
                     STATENEW(i,2) = PHI-3.0*C44*DLAMBDA
                     STATENEW(i,3) = F
                     STATENEW(i,4) = SIGMAY
                     STATENEW(i,5) = DGAMA
                     STATENEW(i,6) = ITER 
                     GOTO 90
                  ELSE ! RMAP has not converged yet
                     IF(ITER.eq.MXITER)THEN
                        write(*,*) 'RMAP has not converged'
                        write(*,*) 'Integration point',i
                        write(*,*) 'Convergence',RESNOR
                        write(*,*) 'dlambda',DLAMBDA,P
                        STOP
                     ENDIF
                  ENDIF          
               ENDDO
!-----------------------------------------------------------------------
!           Elastic point
!-----------------------------------------------------------------------
            ELSE
!-----------------------------------------------------------------------
!              Update the stress tensor
!-----------------------------------------------------------------------
               STRESSNEW(i,1) = STRESS(1)
               STRESSNEW(i,2) = STRESS(2)
               STRESSNEW(i,3) = STRESS(3)
               STRESSNEW(i,4) = STRESS(4)
               STRESSNEW(i,5) = STRESS(5)
               STRESSNEW(i,6) = STRESS(6)
!-----------------------------------------------------------------------
!              Update the history variables
!-----------------------------------------------------------------------
               STATENEW(i,1) = P
               STATENEW(i,2) = PHI
               STATENEW(i,3) = F
               STATENEW(i,4) = SIGMAY
               STATENEW(i,5) = 0.0
               STATENEW(i,6) = 0
            ENDIF
!-----------------------------------------------------------------------
!        End of loop over integration points
!-----------------------------------------------------------------------
         ENDIF           
  90     CONTINUE  
      ENDDO
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END