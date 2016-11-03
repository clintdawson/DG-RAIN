C **************************************************************
C
C     SUBROUTINE CHANNEL_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C
C     1. Calculates the values of the necessary variables at the
C        edge gauss points for INTERNAL edges
C      2. Calls the appropriate subroutine to compute the flux at these
C      points.
C      3. Calls the appropriate subroutine to compute the boundary
C      integrals.

C ***************************************************************

      SUBROUTINE CHANNEL_EDGE_HYDRO(IT)


#ifdef DGSHED
C.....Use appropriate modules

      USE GLOBAL
      USE DG
C      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset

      IMPLICIT NONE

C.....Declare local variables
      INTEGER L, LED, GED, I, K, IT
      Real(SZ) DEN
      REAL(SZ) TX, TY, W_IN 
      Real(SZ) EDG_LEN

      test_el = 0

      DO 1000 L = 1, NCHEDS

C.....Retrieve the global and local edge number
        GED = CHEDS(L)
        LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge
        EL_IN = NEDEL(1,GED)

C        WRITE(*,*) "EL_IN = ", EL_In

        pa = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
      
        test_el = test_el + 1

C.....Retrieve the components of the normal vector to the edge

        NX = COSNX(GED)
        NY = SINNX(GED)

C......Calculate the tangential components
        TX = -NY
        TY = NX
        DEN = 1.D0/(NX*TY - NY*TX)

C.....Retrieve the nodes of the edge
C        N1 = NEDNO(1,GED)
C        N2 = NEDNO(2,GED)

C.....Retrieve the values of water height and normal flow at the exterior state. This will be a constant over an edge.

        EDG_LEN = XLEN(GED)
        ZE_EX = CHVALS(1,L)
        Q_N_EXT = CHVALS(2,L)/EDG_LEN 

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)

            HB_IN = BATHED(I,LED,EL_IN,pa)

C.....Compute the solution at the interior state

            DO K = 2,DOFS(EL_IN)

               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)               

            ENDDO
             

C.....Compute velocity in the tangetial direction
            Q_T_INT = QX_IN*TX + QY_IN*TY
            Q_T_EXT = Q_T_INT

C.....Check to see if SWEVALS has been initialized
            IF (ZE_EX.LT.(-99999) .AND. Q_N_EXT  < -99999 ) THEN
              ZE_EX = ZE_IN 
              Q_N_EXT = QX_IN*NX + QY_IN*NY
            ENDIF

C.....Compute the x and y components of the external state flow
            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

C.....Set the exterior value of the bathymetry equal to the interior
C.....minus the dry bed level

            HB_EX = HB_IN

C            IF (LoadGeoidOffset) then
C               ZE_EX = ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
C            endif

C.....Compute the numerical flux
             CALL NUMERICAL_FLUX(IT, test_el)

C.....Compute the edge integral
            DO K = 1,DOFS(EL_IN)

               W_IN = 2.0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

           ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
#endif

      END SUBROUTINE

