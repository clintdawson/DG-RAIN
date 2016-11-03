C***********************************************************************
C     
C     SUBROUTINE FLOODPLAIN_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for floodplain edges that are connected to channels
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C    
C     
C***********************************************************************

       SUBROUTINE FLOODPLAIN_EDGE_HYDRO(IT)
 
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables
      INTEGER L, LED_IN, GP_IN, GED, I, K, IT, LED_EX, GP_EX
      REAL(SZ) DEN, TX, TY, W_IN, W_EX
      REAL(SZ) EDG_LEN
      
      test_el = 0
      DO 1000 L = 1,NFPEDS

C.....Check to see if CHTOFPVALS has been filled to see if any
C.....calculation needs to be done here
        IF (ABS(CHTOFPVALS(L)).LE.1e-13) GOTO 1000
C.....Retrieve the global and local edge number

         GED = FPEDS(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.....Retrieve the first element to which the edge belongs

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2, GED)

C.....Organize the dofs for the p_adaptivity
         EL = EL_EX
         IF (DOFS(EL_EX) .LT. DOFS(EL_IN)) THEN
            EL = EL_IN
         ENDIF
         PA = PDG_EL(EL)

#ifdef P0
         if (PA.eq.0) then
            PA = 1
         endif
#endif


C.....Retrieve the components of the normal vector to the edge
         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX
         DEN = 1.D0/(NX*TY - NY*TX)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)
          
            ZE_IN = 0.D0
            QX_IN = 0.D0
            QY_IN = 0.D0

            HB_IN = 0.D0

            
C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_IN)

               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED_IN,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED_IN,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED_IN,pa)
               

               HB_IN = HB_IN + HB(K,EL_IN,1)*PHI_EDGE(K,I,LED_IN,pa)
            ENDDO
#ifdef P0
            IF (DOFS(EL_IN).EQ.1) THEN
               DO K = 2,3
                  HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,I,LED_IN,pa)
               ENDDO
            ENDIF
#endif
 
C.....Calculate the tangential velocity at the interior state
            Q_T_INT = QX_IN*TX + QY_IN*TY
            Q_T_EXT = Q_T_INT

           
C.....Set the exterior bed and surface elevation equal to the interior

            ZE_EX = ZE_IN
            HB_EX = HB_IN

C.....Retrieve the values of water height and normal flow at the
C exterior state. This will be constant over an edge.

            EDG_LEN = XLEN(GED)
            Q_N_EXT = CHTOFPVALS(L)/EDG_LEN

C            WRITE(*,*) "L = ", L , "I = ", I, "EDG_LEN = ", EDG_LEN

C            WRITE(*,*) "excessWaterRate = ", CHTOFPVALS(L), " Q_N_EXT = ", Q_N_EXT
C            WRITE(*,*) "El_IN = ", EL_IN, " EL_EX = ", EL_EX

C.....      Compute the x and y components of the external state flow
            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

C.....Com   pute the numerical flux
            CALL NUMERICAL_FLUX(IT, test_el) 
C            WRITE(*,*) " FHAT = ", F_HAT, "GHAT = ", G_HAT, "HHAT = ", H_HAT

C.....Com   pute the edge integral

            DO K = 1,DOFS(EL_IN)

               W_IN = 2.0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED_IN,pa)*WEGP(I,pa)


               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

C               WRITE(*,*) "RHS_ZE(K,EL_IN,IRK) = ",RHS_ZE(K,EL_IN,IRK)  
            ENDDO

         ENDDO


C.....Retrieve the second element to which the edge belongs

         EL_EX = NEDEL(2,GED)
         LED_EX = NEDSD(2, GED)


C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         DO I = 1,NEGP(pa)

            ZE_EX = 0.D0
            QX_EX = 0.D0
            QY_EX = 0.D0
            
            HB_EX = 0.D0


            GP_EX = NEGP(pa) - I + 1
            
C.....Compute the solution at the interior state

            DO K = 1,DOFS(EL_EX)

               ZE_EX = ZE_EX + ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               QX_EX = QX_EX + QX(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               QY_EX = QY_EX + QY(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                          

               HB_EX = HB_EX + HB(K,EL_EX,1)*PHI_EDGE(K,GP_EX,LED_EX,pa)
            ENDDO

#ifdef P0
            IF (DOFS(EL_EX).EQ.1) THEN
               DO K = 2,3
                  HB_EX = HB_EX + HB(K,EL_EX,1  )*PHI_EDGE(K,GP_EX,LED_EX,pa)
               ENDDO
            ENDIF
#endif
 
C.....Calculate the tangential velocity at the exterior state
            Q_T_EXT = QX_EX*TX + QY_EX*TY
            Q_T_INT = Q_T_EXT

           
C.....Set the exterior bed and surface elevation equal to the interior

            ZE_IN = ZE_EX
            HB_IN = HB_EX

C.....Retrieve the values of water height and normal flow at the
C exterior state. This will be constant over an edge.

            Q_N_INT = -CHTOFPVALS(L)/EDG_LEN

C            WRITE(*,*) "L = ", L , "I = ", I, "EDG_LEN = ", EDG_LEN

C            WRITE(*,*) "excessWaterRate = ", CHTOFPVALS(L), " Q_N_EXT = ", Q_N_EXT

C.....Com   pute the x and y components of the external state flow
            QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)*DEN
            QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)*DEN

C.....Com   pute the numerical flux
            CALL NUMERICAL_FLUX(IT, test_el) 
C            WRITE(*,*) " FHAT = ", F_HAT, "GHAT = ", G_HAT, "HHAT = ", H_HAT

C.....Com   pute the edge integral

            DO K = 1,DOFS(EL_EX)

               W_EX = 2.0*M_INV(K,pa)/AREAS(EL_EX)*XLEN(GED)*
     &              PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
C               WRITe(*,*) " W_IN = ", W_IN
C               write(*,*) "F_hat = ", F_HAT

               RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT
               RHS_QX(K,EL_EX,IRK) = RHS_QX(K,EL_EX,IRK) + W_EX*G_HAT
               RHS_QY(K,EL_EX,IRK) = RHS_QY(K,EL_EX,IRK) + W_EX*H_HAT

C               WRITE(*,*) "RHS_ZE(K,EL_IN,IRK) = ",RHS_ZE(K,EL_IN,IRK)  
            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
