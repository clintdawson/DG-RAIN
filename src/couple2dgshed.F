
      MODULE Couple2Dgshed

#ifdef DGSHED
      
      USE SIZES, ONLY: SZ
      USE GLOBAL
      USE DG

      CONTAINS

      SUBROUTINE INIT_DGSHED()

        INTERFACE
          SUBROUTINE DGSHED_SETUP() BIND(C)
            IMPLICIT NONE
          END SUBROUTINE DGSHED_SETUP
        END INTERFACE

        CALL DGSHED_SETUP()
      END SUBROUTINE INIT_DGSHED

      SUBROUTINE POPULATE_SWE_ARRAY()

        REAL(SZ) Q_N_SUM, ZE_AVG

        DO L = 1, NCHEDS

C.....Retrive the global and local edge number
          GED = CHEDS(L)
          LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge
          EL_IN = NEDEL(1,GED)
C          WRITE (*,*) "EL_IN = ", EL_IN 

          pa = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.....Retrieve the components of the normal vector to the edge

          NX = COSNX(GED)
          NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point
          ZE_AVG = 0.D0
          Q_N_SUM = 0.D0

          DO I = 1,NEGP(pa)
            
            ZE_IN = ZE(1,EL_IN,1)
            QX_IN = QX(1,EL_IN,1)
            QY_IN = QY(1,EL_IN,1)
C            HB_IN = BATHED(I, LED, EL_IN,pa)

            DO K = 2,DOFS(EL_IN)
              ZE_IN = ZE_IN + ZE(K,EL_IN,1)*PHI_EDGE(K,I,LED,pa)
              QX_IN = QX_IN + QX(K,EL_IN,1)*PHI_EDGE(K,I,LED,pa)
              QY_IN = QY_IN + QY(K,EL_IN,1)*PHI_EDGE(K,I,LED,pa)
            ENDDO

            Q_N_INT = QX_IN*NX + QY_IN*NY

            ZE_AVG = ZE_AVG + ZE_IN*WEGP(I,pa)

            Q_N_SUM = Q_N_SUM + Q_N_INT*WEGP(I,pa)

          ENDDO

          SWEVALS(1,L) = 0.5*ZE_AVG
          SWEVALS(2,L) = 0.5*XLEN(GED)*Q_N_SUM

C          WRITE(*,*) "ZE_AVG_DGSWEM = ", SWEVALS(1,L), "Q_N_DGSWEM
C     & = ", SWEVALS(2,L)
        ENDDO 
        

      END SUBROUTINE POPULATE_SWE_ARRAY

      SUBROUTINE GET_CHANNELS_VALUES()
        
        INTERFACE
          SUBROUTINE OBTAIN_CHANNELS_VALUES(NCHEDS, CHNUMS,CHVALS) BIND(C)
             USE, INTRINSIC :: ISO_C_BINDING
             IMPLICIT NONE
             INTEGER(C_INT), INTENT(IN) :: NCHEDS 
             INTEGER(C_INT), INTENT(IN), DIMENSION(2,NCHEDS) :: CHNUMS
             REAL(C_DOUBLE), INTENT(OUT), DIMENSION(2,NCHEDS) :: CHVALS
          END SUBROUTINE OBTAIN_CHANNELS_VALUES
        END INTERFACE

        CALL OBTAIN_CHANNELS_VALUES(NCHEDS, CHNUMS, CHVALS)
C        WRITE(*,*) "CHVAL1 = ", CHVALS(1,1), " CHVAL2 = ",
C     &CHVALS(2,1)
       
       END SUBROUTINE GET_CHANNELS_VALUES

       SUBROUTINE SHARE_SWE_VALUES()
        
        INTERFACE
          SUBROUTINE RECEIVE_SWE_VALUES(NCHEDS, CHNUMS, SWEVALS) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: NCHEDS
            INTEGER(C_INT), INTENT(IN), DIMENSION(2,NCHEDS) :: CHNUMS
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(2,NCHEDS) ::SWEVALS 
          END SUBROUTINE RECEIVE_SWE_VALUES
        END INTERFACE
        
        CALL RECEIVE_SWE_VALUES(NCHEDS, CHNUMS, SWEVALS)

       END SUBROUTINE SHARE_SWE_VALUES

      SUBROUTINE ADVANCE_DGSHED(CURRTIME, DT)
        REAL(SZ), INTENT(IN) :: CURRTIME, DT

        INTERFACE
          SUBROUTINE TIME_EVOLUTION(CURRTIME, DT) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: CURRTIME
            REAL(C_DOUBLE), INTENT(IN) :: DT
          END SUBROUTINE
        END INTERFACE

C        WRITE(*,*) "CURRTIME = ", CURRTIME
C        WRITE(*,*) "DT = ", DT
        CALL TIME_EVOLUTION(CURRTIME, DT)
      END SUBROUTINE ADVANCE_DGSHED

      SUBROUTINE COUNT_CONNECTION_WITH_FLOODPLAINS()
        INTERFACE
          SUBROUTINE CALCULATE_NUMBER_OF_BANKS(NFPEDS, FPCHCON) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN):: NFPEDS
            INTEGER(C_INT), INTENT(IN), DIMENSION(2,NFPEDS) :: FPCHCON
          END SUBROUTINE CALCULATE_NUMBER_OF_BANKS
        END INTERFACE

        CALL CALCULATE_NUMBER_OF_BANKS(NFPEDS, FPCHCON)

      END SUBROUTINE COUNT_CONNECTION_WITH_FLOODPLAINS

      SUBROUTINE RECEIVE_FLOODWATER()
        
        INTERFACE
          SUBROUTINE SEND_WATER_TO_FLOODPLAINS(NFPEDS, FPCHCON, 
     &     CHTOFPVALS) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN) :: NFPEDS
            INTEGER(C_INT), INTENT(IN), DIMENSION(2,NFPEDS) :: FPCHCON
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(NFPEDS) ::CHTOFPVALS 
          END SUBROUTINE SEND_WATER_TO_FLOODPLAINS
        END INTERFACE
        
        CALL SEND_WATER_TO_FLOODPLAINS(NFPEDS, FPCHCON,
     &   CHTOFPVALS)

      END SUBROUTINE RECEIVE_FLOODWATER

      SUBROUTINE PASS_AVG_FLOODPLAIN_ELEVATION()
        
        INTEGER L, I, K, LED1, GP1, GED, LED2, GP2, EL1, EL2, EL
        REAL(SZ) ZE1, ZE2, W1, W2, AVGZE1, AVGZE2, AVGZE, QX1, QX2
        REAL(SZ) QY1, QY2, TOTALQN1, TOTALQN2, TOTALQN

        INTERFACE
          SUBROUTINE UPDATE_CHANNEL_DEPTH(NFPEDS, FPCHCON, AVGFPZE,
     &             AVGFPQN) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(C_INT), INTENT(IN):: NFPEDS
            INTEGER(C_INT), INTENT(IN), DIMENSION(2,NFPEDS) :: FPCHCON
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(NFPEDS) :: AVGFPZE
            REAL(C_DOUBLE), INTENT(IN), DIMENSION(NFPEDS) :: AVGFPQN
          END SUBROUTINE UPDATE_CHANNEL_DEPTH
        END INTERFACE

        DO 1000 L = 1,NFPEDS
C.........Retrieve the global and local edge number          
          GED = FPEDS(L)
          LED1 = NEDSD(1,GED)
          LED2 = NEDSD(2,GED)

          NX = COSNX(GED)
          NY = SINNX(GED)

C.........Retreive the two elements to which the edge belongs

          EL1 = NEDEL(1,GED)
          EL2 = NEDEL(2,GED)

C.....Organize the dofs for the p_adaptivity
          EL = EL2
          IF (DOFS(EL2) .LT. DOFS(EL1)) THEN
             EL = EL1
          ENDIF
          PA = PDG_EL(EL)

#ifdef P0 
          if (PA.eq.0) then
             PA = 1
          endif
#endif


C........Compute the edge integral of ZE 
          AVGZE1 = 0.D0
          AVGZE2 = 0.D0
          TOTALQN1 = 0.D0
          TOTALQN2 = 0.D0

          AVGZE = 0.D0
          TOTALQN = 0.D0
          DO I = 1, NEGP(PA)
            
            GP1 = I
            GP2 = NEGP(PA) - I + 1
            
            ZE1 = 0.D0
            ZE2 = 0.D0
            QX1 = 0.D0
            QX2 = 0.D0
            QY1 = 0.D0
            QY2 = 0.D0
            
            DO K = 1, DOFS(EL1)
              ZE1 = ZE1 + ZE(K, EL1, 1)*PHI_EDGE(K,GP1,LED1,PA)
              QX1 = QX1 + QX(K, EL1, 1)*PHI_EDGE(K,GP1,LED1,PA)
              QY1 = QY1 + QY(K, EL1, 1)*PHI_EDGE(K,GP1,LED1,PA)
              
              ZE2 = ZE2 + ZE(K, EL2, 1)*PHI_EDGE(K,GP2,LED2,PA)
              QX2 = QX2 + QX(K, EL2, 1)*PHI_EDGE(K,GP2,LED2,PA)
              QY2 = QY2 + QY(K, EL2, 1)*PHI_EDGE(K,GP2,LED2,PA)
              
            ENDDO
          
            AVGZE1 = AVGZE1 + ZE1*WEGP(GP1, PA)
            AVGZE2 = AVGZE2 + ZE2*WEGP(GP2, PA)
            
            TOTALQN1 = TOTALQN1 + WEGP(GP1, PA)*(QX1*NX + QY1*NY)
            TOTALQN2 = TOTALQN2 - WEGP(GP2, PA)*(QX2*NX + QY2*NY)
    
          ENDDO
          
          AVGZE1 = 0.5*AVGZE1
          AVGZE2 = 0.5*AVGZE2
          AVGZE = 0.5*(AVGZE1+AVGZE2)
          TOTALQN1 = 0.5*TOTALQN1;
          TOTALQN2 = 0.5*TOTALQN2;
          TOTALQN = TOTALQN1+TOTALQN2
          
          AVGFPZE(L) = AVGZE
          TOTALFPQN(L) = TOTALQN

 1000   CONTINUE
        CALL UPDATE_CHANNEL_DEPTH(NFPEDS, FPCHCON, AVGFPZE, TOTALFPQN)

      END SUBROUTINE PASS_AVG_FLOODPLAIN_ELEVATION

#endif

      END MODULE Couple2Dgshed

