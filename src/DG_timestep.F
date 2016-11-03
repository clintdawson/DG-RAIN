
C     
C     SUBROUTINE DG_TIMESTEP(IT)
C     
C     This is the main subroutine for the DG hydro
C     
C     Written by Shintaro Bunya (01-01-2007)
C     
C***********************************************************************

      SUBROUTINE DG_TIMESTEP(IT)

C.....Use appropriate modules
      
      USE SIZES, ONLY : SZ,layers
      USE GLOBAL, ONLY : SEDFLAG, DTDP, STATIM


#ifdef CMPI
      USE MESSENGER_ELEM
#endif

#ifdef DGSHED
      USE Couple2Dgshed
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER IT,OUTIT
      Real(SZ) MITDT, CURRTIME

#ifdef DGSHED
C.....Advance DGSHED to the next time step
C.....Usually time steps for the DGSHED part will be smaller than that
C.....for DGSWEM
      MITDT = (IT-1)*DTDP
      CURRTIME = MITDT + STATIM*86400.D0
      
C      IF (NCHEDS.GT.0) CALL ADVANCE_DGSHED(CURRTIME,DT)
      IF (NCHEDS.GT.0 .OR. NFPEDS.GT.0) CALL ADVANCE_DGSHED(CURRTIME,DT)
#endif

C.....Hydrodynamics

      CALL DG_HYDRO_TIMESTEP(IT)

#ifdef DGSHED

C.....Pass values from shallow water region to channels

C.....Get values from channels

      IF(NCHEDS.GT.0) CALL GET_CHANNELS_VALUES()

      IF (NCHEDS.GT.0) THEN
        CALL POPULATE_SWE_ARRAY()
        CALL SHARE_SWE_VALUES()

      ENDIF

C.....Receive the water flooding out of channels into their banks
      IF (NFPEDS.GT.0) THEN
        CALL RECEIVE_FLOODWATER()
      ENDIF

C.....Update the depth of channels for the next time step
      IF (NFPEDS.GT.0) THEN
        CALL PASS_AVG_FLOODPLAIN_ELEVATION()
      ENDIF

#endif

      
C.....Write out results

#ifdef CMPI
#ifdef BLKOUT
      CALL MSG_BLOCKSYNC_START()
#endif
#endif
      CALL WRITE_RESULTS(IT,.FALSE.)
#ifdef CMPI
#ifdef BLKOUT
      CALL MSG_BLOCKSYNC_FINISH()
#endif
#endif

      CALL SCRUTINIZE_SOLUTION()

      RETURN
      END SUBROUTINE


C***********************************************************************
C     
C     SUBROUTINE SCRUTINIZE_SOLUTION(IT)
C     
C     Detect a negative water depth and a NaN in solution arrays and 
C     terminate if either of them is found.
C     
C     Written by Shintaro Bunya (01-01-2007)
C     
C***********************************************************************

      SUBROUTINE SCRUTINIZE_SOLUTION(IT)

C.....Use appropriate modules
      
#ifdef CMPI
      USE SIZES, ONLY : SZ,MYPROC,layers
      USE MESSENGER_ELEM, ONLY: MESSAGE_FINI, ErrorElevSum
#else
      USE SIZES, ONLY : SZ,layers
#endif
      USE GLOBAL, ONLY : NE,NM,DP,X,Y,SLAM,SFEA,DEG2RAD,tracer_flag,chem_flag,pdg_el
      USE DG, ONLY : DOF,ZE,QX,QY,iota,iota2,dynP,bed,ncheck,dp_node

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER IT,J,K,l
      LOGICAL Detected
      REAL(SZ) DPAVG

      INTEGER ErrorElevExceeded

C.....Detect negative mass
      ErrorElevExceeded=0
      DO J = 1,NE
         DPAVG =  1.D0/ncheck(pdg_el(j)) * (sum(DP_NODE(:,J,pdg_el(j)))) !(DP(NM(J,1))+DP(NM(J,2))+DP(NM(J,3)))/3.D0
         Detected = .FALSE.
         IF((ZE(1,J,1)+DPAVG).LE.0.D0) THEN
            PRINT *, ''
#ifdef CMPI
            PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
            PRINT *, '  ZE(',1,',',J,',1) = ',ZE(1,J,1)
            PRINT *, '  HE(',1,',',J,',1) = ',ZE(1,J,1)+DPAVG
            PRINT *, '  DP = ',DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3))
            PRINT *, '  DPAVG = ',DPAVG
            PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $           SFEA(NM(J,1))/DEG2RAD
            Detected = .TRUE.
            ErrorElevExceeded=1    
         ENDIF
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
      IF(ErrorElevExceeded.NE.0) THEN
         if (myproc.eq.0) then
            PRINT *, ''
            PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
            PRINT *, ''
         endif

         call message_fini ()
         STOP
      endif
#else
c$$$      IF(ErrorElevExceeded.NE.0) THEN
c$$$         PRINT *, ''
c$$$         PRINT *,'  PROGRAM WILL BE TERMINATED in DG_timestep'
c$$$         PRINT *, ''
c$$$         stop
c$$$      ENDIF
#endif


c$$$
c$$$         IF(Detected) THEN
c$$$            PRINT *, ''
c$$$            PRINT *,'  WRITING OUT THE LAST-MOMENT RESULTS...'
c$$$            PRINT *, ''
c$$$            CALL WRITE_RESULTS(IT,.TRUE.)
c$$$            PRINT *, ''
c$$$            PRINT *,'  PROGRAM WILL BE TERMINATED.'
c$$$            PRINT *, ''
c$$$#ifdef CMPI
c$$$            call message_fini ()
c$$$            STOP
c$$$#else
c$$$            stop
c$$$#endif
c$$$         ENDIF
c$$$      ENDDO

C.....Detect NaN
      
      DO J = 1,NE
         DO K = 1,DOF

            Detected = .FALSE.
            IF(ZE(K,J,1).NE.ZE(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  ZE(',K,',',J,',1) IS ',ZE(K,J,1)
               PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $              SFEA(NM(J,1))/DEG2RAD
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
            IF(QX(K,J,1).NE.QX(K,J,1)) THEN
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, ''
               PRINT *, '  QX(',K,',',J,',1) IS ',QX(K,J,1)
               PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $              SFEA(NM(J,1))/DEG2RAD
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
            IF(QY(K,J,1).NE.QY(K,J,1)) THEN
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, ''
               PRINT *, '  QY(',K,',',J,',1) IS ',QY(K,J,1)
               PRINT *, ' x,  y ',SLAM(NM(J,1))/DEG2RAD,
     $              SFEA(NM(J,1))/DEG2RAD
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF

#ifdef TRACE
            IF(iota(K,J,1).NE.iota(K,J,1)) THEN
               PRINT *, ''

#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  iota(',K,',',J,',1) IS ',iota(K,J,1)
               PRINT *, ' x,  y ',X(NM(J,1)),Y(NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif
            
#ifdef CHEM
            IF(iota(K,J,1).NE.iota(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  iota(',K,',',J,',1) IS ',iota(K,J,1)
               PRINT *, ' x,  y ',X(NM(J,1)),Y(NM(J,1))
               Detected = .TRUE.
            ENDIF
            IF(iota2(K,J,1).NE.iota2(K,J,1)) THEN
               PRINT *, ''
#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  iota2(',K,',',J,',1) IS ',iota2(K,J,1)
               PRINT *, ' x,  y ',X(NM(J,1)),Y(NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif

#ifdef DYNP
            IF(dynP(K,J,1).NE.dynP(K,J,1)) THEN
               PRINT *, ''

#ifdef CMPI
               PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
               PRINT *, '  dynP(',K,',',J,',1) IS ',dynP(K,J,1)
               PRINT *, ' x,  y ',X(NM(J,1)),Y(NM(J,1))
               Detected = .TRUE.
               ErrorElevExceeded = 1
            ENDIF
#endif
#ifdef SED_LAY
            
            do l = 1,layers
               IF(bed(K,J,1,l).NE.bed(K,J,1,l)) THEN
                  PRINT *, ''
                  
#ifdef CMPI
                  PRINT *, 'IN SUBDOMAIN ',MYPROC
#endif
                  PRINT *, '  bed(',K,',',J,',1,',l,') IS ', bed(K,J,1,l)
                  PRINT *, ' x,  y ',X(NM(J,1)),Y(NM(J,1))
                  Detected = .TRUE.
                  ErrorElevExceeded = 1
               ENDIF
            enddo
#endif

c$$$            IF(Detected) THEN
c$$$               PRINT *, ''
c$$$               PRINT *,'  WRITING OUT THE LAST-MOMENT RESULTS...'
c$$$               PRINT *, ''
c$$$               CALL WRITE_RESULTS(IT,.TRUE.)
c$$$               PRINT *, ''
c$$$               PRINT *,'  PROGRAM WILL BE TERMINATED.'
c$$$               PRINT *, ''
c$$$#ifdef CMPI
c$$$               call message_fini()
c$$$               STOP
c$$$#else
c$$$               stop
c$$$#endif
c$$$            ENDIF
         ENDDO
      ENDDO
#ifdef CMPI
      CALL ErrorElevSum(ErrorElevExceeded)
c      write(*,*) myproc,ErrorElevExceeded
      IF (ErrorElevExceeded.NE.0) THEN
         CALL MESSAGE_FINI()
         stop
      ENDIF
#else
!      write(*,*)ErrorElevExceeded
c$$$      IF (ErrorElevExceeded.NE.0) THEN
c$$$         STOP
c$$$      ENDIF
#endif

      RETURN
      END SUBROUTINE


