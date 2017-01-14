      MODULE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
C====================================================
C====================================================
      PARAMETER (NOA=8,NCRYSTAL=4)
      PARAMETER (NP=3000,NPA=1000)
      PARAMETER (NE=2000,NEA=1000)
      PARAMETER (NF=10,NFA=10,NU=3000,NUA=10)
      PARAMETER (NQ=10,NT=3000)
      PARAMETER (NCUNIT=15000)
      PARAMETER (IIIUNIT=15000,JJJUNIT=8000)
      PARAMETER (ALG=8.0,PERCENT=0.05)
      PARAMETER (NAA=3000,NCC=3000,NAC=3000)
      PARAMETER (NCT=5000)
      PARAMETER (NATOM=NPA*NOA)
      PARAMETER (master=0)
C====================================================
C====================================================
      INTEGER MP,ME,MPA,MEA,MCUNIT,MG
      INTEGER MATOM
      INTEGER MOA,ICRYSTAL,ICASE
      INTEGER MF,MU,MFA,MUA,MTG,MQ
      INTEGER ISTRESS,MPRINT
      INTEGER MVMDC,MVMDT 
      INTEGER MTIME1,MTIME2,MTIME3,JTIME
      INTEGER IER
      INTEGER IIMAX,JJMAX
      INTEGER numid,ierr,myid,numprocs
      REAL*8  TIME,DTIME,DAMPING 
      REAL*8  CUTOFF,CUTN,CUTMAX,SIGMA,EPSLON
      REAL*8  UVOL
      REAL*8  ALATTICE
      REAL*8  XCRACK1,XCRACK2,ZCRACK1,ZCRACK2
      REAL*8  XVMD1,ZVMD1,XVMD2,ZVMD2
      REAL*8  FMAXC,FMAXA
      REAL*8  TUMASS,TGLOBAL,TDESIRE,TKAI,DTKAI
      REAL*8  BOLTZ,TAU
      REAL*8  TZERO
C====================================================
C====================================================
      DIMENSION XX(3,NP),MARKER(NP),MGROUP(NP)
      DIMENSION XA(3,NPA)
      DIMENSION DISPC(3,NP),POLARC(3,NP),FORCEC(3,NP)
      DIMENSION DISPA(3,NPA),POLARA(3,NPA),FORCEA(3,NPA)
      DIMENSION ODISPC(3,NP),ODISPA(3,NPA)
      DIMENSION XU(3,NCUNIT)
      DIMENSION JELE(NCUNIT)
      DIMENSION TSR(8,NCUNIT)
      DIMENSION ISELF(NP),WEIGHT(NP)
      DIMENSION LISTA(NPA),IAWHICH(NAA,NPA)
      DIMENSION LISTC(NP),ICWHICH(NCC,NP)
      DIMENSION LISTT(NP),ITWHICH(NCT,NP)
      DIMENSION LISTG(0:29),IGWHICH(NCC,0:29)
      DIMENSION ITGROUP(0:29),TSCALE(0:29)
      DIMENSION ALENGTH(NP)
      DIMENSION WTUNIT(IIIUNIT),IIJJ(JJJUNIT,IIIUNIT),IIUNIT(IIIUNIT)
      DIMENSION LISTAC(NPA),IACWHICH(NAC,NPA)
      DIMENSION LISTVMD(NCUNIT)
      DIMENSION AM(0:20),BM(0:20),CM(0:20),DM(0:20)
      DIMENSION AMASS(0:20),CHARGE(0:20),LTYPE(NOA)
      DIMENSION TM(NOA)
      DIMENSION V1(3),V2(3),V3(3),YY(3,NOA)
      DIMENSION ID(16),DELTA(3,3)
      DIMENSION R(NP,NOA,3),F(NP,NOA,3),U(NP,NOA,3),FCORRE(NP,NOA,3)
      DIMENSION VEL(NP,NOA,3),VPC(NP,NOA,3),ACC(NP,NOA,3)
      DIMENSION XATOM(NATOM,3),FATOM(NATOM,3),UATOM(NATOM,3)
      DIMENSION FACORRE(NATOM,3)
      DIMENSION FACLUST(NATOM,3),FCLUST(NP,NOA,3)
      DIMENSION VATOM(NATOM,3),VPA(NATOM,3),ACA(NATOM,3)
      DIMENSION IFNODE(NF),IFALPHA(NF),IFCOMPO(NF),FSCALE(NF)
      DIMENSION IQNODE(NQ),QSCALE(NQ)
      DIMENSION IFATOM(NFA),IFC(NFA),FSCALEA(NFA)
      DIMENSION IUNODE(NU),IUALPHA(NU),IUCOMPO(NU),USCALE(NU)
      DIMENSION IUATOM(NUA),IUC(NUA),USCALEA(NUA)
      DIMENSION ELEC(3), BFIE(3)
      DIMENSION STRESSC(3,3,NP),STRESSA(3,3,NPA)
      DIMENSION IPRINT(80)
      DIMENSION IJKC(8,NE),IJKA(8,NEA)
      DIMENSION TEMPT(NP),QFLUX(NP,3),TDOT(NP)
      DIMENSION TKAIP(NP),DTKAIP(NP)
      DIMENSION FTZERO(NP,NOA,3),HKIJ(3,NP,NP)
      DIMENSION SHAPE(8,9),GRADE(3,8,9),AJCOB(9,NE),CMX(3,8,9,NE)
      DIMENSION ELEE(3,NP),VOLT(NP)
      DIMENSION TKAIG(0:29),DTKAIG(0:29),TEMPG(0:29)
C====================================================
C====================================================
      END MODULE SHAREDATA
C====================================================
C====================================================
C====================================================
      PROGRAM FEAFT
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
      integer mstatus(MPI_STATUS_SIZE)
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
C====================================================
C====================================================
      OPEN (5, FILE ='infile-first',STATUS='UNKNOWN',FORM='FORMATTED') 
      if(myid.eq.master) then
      OPEN (6, FILE ='outfile',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (7, FILE ='vmd.xyz',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (8, FILE ='tecplot.dat',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (9, FILE ='relax',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (10, FILE ='tset',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (11, FILE ='infile-second',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (12, FILE ='infile-relax',STATUS='UNKNOWN',FORM='FORMATTED')
      OPEN (13, FILE ='tcheck',STATUS='UNKNOWN',FORM='FORMATTED')
      end if 
C====================================================
C====================================================
      IER=0
      CALL INPUT
      CALL START
      CALL INITIAL1
      DO 1000 JTIME=1,MTIME1
      CALL CENTRAL1
      CALL OUTRELAX
1000  CONTINUE
      CALL RELAY1
      CALL RELAY1
      TDESIRE=TZERO
      TBEGIN=TDESIRE*2.0D0
      CALL INITIAL2
      CALL RANDOM(TBEGIN)
      DO 2000 JTIME=1,MTIME2
      CALL CENTRAL2
      CALL OUTCHECK
2000  CONTINUE
      CALL RELAY2
      CALL RELAY2
      CALL INITIAL3
      DO 3000 JTIME=1,MTIME3
      CALL CENTRAL3
      CALL OUTPUT
3000  CONTINUE
9000  CONTINUE
C====================================================
C====================================================
      CLOSE(5)
      if(myid.eq.master) then
      CLOSE(6)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      end if
      STOP
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE INPUT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
      CALL CONTROL
      CALL CRACK
      CALL VMDZONE
      CALL DATABASE
      IF(MP.NE.0) CALL NODE
      IF(MPA.NE.0) CALL ATOMS
      IF(MP.NE.0) CALL UNITCELL
      CALL BC
      CALL RUNTIME
      CALL PRINTING
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE START
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
      IF(ID(6).NE.2) CALL ULIST
      IF(ID(6).NE.2) CALL CLIST
      IF(ID(6).EQ.3) CALL ACLIST
      IF(ID(6).NE.1) CALL ALIST
      CALL TLIST
      CALL SHAPEFN
      TZERO=300.0D0
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE RELAY1
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
1     FORMAT(3E15.8)
2     FORMAT(/10X,'I am after Relay1')
C====================================================
C====================================================
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      DO 100 IC=1,3
      if(myid.eq.master) WRITE(12,1) R(IP,IA,IC),
     *VEL(IP,IA,IC),ACC(IP,IA,IC)
100   CONTINUE 
      if(myid.eq.master) WRITE(6,2)
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE RELAY2
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
1     FORMAT(3E15.8)
2     FORMAT(/10X,'I am after Relay2')
C====================================================
C====================================================
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      DO 100 IC=1,3
      if(myid.eq.master) WRITE(11,1) R(IP,IA,IC),
     *VEL(IP,IA,IC),ACC(IP,IA,IC)
100   CONTINUE 
      if(myid.eq.master) WRITE(6,2)
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CONTROL
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
1     FORMAT(10I8)
2     FORMAT(10X,'MP     =',I8/
     *10X,'ME     =',I8/
     *10X,'MCUNIT =',I8/
     *10X,'MPA    =',I8/
     *10X,'MEA    =',I8/
     *10X,'MG    =',I8/)
3     FORMAT(10X,'ICRYSTAL     =',I8/
     *10X,'ICASE     =',I8/)
4     FORMAT(F15.8)
5     FORMAT(10X, 'Damping    =',F15.8)
C====================================================
C====================================================
      READ(5,1) (ID(I),I=1,10)
C=========================================================C
C                                                         C
C     ID(1)=  1 : Tecplot includes atoms                  C
C             2 : Tecplot does not includes atoms         C
C                                                         C
C     ID(2)=  1 : There is no external E field            C
C             2 : Uniform external E field                C
C             3 : any given external E field              C
C                                                         C
C     ID(3)=  1 : There is no external B field            C
C             2 : Uniform external B field                C
C             3 : any given external B field              C
C                                                         C
C     ID(4)=  1 : Doesn't matter whether there is a crack C
C                 It is not considered as a barrier       C
C             2 : There is a crack on x-z plane           C
C                 It is considered as a barrier           C
C                                                         C
C     ID(5)=  1 : there is no damping                     C
C             2 : there is a damping force                C
C                                                         C
C     ID(6)=  1 : Continua only                           C
C             2 : Atoms only                              C
C             3 : Concurrent atoms/continua               C
C                                                         C
C     ID(7)=  1 : simple Nose-Hoover Thermostat           C
C             2 : modified Nose-Hoover Thermostat         C
C             3 : pair-wise thermostat                    C
C                                                         C
C                                                         C
C     ID(8)=                                              C
C                                                         C
C     ID(9)=                                              C
C                                                         C
C     ID(10)=                                             C
C                                                         C
C=========================================================C
      if(myid.eq.master) WRITE(6,1) (ID(I),I=1,10)
      DELTA=0.0D0
      FORALL(I=1:3)
        DELTA(I,I)=1.0D0
      END FORALL
C====================================================
C====================================================
      READ(5,1) MP,ME,MCUNIT,MPA,MEA,MG
      if(myid.eq.master) WRITE(6,2)MP,ME,MCUNIT,MPA,MEA,MG
      READ(5,1) ICRYSTAL,ICASE
      if(myid.eq.master) WRITE(6,3)ICRYSTAL,ICASE 
      READ(5,4) DAMPING
      if(myid.eq.master) WRITE(6,5)DAMPING
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE NODE
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(I8,3F15.8,2I8)
2     FORMAT(10X/10X,'Node coordinate in continuum region'/)
5     FORMAT(10X/10X,'Connectivity of element in continuum region'/)
6     FORMAT(9I8)
C====================================================
C====================================================
      if(myid.eq.master) WRITE(6,2)
      DO 100 I=1,MP
      READ(5,1) L,XX(1,I),XX(2,I),XX(3,I),MARKER(I),MGROUP(I)
      if(myid.eq.master) 
     *WRITE(6,1) I,(XX(K,I),K=1,3),MARKER(I),MGROUP(I)
100   CONTINUE 
      if(myid.eq.master) WRITE(6,5)
      DO 400 I=1,ME
      READ(5,6) L,(IJKC(K,I),K=1,8)
      if(myid.eq.master) WRITE(6,6) I,(IJKC(K,I),K=1,8)
400   CONTINUE 
C====================================================
      LISTG=0
      DO 500 I=1,MP
      IG=MGROUP(I)
      LISTG(IG)=LISTG(IG)+1
      IGWHICH(LISTG(IG),IG)=I
500   CONTINUE
C====================================================
      DO 300 IP=1,MP
      DO 300 IA=1,MOA
      DO 300 IC=1,3
      R(IP,IA,IC)=XX(IC,IP)+YY(IC,IA)
300   CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE UNITCELL
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(10X,'WEIGHT(',I5,')=',F15.8)
2     FORMAT(10X,'TOTAL NUMBER OF UNIT CELLS IN C =',F15.8)
6     FORMAT(I8,3F15.8)
7     FORMAT(10X,'Unitcell coordinates in continuum region')
C====================================================
C====================================================
      DIMENSION RL(3),TSHAPE(8)
C====================================================
C====================================================
      WEIGHT=0.0
      if(myid.eq.master) WRITE(6,7)
      DO 100 I=1,MCUNIT
      READ(5,6) L,XU(1,I),XU(2,I),XU(3,I)
      if(myid.eq.master) WRITE(6,6) I,(XU(K,I),K=1,3)
100   CONTINUE 
C====================================================
C====================================================
      DO 200 IUNIT=1,MCUNIT
      RL=XU(:,IUNIT)
      CALL LOCATION(RL,TSHAPE,JE)
      IF(JE.LT.1.OR.JE.GT.ME) GO TO 200
      JELE(IUNIT)=JE
      DO 40 IG=1,8
      JG=IJKC(IG,JE)
      WEIGHT(JG)=WEIGHT(JG)+TSHAPE(IG)
      TSR(IG,IUNIT)=TSHAPE(IG)
40    CONTINUE
200   CONTINUE
C====================================================
C====================================================
      DO 400 IP=1,MP
      if(myid.eq.master) WRITE(6,1) IP,WEIGHT(IP)
      TWEIGHT=TWEIGHT+WEIGHT(IP)
400   CONTINUE
      if(myid.eq.master) WRITE(6,2) TWEIGHT
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ACLIST
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
5     FORMAT(15X,'WRONG LISTAC(',I5,') =',I5, '>',I5/)
6     FORMAT(I8,3F15.8)
C====================================================
      DO 600 I=1,MPA
      LISTAC(I)=0
      DO 650 IUNIT=1,MCUNIT
      DX=XA(1,I)-XU(1,IUNIT)  
      DY=XA(2,I)-XU(2,IUNIT)  
      DZ=XA(3,I)-XU(3,IUNIT)  
      D=DSQRT(DX*DX+DY*DY+DZ*DZ)
      IF (D.GT.CUTMAX) GO TO 650
C====================================================
C====================================================
      IYES=1
      IF(ID(4).EQ.2) CALL BARRIER(XA(1,I),XA(3,I),
     *XU(1,IUNIT),XU(3,IUNIT),
     *XCRACK1,ZCRACK1,XCRACK2,ZCRACK2,IYES)
      IF(IYES.EQ.0) GO TO 650
C====================================================
C====================================================
      LISTAC(I)=LISTAC(I)+1
      IACWHICH(LISTAC(I),I)=IUNIT
650   CONTINUE
      IF(LISTAC(I).GT.NAC) IER=1
      IF(IER.EQ.1.and.myid.eq.master)
     *WRITE(6,5) I,LISTAC(I),NAC 
C     WRITE(6,5) I,LISTAC(I),NAC 
      IF(IER.EQ.1) GO TO 1000
600   CONTINUE
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ATOMS
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(I8,3F15.8)
2     FORMAT(5X,I8,'  XA =',3F15.8)
4     FORMAT(15X,'Unitcell coordinates in atomic region'/)
5     FORMAT(15X,'Connectivity of element in atomic region'/)
6     FORMAT(9I8)
C====================================================
C====================================================
      if(myid.eq.master) WRITE(6,4)
      DO 100 I=1,MPA
      READ(5,1) K,(XA(L,I),L=1,3)
      if(myid.eq.master) WRITE(6,2) I,(XA(L,I),L=1,3)
100   CONTINUE
      if(myid.eq.master) WRITE(6,5)
      DO 110 I=1,MEA
      READ(5,6) K,(IJKA(L,I),L=1,8)
      if(myid.eq.master) WRITE(6,6) I,(IJKA(L,I),L=1,8)
110   CONTINUE
C====================================================
C====================================================
      MATOM=MPA*MOA
C====================================================
C====================================================
      K=0
      DO 150 I=1,MPA
      DO 150 J=1,MOA
      K=K+1
      XATOM(K,1)=XA(1,I)+YY(1,J)
      XATOM(K,2)=XA(2,I)+YY(2,J)
      XATOM(K,3)=XA(3,I)+YY(3,J)
150   CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ALIST
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(15X,'WRONG LISTA(',I5,') =',I5, '>',I5/)
C====================================================
C====================================================
      DO 200 I=1,MPA
      LISTA(I)=0
      DO 300 J=I,MPA
      DX=XA(1,I)-XA(1,J)
      DY=XA(2,I)-XA(2,J)
      DZ=XA(3,I)-XA(3,J)
      D=DSQRT(DX*DX+DY*DY+DZ*DZ)
      IF(D.GT.CUTMAX) GO TO 300
C====================================================
C====================================================
      IYES=1
      IF(ID(4).EQ.2) CALL BARRIER(XA(1,I),XA(3,I),
     *XA(1,J),XA(3,J),
     *XCRACK1,ZCRACK1,XCRACK2,ZCRACK2,IYES)
      IF(IYES.EQ.0) GO TO 300
C====================================================
C====================================================
      LISTA(I)=LISTA(I)+1
      IAWHICH(LISTA(I),I)=J
300   CONTINUE
      IF(LISTA(I).GT.NAA) IER=1
      IF(IER.EQ.1) THEN
        if(myid.eq.master) WRITE(6,1) I,LISTA(I),NAA
      END IF 
C       WRITE(6,1) I,LISTA(I),NAA
      IF(IER.EQ.1) GO TO 9999
200   CONTINUE
9999  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CCSTRESS
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RKC(3),RLC(3),RRC(3)
C====================================================
C====================================================
C====================================================
      FACTOR=0.5D0/UVOL
C====================================================
C====================================================
      DO 2000 IP=1,MP
      ITEST=0
      DO 1000 II=1,LISTC(IP)
      IUNIT=ICWHICH(II,IP)
      IF(ISELF(IP).EQ.IUNIT) ITEST=1
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RLC=0.0D0
C====================================================
C====================================================
      DO 100 I=1,8
      IIE=JELE(IUNIT)
      JP=IJKC(I,IIE)
      RLC(1)=RLC(1)+TSR(I,IUNIT)*R(JP,IBETA,1)
      RLC(2)=RLC(2)+TSR(I,IUNIT)*R(JP,IBETA,2)
      RLC(3)=RLC(3)+TSR(I,IUNIT)*R(JP,IBETA,3)
100   CONTINUE
C====================================================
C====================================================
      DO 1000 IALPHA=1,MOA
      IF(ITEST.EQ.1.AND.IALPHA.EQ.IBETA) GO TO 1000
      RKC(1)=R(IP,IALPHA,1)
      RKC(2)=R(IP,IALPHA,2)
      RKC(3)=R(IP,IALPHA,3)
      RRC=RKC-RLC
      DC=DSQRT(RRC(1)*RRC(1)+RRC(2)*RRC(2)+RRC(3)*RRC(3))
      IF(DC.GT.CUTN) GO TO 1000
      CALL AKF(AK,DC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DC)
C====================================================
C====================================================
      FAK=FACTOR*AK
      DO 250 IC=1,3
      DO 250 JC=1,3
      STRESSC(IC,JC,IP)=STRESSC(IC,JC,IP)+
     *FAK*RRC(IC)*RRC(JC)
250   CONTINUE
C====================================================
C====================================================
1000  CONTINUE
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AFORCE(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRA(3),FFA(3)
C====================================================
C====================================================
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
C====================================================
C====================================================
      DO 1000 II=1,LISTA(I)
      J=IAWHICH(II,I)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RRA=0.0D0
      IF(I.EQ.J.AND.IALPHA.EQ.IBETA) GO TO 1000
      JATOM=(J-1)*MOA+IBETA
      RRA(1)=XATOM(IATOM,1)-XATOM(JATOM,1)
      RRA(2)=XATOM(IATOM,2)-XATOM(JATOM,2)
      RRA(3)=XATOM(IATOM,3)-XATOM(JATOM,3)
      DA=DSQRT(RRA(1)*RRA(1)+RRA(2)*RRA(2)+RRA(3)*RRA(3))
      IF(DA.GT.CUTN) GO TO 1000
      CALL AKF(AK,DA,IALPHA,IBETA)
C     CALL AKFLJ(AK,DA)
      FFA=-AK*RRA
      FATOM(IATOM,:)=FATOM(IATOM,:)+FFA
      FATOM(JATOM,:)=FATOM(JATOM,:)-FFA
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AASTRESS
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRA(3),FFA(3)
C====================================================
C====================================================
      FACTOR=0.5D0/UVOL
      DO 3000 IATOM=1,MATOM
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
C====================================================
C====================================================
      DO 1000 II=1,LISTA(I)
      J=IAWHICH(II,I)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RRA=0.0D0
      IF(I.EQ.J.AND.IALPHA.EQ.IBETA) GO TO 1000
      JATOM=(J-1)*MOA+IBETA
      RRA(1)=XATOM(IATOM,1)-XATOM(JATOM,1)
      RRA(2)=XATOM(IATOM,2)-XATOM(JATOM,2)
      RRA(3)=XATOM(IATOM,3)-XATOM(JATOM,3)
      DA=DSQRT(RRA(1)*RRA(1)+RRA(2)*RRA(2)+RRA(3)*RRA(3))
      IF(DA.GT.CUTN) GO TO 1000
      CALL AKF(AK,DA,IALPHA,IBETA)
C     CALL AKFLJ(AK,DA)
      FAK=FACTOR*AK
C====================================================
C====================================================
      DO 2000 IC=1,3
      DO 2000 JC=1,3
      STRESSA(IC,JC,I)=STRESSA(IC,JC,I)+FAK*RRA(IC)*RRA(JC)
2000  CONTINUE
1000  CONTINUE
3000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ACFORCE(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RLAC(3),RKAC(3),RRAC(3),FFAC(3)
C====================================================
C====================================================
      RKAC=0.0D0
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
      RKAC=XATOM(IATOM,:)
C====================================================
C====================================================
      DO 2000 II=1,LISTAC(I)
      IUNIT=IACWHICH(II,I)
      IE=JELE(IUNIT)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RLAC=0.0D0
C====================================================
C====================================================
      DO 100 J=1,8
      IP=IJKC(J,IE)
      RLAC(1)=RLAC(1)+TSR(J,IUNIT)*R(IP,IBETA,1)
      RLAC(2)=RLAC(2)+TSR(J,IUNIT)*R(IP,IBETA,2)
      RLAC(3)=RLAC(3)+TSR(J,IUNIT)*R(IP,IBETA,3)
100   CONTINUE
C====================================================
C====================================================
      RRAC=RKAC-RLAC
      DAC=DSQRT(RRAC(1)*RRAC(1)+RRAC(2)*RRAC(2)+RRAC(3)*RRAC(3))
      IF(DAC.GT.CUTN) GO TO 1000
      CALL AKF(AK,DAC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DAC)
      FFAC=-AK*RRAC
C====================================================
C====================================================
      DO 200 IC=1,3
      FATOM(IATOM,IC)=FATOM(IATOM,IC)+FFAC(IC)
      DO 300 LAMDA=1,8
      JP=IJKC(LAMDA,IE)
      F(JP,IBETA,IC)=F(JP,IBETA,IC)-
     *FFAC(IC)*TSR(LAMDA,IUNIT)
300   CONTINUE
200   CONTINUE
C====================================================
C====================================================
1000  CONTINUE
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ACSTRESS
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RLAC(3),RKAC(3),RRAC(3),FFAC(3)
C====================================================
C====================================================
      FACTOR=0.5D0/UVOL
      DO 4000 IATOM=1,MATOM
      RKAC=0.0D0
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
      RKAC=XATOM(IATOM,:)
C====================================================
C====================================================
      DO 2000 II=1,LISTAC(I)
      IUNIT=IACWHICH(II,I)
      IE=JELE(IUNIT)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RLAC=0.0D0
C====================================================
C====================================================
      DO 100 J=1,8
      IP=IJKC(J,IE)
      RLAC(1)=RLAC(1)+TSR(J,IUNIT)*R(IP,IBETA,1)
      RLAC(2)=RLAC(2)+TSR(J,IUNIT)*R(IP,IBETA,2)
      RLAC(3)=RLAC(3)+TSR(J,IUNIT)*R(IP,IBETA,3)
100   CONTINUE
C====================================================
C====================================================
      RRAC=RKAC-RLAC
      DAC=DSQRT(RRAC(1)*RRAC(1)+RRAC(2)*RRAC(2)+RRAC(3)*RRAC(3))
      IF(DAC.GT.CUTN) GO TO 1000
      CALL AKF(AK,DAC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DAC)
C====================================================
C====================================================
      FAK=FACTOR*AK
      DO 3000 IC=1,3
      DO 3000 JC=1,3
      STRESSA(IC,JC,I)=STRESSA(IC,JC,I)+FAK*RRAC(IC)*RRAC(JC)
3000  CONTINUE
1000  CONTINUE
2000  CONTINUE
4000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE DATABASE
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      BOLTZ=3.166815D-6
C====================================================
      V1=0.0D0
      V2=0.0D0
      V3=0.0D0
      YY=0.0D0
      AM=0.0D0
      BM=0.0D0
      CM=0.0D0
      DM=0.0D0
      SIGMA=0.0D0
      EPSLON=0.0D0
      AMASS=0.0D0
      CHARGE=0.0D0
      LTYPE=30
      GO TO (1000,2000,3000,4000),ICRYSTAL
1000  CONTINUE
      GO TO (1001,1002,1003),ICASE
1001  CONTINUE
      AA=4.2*1.889726
      ALATTICE=AA
      UVOL=ALATTICE**3
      V1(1)=AA
      V2(2)=AA
      V3(3)=AA
      YY(1,2)=0.5*AA
      YY(2,2)=0.5*AA
      YY(1,3)=0.5*AA
      YY(3,3)=0.5*AA
      YY(2,4)=0.5*AA
      YY(3,4)=0.5*AA
      YY(1,5)=0.5*AA
      YY(2,6)=0.5*AA
      YY(3,7)=0.5*AA
      YY(1,8)=0.5*AA
      YY(2,8)=0.5*AA
      YY(3,8)=0.5*AA
      CUTOFF=12.0*1.889726
      CUTN=4.0D0*CUTOFF
      CUTMAX=1.2D0*CUTN
      AM(0)=9547.96*0.0367493
      BM(0)=0.21916*1.889726
      CM(0)=32.32*0.0367493*1.889726**6
      DM(0)=32.32*0.0367493*1.889726**12
      AM(1)=1284.38*0.0367493
      BM(1)=0.2997*1.889726
      CM(1)=0.0D0
      DM(1)=0.0D0
      MOA=8
      LTYPE(1)=1
      LTYPE(2)=1
      LTYPE(3)=1
      LTYPE(4)=1
      LTYPE(5)=0
      LTYPE(6)=0
      LTYPE(7)=0
      LTYPE(8)=0
      AMASS(0)=15.9994*1822.8885
      AMASS(1)=24.305*1822.8885
      CHARGE(0)=-2
      CHARGE(1)=2
      GO TO 9000
C====================================================
C====================================================
1002  CONTINUE
      GO TO 9000
C====================================================
C====================================================
1003  CONTINUE
      GO TO 9000
C====================================================
C====================================================
2000  CONTINUE
      GO TO (2001,2002,2003),ICASE
C====================================================
C====================================================
2001  CONTINUE
      AA=3.943*1.889726
      ALATTICE=AA
      UVOL=ALATTICE**3
      V1(1)=AA
      V2(2)=AA
      V3(3)=AA
      YY(1,2)=0.5*AA
      YY(2,2)=0.5*AA
      YY(3,2)=0.5*AA
      YY(1,3)=0.5*AA
      YY(2,3)=0.5*AA
      YY(1,4)=0.5*AA
      YY(3,4)=0.5*AA
      YY(2,5)=0.5*AA
      YY(3,5)=0.5*AA
      CUTOFF=12.0*1.889726
      CUTN=4.0D0*CUTOFF
      CUTMAX=1.2D0*CUTN
      AM(0)=25.410*0.0367493
      BM(0)=0.6937*1.889726
      CM(0)=32.32*0.0367493*1.889726**6
      DM(0)=32.32*0.0367493*1.889726**12
      AM(2)=4818.416*0.0367493
      BM(2)=0.3067*1.889726
      CM(2)=0.0D0
      DM(2)=0.0D0
      AM(3)=4545.823*0.0367493
      BM(3)=0.261*1.889726
      CM(3)=0.0D0
      DM(3)=0.0D0
      MOA=8
      LTYPE(1)=2
      LTYPE(2)=3
      LTYPE(3)=0
      LTYPE(4)=0
      LTYPE(5)=0
      AMASS(0)=15.9994*1822.8885
      AMASS(2)=137.327*1822.8885
      AMASS(3)=47.867*1822.8885
      CHARGE(0)=-2
      CHARGE(2)=2
      CHARGE(3)=4
      GO TO 9000
C====================================================
C====================================================
2002  CONTINUE
      GO TO 9000
C====================================================
C====================================================
2003  CONTINUE
      GO TO 9000
C====================================================
C====================================================
3000  CONTINUE
      GO TO (3001,3002,3003),ICASE
C====================================================
C====================================================
3001  CONTINUE
      AA=3.614*1.889726
      MOA=4
      ALATTICE=AA
      UVOL=ALATTICE**3
      V1(1)=AA
      V2(2)=AA
      V3(3)=AA
      YY(1,2)=0.5*AA
      YY(2,2)=0.5*AA
      YY(1,3)=0.5*AA
      YY(3,3)=0.5*AA
      YY(2,4)=0.5*AA
      YY(3,4)=0.5*AA
      SIGMA=2.277*1.889726
      EPSLON=0.415*0.0367493
      CUTOFF=4*SIGMA
      CUTN=1.0D0*CUTOFF
      CUTMAX=1.5D0*CUTN
      LTYPE(1)=4
      LTYPE(2)=4
      LTYPE(3)=4
      LTYPE(4)=4
      AMASS(4)=63.546*1822.8885
      GO TO 9000
C====================================================
C====================================================
3002  CONTINUE
      GO TO 9000
C====================================================
C====================================================
3003  CONTINUE
      GO TO 9000
C====================================================
C====================================================
4000  CONTINUE
      GO TO 9000
C====================================================
C====================================================
9000  CONTINUE
C====================================================
C====================================================
      TOTAL=0.0D0
      DO 200 IA=1,MOA
      KIA=LTYPE(IA)
      TOTAL=TOTAL+AMASS(KIA)
200   CONTINUE
      DO 250 IA=1,MOA
      KIA=LTYPE(IA)
      TM(IA)=AMASS(KIA)/TOTAL
250   CONTINUE
      TUMASS=TOTAL
C====================================================
C====================================================
C     CALL CENTROID
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AKF(AK,DD,IALPHA,IBETA)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      KIA=LTYPE(IALPHA)
      KIB=LTYPE(IBETA)
      ANAB=CHARGE(KIA)*CHARGE(KIB)
      AK=-ANAB/DD/DD/DD
      IF(DD.GT.CUTOFF) GO TO 1000
      IF(KIA.EQ.0.OR.KIB.EQ.0) THEN
          ITYPE=KIA+KIB
          A=AM(ITYPE)
          B=BM(ITYPE)
          C=CM(ITYPE)
          D=DM(ITYPE)
          FACTOR=DEXP(-DD/B)
C         AK=AK-A*FACTOR/B/DD+6.0D0*C/DD**8-12.0D0*D/DD**13
          AK=AK-A*FACTOR/B/DD+6.0D0*C/DD**8-12.0D0*D/DD**14
      END IF
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AKFLJ(AK,D)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
      IF(D.GT.CUTN) THEN
      AK=0.0D0
      ELSE
      AK=24*EPSLON*SIGMA**6*D**(-7)*(1-2*SIGMA**6*D**(-6))/D
      END IF
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE POTENTIAL(EPAB,DD,IALPHA,IBETA)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      KIA=LTYPE(IALPHA)
      KIB=LTYPE(IBETA)
      ANAB=CHARGE(KIA)*CHARGE(KIB)
      EPAB=ANAB/DD
      IF(DD.GT.CUTOFF) GO TO 1000
      IF(KIA.EQ.0.OR.KIB.EQ.0) THEN
          ITYPE=KIA+KIB
          A=AM(ITYPE)
          B=BM(ITYPE)
          C=CM(ITYPE)
          D=DM(ITYPE)
          FACTOR=DEXP(-DD/B)
          EPAB=EPAB+A*FACTOR-C/DD**6+D/DD**12
      END IF
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CENTRAL1
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      FACTOR=1.0D0+DAMPING*DTIME*0.5
      VPC=VEL+ACC*DTIME*0.5
C====================================================
C====================================================
      R=R+VPC*DTIME
      U=U+VPC*DTIME
C====================================================
C====================================================
      TIME=TIME+DTIME
C====================================================
C====================================================
      CALL GETFORCE1
C====================================================
C====================================================
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      TM2=AMASS(LTYPE(IA))
      DO 100 IC=1,3
      ACC(IP,IA,IC)=F(IP,IA,IC)/TM2/WEIGHT(IP)
100   CONTINUE
      ACC=(ACC-DAMPING*VPC)/FACTOR
C====================================================
C====================================================
C====================================================
      VEL=VPC+ACC*DTIME*0.5
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CENTRAL2
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      VPC=VEL+ACC*DTIME*0.5
C====================================================
C====================================================
      R=R+VPC*DTIME
      U=U+VPC*DTIME
C====================================================
C====================================================
      CALL TEMPC(VPC)
      DTKAI=(TGLOBAL-TDESIRE)/TDESIRE/TAU/TAU
C====================================================
C====================================================
      TIME=TIME+DTIME
      TKAI=TKAI+DTIME*DTKAI
      FACTOR=1.0D0+TKAI*DTIME*0.5D0
C====================================================
C====================================================
      CALL GETFORCE2
C====================================================
C====================================================
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      TM2=AMASS(LTYPE(IA))
      DO 100 IC=1,3
      ACC(IP,IA,IC)=F(IP,IA,IC)/TM2/WEIGHT(IP)
100   CONTINUE
      ACC=(ACC-TKAI*VPC)/FACTOR
C====================================================
C====================================================
      VEL=VPC+ACC*DTIME*0.5
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CENTRAL3
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION FACTOR(NP),VBAR(3)
C====================================================
C====================================================
      VPC=VEL+ACC*DTIME*0.5
C====================================================
      ISTRESS=0
      DO 10 I=1,MPRINT
      J=IPRINT(I)
      IF(J.EQ.JTIME) ISTRESS=1
10    CONTINUE
C====================================================
      L=2
      CALL AMPU(L)
C====================================================
C====================================================
      R=R+VPC*DTIME
      U=U+VPC*DTIME
      L=1
      CALL AMPU(L)
C====================================================
C====================================================
      CALL GETEMP(VPC)
      DO 100 I=1,MTG
      IG=ITGROUP(I)
      DTKAIG(IG)=(TEMPG(IG)-TSCALE(IG))/TSCALE(IG)/TAU/TAU
C============================================
C============================================
100   CONTINUE
C====================================================
C====================================================
      TIME=TIME+DTIME
      TKAIG=TKAIG+DTIME*DTKAIG
C====================================================
C====================================================
      CALL GETFORCE3
C====================================================
C====================================================
      DO 200 IP=1,MP
      DO 200 IA=1,MOA
      TM2=AMASS(LTYPE(IA))
      DO 200 IC=1,3
      ACC(IP,IA,IC)=F(IP,IA,IC)/TM2/WEIGHT(IP)
200   CONTINUE
C====================================================
C====================================================
C====================================================
      CALL NOSE
C====================================================
C====================================================
C====================================================
      VEL=VPC+ACC*DTIME*0.5
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE TEMPC(VTEMP)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
      DIMENSION VKALPHA(3)
      DIMENSION VTEMP(NP,NOA,3)
C====================================================
C====================================================
      TGLOBAL=0.0D0
      DO 300 IU=1,MCUNIT
      IE=JELE(IU)
      DO 300 IA=1,MOA
      VKALPHA=0.0D0
      DO 400 IJ=1,8
      IP=IJKC(IJ,IE)
      VKALPHA(1)=VKALPHA(1)+TSR(IJ,IU)*VTEMP(IP,IA,1)
      VKALPHA(2)=VKALPHA(2)+TSR(IJ,IU)*VTEMP(IP,IA,2)
      VKALPHA(3)=VKALPHA(3)+TSR(IJ,IU)*VTEMP(IP,IA,3)
400   CONTINUE
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      DO 300 IC=1,3
      TGLOBAL=TGLOBAL+TM2*VKALPHA(IC)**2/3
300   CONTINUE
      TGLOBAL=TGLOBAL/MOA/MCUNIT/BOLTZ
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE OUTRELAX 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
      DIMENSION RLC(3),RKC(3),RRC(3)
C====================================================
C     EK=0.0D0
C     EP=0.0D0
C     ET=0.0D0
C     DO 100 IP=1,MP
C     DO 100 IA=1,MOA
C     LL=LTYPE(IA)
C     TM2=AMASS(LL)
C     DO 100 IC=1,3
C     EK=EK+0.5D0*WEIGHT(IP)*TM2*VEL(IP,IA,IC)**2
C100   CONTINUE
C====================================================
C====================================================
C====================================================
C====================================================
C     DO 2000 IP=1,MP
C     ITEST=0
C     DO 1000 II=1,LISTC(IP)
C     IUNIT=ICWHICH(II,IP)
C     IF(ISELF(IP).EQ.IUNIT) ITEST=1
C====================================================
C====================================================
C     DO 1000 IBETA=1,MOA
C     RLC=0.0D0
C     CALL WHEREATOM(IUNIT,IBETA,RLC)
C====================================================
C====================================================
C     DO 1000 IALPHA=1,MOA
C     IF(ITEST.EQ.1.AND.IALPHA.EQ.IBETA) GO TO 1000
C     RKC(1)=R(IP,IALPHA,1)
C     RKC(2)=R(IP,IALPHA,2)
C     RKC(3)=R(IP,IALPHA,3)
C     RRC=RKC-RLC
C     DC=DSQRT(RRC(1)*RRC(1)+RRC(2)*RRC(2)+RRC(3)*RRC(3))
C     IF(DC.GT.CUTN) GO TO 1000
C     CALL POTENTIAL(EPAB,DC,IALPHA,IBETA)
C     EP=EP+EPAB*WEIGHT(IP)*0.5
C1000  CONTINUE
C2000  CONTINUE
C     ET=EK+EP
      FMAX=0.0D0
      DO 3000 IP=1,MP
      DO 3000 IA=1,MOA
      DO 3000 IC=1,3
      FFF=DABS(F(IP,IA,IC))
      IF(FFF.GT.FMAX) FMAX=FFF
3000  CONTINUE
C====================================================
C====================================================
      IF(myid.eq.master) WRITE(9,1) JTIME,FMAX
1     FORMAT(5X,I6,9(2X,E15.8))
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE OUTCHECK 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
      DIMENSION RLC(3),RKC(3),RRC(3)
C====================================================
C     EK=0.0D0
C     EP=0.0D0
C     ET=0.0D0
C     DO 100 IP=1,MP
C     DO 100 IA=1,MOA
C     LL=LTYPE(IA)
C     TM2=AMASS(LL)
C     DO 100 IC=1,3
C     EK=EK+0.5D0*WEIGHT(IP)*TM2*VEL(IP,IA,IC)**2
C100   CONTINUE
C====================================================
C====================================================
C====================================================
C====================================================
C     DO 2000 IP=1,MP
C     ITEST=0
C     DO 1000 II=1,LISTC(IP)
C     IUNIT=ICWHICH(II,IP)
C     IF(ISELF(IP).EQ.IUNIT) ITEST=1
C====================================================
C====================================================
C     DO 1000 IBETA=1,MOA
C     RLC=0.0D0
C     CALL WHEREATOM(IUNIT,IBETA,RLC)
C====================================================
C====================================================
C     DO 1000 IALPHA=1,MOA
C     IF(ITEST.EQ.1.AND.IALPHA.EQ.IBETA) GO TO 1000
C     RKC(1)=R(IP,IALPHA,1)
C     RKC(2)=R(IP,IALPHA,2)
C     RKC(3)=R(IP,IALPHA,3)
C     RRC=RKC-RLC
C     DC=DSQRT(RRC(1)*RRC(1)+RRC(2)*RRC(2)+RRC(3)*RRC(3))
C     IF(DC.GT.CUTN) GO TO 1000
C     CALL POTENTIAL(EPAB,DC,IALPHA,IBETA)
C     EP=EP+EPAB*WEIGHT(IP)*0.5
C1000  CONTINUE
C2000  CONTINUE
C     ET=EK+EP
      IF(myid.eq.master) WRITE(10,1) JTIME,TGLOBAL,DTKAI,TKAI
1     FORMAT(5X,I6,9(2X,E15.8))
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE GETFORCE1
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION BUFF1(1000000)
      DIMENSION BUFF2(1000000)
C====================================================
C====================================================
      F=0.0
      FATOM=0.0
C====================================================
C====================================================
      GO TO (101,102,103),ID(6)
101   CONTINUE
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
      GO TO 150
C====================================================
C====================================================
102   CONTINUE
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
      GO TO 150
C====================================================
C====================================================
103   CONTINUE
      IATOM1=1+myid*MATOM/numprocs
      IATOM2=(myid+1)*MATOM/numprocs
      DO JA=IATOM1,IATOM2
      CALL ACFORCE(JA)
      END DO
C====================================================
C====================================================
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
C====================================================
C====================================================
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
C====================================================
C====================================================
150   CONTINUE
C====================================================
C====================================================
C    sum up the CFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              BUFF1(I)=F(IP,IALPHA,IC)
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              F(IP,IALPHA,IC)=BUFF2(I)
            END DO
          END DO
        END DO
      END IF
C====================================================
C    sum up the AFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          BUFF1(I)=FATOM(IATOM,IC)
        END DO
      END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          FATOM(IATOM,IC)=BUFF2(I)
        END DO
      END DO
      END IF
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE GETFORCE2
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION BUFF1(1000000)
      DIMENSION BUFF2(1000000)
C====================================================
C====================================================
      F=0.0
      FATOM=0.0
C====================================================
C====================================================
      GO TO (101,102,103),ID(6)
101   CONTINUE
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
      GO TO 150
C====================================================
C====================================================
102   CONTINUE
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
      GO TO 150
C====================================================
C====================================================
103   CONTINUE
      IATOM1=1+myid*MATOM/numprocs
      IATOM2=(myid+1)*MATOM/numprocs
      DO JA=IATOM1,IATOM2
      CALL ACFORCE(JA)
      END DO
C====================================================
C====================================================
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
C====================================================
C====================================================
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
C====================================================
C====================================================
150   CONTINUE
C====================================================
C====================================================
C    sum up the CFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              BUFF1(I)=F(IP,IALPHA,IC)
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              F(IP,IALPHA,IC)=BUFF2(I)
            END DO
          END DO
        END DO
      END IF
C====================================================
C    sum up the AFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          BUFF1(I)=FATOM(IATOM,IC)
        END DO
      END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          FATOM(IATOM,IC)=BUFF2(I)
        END DO
      END DO
      END IF
C====================================================
C====================================================
C     F=F+FTZERO
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE GETENERGY
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION VELUC(NP,3)
      DIMENSION TTT(NP)
C====================================================
C====================================================
      VELUC=0.0D0
      TTT=0.0D0
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      DO 100 IC=1,3
      VELUC(IP,IC)=VELUC(IP,IC)+TM(IA)*VEL(IP,IA,IC)
100   CONTINUE
C====================================================
C====================================================
      CALL HFLUX(VEL)
C====================================================
C====================================================
      DO 3000 IE=1,ME
      DO 3000 IG=1,8
      DO 3000 II=1,8
      I=IJKC(II,IE)
      DO 3000 IC=1,3
      TTT(I)=TTT(I)+QFLUX(I,IC)*CMX(IC,II,IG,IE)*AJCOB(IG,IE)
3000  CONTINUE
C====================================================
C====================================================
      DO 4000 IE=1,ME
      DO 4000 IG=1,8
      DO 4000 JJ=1,8
      J=IJKC(JJ,IE)
      DO 4000 KK=1,8
      K=IJKC(KK,IE)
      DO 4000 LL=1,8
      L=IJKC(LL,IE)
      DO 4000 IC=1,3
      TTT(J)=TTT(J)-TEMPT(K)*VELUC(L,IC)*BOLTZ*CMX(IC,LL,IG,IE)*
     *SHAPE(KK,IG)*SHAPE(JJ,IG)*AJCOB(IG,IE)/UVOL
4000  CONTINUE
C====================================================
C====================================================
      CALL AMPQ(TTT)
C====================================================
C====================================================
      DO IP=1,MP
      TDOT(IP)=TTT(IP)/BOLTZ/1.5/WEIGHT(IP)
      END DO
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE GETFORCE3
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION BUFF1(1000000)
      DIMENSION BUFF2(1000000)
C====================================================
C====================================================
      FCORRE=0.0D0
      FACORRE=0.0D0
      FCLUST=0.0
      FACLUST=0.0
C     IF(JTIME.EQ.MTIME2) CALL FCO
C====================================================
C====================================================
      F=0.0
      FATOM=0.0
      STRESSC=0.0D0
      STRESSA=0.0D0 
C====================================================
C====================================================
      GO TO (101,102,103),ID(6)
101   CONTINUE
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
      GO TO 150
C====================================================
C====================================================
102   CONTINUE
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
      GO TO 150
C====================================================
C====================================================
103   CONTINUE
      IATOM1=1+myid*MATOM/numprocs
      IATOM2=(myid+1)*MATOM/numprocs
      DO JA=IATOM1,IATOM2
      CALL ACFORCE(JA)
      END DO
C====================================================
C====================================================
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFORCE(IU)
      END DO
C====================================================
C====================================================
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFORCE(JA)
      END DO
C====================================================
C====================================================
150   CONTINUE
C====================================================
C====================================================
C    sum up the CFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              BUFF1(I)=F(IP,IALPHA,IC)
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              F(IP,IALPHA,IC)=BUFF2(I)
            END DO
          END DO
        END DO
      END IF
C====================================================
C    sum up the AFORCE calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          BUFF1(I)=FATOM(IATOM,IC)
        END DO
      END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          FATOM(IATOM,IC)=BUFF2(I)
        END DO
      END DO
      END IF
C====================================================
C     Correcting force (fpair-fcluster, force due to reference T 
C====================================================
      F=F+FCORRE
      FATOM=FATOM+FACORRE
C====================================================
C====================================================
      IF(ID(2).NE.1.or.ID(3).NE.1) CALL EFORCE
C====================================================
C====================================================
      IF (ISTRESS.EQ.1) THEN
      CALL CCSTRESS
C     CALL CASTRESS
C     CALL ACSTRESS
C     CALL AASTRESS
      END IF
C====================================================
C     Force specified boundary condition
C====================================================
      CALL AMPF
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE EFORCE
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C=========================================================
C=========================================================
      DIMENSION RR(3),EEE(3),BBB(3)
C=========================================================
C=========================================================
      IF(ID(2).EQ.1.AND.ID(3).EQ.1) GO TO 9000
C=========================================================
C=========================================================
      DO 1000 IP=1,MP
      WT=WEIGHT(IP)
      DO 1000 IA=1,MOA
      KIA=LTYPE(IA)
      Q=CHARGE(KIA)
      QQ=Q/CZERO
      RR(1)=R(IP,IA,1)
      RR(2)=R(IP,IA,2)
      RR(3)=R(IP,IA,3)
      CALL EM(TIME,RR,EEE,BBB)
      F(IP,IA,1)=F(IP,IA,1)+WT*(Q*EEE(1)+
     *QQ*(VPC(IP,IA,2)*BBB(3)-VPC(IP,IA,3)*BBB(2)))
      F(IP,IA,2)=F(IP,IA,2)+WT*(Q*EEE(2)+
     *QQ*(VPC(IP,IA,3)*BBB(1)-VPC(IP,IA,1)*BBB(3)))
      F(IP,IA,3)=F(IP,IA,3)+WT*(Q*EEE(3)+
     *QQ*(VPC(IP,IA,1)*BBB(2)-VPC(IP,IA,2)*BBB(1)))
1000  CONTINUE
C=========================================================
C=========================================================
      DO 2000 I=1,MPA
      DO 2000 IA=1,MOA
      IATOM=(I-1)*MOA+IA
      RR(1)=XATOM(IATOM,1)
      RR(2)=XATOM(IATOM,2)
      RR(3)=XATOM(IATOM,3)
C=========================================================
C=========================================================
      CALL EM(TIME,RR,EEE,BBB)
      KIA=LTYPE(IA)
      Q=CHARGE(KIA)
      QQ=Q/CZERO
      FATOM(IATOM,1)=FATOM(IATOM,1)+Q*EEE(1)+
     *QQ*(VPA(IATOM,2)*BBB(3)-VPA(IATOM,3)*BBB(2))
      FATOM(IATOM,2)=FATOM(IATOM,2)+Q*EEE(2)+
     *QQ*(VPA(IATOM,3)*BBB(1)-VPA(IATOM,1)*BBB(3))
      FATOM(IATOM,3)=FATOM(IATOM,3)+Q*EEE(3)+
     *QQ*(VPA(IATOM,1)*BBB(2)-VPA(IATOM,2)*BBB(1))
2000  CONTINUE
C=========================================================
C=========================================================
9000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE EM(T,RRR,EEE,BBB)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRR(3),EEE(3),BBB(3)
C====================================================
C====================================================
      IF(ID(2).EQ.3.OR.ID(3).EQ.3) GO TO 2000
      EEE=ELEC
      BBB=BFIE
      GO TO 3000
2000  CONTINUE
C===================================================
C
C   user specified EM field
C
C===================================================
      EEE=0.0
      BBB=0.0
C===================================================
C===================================================
3000  CONTINUE 
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE INITIAL1
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
2     FORMAT(2I5,6(2X,E15.8))
3     FORMAT(/5X, 'Total initial force and FMAXC =',5E15.8/)
4     FORMAT(I10,6(2X,E15.8))
5     FORMAT(/5X, 'Total initial force and FMAXA =',5E15.8/)
C====================================================
C====================================================
      JTIME=0
      TIME=0.0D0
C====================================================
C====================================================
      F=0.0D0
      U=0.0D0
      VEL=0.0D0
      FATOM=0.0D0
      UATOM=0.0D0
      VATOM=0.0D0
C====================================================
C====================================================
      CALL GETFORCE1
C====================================================
C====================================================
      DO 300 IP=1,MP
      DO 300 IA=1,MOA
      KIA=LTYPE(IA)
      TM2=AMASS(KIA)
      DO 300 IC=1,3
      ACC(IP,IA,IC)=F(IP,IA,IC)/TM2/WEIGHT(IP)
300   CONTINUE
C====================================================
C====================================================
      DO 400 IP=1,MPA
      DO 400 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      KIA=LTYPE(IA)
      TM2=AMASS(KIA)
      DO 400 IC=1,3
      ACA(IATOM,IC)=FATOM(IATOM,IC)/TM2
400   CONTINUE
C====================================================
C====================================================
      CALL RIGHT
C====================================================
C====================================================
C     CALL VMDSTART
C     CALL VMDPLOTC
C     CALL VMDPLOTA
      CALL TECPSTART
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE INITIAL2
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C====================================================
C====================================================
      JTIME=0
      TIME=0.0D0
      U=0.0D0
      UATOM=0.0D0
      TKAI=0.0
C====================================================
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE INITIAL3
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(10/10X,'I am in INITIAL3'/)
C====================================================
C====================================================
      JTIME=0
      TIME=0.0D0
      if(myid.eq.master) WRITE(6,1)
C====================================================
C====================================================
      F=0.0D0
      U=0.0D0
      FATOM=0.0D0
      UATOM=0.0D0
      CALL GETEMP(VEL)
      CALL GETFORCE3
C====================================================
C====================================================
      DO 300 IP=1,MP
      DO 300 IA=1,MOA
      KIA=LTYPE(IA)
      TM2=AMASS(KIA)
      DO 310 IC=1,3
      ACC(IP,IA,IC)=F(IP,IA,IC)/TM2/WEIGHT(IP)
310   CONTINUE
300   CONTINUE
C====================================================
C====================================================
      DO 400 IP=1,MPA
      DO 400 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      KIA=LTYPE(IA)
      TM2=AMASS(KIA)
      DO 400 IC=1,3
      ACA(IATOM,IC)=FATOM(IATOM,IC)/TM2
400   CONTINUE
C====================================================
C====================================================
      TKAIG=0.0D0
      DTKAIG=0.0D0
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AMPF
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C
C     user specified subroutine
C
C====================================================
      FT=1.0D0
C====================================================
      DO 200 I=1,MF
      IP=IFNODE(I)
      IA=IFALPHA(I)
      IC=IFCOMPO(I)
      F(IP,IA,IC)=F(IP,IA,IC)+FSCALE(I)*FT
200   CONTINUE
      DO 300 I=1,MFA
      IA=IFATOM(I)
      IC=IFC(I)
      FATOM(IA,IC)=FATOM(IA,IC)+FSCALEA(I)*FT
300   CONTINUE
      RETURN
      END      
C====================================================
C====================================================
C====================================================
      SUBROUTINE AMPT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C     TFACTOR=1.0D0
C====================================================
C====================================================
      DO 50 I=1,MTG
      IG=ITGROUP(I)
      TEMPG(IG)=TSCALE(I)
50    CONTINUE
C====================================================
C====================================================
      RETURN
      END      
C====================================================
C====================================================
C====================================================
      SUBROUTINE AMPQ(TT)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION TT(NP)
C====================================================
C====================================================
      QFACTOR=1.0D0
C====================================================
C====================================================
      DO 50 I=1,MQ
      IP=IQNODE(I)
      QQ=QSCALE(I)*QFACTOR
      TT(IP)=TT(IP)-QQ
50    CONTINUE
C====================================================
C====================================================
C====================================================
      RETURN
      END      
C====================================================
C====================================================
C====================================================
      SUBROUTINE AMPU(IC)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
      TFINAL=MTIME3*DTIME
      TAPPLY=40000.0
      IF(TIME.GT.TAPPLY) GO TO 2000
      NH=2
      FACT=8.0D0
      UT=FACT*DSIN(0.5*3.14159268*TIME*NH/TAPPLY)
      VT=FACT*0.5*3.14159265*NH/TAPPLY*
     *DCOS(0.5*3.14159265*TIME*NH/TAPPLY)
      GO TO 1000
2000  CONTINUE
      UT=0.0D0
      VT=0.0D0
1000  CONTINUE
C====================================================
C====================================================
      GO TO (100,120),IC
100   CONTINUE
      DO 50 I=1,MU
      IP=IUNODE(I)
      IA=IUALPHA(I)
      IC=IUCOMPO(I)
      TTT=USCALE(I)*UT
      R(IP,IA,IC)=R(IP,IA,IC)+TTT-U(IP,IA,IC)
      U(IP,IA,IC)=TTT
50    CONTINUE
      DO 60 I=1,MUA
      IA=IUATOM(I)
      IC=IUC(I)
      TTT=USCALEA(I)*UT
      XATOM(IA,IC)=XATOM(IA,IC)+TTT-UATOM(IA,IC)
      UATOM(IA,IC)=TTT
60    CONTINUE
      GO TO 9000
C====================================================
C====================================================
120   CONTINUE
      DO 150 I=1,MU
      IP=IUNODE(I)
      IA=IUALPHA(I)
      IC=IUCOMPO(I)
      TTT=USCALE(I)*VT
      VPC(IP,IA,IC)=TTT
150   CONTINUE
      DO 160 I=1,MUA
      IA=IUATOM(I)
      IC=IUC(I)
      TTT=USCALEA(I)*VT
      VPA(IA,IC)=TTT
160   CONTINUE
9000  CONTINUE
C====================================================
      RETURN
      END      
C====================================================
C====================================================
C====================================================
      SUBROUTINE OUTPUT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C     CALL CAOUT
      CALL OUTEMP
C     CALL TDEAL
      IF(ISTRESS.EQ.1) CALL CAPRINT
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CAOUT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
1     FORMAT(5X,I8,12(2X,E15.8))
C====================================================
      CALL POLAROUT
      CALL DISPOUT
      IC1=1043
      IC2=1098
      IC3=1148
      IC4=1198
      IC5=1248
      IC6=1298
      IC7=1348
      IC8=1398
      IC9=1448
      IC10=1498
      IC11=1548
      if(myid.eq.master) WRITE(15,1) JTIME,DISPC(3,IC1),
     *DISPC(3,IC2),DISPC(3,IC3),DISPC(3,IC4),DISPC(3,IC5),
     *DISPC(3,IC6),DISPC(3,IC7),DISPC(3,IC8),DISPC(3,IC9),
     *DISPC(3,IC10),DISPC(3,IC11)
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE OUTEMP
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
1     FORMAT(I6,23(2X,F12.6))
C====================================================
      if(myid.eq.master) write(13,1) JTIME,(TEMPG(L),L=1,MG)
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CAPRINT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(10X/)
2     FORMAT(5X,'==============================================')
3     FORMAT(10X/5X,
     *'JTIME =',I5,'  Time =',E15.8,' t-bar (Me Bohr^2/h-bar)'/)
4     FORMAT(I5,2X,3(2X,E15.8),'  Bohr   (U)',2X,E15.8,'     K.')
41    FORMAT(7X,3(2X,E15.8),'  Bohr   (OU)')
5     FORMAT(10X/5X,
     *'The followings are informations at lattice point :'/5X,
     *'   (Note: Hartree = h-bar^2/Me/Bohr^2)'/)
6     FORMAT(/10X,'Total Forces =',3E15.8/)
7     FORMAT(/10X,'FMAX =',F15.8,'  UMAX =',F15.8,'  PMAX =',F15.8/)
9     FORMAT(/10X,'Average Polarization =',3F15.8/)
8     FORMAT(7X,3(2X,E15.8),'  Hartree/Bohr')
11    FORMAT(7X,3(2X,E15.8),'  e Bohr')
12    FORMAT(7X,3(2X,E15.8),'  Hartree/Bohr^3')
14    FORMAT(7X,3(2X,E15.8),'  Bohr   (O-U)')
15    FORMAT(7X,3(2X,E15.8),'  Hartree/e/Bohr (electric field)')
16    FORMAT(7X,3(2X,E15.8),'  Hartree/e (voltage)')
33    FORMAT('  A',I6,1X,3(2X,E15.8),'  Bohr   (U)')
333   FORMAT('  A',I6,1X,3(2X,E15.8),'  Bohr   (X)')
44    FORMAT(3X,I6,1X,3(2X,E15.8),'  Bohr   (U)')
45    FORMAT(I5,2X,3(2X,E15.8),'  Bohr   (X)')
88    FORMAT(3X,I6,1X,3(2X,E15.8),'  Hartree/Bohr')
99    FORMAT(10X,'Something is wrong in OUTPUT1,'/)
C==========================================
C==========================================
      CALL POLAROUT
      CALL DISPOUT
      CALL VOLTAGE
      CALL TECPLOT
C==========================================
C==========================================
      if(myid.eq.master) WRITE(6,1)
      if(myid.eq.master) WRITE(6,2)
      if(myid.eq.master) WRITE(6,2)
      if(myid.eq.master) WRITE(6,1)
      if(myid.eq.master) WRITE(6,3) JTIME,TIME
C======================================================
C======================================================
      if(myid.eq.master) WRITE(6,5)
      TFX=0.0D0
      TFY=0.0D0
      TFZ=0.0D0
      DO 100 IP=1,MP
      if(myid.eq.master) then
      WRITE(6,4) IP,(DISPC(K,IP),K=1,3),TEMPT(IP)
      WRITE(6,41) (ODISPC(K,IP),K=1,3)
      WRITE(6,8) (FORCEC(K,IP),K=1,3)
      WRITE(6,11) (POLARC(K,IP),K=1,3)
      WRITE(6,12) (STRESSC(1,K,IP),K=1,3)
      WRITE(6,12) (STRESSC(2,K,IP),K=1,3)
      WRITE(6,12) (STRESSC(3,K,IP),K=1,3)
      WRITE(6,15) (ELEE(K,IP),K=1,3)
      WRITE(6,16)  VOLT(IP)
      end if
      DO 150 IA=1,MOA
C     if(myid.eq.master) WRITE(6,44) IA,(U(IP,IA,IC),IC=1,3)
150   CONTINUE
      DO 160 IA=1,MOA
C     if(myid.eq.master) WRITE(6,88) IA,(F(IP,IA,IC),IC=1,3)
160   CONTINUE
      TFX=TFX+FORCEC(1,IP)
      TFY=TFY+FORCEC(2,IP)
      TFZ=TFZ+FORCEC(3,IP)
100   CONTINUE
C======================================================
C======================================================
      DO 500 IP=1,MPA
      if(myid.eq.master) then
      WRITE(6,4) IP,(DISPA(K,IP),K=1,3)
      WRITE(6,41) (ODISPA(K,IP),K=1,3)
      WRITE(6,8) (FORCEA(K,IP),K=1,3)
      WRITE(6,11) (POLARA(K,IP),K=1,3)
      WRITE(6,12) (STRESSA(1,K,IP),K=1,3)
      WRITE(6,12) (STRESSA(2,K,IP),K=1,3)
      WRITE(6,12) (STRESSA(3,K,IP),K=1,3)
      end if
      DO 550 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      if(myid.eq.master) WRITE(6,44) IA,(UATOM(IATOM,IC),IC=1,3)
550   CONTINUE
      DO 560 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      if(myid.eq.master) WRITE(6,88) IA,(FATOM(IATOM,IC),IC=1,3)
      TFX=TFX+FATOM(IATOM,1)
      TFY=TFY+FATOM(IATOM,2)
      TFZ=TFZ+FATOM(IATOM,3)
560   CONTINUE
500   CONTINUE
C======================================================
C======================================================
600   CONTINUE
      if(myid.eq.master) WRITE(6,6) TFX,TFY,TFZ
C======================================================
C======================================================
C     CALL VMDPLOTC
C     CALL VMDPLOTA
C======================================================
C======================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE POLAROUT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      POLARC=0.0D0
      POLARA=0.0D0
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      KIA=LTYPE(IA)
      DO 100 IC=1,3
      PQ=CHARGE(KIA)*R(IP,IA,IC)/UVOL
      POLARC(IC,IP)=POLARC(IC,IP)+PQ
100   CONTINUE
C====================================================
C====================================================
      DO 200 IP=1,MPA
      DO 200 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      KIA=LTYPE(IA)
      DO 200 IC=1,3
      PQ=CHARGE(KIA)*XATOM(IATOM,IC)/UVOL
      POLARA(IC,IP)=POLARA(IC,IP)+PQ
200   CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE DISPOUT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DISPC=0.0D0
      DISPA=0.0D0
      ODISPC=0.0D0
      ODISPA=0.0D0
      FORCEC=0.0D0
      FORCEA=0.0D0
C====================================================
C====================================================
      FMAXC=0.0
      FMAXA=0.0
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      DO 100 IC=1,3
      DISPC(IC,IP)=DISPC(IC,IP)+TM(IA)*U(IP,IA,IC)
      FORCEC(IC,IP)=FORCEC(IC,IP)+F(IP,IA,IC)
      FM=DABS(F(IP,IA,IC))
      IF(FM.GT.FMAXC) FMAXC=FM
100   CONTINUE
      DO 150 IP=1,MP
      DO 150 IC=1,3
      DO 160 IA=1,MOA
      ODISPC(IC,IP)=ODISPC(IC,IP)+TM(IA)*(U(IP,IA,IC)-
     *DISPC(IC,IP))**2
160   CONTINUE
      ODISPC(IC,IP)=DSQRT(ODISPC(IC,IP))
150   CONTINUE
C====================================================
C====================================================
      DO 200 IP=1,MPA
      DO 200 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      DO 200 IC=1,3
      DISPA(IC,IP)=DISPA(IC,IP)+TM(IA)*UATOM(IATOM,IC)
      FORCEA(IC,IP)=FORCEA(IC,IP)+FATOM(IATOM,IC)
      FM=DABS(FATOM(IATOM,IC))
      IF(FM.GT.FMAXA) FMAXA=FM
200   CONTINUE
C====================================================
      DO 250 IP=1,MPA
      DO 250 IC=1,3
      DO 260 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      ODISPA(IC,IP)=ODISPA(IC,IP)+TM(IA)*(UATOM(IATOM,IC)-
     *DISPA(IC,IP))**2
260   CONTINUE
      ODISPA(IC,IP)=DSQRT(ODISPA(IC,IP))
250   CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CASTRESS
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RK(3),RL(3),RRR(3)
C====================================================
C====================================================
C====================================================
      FACTOR=0.5D0/UVOL
C====================================================
C====================================================
      DO 2000 JP=1,MPA
      RL=0.0D0
      DO 2000 IBETA=1,MOA
      IATOM=(JP-1)*MOA+IBETA
      RL(1)=XATOM(IATOM,1)
      RL(2)=XATOM(IATOM,2)
      RL(3)=XATOM(IATOM,3)
C====================================================
C====================================================
      DO 1000 IP=1,MP
      DO 1000 IALPHA=1,MOA
      RK=0.0D0
C====================================================
C====================================================
      RK(1)=R(IP,IALPHA,1)
      RK(2)=R(IP,IALPHA,2)
      RK(3)=R(IP,IALPHA,3)
C====================================================
C====================================================
      RRR=RK-RL
      D=DSQRT(RRR(1)*RRR(1)+RRR(2)*RRR(2)+RRR(3)*RRR(3))
      IF(D.GT.CUTN) GO TO 1000
C====================================================
C====================================================
      IYES=1
      IF(ID(4).EQ.2) CALL BARRIER(RK(1),RK(3),
     *RL(1),RL(3),
     *XCRACK1,ZCRACK1,XCRACK2,ZCRACK2,IYES)
      IF(IYES.EQ.0) GO TO 1000
C====================================================
C====================================================
      CALL AKF(AK,D,IALPHA,IBETA)
C     CALL AKFLJ(AK,D)
      FAK=FACTOR*AK
      DO 200 IC=1,3
      DO 200 JC=1,3
      STRESSC(IC,JC,IP)=STRESSC(IC,JC,IP)+
     *FAK*RRR(IC)*RRR(JC)
200   CONTINUE
C====================================================
C====================================================
1000  CONTINUE
C====================================================
C====================================================
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE BARRIER(X1,Y1,X2,Y2,X3,Y3,X4,Y4,IYES)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C============================================
C============================================
      IYES=1
      YC=Y3
      IF(Y1*Y2.GE.0.0D0) GO TO 3000
      IF(DABS(X1-X2).LT.1.0D-10) GO TO 100
      S=(Y1-Y2)/(X1-X2)
      C=Y1-S1*X1
      XC=-C/S
      GO TO 1000
100   CONTINUE
      XC=X1
1000  CONTINUE
      IF(XC.LE.X4.AND.XC.GE.X3) GO TO 2000
      GO TO 3000
2000  CONTINUE
      IYES=0
3000  CONTINUE
C============================================
C============================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CENTROID
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION XCENTER(3)
C====================================================
C====================================================
      XCENTER=0.0D0
      T1=0.0D0
      DO 100 I=1,MOA
      K=LTYPE(I)
      T1=T1+AMASS(K)
      XCENTER(1)=XCENTER(1)+AMASS(K)*YY(1,I) 
      XCENTER(2)=XCENTER(2)+AMASS(K)*YY(2,I) 
      XCENTER(3)=XCENTER(3)+AMASS(K)*YY(3,I)
100   CONTINUE
      XCENTER=XCENTER/T1
C====================================================
C====================================================
      DO 200 I=1,MOA
      YY(1,I)=YY(1,I)-XCENTER(1)
      YY(2,I)=YY(2,I)-XCENTER(2)
      YY(3,I)=YY(3,I)-XCENTER(3)
200   CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE BC
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(6I8,F15.8)
2     FORMAT(5X,I8,'  INODE, IA, IC =',3I8,'  USCALE=',F15.8)
3     FORMAT(5X,I8,'  INODE, IA, IC =',3I8,'  FSCALE=',F15.8)
4     FORMAT(10X/10X,'WRONG: MU,NU,MF,NF,MUA,NUA,MFA,NFA =',8I8/)
5     FORMAT(10X/10X,'The displacement-specified BC are:'/)
6     FORMAT(10X/10X,'The force-specified BC are:'/)
7     FORMAT(3I8,F15.8)
8     FORMAT(5X,I8,'  IATOM, IC =',2I8,'  USCALEA=',F15.8)
9     FORMAT(5X,I8,'  IATOM, IC =',2I8,'  FSCALEA=',F15.8)
10    FORMAT(10X/10X,'Boundary conditions:'/)
11    FORMAT(10X,'MU,MF,MUA,MFA,MTP,MQ:',6I8)
12    FORMAT(10X/10X,'The temperature-specified BC are:'/)
13    FORMAT(2I8,F15.8)
14    FORMAT(5X,I8,'  TNODE=',I8,'  TSCALE=',F15.8)
15    FORMAT(10X/10X,'The heatflux-specified BC are:'/)
16    FORMAT(5X,I8,'  QNODE=',I8,'  QSCALE=',F15.8)
17    FORMAT(4I8,F15.8)
C====================================================
C====================================================
      READ(5,1) MU,MF,MUA,MFA,MTG,MQ
      if(myid.eq.master) WRITE(6,10)
      if(myid.eq.master) WRITE(6,11) MU,MF,MUA,MFA,MTG,MQ
C====================================================
C====================================================
      IF(MU.NE.0.OR.MUA.NE.0) Then
      if(myid.eq.master) WRITE(6,5)
      END IF
      DO 100 I=1,MU
      READ(5,17) K,IUNODE(I),IUALPHA(I),IUCOMPO(I),USCALE(I)
      if(myid.eq.master) WRITE(6,2) I,IUNODE(I),IUALPHA(I),
     *IUCOMPO(I),USCALE(I)
100   CONTINUE
      DO 150 I=1,MUA
      READ(5,7) K,IUATOM(I),IUC(I),USCALEA(I)
      if(myid.eq.master) WRITE(6,8) I,IUATOM(I),IUC(I),USCALEA(I)
150   CONTINUE
C====================================================
C====================================================
      IF(MF.NE.0.OR.MFA.NE.0) THEN
      if(myid.eq.master) WRITE(6,6)
      END IF
      DO 200 I=1,MF
      READ(5,1) K,IFNODE(I),IFALPHA(I),IFCOMPO(I),FSCALE(I)
      if(myid.eq.master) WRITE(6,3) I,IFNODE(I),IFALPHA(I),
     *IFCOMPO(I),FSCALE(I)
200   CONTINUE
      DO 250 I=1,MFA
      READ(5,7) K,IFATOM(I),IFC(I),FSCALEA(I)
      if(myid.eq.master) WRITE(6,9) I,IFATOM(I),IFC(I),FSCALEA(I)
250   CONTINUE
C====================================================
      IF(MTG.NE.0) THEN
      if(myid.eq.master) WRITE(6,12)
      END IF
      DO 300 I=1,MTG
      READ(5,13) K,L,TEMP
      ITGROUP(I)=L
      TSCALE(L)=TEMP
      if(myid.eq.master) WRITE(6,14) I,ITGROUP(I),TSCALE(L)
300   CONTINUE
C====================================================
C====================================================
      IF(MQ.NE.0) THEN
      if(myid.eq.master) WRITE(6,15)
      END IF
      DO 350 I=1,MQ
      READ(5,13) K,IQNODE(I),QSCALE(I)
      if(myid.eq.master) WRITE(6,16) I,IQNODE(I),QSCALE(I)
350   CONTINUE
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CRACK
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(2F15.8)
2     FORMAT(10X/10X,'The crack is between '/
     *10X,2F15.8, 5X,2F15.8/)
C====================================================
C====================================================
      ZCRACK1=0.0D0
      ZCRACK2=0.0D0
      READ(5,1) XCRACK1,XCRACK2
      if(myid.eq.master) WRITE(6,2) XCRACK1,ZCRACK1,XCRACK2,ZCRACK2
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE VMDZONE
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(4F15.8)
2     FORMAT(10X/10X,'The VMDZONE is between '/
     *10X,2F15.8, 5X,2F15.8/)
C====================================================
C====================================================
      READ(5,1) XVMD1,ZVMD1,XVMD2,ZVMD2
      if(myid.eq.master) WRITE(6,2) XVMD1,ZVMD1,XVMD2,ZVMD2
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE VMDSTART
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      LISTVMD=0
      I=0
      DO 1000 IUNIT=1,MCUNIT
      X1=XU(1,IUNIT)
      Z1=XU(3,IUNIT)
      IF(X1.GT.XVMD2) GO TO 1000
      IF(X1.LT.XVMD1) GO TO 1000
      IF(Z1.GT.ZVMD2) GO TO 1000
      IF(Z1.LT.ZVMD1) GO TO 1000
      I=I+1
      LISTVMD(IUNIT)=1
1000  CONTINUE
      MVMDC=I*MOA
      MVMDT=MVMDC+MATOM
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE VMDPLOTC
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
C     DIMENSION TTT(3),SHAPE(8)
C================================================
C================================================
1     FORMAT(A8,3(2X,E15.8))
2     FORMAT(I5/5X,'atoms in single crystal')
C================================================
      DIMENSION RL(3)
C================================================
      if(myid.eq.master) WRITE(7,2) MVMDT
C================================================
C================================================
      DO 1000 IUNIT=1,MCUNIT
      IF(LISTVMD(IUNIT).EQ.0) GO TO 1000
      DO 2000 IA=1,MOA
      RL=0.0D0
      J=LTYPE(IA)
      DO 3000 K=1,8
      IIE=JELE(IUNIT)
      IP=IJKC(K,IIE)
      RL(1)=RL(1)+TSR(K,IUNIT)*R(IP,IA,1)
      RL(2)=RL(2)+TSR(K,IUNIT)*R(IP,IA,2)
      RL(3)=RL(3)+TSR(K,IUNIT)*R(IP,IA,3)
3000  CONTINUE
      if(myid.eq.master) then
      IF(J.EQ.1) WRITE(7,1) 'Mg',(RL(L),L=1,3)
      IF(J.EQ.2) WRITE(7,1) 'Ba',(RL(L),L=1,3)
      IF(J.EQ.3) WRITE(7,1) 'Ti',(RL(L),L=1,3)
      IF(J.EQ.0) WRITE(7,1) 'O', (RL(L),L=1,3)
      end if
2000  CONTINUE
1000  CONTINUE
C================================================
C================================================
      RETURN
      END
C================================================
C================================================
C================================================
      SUBROUTINE VMDPLOTA
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
1     FORMAT(A8,3(2X,E15.8))
C================================================
C================================================
      DIMENSION TTT(3)
C================================================
C================================================
      DO 1000 I=1,MPA
      DO 1000 IA=1,MOA
      IATOM=(I-1)*MOA+IA
      TTT=XATOM(IATOM,:)
      J=LTYPE(IA)
      if(myid.eq.master) then
      IF(J.EQ.1) WRITE(7,1) 'Mg',(TTT(K),K=1,3)
      IF(J.EQ.2) WRITE(7,1) 'Ba',(TTT(K),K=1,3)
      IF(J.EQ.3) WRITE(7,1) 'Ti',(TTT(K),K=1,3)
      IF(J.EQ.0) WRITE(7,1) 'O', (TTT(K),K=1,3)
      end if
1000  CONTINUE
C================================================
C================================================
      RETURN
      END
C================================================
C================================================
C================================================
      SUBROUTINE TECPSTART
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION XXX(3,NP+NPA),POLAR(3,NP+NPA),LMN(8,ME+MEA)
C====================================================
C====================================================
10    FORMAT('Title="Tecplot Results"')
11    FORMAT('Variables="X","Y","Z","DISP1","DISP2","DISP3","POLAR1",
     *"POLAR2","POLAR3","T11","T22","T33","T23","T31","T12","TEMP",
     *"ELEC1","ELEC2","ELEC3","VOLT"')
12    FORMAT('Zone T="Load step 0",N=',I5,' E=',I5,' datapacking=point,
     *zonetype=febrick')
13    FORMAT(20(2X,E15.8))
15    FORMAT(8I10)
C====================================================
C  write node information
C====================================================
      MPT=MP
      MET=ME
      IF(ID(1).EQ.1) MPT=MP+MPA
      IF(ID(1).EQ.1) MET=ME+MEA
      if(myid.eq.master) WRITE(8,10) 
      if(myid.eq.master) WRITE(8,11)
      if(myid.eq.master) WRITE(8,12) MPT,MET
      DO 100 I=1, MPT
      IF(I.LE.MP) THEN
       XXX(:,I)=XX(:,I)
       POLAR(:,I)=POLARC(:,I)
      ELSE
       XXX(:,I)=XA(:,I-MP)
       POLAR(:,I)=POLARA(:,I-MP)
      END IF
100   CONTINUE
      DO 150 I=1,MPT
      if(myid.eq.master) WRITE(8,13) XXX(1,I),XXX(2,I),XXX(3,I),
     *0.0,0.0,0.0,
     *POLAR(1,I),POLAR(2,I),POLAR(3,I),0.0,0.0,0.0,
     *0.0,0.0,0.0,TZERO,0.0,0.0,0.0,0.0
 150  CONTINUE
C====================================================
C  write element information (connectivity)
C====================================================
      DO 170 I=1,MET
      IF(I.LE.ME) THEN
        LMN(:,I)=IJKC(:,I)
      ELSE
        LMN(:,I)=IJKA(:,I-ME)+MP
      END IF
170   CONTINUE
      DO 175 I=1,MET
      if(myid.eq.master) WRITE(8,15) (LMN(J,I),J=5,8),
     *(LMN(K,I),K=1,4)
175   CONTINUE
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE TECPLOT 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
      DIMENSION POLAR(3,NP+NPA),DISP(3,NP+NPA),STRESS(3,3,NP+NPA)
C================================================
C================================================
1     FORMAT('Zone T="',I5,'-th LOAD STEP"',' N=',I5,' E=',I5,
     *'  Datapacking=point,
     *zonetype=febrick,varsharelist=([1-3]=1),
     *connectivitysharezone=1') 
2     FORMAT(17(2X,E15.8))      
C================================================
C================================================
      DO 300 IG=1,MG
      DO 300 II=1,LISTG(IG)
      IP=IGWHICH(II,IG)
      TEMPT(IP)=TEMPG(IG)
300   CONTINUE
C================================================
C================================================
      MPT=MP
      MET=ME
      IF(ID(1).EQ.1) MPT=MP+MPA
      IF(ID(1).EQ.1) MET=ME+MEA
      if(myid.eq.master) WRITE(8,1) JTIME,MPT,MET
      DO 100 I=1, MPT
      IF(I.LE.MP) THEN
       DISP(:,I)=DISPC(:,I)
       POLAR(:,I)=POLARC(:,I)
       STRESS(:,:,I)=STRESSC(:,:,I)
      ELSE
       DISP(:,I)=DISPA(:,I-MP)
       POLAR(:,I)=POLARA(:,I-MP)
       STRESS(:,:,I)=STRESSA(:,:,I-MP)
      END IF
100   CONTINUE
      DO 200 IP=1,MPT
      if(myid.eq.master) WRITE(8,2) (DISP(I,IP),I=1,3),
     *(POLAR(J,IP),J=1,3),
     *STRESS(1,1,IP),STRESS(2,2,IP),STRESS(3,3,IP),STRESS(2,3,IP),
     *STRESS(3,1,IP),STRESS(1,2,IP),TEMPT(IP),
     *ELEE(1,IP),ELEE(2,IP),ELEE(3,IP),VOLT(IP)
200   CONTINUE
C================================================
C================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE RUNTIME
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
1     FORMAT(3I8,F15.8)
2     FORMAT(5X/,
     *5X,'Time steps for relaxation=',I8/
     *5X,'Time steps for equlibrium=',I8/
     *5X,'Time steps for real problem=',I8/
     *5X,'Time step=',F15.8/
     *5X,'Tau=',F15.8/)
C================================================
C================================================
      READ(5,1) MTIME1,MTIME2,MTIME3,DTIME
      TAU=2000*DTIME
      if(myid.eq.master) WRITE(6,2) MTIME1,MTIME2,MTIME3,DTIME,TAU
C================================================
C================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE PRINTING
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
1     FORMAT(10I8)
2     FORMAT(5X,'IPRINT(',I2,')',I8)
3     FORMAT(10X/10X,'PRINTING information:'/)
C================================================
C================================================
      READ(5,1) MPRINT
      READ(5,1) (IPRINT(K),K=1,MPRINT)
      if(myid.eq.master) WRITE(6,3)
      DO 100 I=1,MPRINT
      if(myid.eq.master) WRITE(6,2) I,IPRINT(I)
100   CONTINUE
C================================================
C================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE TLIST
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
      DIMENSION RK(3)
C================================================
C================================================
      LISTT=0
C========================================
C========================================
      DO 2000 IUNIT=1,MCUNIT
      ALMIN=1.0D10
      IE=JELE(IUNIT)
      RK=XU(:,IUNIT)
C========================================
C========================================
      DO 100 I=1,8
      IP=IJKC(I,IE)
      AL=DSQRT((RK(1)-XX(1,IP))**2
     *+(RK(2)-XX(2,IP))**2
     *+(RK(3)-XX(3,IP))**2)
      IF(AL.LT.ALMIN) IMIN=IP
      IF(AL.LT.ALMIN) ALMIN=AL
100   CONTINUE
C========================================
C========================================
      LISTT(IMIN)=LISTT(IMIN)+1
      ITWHICH(LISTT(IMIN),IMIN)=IUNIT
2000  CONTINUE
C========================================
C========================================
C========================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ULIST
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
      DIMENSION RK(3),RL(3),RRR(3)
      DIMENSION KP(NCUNIT)
      DIMENSION UP(NP)
      DIMENSION TEMP(NP)
C================================================
C================================================
1     FORMAT(10X,'WRONG in ULIST!'/)
2     FORMAT(10X,I5,'    LOOK in ULIST  number of unit cell, 
     *no. of point, per=',3F10.4)
C========================================
C========================================
      KP=0
      UP=0
      TEMP=0
C========================================
C========================================
      ALENGTH=0.0D0
      IIUNIT=0
      IIJJ=0
      IIMAX=0
      JJMAX=0
      WTUNIT=0.0D0
C========================================
C========================================
C========================================
      DO 1000 IE=1,ME
      AVGLEN=0.0D0
      IP=IJKC(1,IE)
      JP=IJKC(7,IE)
      AVGLEN=AVGLEN+0.25*DSQRT(
     *(XX(1,IP)-XX(1,JP))**2+
     *(XX(2,IP)-XX(2,JP))**2+
     *(XX(3,IP)-XX(3,JP))**2)
      IP=IJKC(2,IE)
      JP=IJKC(8,IE)
      AVGLEN=AVGLEN+0.25*DSQRT(
     *(XX(1,IP)-XX(1,JP))**2+
     *(XX(2,IP)-XX(2,JP))**2+
     *(XX(3,IP)-XX(3,JP))**2)
      IP=IJKC(3,IE)
      JP=IJKC(5,IE)
      AVGLEN=AVGLEN+0.25*DSQRT(
     *(XX(1,IP)-XX(1,JP))**2+
     *(XX(2,IP)-XX(2,JP))**2+
     *(XX(3,IP)-XX(3,JP))**2)
      IP=IJKC(4,IE)
      JP=IJKC(6,IE)
      AVGLEN=AVGLEN+0.25*DSQRT(
     *(XX(1,IP)-XX(1,JP))**2+
     *(XX(2,IP)-XX(2,JP))**2+
     *(XX(3,IP)-XX(3,JP))**2)
      AVGLEN=PERCENT*AVGLEN
      ALENGTH(IE)=ALG
      IF(AVGLEN.GT.ALG) ALENGTH(IE)=AVGLEN
1000  CONTINUE
C========================================
C========================================
C========================================
      ICOUNT=0
      DO 2000 IUNIT=1,MCUNIT
      ALMIN=1.0D10
      IMIN=0
      IE=JELE(IUNIT)
      RK=XU(:,IUNIT)
C========================================
C========================================
      DO 100 I=1,8
      IP=IJKC(I,IE)
      AL=DSQRT((RK(1)-XX(1,IP))**2
     *+(RK(2)-XX(2,IP))**2
     *+(RK(3)-XX(3,IP))**2)
      IF(AL.LT.ALMIN) IMIN=IP
      IF(AL.LT.ALMIN) ALMIN=AL
100   CONTINUE
C========================================
C========================================
      IF(ALMIN.GT.ALENGTH(IE)) GO TO 2000
C========================================
C========================================
      ICOUNT=ICOUNT+1
C========================================
C========================================
      DO I=1,8
      JP=IJKC(I,IE)
      UP(JP)=UP(JP)+TSR(I,IUNIT)
      END DO
C========================================
C========================================
      IIUNIT(ICOUNT)=IUNIT
      KP(IUNIT)=IMIN
      JCOUNT=0
      DO 3000 JUNIT=1,MCUNIT
      RL=XU(:,JUNIT)
      RRR=RK-RL
      D=DSQRT(RRR(1)*RRR(1)+RRR(2)*RRR(2)+RRR(3)*RRR(3))
      IF(D.GT.CUTMAX) GO TO 3000
      IF(ID(4).EQ.1) GO TO 3500
      CALL BARRIER(RK(1),RK(3),RL(1),RL(3),
     *XCRACK1,ZCRACK1,XCRACK2,ZCRACK2,IYES)
      IF(IYES.EQ.0) GO TO 3000
3500  CONTINUE
      JCOUNT=JCOUNT+1
      IIJJ(JCOUNT,ICOUNT)=JUNIT
3000  CONTINUE
      IF(JCOUNT.GT.JJMAX) JJMAX=JCOUNT
2000  CONTINUE
      IIMAX=ICOUNT
      IF(IIMAX.GT.IIIUNIT) IER=1
      IF(JJMAX.GT.JJJUNIT) IER=1
      IF(IER.NE.0.and.myid.eq.master) WRITE(6,1)
      IF(IER.NE.0) GO TO 9999
C========================================
C========================================
C========================================
      DO 4000 IP=1,MP
      TEMP(IP)=WEIGHT(IP)/UP(IP)
      if(myid.eq.master) WRITE(6,2) IP,WEIGHT(IP),UP(IP),TEMP(IP)
4000  CONTINUE
C========================================
C========================================
C========================================
      DO 5000 I=1,IIMAX
      IUNIT=IIUNIT(I)
      IP=KP(IUNIT)
      IF(IP.EQ.0) GO TO 5000
      WTUNIT(I)=TEMP(IP)
5000  CONTINUE
9999  CONTINUE
C========================================
C========================================
C========================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CFORCE(I)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C================================================
C================================================
      DIMENSION RRR(3),FFF(3)
      DIMENSION RK(3),RL(3)
C============================================
C============================================
C============================================
      IUNIT=IIUNIT(I)
      IIE=JELE(IUNIT)
C============================================
C============================================
C============================================
      DO 2000 J=1,JJMAX
      JUNIT=IIJJ(J,I)
      IF(JUNIT.EQ.0) GO TO 2000
      JJE=JELE(JUNIT)
C============================================
C============================================
C============================================
      DO 1000 IALPHA=1,MOA
      RK=0.0D0
      DO 110 L=1,8
      IP=IJKC(L,IIE)
      RK=RK+TSR(L,IUNIT)*R(IP,IALPHA,:)
110   CONTINUE
C============================================
C============================================
C============================================
      DO 1000 IBETA=1,MOA
      RL=0.0D0
      IF(JUNIT.EQ.IUNIT.AND.IBETA.EQ.IALPHA) GO TO 1000
      DO 120 L=1,8
      IP=IJKC(L,JJE)
      RL=RL+TSR(L,JUNIT)*R(IP,IBETA,:)
120   CONTINUE
      RRR=RK-RL
      D=DSQRT(RRR(1)*RRR(1)+RRR(2)*RRR(2)+RRR(3)*RRR(3))
      IF(D.GT.CUTN) GO TO 1000
C============================================
C============================================
      CALL AKF(AK,D,IALPHA,IBETA)
C     CALL AKFLJ(AK,D)
C============================================
C============================================
      FFF=-AK*RRR/2.0D0*WTUNIT(I)
C============================================
C============================================
C============================================
      DO 100 IC=1,3
      DO 100 LAMDA=1,8
      IP=IJKC(LAMDA,IIE)
      JP=IJKC(LAMDA,JJE)
      F(IP,IALPHA,IC)=F(IP,IALPHA,IC)+
     *FFF(IC)*TSR(LAMDA,IUNIT)
      F(JP,IBETA,IC)=F(JP,IBETA,IC)-
     *FFF(IC)*TSR(LAMDA,JUNIT)
100   CONTINUE
C============================================
C============================================
C============================================
1000  CONTINUE
2000  CONTINUE
C============================================
C============================================
C============================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CLIST
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
4     FORMAT(15X,'WRONG LISTC(',I5,') =',I5, '>',I5/)
C====================================================
C====================================================
      DO 500 IP=1,MP
      LISTC(IP)=0
      DO 550 IUNIT=1,MCUNIT
      DX=XX(1,IP)-XU(1,IUNIT)
      DY=XX(2,IP)-XU(2,IUNIT)
      DZ=XX(3,IP)-XU(3,IUNIT)
      D=DSQRT(DX*DX+DY*DY+DZ*DZ)
      CUTAA=CUTMAX
      IF (D.GT.CUTAA) GO TO 550
      IF (D.LT.1.0D-6) ISELF(IP)=IUNIT
C====================================================
C====================================================
      IYES=1
      IF(ID(4).EQ.2) CALL BARRIER(XX(1,IP),XX(3,IP),
     *XU(1,IUNIT),XU(3,IUNIT),
     *XCRACK1,ZCRACK1,XCRACK2,ZCRACK2,IYES)
      IF(IYES.EQ.0) GO TO 550
C====================================================
C====================================================
      LISTC(IP)=LISTC(IP)+1
      ICWHICH(LISTC(IP),IP)=IUNIT
550   CONTINUE
      IF(LISTC(IP).GT.NCC) IER=1
      IF(IER.EQ.1) Then
      if(myid.eq.master) WRITE(6,4) IP,LISTC(IP),NCC
      end if 
      IF(IER.EQ.1) GO TO 9999
500   CONTINUE 
C====================================================
C====================================================
9999  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C============================================
C============================================
C============================================
      SUBROUTINE LOCATION(RRR,SHAPEF,JE)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRR(3),SHAPEF(8),XPOLY(4),ZPOLY(4)
C======================================
C======================================
1     FORMAT(/10X,'Something is wrong. DET =',4F15.8)
C======================================
C======================================
      EIPSILON=1.0D-5
      X=RRR(1)
      Y=RRR(2)
      Z=RRR(3)
      IEI=1
      IEF=ME
C======================================
C======================================
      DO 100 IE=IEI,IEF
      Y1=XX(2,IJKC(1,IE))-EIPSILON
      Y2=XX(2,IJKC(4,IE))+EIPSILON
      IF(Y.GT.Y2) GO TO 100
      IF(Y.LT.Y1) GO TO 100
      SS=-1.0D0+2.0D0*(Y-Y1)/(Y2-Y1)
      XPOLY(1)=XX(1,IJKC(5,IE))
      XPOLY(2)=XX(1,IJKC(6,IE))
      XPOLY(3)=XX(1,IJKC(2,IE))
      XPOLY(4)=XX(1,IJKC(1,IE))
      ZPOLY(1)=XX(3,IJKC(5,IE))
      ZPOLY(2)=XX(3,IJKC(6,IE))
      ZPOLY(3)=XX(3,IJKC(2,IE))
      ZPOLY(4)=XX(3,IJKC(1,IE))
      DO 200 I=1,4
      X1=XPOLY(I)
      Z1=ZPOLY(I)
      J=I+1
      IF(I.EQ.4) J=1
      X2=XPOLY(J)
      Z2=ZPOLY(J)
      TEST=(X1-X)*(Z2-Z)-(X2-X)*(Z1-Z)
      DTEST=DABS(TEST)
      IF(DTEST.GT.5.0D-5.AND.TEST.LT.0.0D0) GO TO 100
200   CONTINUE
      JE=IE
      GO TO 300
100   CONTINUE
      JE=0
      GO TO 9000
300   CONTINUE
C=================================
C=================================
      X1=XX(1,IJKC(5,JE))
      X2=XX(1,IJKC(6,JE))
      X3=XX(1,IJKC(2,JE))
      X4=XX(1,IJKC(1,JE))
      Z1=XX(3,IJKC(5,JE))
      Z2=XX(3,IJKC(6,JE))
      Z3=XX(3,IJKC(2,JE))
      Z4=XX(3,IJKC(1,JE))
C=================================
C=================================
      A1=X2+X3-X1-X4
      B1=X3+X4-X1-X2
      D1=X1+X3-X2-X4
      C1=4.0D0*X-(X1+X2+X3+X4)
      A2=Z2+Z3-Z1-Z4
      B2=Z3+Z4-Z1-Z2
      D2=Z1+Z3-Z2-Z4
      C2=4.0D0*Z-(Z1+Z2+Z3+Z4)
C
C
C     The governing equations are:
C
C     D1*RS+A1*R+B1*S=C1
C
C     D2*RS+A2*R+B2*S=C2
C
C
C=================================
      IF(DABS(D1).LT.1.0D-9.AND.DABS(D2).LT.1.0D-9) 
     *GO TO 1000
      IF(DABS(D1).LT.1.0D-9) GO TO 2000
      IF(DABS(D2).LT.1.0D-9) GO TO 3000
      GO TO 4000
1000  CONTINUE
      DET=A1*B2-A2*B1
      RR=(B2*C1-B1*C2)/DET
      TT=(A1*C2-A2*C1)/DET
      GO TO 5000
2000  CONTINUE
      IF(DABS(A1).LT.1.0D-9) GO TO 2500
      AA=-D2*B1
      BB=D2*C1+A1*B2-A2*B1
      CC=A2*C1-A1*C2
      CALL TWICE(AA,BB,CC,TT)
      RR=(C1-B1*TT)/A1
      GO TO 5000
2500  CONTINUE
      RR=(C2*B1-B2*C1)/(D2*C1+B1*A2)
      TT=C1/B1
      GO TO 5000
3000  CONTINUE
      IF(DABS(A2).LT.1.0D-9) GO TO 3500
      AA=-D1*B2
      BB=D1*C2+A2*B1-A1*B2
      CC=A1*C2-A2*C1
      CALL TWICE(AA,BB,CC,TT)
      RR=(C2-B2*TT)/A2
      GO TO 5000
3500  CONTINUE
      RR=(C1*B2-B1*C2)/(D1*C2+B2*A1)
      TT=C2/B2
      GO TO 5000
4000  CONTINUE
      ALPHA=A1*D2-A2*D1
      BETA=B1*D2-B2*D1
      GAMA=C1*D2-C2*D1
      IF(DABS(ALPHA).LT.1.0D-9) GO TO 4500
      AA=-D1*BETA
      BB=D1*GAMA+B1*ALPHA-A1*BETA
      CC=A1*GAMA-ALPHA*C1
      CALL TWICE(AA,BB,CC,TT)
      RR=(GAMA-BETA*TT)/ALPHA
      GO TO 5000
4500  CONTINUE
      RR=(C1*BETA-B1*GAMA)/(D1*GAMA+A1*BETA)
      TT=GAMA/BETA
5000  CONTINUE
C=================================
C=================================
      SHAPEF(1)=0.125*(1.0-RR)*(1.0-SS)*(1.0+TT)
      SHAPEF(2)=0.125*(1.0+RR)*(1.0-SS)*(1.0+TT)
      SHAPEF(3)=0.125*(1.0+RR)*(1.0+SS)*(1.0+TT)
      SHAPEF(4)=0.125*(1.0-RR)*(1.0+SS)*(1.0+TT)
      SHAPEF(5)=0.125*(1.0-RR)*(1.0-SS)*(1.0-TT)
      SHAPEF(6)=0.125*(1.0+RR)*(1.0-SS)*(1.0-TT)
      SHAPEF(7)=0.125*(1.0+RR)*(1.0+SS)*(1.0-TT)
      SHAPEF(8)=0.125*(1.0-RR)*(1.0+SS)*(1.0-TT)
C=================================
C=================================
9000  CONTINUE
      RETURN
      END
C============================================
C============================================
C============================================
      SUBROUTINE TWICE(A,B,C,X)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
      IF(DABS(A).LT.1.0D-9) GO TO 1000
      TEST=B*B-4.0*A*C
      X1=(-B+DSQRT(TEST))*0.5/A
      X2=(-B-DSQRT(TEST))*0.5/A
      TEST1=(X1-1.0)*(X1+1.0)
      TEST2=(X2-1.0)*(X2+1.0)
      X=X1
      IF(TEST1.GT.1.0D-5) X=X2
      GO TO 2000
1000  CONTINUE
      X=-C/B
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C============================================
C============================================
C============================================
      SUBROUTINE FCO
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
1     FORMAT(10X,'FCORRE(',2I5,' IC)=',3(2X,F15.8))
2     FORMAT(10X,'FACORRE(',I5,' IC)=',3(2X,F15.8))
3     FORMAT(10X,/10X,'Total FCORRE in continuum=',3(2X,F15.8))
4     FORMAT(10X,/10X,'Total FCORRE after atoms=',3(2X,F15.8))
C====================================================
C====================================================
      CALL FPAIRS
      CALL FCLUSTER
C====================================================
C====================================================
      TFX=0.0D0
      TFY=0.0D0
      TFZ=0.0D0
      DO 101 IP=1,MP
      DO 101 IA=1,MOA
      if(myid.eq.master) WRITE(6,1) IP,IA,(FCORRE(IP,IA,IC),IC=1,3)
      TFX=TFX+FCORRE(IP,IA,1)
      TFY=TFY+FCORRE(IP,IA,2)
      TFZ=TFZ+FCORRE(IP,IA,3)
101   CONTINUE
      WRITE(6,3)TFX, TFY, TFZ
      DO 102 IP=1,MPA
      DO 102 IA=1,MOA
      IATOM=(IP-1)*MOA+IA
      if(myid.eq.master) WRITE(6,2) IATOM,(FACORRE(IATOM,IC),IC=1,3)
      TFX=TFX+FACORRE(IATOM,1)
      TFY=TFY+FACORRE(IATOM,2)
      TFZ=TFZ+FACORRE(IATOM,3)
102   CONTINUE
      if(myid.eq.master) WRITE(6,4)TFX, TFY, TFZ
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE FPAIRS 
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C============================================
C============================================
      DIMENSION BUFF1(1000000)
      DIMENSION BUFF2(1000000)
C============================================
C============================================
      MU1=1+myid*MCUNIT/numprocs
      MU2=(myid+1)*MCUNIT/numprocs
      DO JA=MU1,MU2
      CALL CFPAIRS(JA)
      END DO
C====================================================
C====================================================
      IATOM1=1+myid*MATOM/numprocs
      IATOM2=(myid+1)*MATOM/numprocs
      DO IA=IATOM1,IATOM2
      CALL ACFPAIRS(IA)
      END DO
C====================================================
C====================================================
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO LA=JATOM1,JATOM2
      CALL AFPAIRS(LA)
      END DO
C====================================================
C====================================================
C    sum up the CFPAIRS calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              BUFF1(I)=FCORRE(IP,IALPHA,IC)
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              FCORRE(IP,IALPHA,IC)=BUFF2(I)
            END DO
          END DO
        END DO
      END IF
C====================================================
C    sum up the AFPAIRS calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          BUFF1(I)=FACORRE(IATOM,IC)
        END DO
      END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          FACORRE(IATOM,IC)=BUFF2(I)
        END DO
      END DO
      END IF
C====================================================
C====================================================
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE FCLUSTER
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C============================================
C============================================
      DIMENSION BUFF1(1000000)
      DIMENSION BUFF2(1000000)
C============================================
C============================================
      MU1=1+myid*IIMAX/numprocs
      MU2=(myid+1)*IIMAX/numprocs
      DO IU=MU1,MU2
      CALL CFCLUSTER(IU)
      END DO
C====================================================
C====================================================
      IATOM1=1+myid*MATOM/numprocs
      IATOM2=(myid+1)*MATOM/numprocs
      DO JA=IATOM1,IATOM2
      CALL ACFCLUSTER(JA)
      END DO
C====================================================
C====================================================
      JATOM1=1+myid*MATOM/numprocs
      JATOM2=(myid+1)*MATOM/numprocs
      DO JA=JATOM1,JATOM2
      CALL AFCLUSTER(JA)
      END DO
C====================================================
C====================================================
C    sum up the CFCLUSTER calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              BUFF1(I)=FCLUST(IP,IALPHA,IC)
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
        DO IP=1,MP
          DO IALPHA=1,MOA
            DO IC=1,3
              I=I+1
              FCLUST(IP,IALPHA,IC)=BUFF2(I)
            END DO
          END DO
        END DO
      END IF
C====================================================
C    sum up the AFCLUSTER calculated by all PEs
C====================================================
      IF(numprocs.GT.1) THEN
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          BUFF1(I)=FACLUST(IATOM,IC)
        END DO
      END DO
        CALL MPI_ALLREDUCE(BUFF1,BUFF2,I,MPI_DOUBLE_PRECISION,
     *  MPI_SUM,MPI_COMM_WORLD,ierr)
C====================================================
C====================================================
      I=0
      DO IATOM=1,MATOM
        DO IC=1,3
          I=I+1
          FACLUST(IATOM,IC)=BUFF2(I)
        END DO
      END DO
      END IF
C====================================================
C====================================================
      FCORRE=FCORRE-FCLUST
      FACORRE=FACORRE-FACLUST
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE CFPAIRS(IUNIT) 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C============================================
C============================================
1     FORMAT(10X,'FCORRE( ',2I5,' IC)='3(2X,F15.8))
2     FORMAT(10X,'FCORRE in FPAIRS')
3     FORMAT(10X,/10X,'TOTAL FCORRES =',3E15.8)
4     FORMAT(10X,/10X,'Standard Deviation  =',E15.8)
C============================================
      DIMENSION RRR(3),FF(3)
      DIMENSION RK(3),RL(3)
C============================================
C============================================
      IIE=JELE(IUNIT)
C============================================
C============================================
C============================================
      DO 2000 JUNIT=IUNIT,MCUNIT
      JJE=JELE(JUNIT)
C============================================
C============================================
C============================================
      DO 1000 IALPHA=1,MOA
      RK=0.0D0 
      CALL WHEREATOM(IUNIT,IALPHA,RK)
C============================================
C============================================
C============================================
      DO 1000 IBETA=1,MOA
      RL=0.0D0
      IF(JUNIT.EQ.IUNIT.AND.IBETA.LE.IALPHA) GO TO 1000
      CALL WHEREATOM(JUNIT,IBETA,RL)
      CALL DISTANCE(RK,RL,DCC)
C============================================
C============================================
      CALL AKF(AK,DCC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DCC)
C============================================
C============================================
      RRR=RK-RL
      FF(1)=-AK*RRR(1)
      FF(2)=-AK*RRR(2)
      FF(3)=-AK*RRR(3)
C============================================
C============================================
C============================================
      DO 100 IC=1,3
      DO 100 LAMDA=1,8
      IP=IJKC(LAMDA,IIE)
      JP=IJKC(LAMDA,JJE)
      FCORRE(IP,IALPHA,IC)=FCORRE(IP,IALPHA,IC)+
     *FF(IC)*TSR(LAMDA,IUNIT)
      FCORRE(JP,IBETA,IC)=FCORRE(JP,IBETA,IC)-
     *FF(IC)*TSR(LAMDA,JUNIT)
100   CONTINUE
C============================================
C============================================
C============================================
1000  CONTINUE
2000  CONTINUE
C============================================
C============================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ACFPAIRS(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RLAC(3),RKAC(3),RRAC(3),FFAC(3)
C====================================================
C====================================================
      RKAC=0.0D0
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
      RKAC=XATOM(IATOM,:)
C====================================================
C====================================================
      DO 2000 IUNIT=1,MCUNIT
      IE=JELE(IUNIT)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      CALL WHEREATOM(IUNIT,IBETA,RLAC)
      CALL DISTANCE(RKAC,RLAC,DAC)
C====================================================
C====================================================
      CALL AKF(AK,DAC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DAC)
      RRAC=RKAC-RLAC
      FFAC=-AK*RRAC
C====================================================
C====================================================
      DO 200 IC=1,3
      FACORRE(IATOM,IC)=FACORRE(IATOM,IC)+FFAC(IC)
      DO 300 LAMDA=1,8
      JP=IJKC(LAMDA,IE)
      FCORRE(JP,IBETA,IC)=FCORRE(JP,IBETA,IC)-
     *FFAC(IC)*TSR(LAMDA,IUNIT)
300   CONTINUE
200   CONTINUE
C====================================================
C====================================================
1000  CONTINUE
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AFPAIRS(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRA(3),FFA(3)
C====================================================
C====================================================
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
C====================================================
C====================================================
      DO 1000 J=1,MPA
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RRA=0.0D0
      IF(I.EQ.J.AND.IALPHA.EQ.IBETA) GO TO 1000
      JATOM=(J-1)*MOA+IBETA
      RRA(1)=XATOM(IATOM,1)-XATOM(JATOM,1)
      RRA(2)=XATOM(IATOM,2)-XATOM(JATOM,2)
      RRA(3)=XATOM(IATOM,3)-XATOM(JATOM,3)
      DA=DSQRT(RRA(1)*RRA(1)+RRA(2)*RRA(2)+RRA(3)*RRA(3))
      CALL AKF(AK,DA,IALPHA,IBETA)
C     CALL AKFLJ(AK,DA)
      FFA=-AK*RRA
      FACORRE(IATOM,:)=FACORRE(IATOM,:)+FFA
      FACORRE(JATOM,:)=FACORRE(JATOM,:)-FFA
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE CFCLUSTER(I) 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C============================================
C============================================
C============================================
      DIMENSION RRR(3),FF(3)
      DIMENSION RK(3),RL(3)
C============================================
C============================================
      IUNIT=IIUNIT(I)
      IIE=JELE(IUNIT)
C============================================
C============================================
C============================================
      DO 2000 J=1,JJMAX
      JUNIT=IIJJ(J,I)
      IF(JUNIT.EQ.0) GO TO 2000
      JJE=JELE(JUNIT)
C============================================
C============================================
C============================================
      DO 1000 IALPHA=1,MOA
      CALL WHEREATOM(IUNIT,IALPHA,RK)
C============================================
C============================================
C============================================
      DO 1000 IBETA=1,MOA
      IF(JUNIT.EQ.IUNIT.AND.IBETA.EQ.IALPHA) GO TO 1000
      CALL WHEREATOM(JUNIT,IBETA,RL)
      CALL DISTANCE(RK,RL,D)
C     IF(D.GT.CUTN) GO TO 1000
C============================================
C============================================
      CALL AKF(AK,D,IALPHA,IBETA)
C     CALL AKFLJ(AK,D)
C============================================
C============================================
      RRR=RK-RL
      FF(1)=-AK*RRR(1)/2.0D0*WTUNIT(I)
      FF(2)=-AK*RRR(2)/2.0D0*WTUNIT(I)
      FF(3)=-AK*RRR(3)/2.0D0*WTUNIT(I)
C============================================
C============================================
C============================================
      DO 100 IC=1,3
      DO 100 LAMDA=1,8
      IP=IJKC(LAMDA,IIE)
      JP=IJKC(LAMDA,JJE)
      FCLUST(IP,IALPHA,IC)=FCLUST(IP,IALPHA,IC)+
     *FF(IC)*TSR(LAMDA,IUNIT)
      FCLUST(JP,IBETA,IC)=FCLUST(JP,IBETA,IC)-
     *FF(IC)*TSR(LAMDA,JUNIT)
100   CONTINUE
C============================================
C============================================
C============================================
1000  CONTINUE
2000  CONTINUE
C============================================
C============================================
C============================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE ACFCLUSTER(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RLAC(3),RKAC(3),RRAC(3),FFAC(3)
C====================================================
C====================================================
      RKAC=0.0D0
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
      RKAC=XATOM(IATOM,:)
C====================================================
C====================================================
      DO 2000 II=1,LISTAC(I)
      IUNIT=IACWHICH(II,I)
      IE=JELE(IUNIT)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      CALL WHEREATOM(IUNIT,IBETA,RLAC)
      CALL DISTANCE(RKAC,RLAC,DAC)
      IF(DAC.GT.CUTN) GO TO 1000
C====================================================
C====================================================
      CALL AKF(AK,DAC,IALPHA,IBETA)
C     CALL AKFLJ(AK,DAC)
      RRAC=RKAC-RLAC
      FFAC=-AK*RRAC
C====================================================
C====================================================
      DO 200 IC=1,3
      FACLUST(IATOM,IC)=FACLUST(IATOM,IC)+FFAC(IC)
      DO 300 LAMDA=1,8
      JP=IJKC(LAMDA,IE)
      FCLUST(JP,IBETA,IC)=FCLUST(JP,IBETA,IC)-
     *FFAC(IC)*TSR(LAMDA,IUNIT)
300   CONTINUE
200   CONTINUE
C====================================================
C====================================================
1000  CONTINUE
2000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE AFCLUSTER(IATOM)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
      DIMENSION RRA(3),FFA(3)
C====================================================
C====================================================
      I=(IATOM-1)/MOA+1
      IALPHA=IATOM-(I-1)*MOA
C====================================================
C====================================================
      DO 1000 II=1,LISTA(I)
      J=IAWHICH(II,I)
C====================================================
C====================================================
      DO 1000 IBETA=1,MOA
      RRA=0.0D0
      IF(I.EQ.J.AND.IALPHA.EQ.IBETA) GO TO 1000
      JATOM=(J-1)*MOA+IBETA
      RRA(1)=XATOM(IATOM,1)-XATOM(JATOM,1)
      RRA(2)=XATOM(IATOM,2)-XATOM(JATOM,2)
      RRA(3)=XATOM(IATOM,3)-XATOM(JATOM,3)
      DA=DSQRT(RRA(1)*RRA(1)+RRA(2)*RRA(2)+RRA(3)*RRA(3))
      IF(DA.GT.CUTN) GO TO 1000
      CALL AKF(AK,DA,IALPHA,IBETA)
C     CALL AKFLJ(AK,DA)
      FFA=-AK*RRA
      FACLUST(IATOM,:)=FACLUST(IATOM,:)+FFA
      FACLUST(JATOM,:)=FACLUST(JATOM,:)-FFA
1000  CONTINUE
C====================================================
C====================================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE WHEREUNIT(IUNIT,RK)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C========================================================================
C========================================================================
      DIMENSION RK(3)
C========================================================================
C========================================================================
      RK=0.0D0
      IE=JELE(IUNIT)
      DO 100 I=1,8
      IP=IJKC(I,IE)
      RK(1)=RK(1)+TSR(I,IUNIT)*XX(1,IP)
      RK(2)=RK(2)+TSR(I,IUNIT)*XX(2,IP)
      RK(3)=RK(3)+TSR(I,IUNIT)*XX(3,IP)
100   CONTINUE
C========================================================================
C========================================================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE WHEREATOM(IUNIT,IALPHA,RK)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C========================================================================
C========================================================================
      DIMENSION RK(3)
C========================================================================
C========================================================================
      RK=0.0D0
      IE=JELE(IUNIT)
      DO 100 I=1,8
      IP=IJKC(I,IE)
      RK(1)=RK(1)+TSR(I,IUNIT)*R(IP,IALPHA,1)
      RK(2)=RK(2)+TSR(I,IUNIT)*R(IP,IALPHA,2)
      RK(3)=RK(3)+TSR(I,IUNIT)*R(IP,IALPHA,3)
100   CONTINUE
C========================================================================
C========================================================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE DISTANCE(RK,RL,D)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C========================================================================
C========================================================================
      DIMENSION RK(3),RL(3)
C========================================================================
C========================================================================
      D=0.0D0
      DO 100 I=1,3
      TT=RK(I)-RL(I)
      D=D+TT*TT
100   CONTINUE
      D=DSQRT(D)
C========================================================================
C========================================================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE RANDOM(TBEGIN)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C==================================================
C==================================================
C======================================
C======================================
C======================================
      DIMENSION JSEED(1)
      DIMENSION ITODAY(3),INOW(3)
      DIMENSION VKALPHA(3)
      EXTERNAL IDATE
C     DIMENSION RR(MATOM)
      DIMENSION RR(NP*NOA)
C     CALL IDATE(ITODAY)
      CALL ITIME(INOW)
C     JSEED(1)=1234
C     JSEED(1)=ITODAY(1)*1000000+INOW(1)*10000+INOW(2)*100+INOW(3)
      JSEED(1)=INOW(1)*10000+INOW(2)*100+INOW(3)
      CALL RANDOM_SEED(SIZE=K)
      CALL RANDOM_SEED(PUT=JSEED(1:K))
C======================================
C======================================
      BOLTZ=3.166815D-6
      DO 1000 IC=1,3
      CALL RANDOM_NUMBER(RR)
      TT1=0.0
      TT2=0.0
      I=0
      DO 100 IP=1,MP
      DO 100 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      I=I+1
      TT1=TT1+TM2*RR(I)
      TT2=TT2+TM2
100   CONTINUE
      AVG=TT1/TT2
      I=0
      DO 150 IP=1,MP
      DO 150 IA=1,MOA
      I=I+1
      RR(I)=RR(I)-AVG
150   CONTINUE
      SIGMA=0.0
      I=0
      DO 200 IP=1,MP
      DO 200 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      I=I+1
      SIGMA=SIGMA+TM2*RR(I)**2
200   CONTINUE
      FACTOR=BOLTZ*TBEGIN*MP*MOA/SIGMA
      FACTOR=DSQRT(FACTOR)
      I=0
      DO 300 IP=1,MP
      DO 300 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      I=I+1
      VEL(IP,IA,IC)=FACTOR*RR(I)
300   CONTINUE
1000  CONTINUE
C======================================
C======================================
      CALL CHECKT(TBEGIN)
C======================================
C======================================
      CALL TEMPC(VEL)
C======================================
C======================================
      IF(myid.eq.master) WRITE(6,1) TGLOBAL
1     FORMAT(10X/10X,'TGLOBAL =',F15.8/)
      FACTOR=DSQRT(TBEGIN/TGLOBAL)
      VEL=VEL*FACTOR
C======================================
C======================================
      CALL TEMPC(VEL)
      IF(myid.eq.master) WRITE(6,1) TGLOBAL
C======================================
C======================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE CHECKT(TBEGIN)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C==================================================
C==================================================
      TOTAL=0.0
      DO 100 IC=1,3
      AVG=0.0
      I=0
      DO 200 IP=1,MP
      DO 200 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      I=I+1
      TOTAL=TOTAL+TM2*VEL(IP,IA,IC)**2
      AVG=AVG+TM2*VEL(IP,IA,IC)
200   CONTINUE
      if(myid.eq.master) WRITE(6,2) IC,AVG
100   CONTINUE
      TOTAL=TOTAL/MP/MOA/3.0/BOLTZ
      if(myid.eq.master) WRITE(6,1) TOTAL,TBEGIN
1     FORMAT(//5X,'TOTAL,TDESIRE =',2E15.8/)
2     FORMAT(/5X,'AVERAGE(',I1,') =',E15.8/)
C======================================
C======================================
C======================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE HFLUX(VTEMP)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C==================================================
C==================================================
      DIMENSION VTEMP(NP,NOA,3),VKAL(3)
      DIMENSION VBAR(3),TT2(3)
C==================================================
C==================================================
C==================================================
      DO 1000 IP=1,MP
C==================================================
C==================================================
C==================================================
      VBAR=0.0D0
      TT1=0.0D0
      TT2=0.0D0
      TTT=0.0D0
C==================================================
C==================================================
      DO 100 II=1,LISTT(IP)
      IUNIT=ITWHICH(II,IP)
      DO 100 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      TT1=TT1+TM2
      VKAL=0.0D0
C====================================================
C====================================================
      DO 200 IJ=1,8
      IIE=JELE(IUNIT)
      JP=IJKC(IJ,IIE)
      VKAL(1)=VKAL(1)+TSR(IJ,IUNIT)*VTEMP(JP,IA,1)
      VKAL(2)=VKAL(2)+TSR(IJ,IUNIT)*VTEMP(JP,IA,2)
      VKAL(3)=VKAL(3)+TSR(IJ,IUNIT)*VTEMP(JP,IA,3)
200   CONTINUE
      TT2(1)=TT2(1)+TM2*VKAL(1)
      TT2(2)=TT2(2)+TM2*VKAL(2)
      TT2(3)=TT2(3)+TM2*VKAL(3)
100   CONTINUE
C==================================================
C==================================================
      VBAR=TT2/TT1
C====================================================
C====================================================
C====================================================
      DO 300 II=1,LISTT(IP)
      IUNIT=ITWHICH(II,IP)
      DO 300 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      VKAL=0.0D0
C====================================================
C====================================================
      DO 400 IJ=1,8
      IIE=JELE(IUNIT)
      JP=IJKC(IJ,IIE)
      VKAL(1)=VKAL(1)+TSR(IJ,IUNIT)*VTEMP(JP,IA,1)
      VKAL(2)=VKAL(2)+TSR(IJ,IUNIT)*VTEMP(JP,IA,2)
      VKAL(3)=VKAL(3)+TSR(IJ,IUNIT)*VTEMP(JP,IA,3)
400   CONTINUE
      TEMP=(VKAL(1)-VBAR(1))**2+(VKAL(2)-VBAR(2))**2+
     *(VKAL(3)-VBAR(3))**2
      QFLUX(IP,1)=QFLUX(IP,1)+TM2*(VKAL(1)-VBAR(1))*TEMP
      QFLUX(IP,2)=QFLUX(IP,2)+TM2*(VKAL(2)-VBAR(2))*TEMP
      QFLUX(IP,3)=QFLUX(IP,3)+TM2*(VKAL(3)-VBAR(3))*TEMP
300   CONTINUE
C==================================================
C==================================================
C==================================================
      QFLUX(IP,:)=QFLUX(IP,:)*0.5/MOA/LISTT(IP)/UVOL
1000  CONTINUE
C==================================================
C==================================================
C==================================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE FTINITIAL
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C========================================================================
C========================================================================
      DIMENSION TB(3,NP)
C============================================
C============================================
C============================================
      TB=0.0D0
C============================================
C============================================
C============================================
      DO 2000 IE=1,ME
      DO 2000 IG=1,8
      DO 2000 JJ=1,8
      J=IJKC(JJ,IE)
      DO 2000 IC=1,3
      TB(IC,J)=TB(IC,J)+CMX(IC,JJ,IG,IE)*AJCOB(IG,IE)*TZERO
2000  CONTINUE
C============================================
C============================================
C============================================
      DO 3000 IA=1,MOA
      LL=LTYPE(IA)
      FACT=TM(IA)*BOLTZ
      DO 3000 IP=1,MP
      DO 3000 IC=1,3
      FTZERO(IP,IA,IC)=FACT*TB(IC,IP)*WEIGHT(IP)
3000  CONTINUE
C============================================
C============================================
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE HMATRIX
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C============================================
C============================================
C============================================
      DO 2000 IE=1,ME
      DO 2000 IG=1,8
      DO 2000 II=1,8
      I=IJKC(II,IE)
      DO 2000 JJ=1,8
      J=IJKC(JJ,IE)
      DO 2000 IC=1,3
      HKIJ(IC,I,J)=HKIJ(IC,I,J)+
     *CMX(IC,JJ,IG,IE)*SHAPE(II,IG)*AJCOB(IG,IE)
2000  CONTINUE
C============================================
C============================================
C============================================
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE CMATRIX
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C=======================================================================
C=======================================================================
C================================================
      DIMENSION XXX(8,3),AJJ(3,3),AJI(3,3)
C================================================
C================================================
      TV=0.0D0
C================================================
C================================================
      DO 1000 IE=1,ME
C================================================
      DO 10 J=1,8
      DO 10 K=1,3
      XXX(J,K) = XX(K,IJKC(J,IE))
10    CONTINUE 
C================================================
      DO 2000 IG=1,9
C================================================
C================================================
      DO 20 I=1,3
      DO 20 J=1,3
      AJJ(I,J)=0.0D0
      DO 20 IALPHA=1,8
20    AJJ(I,J)=AJJ(I,J)+GRADE(J,IALPHA,IG)*XXX(IALPHA,I)
C================================================
      CALL INV(AJJ,AJI,D)
      AJCOB(IG,IE)=D
C================================================
C================================================
      DO 500 I=1,3
      DO 500 IALPHA=1,8
      CMX(I,IALPHA,IG,IE)=0.0D0
      DO 500 K=1,3
      CMX(I,IALPHA,IG,IE)=CMX(I,IALPHA,IG,IE)+GRADE(K,IALPHA,IG)*
     *                    AJI(K,I)
500   CONTINUE
C================================================
C================================================
2000  CONTINUE
C================================================
C================================================
1000  CONTINUE
C================================================
C================================================
C     CALL CHECK2
C================================================
C================================================
      RETURN
      END
C=========================================================================
C=========================================================================
C=========================================================================
      SUBROUTINE INV(A,AI,DET)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C=========================================================================
C=========================================================================
C-------------------------------------------------------------------------C
C       Program to compute the inverse matrix of the matrix A             C
C                                                                         C
C-------------------------------------------------------------------------C
C=========================================================================
C=========================================================================
      DIMENSION A(3,3),AI(3,3),AAI(3,3)
C=========================================================================
C=========================================================================
      DET = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
     &      + A(2,1)*(A(1,3)*A(3,2)-A(1,2)*A(3,3))
     &      + A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      AI(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      AI(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
      AI(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      AI(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
      AI(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
      AI(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
      AI(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
      AI(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
      AI(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      DO 10 I=1,3
      DO 10 J=1,3
      AI(I,J)=AI(I,J)/DET
  10  CONTINUE
      TT=DET*1.0D-8
      DO 100 I=1,3
      DO 100 J=1,3
      AAI(I,J)=0.0D0
      DO 200 K=1,3
      AAI(I,J)=AAI(I,J)+A(I,K)*AI(K,J)
200   CONTINUE
      TEST=DABS(AAI(I,J)-DELTA(I,J))
      IF(TEST.GT.TT) IER=1 
100   CONTINUE
      IF(IER.EQ.0) GO TO 300
      if(myid.eq.master) WRITE(6,1) DET,TT,JTIME
      if(myid.eq.master) WRITE(6,3) 
      DO 401 I=1,3
      if(myid.eq.master) WRITE(6,2) A(I,1),A(I,2),A(I,3)
401   CONTINUE
      if(myid.eq.master) WRITE(6,3) 
      DO 402 I=1,3
      if(myid.eq.master) WRITE(6,2) AI(I,1),AI(I,2),AI(I,3)
402   CONTINUE
      if(myid.eq.master) WRITE(6,3) 
      DO 403 I=1,3
      if(myid.eq.master) WRITE(6,2) AAI(I,1),AAI(I,2),AAI(I,3)
403   CONTINUE
      if(myid.eq.master) WRITE(6,3) 
C=========================================================================
C=========================================================================
1     FORMAT(10X/5X,
     *'Very wrong in INV! DET, TT=',2E15.8,' ITIME =',I5/)
2     FORMAT(3(5X,E15.8))
3     FORMAT(/)
C=========================================================================
C=========================================================================
300   CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE CHECK2
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C=======================================================================
C=======================================================================
C----------------------------------------------------------------------C      
C                                                                      C
C     Date checked: Jan. 28, 2007                                      C
C                                                                      C
C----------------------------------------------------------------------C      
C=======================================================================
C=======================================================================
1     FORMAT(/4(5X,F15.8))
2     FORMAT(/I5,5X,E15.8)
3     FORMAT(/' S',8(1X,F11.4))
4     FORMAT(/' C',8(1X,E11.4))
5     FORMAT(/'    Total =',E15.8)
6     FORMAT(/5X,
     *'It shows the sharing of unity and condition of consistency.'/)
C==============================================
C==============================================
      IE=(1+ME)/2
      TOTAL=0.0D0
      TEST=0.0
      DO 200 IG=1,8
      TOTAL=TOTAL+AJCOB(IG,IE)
      if(myid.eq.master) WRITE(6,2) IG,AJCOB(IG,IE)
      if(myid.eq.master) WRITE(6,3) (SHAPE(K,IG),K=1,8)
      if(myid.eq.master) WRITE(6,4) (CMX(1,K,IG,IE),K=1,8)
      if(myid.eq.master) WRITE(6,4) (CMX(2,K,IG,IE),K=1,8)
      if(myid.eq.master) WRITE(6,4) (CMX(3,K,IG,IE),K=1,8)
      T=-1.0D0
      TX=0.0D0
      TY=0.0D0
      TZ=0.0D0
      DO 100 IS=1,8
      T=T+SHAPE(IS,IG)
      TX=TX+CMX(1,IS,IG,IE)
      TY=TY+CMX(2,IS,IG,IE)
      TZ=TZ+CMX(3,IS,IG,IE)
100   CONTINUE
      TEST=TEST+T+TX+TY+TZ
      if(myid.eq.master) WRITE(6,1) T,TX,TY,TZ
200   CONTINUE
      WRITE(6,5) TOTAL
      if(myid.eq.master) WRITE(6,5) TOTAL
      if(myid.eq.master.AND.DABS(TEST).LT.1.0D-9) WRITE(6,6) 
C==============================================
C==============================================
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE SHAPEFN
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C=======================================================================
C=======================================================================
C----------------------------------------------------------------------C      
C                                                                      C
C     Date checked: Jan. 28, 2007                                      C
C                                                                      C
C----------------------------------------------------------------------C      
C----------------------------------------------------------------------C      
C     Subroutine to calculate the shape funtions and their derivatives C
C     with respect to the natural coordinates of the 8-nodes solid     C
C     element.                                                         C
C                                                                      C
C----------------------------------------------------------------------C
C     SHAPE(I,J)  :  shape functions                                   C
C     SHAPER=D(shapef)/Dr                                              C
C     SHAPES=D(shapef)/Ds                                              C
C     SHAPET=D(shapef)/Dt                                              C
C     R,S,T        :  the natural coordinates                          C
C----------------------------------------------------------------------C
C==========================================================
C==========================================================
      C = 0.57735
      C=0.5
      C=1.0D0
      DO 100 IG=1,9
      GO TO (101,102,103,104,105,106,107,108,109),IG
101   CONTINUE
      RR =  -C
      SS =  -C
      TT =   C
      GO TO 150
102   CONTINUE
      RR =   C
      SS =  -C
      TT =   C
      GO TO 150
103   CONTINUE
      RR =   C
      SS =   C
      TT =   C
      GO TO 150
104   CONTINUE
      RR =  -C
      SS =   C
      TT =   C
      GO TO 150
105   CONTINUE
      RR =  -C
      SS =  -C
      TT =  -C
      GO TO 150
106   CONTINUE
      RR =   C
      SS =  -C
      TT =  -C
      GO TO 150
107   CONTINUE
      RR =   C
      SS =   C
      TT =  -C
      GO TO 150
108   CONTINUE
      RR =  -C
      SS =   C
      TT =  -C
      GO TO 150
109   CONTINUE
      RR =   0
      SS =   0
      TT =   0
150   CONTINUE
      DO 200 I=1,8
      GO TO (201,202,203,204,205,206,207,208),I
201   CONTINUE
      II=-1
      JJ=-1
      KK=+1
      GO TO 250
202   CONTINUE
      II=+1
      JJ=-1
      KK=+1
      GO TO 250
203   CONTINUE
      II=+1
      JJ=+1
      KK=+1
      GO TO 250
204   CONTINUE
      II=-1
      JJ=+1
      KK=+1
      GO TO 250
205   CONTINUE
      II=-1
      JJ=-1
      KK=-1
      GO TO 250
206   CONTINUE
      II=+1
      JJ=-1
      KK=-1
      GO TO 250
207   CONTINUE
      II=+1
      JJ=+1
      KK=-1
      GO TO 250
208   CONTINUE
      II=-1
      JJ=+1
      KK=-1
250   CONTINUE
      SHAPE(I,IG) = 0.125*(1.+II*RR)*(1.+JJ*SS)*(1.+KK*TT)
      GRADE(1,I,IG) = 0.125*II*(1.+JJ*SS)*(1.+KK*TT)
      GRADE(2,I,IG) = 0.125*(1.+II*RR)*JJ*(1.+KK*TT)
      GRADE(3,I,IG) = 0.125*(1.+II*RR)*(1.+JJ*SS)*KK
200   CONTINUE
100   CONTINUE
C==========================================================
C==========================================================
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE RIGHT
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C=======================================================================
C=======================================================================
2     FORMAT(2I5,6(2X,E15.8))
3     FORMAT(/5X, 'Total initial force and FMAXC =',5E15.8/)
4     FORMAT(I10,6(2X,E15.8))
5     FORMAT(/5X, 'Total initial force and FMAXA =',5E15.8/)
C====================================================
C====================================================
      TX=0.0D0
      TY=0.0D0
      TZ=0.0D0
      FMAXC=0.0D0
      DO 600 IP=1,MP
      DO 600 IA=1,MOA
      DO 610 IC=1,3
      FKIA=DABS(F(IP,IA,IC))
      IF(FKIA.GT.FMAXC) FMAXC=FKIA
610   CONTINUE
      TX=TX+F(IP,IA,1)
      TY=TY+F(IP,IA,2)
      TZ=TZ+F(IP,IA,3)
      if(myid.eq.master.AND.IA.eq.1) WRITE(6,2) IP,IA,(F(IP,IA,L),L=1,3)
     *,TEMPT(IP),TDOT(IP)
      if(myid.eq.master.AND.IA.ne.1) WRITE(6,2) IP,IA,(F(IP,IA,L),L=1,3)
600   CONTINUE
      if(myid.eq.master) WRITE(6,3) TX,TY,TZ,FMAXC
C====================================================
C====================================================
      FMAXA=0.0D0
      DO 700 IATOM=1,MATOM
      DO 710 IC=1,3
      FKIA=DABS(FATOM(IATOM,IC))
      IF(FKIA.GT.FMAXA) FMAXA=FKIA
710   CONTINUE
      TX=TX+FATOM(IATOM,1)
      TY=TY+FATOM(IATOM,2)
      TZ=TZ+FATOM(IATOM,3)
      if(myid.eq.master) WRITE(6,4) IATOM,(FATOM(IATOM,L),L=1,3)
700   CONTINUE
      if(myid.eq.master) WRITE(6,5) TX,TY,TZ,FMAXA 
C====================================================
C====================================================
      RETURN
      END
C====================================================
C====================================================
C====================================================
      SUBROUTINE TEST
      USE MPI
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C====================================================
C====================================================
C====================================================
C====================================================
C     F=FTZERO
C====================================================
C     Force due to temperature variation
C====================================================
      DO 3000 JP=1,MP
      DO 3000 IA=1,MOA
      FACTOR=TM(IA)*BOLTZ
      DO 3000 IC=1,3
      DO 3000 IP=1,MP
      F(JP,IA,IC)=F(JP,IA,IC)-FACTOR*HKIJ(IC,JP,IP)*TZERO
3000  CONTINUE
C====================================================
C====================================================a
      CALL RIGHT
C====================================================
C====================================================a
      RETURN
      END
C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE VOLTAGE 
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C============================================
C============================================
      DIMENSION RK(3),RL(3),RRR(3)
      DIMENSION UP(3,MCUNIT)
C============================================
C============================================
C============================================
      UP=0.0D0
C============================================
C============================================
C     Compute Polarization in each unit cell
C============================================
C============================================
      DO 5000 IUNIT=1,MCUNIT
      IIE=JELE(IUNIT)
      DO 5000 IA=1,MOA
      KIA=LTYPE(IA)
      Q=CHARGE(KIA)
C============================================
C============================================
      RK=0.0D0
      DO 5100 I=1,8
      II=IJKC(I,IIE)
      RK(1)=RK(1)+TSR(I,IUNIT)*R(II,IA,1)
      RK(2)=RK(2)+TSR(I,IUNIT)*R(II,IA,2)
      RK(3)=RK(3)+TSR(I,IUNIT)*R(II,IA,3)
5100  CONTINUE
      UP(1,IUNIT)=UP(1,IUNIT)+Q*RK(1)
      UP(2,IUNIT)=UP(2,IUNIT)+Q*RK(2)
      UP(3,IUNIT)=UP(3,IUNIT)+Q*RK(3)
5000  CONTINUE
C============================================
C============================================
C============================================
C============================================
C============================================
      VOLT=0.0D0
      ELEE=0.0D0
      DO 2000 IP=1,MP
C============================================
C============================================
C============================================
      RK=0.0D0
      DO 2500 IA=1,MOA
      RK(1)=RK(1)+R(IP,IA,1)*TM(IA)
      RK(2)=RK(2)+R(IP,IA,2)*TM(IA)
      RK(3)=RK(3)+R(IP,IA,3)*TM(IA)
2500  CONTINUE
C============================================
C============================================
      DO 3000 JUNIT=1,MCUNIT
      JJE=JELE(JUNIT)
C============================================
C============================================
C============================================
      RL=0.0D0
      DO 3500 IALPHA=1,MOA
      DO 3500 I=1,8
      II=IJKC(I,JJE)
      RL(1)=RL(1)+TSR(I,JUNIT)*R(II,IALPHA,1)*TM(IALPHA)
      RL(2)=RL(2)+TSR(I,JUNIT)*R(II,IALPHA,2)*TM(IALPHA)
      RL(3)=RL(3)+TSR(I,JUNIT)*R(II,IALPHA,3)*TM(IALPHA)
3500  CONTINUE
C============================================
C============================================
C============================================
      RRR=RK-RL
      DIST=DSQRT(RRR(1)**2+RRR(2)**2+RRR(3)**2)
      IF(DIST.LT.1.0D0) GO TO 3000
C============================================
C============================================
C============================================
      TTE=0.0D0
      DO 3600 I=1,3
      TTE=TTE+UP(I,JUNIT)*RRR(I)
3600  CONTINUE
C============================================
      DO 3700 IC=1,3
      VOLT(IP)=VOLT(IP)+UP(IC,JUNIT)*RRR(IC)/DIST**3
      ELEE(IC,IP)=ELEE(IC,IP)-UP(IC,JUNIT)/DIST**3
      ELEE(IC,IP)=ELEE(IC,IP)+3*TTE*RRR(IC)/DIST**5
3700  CONTINUE
C============================================
C============================================
3000  CONTINUE
C============================================
C============================================
2000  CONTINUE
C============================================
C============================================
C============================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE GETEMP(VTEMP)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C==================================================
C==================================================
      DIMENSION VTEMP(NP,NOA,3),VKAL(3)
      DIMENSION VBAR(3),TT2(3)
C==================================================
C==================================================
C==================================================
      TEMPG=0.0D0
C==================================================
C==================================================
C==================================================
      DO 5000 IG=1,MG
      VBAR=0.0D0
      TT1=0.0D0
      TT2=0.0D0
      ICOUNT=0
C==================================================
C==================================================
C==================================================
      DO 3000 JJ=1,LISTG(IG)
      IP=IGWHICH(JJ,IG)
C==================================================
C==================================================
      DO 100 II=1,LISTT(IP)
      IUNIT=ITWHICH(II,IP)
      JJE=JELE(IUNIT)
C====================================================
C====================================================
      DO 10 IK=1,8
      KP=IJKC(IK,JJE)
C====================================================
C====================================================
      IF(MARKER(KP).EQ.2) GO TO 100
C====================================================
C====================================================
10    CONTINUE    
C====================================================
C====================================================
      ICOUNT=ICOUNT+1
C====================================================
C====================================================
      DO 200 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      TT1=TT1+TM2
      VKAL=0.0D0
C====================================================
C====================================================
      DO 300 IJ=1,8
      JP=IJKC(IJ,JJE)
      VKAL(1)=VKAL(1)+TSR(IJ,IUNIT)*VTEMP(JP,IA,1)
      VKAL(2)=VKAL(2)+TSR(IJ,IUNIT)*VTEMP(JP,IA,2)
      VKAL(3)=VKAL(3)+TSR(IJ,IUNIT)*VTEMP(JP,IA,3)
300   CONTINUE
C====================================================
C====================================================
C====================================================
      TT2(1)=TT2(1)+TM2*VKAL(1)
      TT2(2)=TT2(2)+TM2*VKAL(2)
      TT2(3)=TT2(3)+TM2*VKAL(3)
200   CONTINUE
C====================================================
C====================================================
C====================================================
100   CONTINUE
C====================================================
C====================================================
C====================================================
3000  CONTINUE
C==================================================
C==================================================
      VBAR=TT2/TT1
C====================================================
C====================================================
C====================================================
C====================================================
      DO 4000 JJ=1,LISTG(IG)
      IP=IGWHICH(JJ,IG)
C====================================================
C====================================================
C====================================================
      DO 500 II=1,LISTT(IP)
      IUNIT=ITWHICH(II,IP)
      JJE=JELE(IUNIT)
C====================================================
C====================================================
      DO 20 IK=1,8
      KP=IJKC(IK,JJE)
C====================================================
C====================================================
      IF(MARKER(KP).EQ.2) GO TO 500
C====================================================
C====================================================
20    CONTINUE    
C====================================================
C====================================================
C====================================================
      DO 600 IA=1,MOA
      LL=LTYPE(IA)
      TM2=AMASS(LL)
      VKAL=0.0D0
C====================================================
C====================================================
      DO 700 IJ=1,8
      JP=IJKC(IJ,JJE)
      VKAL(1)=VKAL(1)+TSR(IJ,IUNIT)*VTEMP(JP,IA,1)
      VKAL(2)=VKAL(2)+TSR(IJ,IUNIT)*VTEMP(JP,IA,2)
      VKAL(3)=VKAL(3)+TSR(IJ,IUNIT)*VTEMP(JP,IA,3)
700   CONTINUE
C====================================================
C====================================================
C====================================================
      TEMP=(VKAL(1)-VBAR(1))**2+(VKAL(2)-VBAR(2))**2+
     *(VKAL(3)-VBAR(3))**2
      TEMPG(IG)=TEMPG(IG)+TEMP*TM2
600   CONTINUE
C====================================================
C====================================================
C====================================================
500   CONTINUE
C====================================================
C====================================================
C====================================================
4000  CONTINUE
C==================================================
C==================================================
      TEMPG(IG)=TEMPG(IG)/3.0/BOLTZ/MOA/ICOUNT
C==================================================
C==================================================
5000  CONTINUE
C==================================================
C==================================================
C==================================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE NOSE
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C==================================================
C==================================================
C==================================================
      DIMENSION VBAR(3)
C==================================================
C==================================================
C==================================================
      DO 1000 I=1,MTG 
      VBAR=0.0
      TTMM=0.0
      IG=ITGROUP(I)
      MCELL=LISTG(IG)
C==================================================
C==================================================
C==================================================
      DO 200 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 200 IA=1,MOA
      TTMM=TTMM+AMASS(LTYPE(IA))
      VBAR(1)=VBAR(1)+AMASS(LTYPE(IA))*VPC(IP,IA,1)
      VBAR(2)=VBAR(2)+AMASS(LTYPE(IA))*VPC(IP,IA,2)
      VBAR(3)=VBAR(3)+AMASS(LTYPE(IA))*VPC(IP,IA,3)
200   CONTINUE
      VBAR=VBAR/TTMM
C====================================================
C====================================================
C====================================================
      GO TO (5000,6000,7000),ID(7)
C==================================================
C==================================================
C==================================================
6000  CONTINUE
      CALL HOOVER(VBAR,TTMM,IG,MCELL)
      GO TO 1000
C==================================================
C==================================================
7000  CONTINUE
C     CALL HOOVER(VBAR,TTMM,IG,MCELL)
C     need a new subroutine
      GO TO 1000
C==================================================
C==================================================
5000  CONTINUE
C==================================================
C==================================================
      DO 400 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 400 IA=1,MOA
      DO 400 IC=1,3
      ACC(IP,IA,IC)=ACC(IP,IA,IC)-
     *TKAIG(IG)*(VPC(IP,IA,IC)-VBAR(IC))
400   CONTINUE
C==================================================
C==================================================
C==================================================
1000  CONTINUE
C==================================================
C==================================================
C==================================================
      RETURN
      END
C==================================================
C==================================================
C==================================================
      SUBROUTINE HOOVER(VBAR,TTMM,IG,MCELL)
      USE SHAREDATA
      IMPLICIT REAL*8 (A-H,O-Z)
C==================================================
C==================================================
C==================================================
      DIMENSION VBAR(3),VT(MCELL,NOA,3)
      DIMENSION RC(3),RR(MCELL,NOA,3)
      DIMENSION PPP(3),HHH(3)
C==================================================
C==================================================
C==================================================
1     FORMAT(2X,I5,'-th group: Moment of Inertia =',E15.8/)
2     FORMAT(10X,'Linear  Momentum =',3(2X,E15.8)/)
3     FORMAT(10X,'Angular Momentum =',3(2X,E15.8)/)
4     FORMAT(10X,'Average Velocity =',3(2X,E15.8)/)
5     FORMAT(10X,'R of MASS CENTER =',3(2X,E15.8)/)
6     FORMAT(10X,'HH1,HH2,ERROR    =',3(2X,E15.8)/)
7     FORMAT(10X,'Number of Iterate=',I5/)
C==================================================
C==================================================
C==================================================
      RC=0.0D0
      PPP=0.0D0
      HHH=0.0D0
      AINERTIA=0.0D0
C==================================================
C==================================================
      DO 100 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 100 IA=1,MOA
      AMIA=AMASS(LTYPE(IA))
      DO 110 IC=1,3
      VT(II,IA,IC)=VPC(IP,IA,IC)-VBAR(IC)
      RC(IC)=RC(IC)+AMIA*R(IP,IA,IC)
      PPP=PPP+AMIA*VT(II,IA,IC)
110   CONTINUE
      HHH(1)=HHH(1)+AMIA*(R(IP,IA,2)*VT(II,IA,3)-
     *R(IP,IA,3)*VT(II,IA,2))
      HHH(2)=HHH(2)+AMIA*(R(IP,IA,3)*VT(II,IA,1)-
     *R(IP,IA,1)*VT(II,IA,3))
      HHH(3)=HHH(3)+AMIA*(R(IP,IA,1)*VT(II,IA,2)-
     *R(IP,IA,2)*VT(II,IA,1))
100   CONTINUE
      RC=RC/TTMM
C==================================================
C==================================================
      DO 200 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 200 IA=1,MOA
      AMIA=AMASS(LTYPE(IA))
      DO 200 IC=1,3
      RR(II,IA,IC)=R(IP,IA,IC)-RC(IC)
      AINERTIA=AINERTIA+AMIA*
     *RR(II,IA,IC)*RR(II,IA,IC)
200   CONTINUE
C==================================================
C==================================================
C     WRITE(6,1) IG,AINERTIA
C     WRITE(6,2) (PPP(N),N=1,3)
C     WRITE(6,3) (HHH(N),N=1,3)
C     WRITE(6,4) (VBAR(N),N=1,3)
C     WRITE(6,5) (RC(N),N=1,3)
C==================================================
C==================================================
      HH1=DSQRT(HHH(1)*HHH(1)+HHH(2)*HHH(2)+HHH(3)*HHH(3))
C==================================================
C==================================================
C==================================================
      ITERATE=0
1000  CONTINUE
      ITERATE=ITERATE+1
C==================================================
C==================================================
C==================================================
      DO 300 II=1,MCELL
      DO 300 IA=1,MOA
C==================================================
      VT(II,IA,1)=VT(II,IA,1)+(RR(II,IA,2)*HHH(3)-
     *RR(II,IA,3)*HHH(2))/AINERTIA
C==================================================
      VT(II,IA,2)=VT(II,IA,2)+(RR(II,IA,3)*HHH(1)-
     *RR(II,IA,1)*HHH(3))/AINERTIA
C==================================================
      VT(II,IA,3)=VT(II,IA,3)+(RR(II,IA,1)*HHH(2)-
     *RR(II,IA,2)*HHH(1))/AINERTIA
C==================================================
300   CONTINUE
C==================================================
C==================================================
C==================================================
      PPP=0.0D0
      HHH=0.0D0
      DO 400 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 400 IA=1,MOA
      AMIA=AMASS(LTYPE(IA))
      PPP(1)=PPP(1)+AMIA*VT(II,IA,1)
      PPP(2)=PPP(2)+AMIA*VT(II,IA,2)
      PPP(3)=PPP(3)+AMIA*VT(II,IA,3)
      HHH(1)=HHH(1)+AMIA*(R(IP,IA,2)*VT(II,IA,3)-
     *R(IP,IA,3)*VT(II,IA,2))
      HHH(2)=HHH(2)+AMIA*(R(IP,IA,3)*VT(II,IA,1)-
     *R(IP,IA,1)*VT(II,IA,3))
      HHH(3)=HHH(3)+AMIA*(R(IP,IA,1)*VT(II,IA,2)-
     *R(IP,IA,2)*VT(II,IA,1))
400   CONTINUE
C==================================================
C==================================================
C==================================================
      HH2=DSQRT(HHH(1)*HHH(1)+HHH(2)*HHH(2)+HHH(3)*HHH(3))
      ERROR=HH2/HH1
C     WRITE(6,2) (PPP(N),N=1,3)
C     WRITE(6,3) (HHH(N),N=1,3)
C     WRITE(6,6) HH1,HH2,ERROR
      IF(ERROR.GT.1.0D-3) GO TO 1000
C     WRITE(6,7) ITERATE
C==================================================
C==================================================
C==================================================
      DO 500 II=1,MCELL
      IP=IGWHICH(II,IG)
      DO 500 IA=1,MOA
      DO 500 IC=1,3
      ACC(IP,IA,IC)=ACC(IP,IA,IC)-
     *TKAIG(IG)*VT(II,IA,IC)
500   CONTINUE
C==================================================
C==================================================
C==================================================
      RETURN
      END
