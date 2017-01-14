      PROGRAM MESH
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     This program generates an input file for MAFT.f .                C
C                                                                      C
C                                                                      C
C     The discrete atoms come in unit cell by unit cell.               C
C                                                                      C
C     Version: March 13 2010                                           C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C======================================
      DIMENSION ID(20),LMN(80),IJK(8,2000),IJKA(8,2000)
      DIMENSION XA(3000),YA(3000),ZA(3000)
      DIMENSION X(3000),Y(3000),Z(3000),MARKER(3000),MGROUP(3000)
      DIMENSION XU(15000),YU(15000),ZU(15000)
      DIMENSION XPP(3,6000),POLAR(3,6000),LM(8,4000)
C======================================
C======================================
      OPEN (6, FILE ='infile-first',STATUS='UNKNOWN',FORM='FORMATTED') 
      OPEN (8, FILE ='tp.dat',STATUS='UNKNOWN',FORM='FORMATTED') 
C======================================
C======================================
C======================================
1     FORMAT(12I8)
2     FORMAT(4F15.8)
3     FORMAT(I8,3F15.8,2I8)
4     FORMAT(4I8,F15.8)
5     FORMAT(3I8,F15.8)
6     FORMAT(2I8,F15.8)
C======================================
C    CONTROL Input
C======================================
      ID(1)=1     
      ID(2)=1
      ID(3)=1
      ID(4)=1
      ID(5)=2
      ID(6)=1
      ID(7)=0
      ID(8)=0
      ID(9)=0
      ID(10)=0
      WRITE(6,1) (ID(K),K=1,10) 
C======================================
C======================================
C   Continue region MP  ME 
C   Continue region number of unitcell MCUNIT
C   Atomic region   MPA MEA 
C======================================
C======================================
      MP=3*3*34
      ME=2*2*33
      MCUNIT=5*5*65
      MPA=0
      MEA=0
      MG=10
      WRITE(6,1) MP,ME,MCUNIT,MPA,MEA,MG
C======================================
C======================================
C   ICYSTAL=1,2,3,4, ICASE=1,2,3,4.... 
C   1==Rocksalt  ----1==MgO
C   2==Perveskite ---1==BaTiO3
C   AA===length of crystal structure
C======================================
C======================================
      AA=4.2*1.889726
      MOA=8
      ICRYSTAL=1
      ICASE=1
      WRITE(6,1) ICRYSTAL,ICASE 
C========================================
C========================================
      DAMPING=1.5D-4
      WRITE(6,2) DAMPING
C======================================
C======================================
C   CRACK location  2 points x-direction
C======================================
C======================================
      XCRACK1=0.0
      XCRACK2=0*AA
      WRITE(6,2) XCRACK1,XCRACK2 
C======================================
C======================================
C   VMDZONE 2 points x and z directions
C======================================
C======================================
      XVMD1=0.0
      ZVMD1=0*AA
      XVMD2=2*AA
      ZVMD2=33*AA
      WRITE(6,2) XVMD1,ZVMD1,XVMD2,ZVMD2 
C======================================
C======================================
C     NODE Continue region NODE Information
C======================================
C======================================
      NX=3
      NY=3
      NZ=34
      L=0
      DO 10 IY=1,NY
      YYY=(IY-1)*AA*1
      DO 10 IZ=1,NZ
      IMARK=0
      IF(IZ.EQ.1.OR.IZ.EQ.34) IMARK=2
      IF(IZ.EQ.2.OR.IZ.EQ.33) IMARK=1
      IG=0
      IF(IZ.GE.3.AND.IZ.LE.32) IG=(IZ-3)/3+1 
      IF(IZ.LE.2) ZZZ=(IZ-1)*AA
      IF(IZ.GE.3.and.IZ.LE.33) ZZZ=(IZ-2)*AA*1+AA
      IF(IZ.EQ.34) ZZZ=33*AA
      DO 10 IX=1,NX
      L=L+1
      XXX=(IX-1)*AA*1
      X(L)=XXX
      Y(L)=YYY
      Z(L)=ZZZ
      MARKER(L)=IMARK
      MGROUP(L)=IG
10    CONTINUE
      MP=L
      DO 100 I=1,MP
      WRITE(6,3) I, X(I), Y(I), Z(I),MARKER(I),MGROUP(I)
100   CONTINUE
C======================================
C======================================
      L=0
      NFAC=102
      DO 11 IY=1,NY-1
      J=(IY-1)*NFAC
      DO 11 IZ=1,NZ-1
      DO 11 IX=1,NX-1
      L=L+1
      IJK(1,L)=J+IZ*NX+IX
      IJK(2,L)=IJK(1,L)+1
      IJK(5,L)=IJK(1,L)-NX
      IJK(6,L)=IJK(5,L)+1
      IJK(4,L)=IJK(1,L)+NFAC
      IJK(3,L)=IJK(2,L)+NFAC
      IJK(8,L)=IJK(5,L)+NFAC
      IJK(7,L)=IJK(6,L)+NFAC
11    CONTINUE
      ME=L
      DO 110 I=1,ME
      WRITE(6,1) I,(IJK(K,I),K=1,8)
110   CONTINUE
C======================================
C======================================
C     ATOMS Atomic region unitcell Information
C======================================
C======================================
      IF(MPA.EQ.0) GO TO 9999
      NXA=24
      NYA=3
      NZA=11
      L=0
      DO 20 IY=1,NYA
      YYY=(IY-1)*AA
      DO 20 IZ=1,NZA
      ZZZ=-5*AA+(IZ-1)*AA
      II=1
      IF(IZ.EQ.6) II=6
      DO 20 IX=II,NXA
      L=L+1
      XXX=(IX-1)*AA
      XA(L)=XXX
      YA(L)=YYY
      ZA(L)=ZZZ
20    CONTINUE
      MPA=L
      DO 200 I=1,MPA
      WRITE(6,3) I, XA(I), YA(I), ZA(I)
200   CONTINUE
      L=0
      NAFAC=259
      DO 21 IY=1,NYA-1
      J=(IY-1)*NAFAC
      DO 21 IZ=1,4
      DO 21 IX=1,NXA-1
      L=L+1
      IJKA(1,L)=J+IZ*NXA+IX
      IJKA(2,L)=IJKA(1,L)+1
      IJKA(5,L)=IJKA(1,L)-NXA
      IJKA(6,L)=IJKA(5,L)+1
      IJKA(4,L)=IJKA(1,L)+NAFAC
      IJKA(3,L)=IJKA(2,L)+NAFAC
      IJKA(8,L)=IJKA(5,L)+NAFAC
      IJKA(7,L)=IJKA(6,L)+NAFAC
21    CONTINUE
      DO 22 IY=1,NYA-1
      J=(IY-1)*NAFAC
      DO 22 IZ=5
      DO 22 IX=6,NXA-1
      L=L+1
      IJKA(1,L)=J+IZ*NXA+IX-5
      IJKA(2,L)=IJKA(1,L)+1
      IJKA(5,L)=IJKA(1,L)-NXA+5
      IJKA(6,L)=IJKA(5,L)+1
      IJKA(4,L)=IJKA(1,L)+NAFAC
      IJKA(3,L)=IJKA(2,L)+NAFAC
      IJKA(8,L)=IJKA(5,L)+NAFAC
      IJKA(7,L)=IJKA(6,L)+NAFAC
22    CONTINUE
      DO 23 IY=1,NYA-1
      J=(IY-1)*NAFAC
      DO 23 IZ=6
      DO 23 IX=6,NXA-1
      L=L+1
      IJKA(1,L)=J+IZ*NXA+19+IX
      IJKA(2,L)=IJKA(1,L)+1
      IJKA(5,L)=IJKA(1,L)-NXA
      IJKA(6,L)=IJKA(5,L)+1
      IJKA(4,L)=IJKA(1,L)+NAFAC
      IJKA(3,L)=IJKA(2,L)+NAFAC
      IJKA(8,L)=IJKA(5,L)+NAFAC
      IJKA(7,L)=IJKA(6,L)+NAFAC
23    CONTINUE
      DO 24 IY=1,NYA-1
      J=(IY-1)*NAFAC
      DO 24 IZ=7,10
      DO 24 IX=1,NXA-1
      L=L+1
      IJKA(1,L)=J+(IZ-7)*NXA+163+IX
      IJKA(2,L)=IJKA(1,L)+1
      IJKA(5,L)=IJKA(1,L)-NXA
      IJKA(6,L)=IJKA(5,L)+1
      IJKA(4,L)=IJKA(1,L)+NAFAC
      IJKA(3,L)=IJKA(2,L)+NAFAC
      IJKA(8,L)=IJKA(5,L)+NAFAC
      IJKA(7,L)=IJKA(6,L)+NAFAC
24    CONTINUE
      MEA=L
      DO 210 I=1,MEA
      WRITE(6,1) I,(IJKA(K,I),K=1,8)
210   CONTINUE
9999  CONTINUE
C======================================
C======================================
C     UNITCELL Continuum region unitcell Information
C======================================
C======================================
      NCUX=3
      NCUY=3
      NCUZ=34
      L=0
      DO 30 IY=1,NCUY
      YYY=(IY-1)*AA
      DO 30 IZ=1,NCUZ
      ZZZ=(IZ-1)*AA
      DO 30 IX=1,NCUX
      L=L+1
      XXX=(IX-1)*AA
      XU(L)=XXX
      YU(L)=YYY
      ZU(L)=ZZZ
30    CONTINUE
      MCUNIT=L
      DO 300 I=1,MCUNIT
      WRITE(6,3) I, XU(I), YU(I), ZU(I)
300   CONTINUE
C======================================
C======================================
C     BC boundary condition 
C     MU----No. of Continue atoms displace
C     MF----No. of Continue atoms force
C     MUA---No. of Continue atoms displace
C     MFA---No. of Continue atoms force 
C======================================
C======================================
      MU=36*24
      MF=0
      MUA=0
      MFA=0
      MTG=2
      MQ=0
      WRITE(6,1) MU, MF, MUA, MFA,MTG,MQ
C======================================
C======================================
      IF(MU.NE.0) THEN
      L=0
      DO 400 IY=1,NY
      DO 400 IZ=1,2
      DO 400 IX=1,NX
      IU=(IY-1)*NFAC+(IZ-1)*3+IX
      DO 400 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,1,0.0
      L=L+1
      WRITE(6,4) L,IU,IA,2,0.0
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
400   CONTINUE
      DO 410 IY=1,NY
      DO 410 IZ=33,34
      DO 410 IX=1,NX
      IU=(IY-1)*NFAC+(IZ-1)*3+IX
      DO 410 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,1,0.0
      L=L+1
      WRITE(6,4) L,IU,IA,2,0.0
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
410   CONTINUE
      END IF
      IF(MTG.NE.0) THEN
      WRITE(6,6) 1,1,300.0
      WRITE(6,6) 2,10,350.0
      END IF
C======================================
C======================================
      IF(MUA.NE.0) THEN
      L=0
      NLOW=1
      NHIGH=4
      DO 500 IY=1,NY
      DO 500 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 500 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
500   CONTINUE
      NLOW=1
      NHIGH=4
      DO 510 IY=1,NY
      DO 510 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 510 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
510   CONTINUE
      END IF
C======================================
C======================================
      IF(MF.NE.0) THEN
      L=0
      NLOW=4
      NHIGH=5
      DO 600 IY=1,NY
      DO 600 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 600 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
600   CONTINUE
      NLOW=1
      NHIGH=4
      DO 610 IY=1,NY
      DO 610 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 610 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
610   CONTINUE
      END IF
C======================================
C======================================
      IF(MFA.NE.0) THEN
      L=0
      NLOW=1
      NHIGH=4
      DO 700 IY=1,NY
      DO 700 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 700 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
700   CONTINUE
      NLOW=1
      NHIGH=4
      DO 710 IY=1,NY
      DO 710 IX=NLOW,NHIGH
      IU=(IY-1)*NSURF+IX
      DO 710 IA=1,MOA
      L=L+1
      WRITE(6,4) L,IU,IA,3,0.0
710   CONTINUE
      END IF
C======================================
C======================================
C    PRINT information
C    MPRINT  no of print points 
C    MTIME  total run time
C    MTIME1 damping time
C    DTIME  time step
C======================================
C======================================
      MPRINT=40
      MTIME1=60000
      MTIME2=4000
      MTIME3=80000
      DTIME=40.0
      WRITE(6,5) MTIME1,MTIME2,MTIME3,DTIME
      WRITE(6,1) MPRINT
      DO 125 L=1,4
      LL=(L-1)*20000+2000
      DO 126 K=1,10
      LMN(K)=(K-1)*2000+LL
126   CONTINUE
      WRITE(6,1) (LMN(K),K=1,10)
125   CONTINUE
C==============================================
C==============================================
C==============================================
91    FORMAT('Title="Tecplot Results"')
92    FORMAT('Variables="X","Y","Z","DISP1","DISP2","DISP3","POLAR1",
     *"POLAR2","POLAR3","T11","T22","T33","T23","T31",
     *"T12"')
93    FORMAT('Zone T="Load step 0",N=',I5,' E=',I5,' datapacking=point,
     *zonetype=febrick')
94    FORMAT(15F15.6)
95    FORMAT(8I10)
C====================================================
C  write node information
C====================================================
      MPT=MP
      MET=ME
      POLAR=0.0D0
      IF(ID(1).EQ.1) MPT=MP+MPA
      IF(ID(1).EQ.1) MET=ME+MEA
      WRITE(8,91) 
      WRITE(8,92)
      WRITE(8,93) MPT,MET
      DO 109 I=1, MPT
      IF(I.LE.MP) THEN
       XPP(1,I)=X(I)
       XPP(2,I)=Y(I)
       XPP(3,I)=Z(I)
      ELSE
       XPP(1,I)=XA(I-MP)
       XPP(2,I)=YA(I-MP)
       XPP(3,I)=ZA(I-MP)
      END IF
109   CONTINUE
      DO 155 I=1, MPT
      WRITE(8,94) XPP(1,I),XPP(2,I),XPP(3,I),0.0,0.0,0.0,
     *POLAR(1,I),POLAR(2,I),POLAR(3,I),0.0,0.0,0.0,
     *0.0,0.0,0.0
155   CONTINUE
C====================================================
C  write element information (connectivity)
C====================================================
      DO 170 I=1,MET
      IF(I.LE.ME) THEN
        LM(:,I)=IJK(:,I)
      ELSE
        LM(:,I)=IJKA(:,I-ME)+MP
      END IF
170   CONTINUE
      DO 175 I=1,MET
      WRITE(8,95) (LM(J,I),J=5,8),
     *(LM(K,I),K=1,4)
175   CONTINUE
C====================================================
C====================================================
      WRITE(6,1) MP,ME,MCUNIT,MPA,MEA
C====================================================
      STOP
      END
