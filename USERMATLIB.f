            SUBROUTINE obrat(N,A,C)
            ! input ...
            ! a(n,n) - array of coefficients for matrix A
            ! n      - dimension
            ! output ...
            ! c(n,n) - inverse matrix of A
            ! comments ...
            ! the original matrix a(n,n) will be destroyed 
            ! during the calculation
            !===========================================================
            implicit none 
            integer n
            double precision A(n,n), C(n,n), OOO(n,2*n), lkl(n,2*n) 
            double precision ll,tt, gaus(n,2*n)
            integer i, j, pp, kk, kl
                  gaus = 0.0d0
                  do i = 1, N
                      do j = 1, N
                          gaus(i,j) = A(i,j)
                      end do
                  end do
                  do i=1, N
                      gaus(i,i+N) = 1.0d0
                  end do
                  OOO = gaus
                  do i = 1,N
                      lkl = OOO
                      kl = i+1
                      If(OOO(i,i).EQ.0.0d0) then
                          do while (OOO(i,i).EQ.0.0d0)
                              kl = kl
                          do j = 1, 2*N
                              OOO(i,j) = lkl(kl,j)
                              OOO(kl,j) = lkl(i,j)
                          end do
                          If(OOO(i,i).EQ.0.0d0) then
                              OOO = lkl
                          end if    
                          kl = kl +1
                          end do
                      end if    
                      tt = OOO(i,i)
                      do pp = 1, 2*N
                         OOO(i,pp) = OOO(i,pp)/tt 
                      end do   
                      do kk = 1,N-i
                          ll = OOO(kk+i,i)
                          do j = 1,2*N
                              OOO(kk+i,j) = OOO(kk+i,j)-OOO(i,j)*ll
                          end do
                      end do
                  end do   
c          обратный ход
                  do i = -N,-1
                      tt = OOO(-i,-i)
                      do pp = 1, 2*N
                         OOO(-i,pp) = OOO(-i,pp)/tt 
                      end do   
                      do kk = i+1, -1
                          ll = OOO(-kk,-i)
                          do j = 1,2*N
                              OOO(-kk,j) = OOO(-kk,j)-OOO(-i,j)*ll
                          end do
                      end do
                  end do   
                  do i=1,N
                      do j=1,N
                          C(i,j) = OOO(i,j+N)
                      end do
                  end do    
            end subroutine obrat
*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   cutFactor, pVolDer, hrmflg, var3, var4,
     &                   var5, var6, var7)
     
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT"::usermat  

#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 pVolDer (3),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
      DOUBLE PRECISION hrmflg
c
      EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d
c      EXTERNAL         usermat_harm

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
      data             var1/0.0d0/
      data             var2/0.0d0/

    
c ***    time domain analysis

      IF(ncomp .GE. 4) THEN
c ***    3d, plane strain and axisymmetric example
         call usermat3d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)

      ELSE IF(nDirect.eq. 2 .and. ncomp .EQ. 3) THEN
c ***    plane stress example
         call usermatps (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)

      ELSE IF(ncomp .EQ. 3) THEN
c ***    3d beam example
         call usermatbm (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)

      ELSE IF(ncomp .EQ. 1) THEN
c ***    1d beam example
         call usermat1d (
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords,
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)

      END IF
      return
      end
*deck,usermat1d    USERDISTRIB  parallel                                gal
      subroutine usermat1d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT1D"::usermat1d 
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, SMALL
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 mcomp      = 1
     &                 )
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c
c
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DOUBLE PRECISION qEl,   pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq,    signTens,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 twoG,  fratio, prpr
      return
      end
*deck,usermat3d    USERDISTRIB  parallel                                gal
      subroutine usermat3d(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad, 
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
          
     
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMAT3D"::usermat3d 
c*************************************************************************
c     *** primary function ***
c     Add to preprocessor /DNOSTDCALL /DARGTRAIL /DPCWIN64_SYS /DPCWINX64_SYS /DPCWINNT_SYS /DCADOE_ANSYS /D__EFL /DFORTRAN /fpp /4Yportlib /auto /c /Fo.\ /MD /W0 и YES
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and beam.
c
c           a 3D material constitutive model can use for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be use.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements or plane elements
c       in plane strain or axisymmetric stress state. 
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variable
c            ustatev(1)                     - equivalent plastic strain
c            ustatev(2) - statev(1+ncomp)   - plastic strain vector
c            ustatev(nStatev)               - von-Mises stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      tsstif   (dp,ar(2),io)             transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change 
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor 
c                                         define it if a smaller step size is wished
c                                         recommended value is 0~1
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresss and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
cc
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut, 
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), kmnp(nStatev),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3), D   (ncomp  ), klk(nStatev),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2), llll(ncomp),
     &                 dsdePl1  (2*ncomp+2,2*ncomp+2)         
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6
     &                 )
     
c *** primary function ***
c
c user defined material constitutive model
c
c Attention:
c     User must define material constitutive law properly
c     according to the stress state such as 3D, plane strain 
c     and axisymmetry, plane stress and 3D/1D beam.
c
c     A 3D material constitutive model can be used for
c     plane strain and axisymmetry cases.
c
c     When using shell elements, a plane stress algorithm
c     must be used. At this time is not realized
c
c     Based on the algorithm implemented in the FORTRAN program
c     developed by De Souza Neto et al. (2008)[5]
c
c --- local variables
c
c      sigElp   (dp,ar(6  ),l)            trial stress
c      dsdeEl   (dp,ar(6,6),l)            elastic moduli
c      sigDev   (dp,ar(6  ),l)            deviatoric stress tensor
c      dfds     (dp,ar(6  ),l)            derivative of the yield function 
c      JM       (dp,ar(6,6),l)            2D matrix for a 4 order tensor
c      pEl      (dp,sc     ,l)            hydrostatic pressure stress
c      qEl      (dp,sc     ,l)            von-mises stress
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      sigy     (dp,sc     ,l)            current equivalent stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,l)            three time of shear moduli
c      rr       (dp,sc     ,l)            hardening constant
c      s        (dp,sc     ,l)            hardening constant
c      gamma    (dp,sc     ,l)            hardening constant
c      Rinf     (dp,sc     ,l)            hardening constant
c      K        (dp,sc     ,l)            bulk modulus
c      GG       (dp,sc     ,l)            shear modulus
c      w        (dp,sc     ,l)            Damage
c      R        (dp,sc     ,l)            equivalent plastic strain
c      sigeqv   (dp,sc     ,l)            equivalent elastic stress
c      epseqv   (dp,sc     ,l)            equivalent strian
c      Int      (dp,sc     ,l)            integrity
c      f1, f2, x, x1, x2, sch (dp,sc     ,l)            variables used to determine actual stress
c      eps      (dp,sc     ,l)            iteration precision
c      dsdePl   (dp,ar(ncomp,ncomp),io)                            material jacobian matrix
c      
c      
c --- temperary variables for solution purpose
c      i, j
c      threeOv2qEl, oneOv3G, qElOv3G, con1, con2, fratio
      EXTERNAL         vzero, vmove, get_ElmData, egen
      double precision egen
      DOUBLE PRECISION sigElp(mcomp), Strainob(mcomp), epsDev(mcomp),
     &                 dsdeEl(mcomp,mcomp), G(mcomp), sig_dev(mcomp),
     &                 sigDev(mcomp), JM(ncomp,ncomp),IxI(ncomp,ncomp),
     &                 sigi(mcomp), strainEl(mcomp), ak(ncomp),   
     &                 diagi(ncomp,ncomp), Ii(ncomp), eye(ncomp,ncomp), 
     &                 Id(ncomp,ncomp), sig5(mcomp), Is(ncomp,ncomp),
     &                 Di(ncomp), Ii2(ncomp), enToPhys(mcomp), 
     &                 physToEn(ncomp), depsPl(ncomp), temp1(mcomp), 
     &                 sig_Dev1(mcomp), strainEl1(mcomp), temp2(mcomp),
     &                 sigElp_tr(mcomp), PhysToEng(mcomp),StrDev(mcomp),
     &                 sigg(mcomp), sigElp1(ncomp), siggT(ncomp,1),
     &                 beta(ncomp), beta0(ncomp), rel_tr(mcomp), 
     &                 funA(ncomp), alfa_k(2*ncomp+2), tipp(ncomp+1),
     &                 invCe(ncomp,ncomp), rel(mcomp), rel2(mcomp),
     &                 sigg_Dev(mcomp), NN(ncomp), DP(ncomp,ncomp),
     &                 N2(ncomp), DCe(ncomp,ncomp),dq_dSig(ncomp),
     &                 dSig_dD(ncomp), dRel_dD(ncomp), B(1,ncomp),
     &                 A(2*ncomp+2,1), dN_dSig(ncomp,ncomp),
     &                 dA1_dSig(ncomp), dN_dBeta(ncomp,ncomp), 
     &                 dA1_dBeta(ncomp),  dA2_dD(ncomp), siggg(1,ncomp), 
     &                 dA2_dSig(ncomp,ncomp), dA3_dD(ncomp), 
     &                 dA3_dSig(ncomp,ncomp), dA3_dP(ncomp),
     &                 dA3_dBetta(ncomp), dY_dSig(ncomp),deps_pl(ncomp),
     &                 dA4_dSig(ncomp),  jac(ncomp), deps_pl_Dev(ncomp),
     &                 tr_dSig_dD(ncomp), tr_dRel_dD(ncomp),
     &                 alpha_k(2*ncomp+2,1), alpha(2*ncomp+2,1),
     &                 dA2_dP(ncomp), dA3_dBeta(ncomp,ncomp),
     &                 dN_dD(ncomp), dA2_dBeta(ncomp,ncomp),
     &                 dsdeEl1(ncomp,ncomp), ijac(2*ncomp+2,2*ncomp+2),
     &                 eyeBig(2*ncomp+2,2*ncomp+2), tippp(ncomp), 
     &                 ffff(ncomp,1), gggg(ncomp,1), ggggg(ncomp),
     &                 dsdePl2(ncomp,ncomp), J222(ncomp)
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
      DATA PhysToEng/1.0D0,1.0D0,1.0D0,0.5D0,0.5D0,0.5D0/
c
      INTEGER          i, j, sch, urt, iter, tutu, incRev
      DOUBLE PRECISION pEl, qEl, pleq_t, sigy_t, sigy, eps_e, a1, a2,
     &                 dpleq, pleq, epEl, w0, R, R0, sig_hyd, plCor, a3,
     &                 young, posn, sigy0, dsigdep, q, dY, F, snorm,
     &                 elast1, elast2, s, Phi, Int0, Intt, plasticMult,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,  con2, dperr(3), eps, lamme, emd,
     &                 f3, K, fsig, f1, f2, young1, posn1, sigeqv, a4,
     &                 epseqv, dd, e, Dam, fsig0, dfsig, rr, fsigy,
     &                 Rinf, gamma, eel, eps_e_tr, norm_f, sigeqv_tr,
     &                 dintt, pr, pr1, GG, Dam0, sig, ff, eps_e_hyd_tr,
     &                 eEl1, qEl1, toch, qtr, qtr2, j2, eEl5, qtr5, w,
     &                 epspoln, epseqv5, qEl2, Str, eEL2, aaa, bbb,
     &                 beta_eqv, q_tr, Y, f4, f6, dJ2_dD, Determine,
     &                 sigg_hyd, J2_tr, sum_st_dSig_dD, norm_A, dq_dD,
     &                 maxiter, maxiter1, maxiter2, f5, f7, dA4_dP,
     &                 sum_diag_dSig_dD, dA1_dD, dA1_dP, dA4_dD, Det3,
     &                 deps_pl_hyd, deps_pl_eqv, M1, M2, M3, Det1, Det2,
     &                 J22, qEl11,v,x1,x2,x,alfa, rrrr, Rconst, 
     &                 sigg_hyd_beta, sigg_Dev_beta(ncomp), kkk,
     &                 princstrs(11), princSm, princSdev(3), cn1, cn2,
     &                 Lodeangle, LimSurfVal, NLCnst, NLPwr      
c*************************************************************************
c
      keycut   = 0
c *** get Young's modulus in MPa and Poisson's ratio, initial yield stress and the hardening parameters
c *** material parameters
      young    = prop(1)
      posn     = prop(2)
c *** initial yield stress     
      sigy0    = prop(3)
c *** isotropic hardening     
      rr       = prop(4)   
c *** isotropic hardening     
      s        = prop(5)   
c *** isotropic hardening     
      Rinf     = prop(6)   
c *** isotropic hardening     
      gamma    = prop(7)
c *** kinematic hardening
      aaa      = prop(8)
c *** kinematic hardening
      bbb      = prop(9) 
      Rconst   = prop(10)
c *** Multiaxial flag, if zero Rv=1
      incRev    = 1
c *** Multiaxial function const 
      NLCnst    = prop(11)
c *** Multiaxial function power
      NLPwr    = prop(12)
c
c     Limit surface, default eq to 1  
      LimSurfVal=1

      

      kkk = CEILING(Time)
      if (MOD(CEILING(TIME),2).EQ.0) then
c          aaa = -aaa
c          bbb = -bbb
c          sigy0 = 0.9*sigy0
      endif    
c *** bulk modulus    
      K  = third*young/(1-2*posn)
c *** shear modulus
      GG = young/(TWO*(ONE+posn))
      R0 = zero
      w0 = zero
      pleq_t = zero
      plasticMult = ZERO 
c *** initial damage
      w0       = ustatev(8)
c *** initial hardering internal variable
      R0       = ustatev(1)
c *** initial integrity
      Int0 = 1 - w0
c *** initial true back stress
        do i=1,ncomp
           beta0(i) = ustatev(i+10)
        end do
c *** array for hydraulic component retrieval
        eyeBig = zero
        Ii = ONE
        physToEn = ONE
        do i=1,(2*ncomp+2)
           eyeBig(i, i) = ONE
        end do
        do i = 4,ncomp
           Ii(i) = zero
           physToEn(i) = TWO
        end do         
c *** conversion arrays      
        diagI = zero
        eye = zero
        do i=1,ncomp
           diagI(i, i) = Ii(i)
           eye(i, i) = 1
        end do    
c          deviatoric projection tensor in Voigt notation
         do i=1,ncomp
            do j=1,ncomp
              Is(i,j) = (eye(i,j) + diagI(i,j))*0.5
              IxI(i,j) = Ii(i)*Ii(j)
              Id(i,j) = Is(i,j) - third*IxI(i,j) 
            end do
         end do
c         
         do i=4,ncomp
              Is(i,i) = ONE
              Id(i,i) = ONE
         end do
      kmnp = zero 
      klk = zero

c
c *** plastic strain tensor
      call vmove(ustatev(2), epsPl(1), ncomp)
      twoG     = young / (ONE+posn)
      threeG   = ONEHALF * twoG
      elast1=young*posn/((1.0D0+posn)*(1.0D0-TWO*posn))
      elast2=HALF*twoG
c
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = elast2
      
c
c *** calculate elastic stiffness matrix (3d)
      dsdeEl(1,1)=(elast1+TWO*elast2)*G(1)*G(1)
      dsdeEl(1,2)=elast1*G(1)*G(2)+elast2*TWO*G(4)*G(4)
      dsdeEl(1,3)=elast1*G(1)*G(3)+elast2*TWO*G(5)*G(5)
      dsdeEl(1,4)=elast1*G(1)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(1,5)=elast1*G(1)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(1,6)=elast1*G(1)*G(6)+elast2*TWO*G(4)*G(5)
      dsdeEl(2,2)=(elast1+TWO*elast2)*G(2)*G(2)
      dsdeEl(2,3)=elast1*G(2)*G(3)+elast2*TWO*G(6)*G(6)
      dsdeEl(2,4)=elast1*G(2)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(2,5)=elast1*G(2)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(2,6)=elast1*G(2)*G(6)+elast2*TWO*G(2)*G(6)
      dsdeEl(3,3)=(elast1+TWO*elast2)*G(3)*G(3)
      dsdeEl(3,4)=elast1*G(3)*G(4)+elast2*TWO*G(5)*G(6)
      dsdeEl(3,5)=elast1*G(3)*G(5)+elast2*TWO*G(5)*G(3)
      dsdeEl(3,6)=elast1*G(3)*G(6)+elast2*TWO*G(6)*G(3)
      dsdeEl(4,4)=elast1*G(4)*G(4)+elast2*(G(1)*G(2)+G(4)*G(4))
      dsdeEl(4,5)=elast1*G(4)*G(5)+elast2*(G(1)*G(6)+G(5)*G(4))
      dsdeEl(4,6)=elast1*G(4)*G(6)+elast2*(G(4)*G(6)+G(5)*G(2))
      dsdeEl(5,5)=elast1*G(5)*G(5)+elast2*(G(1)*G(3)+G(5)*G(5))
      dsdeEl(5,6)=elast1*G(5)*G(6)+elast2*(G(4)*G(3)+G(5)*G(6))
      dsdeEl(6,6)=elast1*G(6)*G(6)+elast2*(G(2)*G(3)+G(6)*G(6))
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
c
c *** get initial stress
      call vzero(sigi(1),ncomp)      
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)
  
c
c *** calculate the trial stress and
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         strainEl(i) = Strain(i) + dStrain(i) - epsPl(i)
      end do
      call vzero(sigElp, 6)
      call vzero(sigElp1, 6)
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i)  = sigElp(i)+dsdeEl(j,i)*strainEl(j)*Int0
         end do
         sigElp(i)  = sigElp(i) + sigi(i)
      end do
      dsdeEl1 = dsdePl   
      IF (incRev==1) THEN
c     calculate principal stresses 
              princstrs=0
              princstrs(1:ncomp) = sigElp(1:ncomp)
              call prinst(princstrs) 
c     Lode angle
              princSm=(princstrs(7)+princstrs(8)+princstrs(9))/3
              princSdev(1:3)=princstrs(7:9)-princSm
              cn1=princSdev(1)*princSdev(2)*princSdev(3)
              Lodeangle=27/2*cn1/princstrs(11)**3
              cn2=1-2/3.14*acos(Lodeangle)
c     multuaxial state
           LimSurfVal=NLCnst+(1-NLCnst)*(cn2*6/3.14)**NLPwr     
       ENDIF
      sigy0 = sigy0/LimSurfVal
c
c *** elastic trial true stress
      sigElp_tr  = sigElp
c      
      If  (ncomp .EQ. 4) then
          strainEl(5) = 0
          strainEl(6) = 0 
      end if
c *** true hydrostatic pressure trial stress
      pEl = THIRD * (sigElp_tr(1) + sigElp_tr(2) + sigElp_tr(3))
c *** compute the true deviatoric trial stress tensor
      sigDev = sigElp_tr - pEl*Ii
c *** compute relative trial stress
      rel_tr = sigDev - beta0
      If  (ncomp .EQ. 4) then
          rel_tr(5) = 0
          rel_tr(6) = 0 
      end if
c *** compute trial von-mises and equivalent elastic true stress    
      J2_tr = sum(rel_tr**2*PhysToEn) / 2  
      q_tr = sqrt(3.0d0*J2_tr)   
c *** compute yield stress
      fsig = sigy0 + Rinf * (1 - exp(-gamma * R0)) + Rconst*R0
      Phi = q_tr/Int0 - fsig
      iter = zero
c *** Check if Yield Criterion is met
      IF (Phi.GE.0) THEN   
c ******* inital guess for the system variables
          sigg = sigElp_tr
          w = 1 - Int0
          If  (ncomp .EQ. 4) then
              sigg(5) = 0
              sigg(6) = 0                       
          end if         
          beta = beta0
c ******* inverse tangent
          maxiter2 = 500
          plasticMult = zero
          norm_A = 1
          alpha_k = zero
          alpha = alpha_k
          CALL obrat(ncomp, dsdePl, dsdePl2)
          invCe = matmul(dsdePl2, eye)
c          
          do while (norm_A.GE.1e-10.and.iter.LE.(maxiter2))
c
                IF (iter.EQ.maxiter2) GO TO 800
               
c ****** inital guess for the hardening variable
                R = R0 + abs(plasticMult)
                Intt = 1 - w
c                
                sigg_hyd = THIRD * (sigg(1) + sigg(2) + sigg(3))
                sigg_Dev = sigg - sigg_hyd*Ii
                sigg_hyd_beta = THIRD * (beta(1) + beta(2) + beta(3))
                sigg_Dev_beta = beta - sigg_hyd_beta*Ii
                rel = sigg_Dev - sigg_Dev_beta
                rel2 = rel * physToEn
                  If  (ncomp .EQ. 4) then
                      rel(5) = 0
                      rel(6) = 0
                      rel2(5) = 0
                      rel2(6) = 0                       
                  end if
                do i=1, ncomp
                    J222(i) = rel(i)*rel(i)*physToEn(i)
                end do
                J22 = sum(rel**2*PhysToEn) / TWO 
                qEl = sqrt(3.0d0*J22)
c
                fsig = sigy0 + Rinf * (1 - exp(-gamma*R)) + Rconst*R 
c ****** derivatives
                dfsig = Rinf * gamma * exp ( -gamma * R ) + Rconst
                NN = ONEHALF * rel / (Intt * qEl)
                N2 = ONEHALF * rel2 / (Intt * qEl)
c 
                do i=1, ncomp
                  siggg(1,i) = sigg(i)
                end do
                do i=1, ncomp
                  siggT(i,1) = sigg(i)
                end do
                B = matmul(siggg,invCe)
                IF (ncomp .EQ. 4) then
                  Y = 
     &            -(1/(2*Intt**2))*(siggT(1,1)*B(1,1)+
     &            siggT(2,1)*B(1,2)+siggT(3,1)*B(1,3)+
     &            siggT(4,1)*B(1,4)) 
                end if
                IF (ncomp .EQ. 6) then
                  Y = 
     &            -(1/(2*Intt**2))*(siggT(1,1)*B(1,1)+
     &            siggT(2,1)*B(1,2)+siggT(3,1)*B(1,3)+
     &            siggT(4,1)*B(1,4)+siggT(5,1)*B(1,5)+
     &            siggT(6,1)*B(1,6))
                end if                       
                DCe = Intt *  dsdeEl1                                  
c ****** some recurring factors
                f1 = qEl ** 2                                         
                f3 = aaa * plasticMult                                  
                f4 = (-Y/rr) ** s                                   
                f5 = Intt ** 2                                      
                f6 = (-Y/rr) ** (s-1)                              
                f7 = plasticMult * s / (rr * Intt)                     
c ****** some recurring derivatives with respect to true stress
                dq_dsig = NN * Intt                                    
c ****** some recurring derivatives with respect to D                                                                                         
                dSig_dD = - sigg / Intt     
                tr_dSig_dD = dSig_dD * Ii                              
                sum_st_dSig_dD = sum(tr_dSig_dD)                        
                 dRel_dD =                                            
     &              dSig_dD - THIRD * sum_st_dSig_dD * Ii
                tr_dRel_dD = dRel_dD * rel * physToEn                 
                dJ2_dD = sum(tr_dRel_dD)                              
                dq_dD = HALF * dJ2_dD * ((3/J22)**(HALF))              
c ****** compute all functions                                                                     
                A(1,1) = qEl / Intt - fsig
                  do i=1,ncomp
                      ffff(i,1) = StrainEl(i)-plasticMult*N2(i) 
                  end do
                  gggg = matmul(DCe,ffff)
                  do i=1, ncomp
                  ggggg(i) = gggg(i,1)
                  end do
                  DO i = 1 , ncomp
                       A(i+1,1)=sigg(i)-ggggg(i)
                  END DO 
                  DO i = 1 , ncomp
                   A((i+ncomp+1),1) = 
     &            beta(i)-beta0(i)-plasticMult*
     &            (aaa*2/3*NN(i)-bbb*beta(i))                 
                   
                  END DO
                  A((2*ncomp+2),1) = w - w0 - (1/Intt)*f4*plasticMult   

c ****** set the norm
                Norm_A = norm2(A)                                      
c ****** compute all necessary derivatives
                  dA1_dSig = N2                                        
                  dA1_dD = qEl/f5 + dq_dD/Intt                        
                  dA1_dP = -dfsig                                     
                  dA1_dBeta = -N2                                    
                do i=1,ncomp
                 do j=1,ncomp  
                   dN_dSig(i,j) = rel(i)*dq_dSig(j)                   
                 end do
                end do  
                dN_dSig = 
     &             -1.5d0*(dN_dSig - qEl*Id)/(Intt*f1)                 
                dN_dD=                                                
     &           1.5d0*(qEl*(rel+Intt*dRel_dD)-
     &           Intt*dq_dD*rel)/(f5*f1)
                do i=1,ncomp
                 do j=1,ncomp
                dN_dBeta(i,j) = rel(i)*dq_dSig(j)
                 end do
                end do 
                dN_dBeta = 
     &                1.5d0*(dN_dBeta- qEl*Is)/(Intt*f1)
                dA2_dSig=Is+plasticMult*Intt*2*GG*dN_dSig            
                  do i=1,ncomp
                      ffff(i,1) = StrainEl(i)
                  end do
                  gggg = matmul(dsdeEl1,ffff)
                  do i=1, ncomp
                  ggggg(i) = gggg(i,1)
                  end do                
                dA2_dD =                                               
     &            dSig_dD+ggggg+2*GG*(Intt*plasticMult*
     &            dN_dD-plasticMult*NN)
                  dA2_dP = Intt*2*GG*NN                               
                  dA2_dBeta = plasticMult*Intt*2*GG*dN_dBeta          
                  dA3_dSig = -f3 * dN_dSig                             
                  dA3_dD = -f3 * dN_dD
                  dA3_dP = bbb * beta - aaa*2/3 * NN                 
                dA3_dBeta=
     &            (ONE+bbb*plasticMult)*Is-f3*dN_dBeta                
                  do i=1,ncomp
                      ffff(i,1) = sigg(i)
                  end do
                  gggg = matmul(invCe,ffff)
                  do i=1, ncomp
                  ggggg(i) = gggg(i,1)
                  end do                     
                  dY_dSig = - (1/(Intt**2))*ggggg                     
                  dA4_dSig = f7 * dY_dSig * f6                      
                dA4_dD = 1 - plasticMult * f4 / f5
                dA4_dP = -(1 / Intt)*f4
                dsdePl1 = zero 
c ****** fill Jacobian
                do i=1,ncomp    
                    dsdePl1(1,i) = dA1_dSig(i)
                    dsdePl1(1, ncomp+1) = dA1_dD
                    dsdePl1(1, ncomp+2) = dA1_dP
                    dsdePl1(1, i+ncomp+2) = dA1_dBeta(i)
                end do    
                do i=1,ncomp
                    do j=1,ncomp
                    dsdePl1(i+1, j) = dA2_dSig(i,j)
                    dsdePl1(i+1, ncomp+1) = dA2_dD(i)
                    dsdePl1(i+1, ncomp+2) = dA2_dP(i)
                    dsdePl1(i+1, j+ncomp+2) = dA2_dBeta(i,j)

                    dsdePl1(i+ncomp+1, j) = dA3_dSig(i,j)
                    dsdePl1(i+ncomp+1, ncomp+1) = dA3_dD(i)
                    dsdePl1(i+ncomp+1, ncomp+2) = dA3_dP(i)
                    dsdePl1(i+ncomp+1, j+ncomp+2) = dA3_dBeta(i,j)
                
                    dsdePl1(2*ncomp+2, j) = dA4_dSig(j)
                    dsdePl1(2*ncomp+2, ncomp+1) = dA4_dD
                    dsdePl1(2*ncomp+2, ncomp+2) = dA4_dP
                    dsdePl1(2*ncomp+2, j+ncomp+2) = 0
                    end do
                end do
c ****** make a vector of the system variables                             
                do i=1,ncomp
                    alpha_k(i,1) = sigg(i)
                    alpha_k(ncomp+1,1) = w
                    alpha_k(ncomp+2,1) = plasticMult
                    alpha_k(i+ncomp+2,1) = beta(i)
                end do  
c ****** solve the linear system
                    alpha = 0
                  CALL obrat(2*ncomp+2, dsdePl1, dsdePl1)
                    alpha = alpha_k - matmul(dsdePl1, A )
                    iter = iter + 1
c ****** update internals               
                do i=1,ncomp
                    sigg(i) = alpha(i,1)
                    w = alpha(ncomp+1,1)
                    plasticMult = alpha(ncomp+2,1)
                    beta(i) = alpha(i+ncomp+2,1)
                end do 
          end do          
c *** secondary system variables
                R = R0 + abs(plasticMult)
                Intt = 1 - w
c                
                sigg_hyd = THIRD * (sigg(1) + sigg(2) + sigg(3))
                sigg_Dev = sigg - sigg_hyd*Ii
                sigg_hyd_beta = THIRD * (beta(1) + beta(2) + beta(3))
                sigg_Dev_beta = beta - sigg_hyd_beta*Ii
                rel = sigg_Dev - sigg_Dev_beta
                rel2 = rel * physToEn
                  If  (ncomp .EQ. 4) then
                      rel(5) = 0
                      rel(6) = 0
                      rel2(5) = 0
                      rel2(6) = 0                       
                  end if
                do i=1, ncomp
                    J222(i) = rel(i)*rel(i)*physToEn(i)
                end do
                J22 = sum(rel**2*PhysToEn) / TWO 
                qEl = sqrt(3.0d0*J22)
c
                fsig = sigy0 + Rinf * (1 - exp(-gamma*R)) + Rconst * R
c ****** derivatives
                dfsig = Rinf * gamma * exp ( -gamma * R ) + Rconst
                NN = ONEHALF * rel / (Intt * qEl)
                N2 = ONEHALF * rel2 / (Intt * qEl)
c 
                do i=1, ncomp
                  siggg(1,i) = sigg(i)
                end do
                do i=1, ncomp
                  siggT(i,1) = sigg(i)
                end do
                B = matmul(siggg,invCe)
                IF (ncomp .EQ. 4) then
                  Y = 
     &            -(1/(2*Intt**2))*(siggT(1,1)*B(1,1)+
     &            siggT(2,1)*B(1,2)+siggT(3,1)*B(1,3)+
     &            siggT(4,1)*B(1,4)) 
                end if
                IF (ncomp .EQ. 6) then
                  Y = 
     &            -(1/(2*Intt**2))*(siggT(1,1)*B(1,1)+
     &            siggT(2,1)*B(1,2)+siggT(3,1)*B(1,3)+
     &            siggT(4,1)*B(1,4)+siggT(5,1)*B(1,5)+
     &            siggT(6,1)*B(1,6))
                end if                       
                DCe = Intt *  dsdeEl1                                 
c ****** some recurring factors
                f1 = qEl ** 2                                         
                f3 = aaa * plasticMult                               
                f4 = (-Y/rr) ** s                                      
                f5 = Intt ** 2                                         
                f6 = (-Y/rr) ** (s-1)                                 
                f7 = plasticMult * s / (rr * Intt)                     
c ****** some recurring derivatives with respect to true stress
                dq_dsig = NN * Intt                                   
c ****** some recurring derivatives with respect to D                                                                                           
                dSig_dD = - sigg / Intt     
                tr_dSig_dD = dSig_dD * Ii                             
                sum_st_dSig_dD = sum(tr_dSig_dD)                        
                 dRel_dD =                                             
     &              dSig_dD - THIRD * sum_st_dSig_dD * Ii
                tr_dRel_dD = dRel_dD * rel * physToEn                  
                dJ2_dD = sum(tr_dRel_dD)                               
                dq_dD = HALF * dJ2_dD * ((3/J22)**(HALF))             
c ****** compute all necessary derivatives
                  dA1_dSig = N2                                       
                  dA1_dD = qEl/f5 + dq_dD/Intt                        
                  dA1_dP = -dfsig                                     
                  dA1_dBeta = -N2                                      
                do i=1,ncomp
                 do j=1,ncomp  
                   dN_dSig(i,j) = rel(i)*dq_dSig(j)                  
                 end do
                end do  
                dN_dSig = 
     &             -1.5d0*(dN_dSig - qEl*Id)/(Intt*f1)                
                dN_dD=                                            
     &           1.5d0*(qEl*(rel+Intt*dRel_dD)-
     &           Intt*dq_dD*rel)/(f5*f1)
                do i=1,ncomp
                 do j=1,ncomp
                dN_dBeta(i,j) = rel(i)*dq_dSig(j)
                 end do
                end do 
                dN_dBeta = 
     &                1.5d0*(dN_dBeta- qEl*Is)/(Intt*f1)
                dA2_dSig=Is+plasticMult*Intt*2*GG*dN_dSig             
                  do i=1,ncomp
                      ffff(i,1) = StrainEl(i)
                  end do
                  gggg = matmul(dsdeEl1,ffff)
                  do i=1, ncomp
                  ggggg(i) = gggg(i,1)
                  end do                
                dA2_dD =                                              
     &            dSig_dD+ggggg+2*GG*(Intt*plasticMult*
     &            dN_dD-plasticMult*NN)
                  dA2_dP = Intt*2*GG*NN                             
                  dA2_dBeta = plasticMult*Intt*2*GG*dN_dBeta           
                  dA3_dSig = -f3 * dN_dSig                            
                  dA3_dD = -f3 * dN_dD
                  dA3_dP = bbb * beta - aaa*2/3 * NN                   
                dA3_dBeta=
     &            (ONE+bbb*plasticMult)*Is-f3*dN_dBeta                
                  do i=1,ncomp
                      ffff(i,1) = sigg(i)
                  end do
                  gggg = matmul(invCe,ffff)
                  do i=1, ncomp
                  ggggg(i) = gggg(i,1)
                  end do                     
                  dY_dSig = - (1/(Intt**2))*ggggg                     
                  dA4_dSig = f7 * dY_dSig * f6                        
                dA4_dD = 1 - plasticMult * f4 / f5
                dA4_dP = -(1 / Intt)*f4
                dsdePl1 = zero 
c ****** fill Jacobian
                do i=1,ncomp    
                    dsdePl1(1,i) = dA1_dSig(i)
                    dsdePl1(1, ncomp+1) = dA1_dD
                    dsdePl1(1, ncomp+2) = dA1_dP
                    dsdePl1(1, i+ncomp+2) = dA1_dBeta(i)
                end do    
                do i=1,ncomp
                    do j=1,ncomp
                    dsdePl1(i+1, j) = dA2_dSig(i,j)
                    dsdePl1(i+1, ncomp+1) = dA2_dD(i)
                    dsdePl1(i+1, ncomp+2) = dA2_dP(i)
                    dsdePl1(i+1, j+ncomp+2) = dA2_dBeta(i,j)

                    dsdePl1(i+ncomp+1, j) = dA3_dSig(i,j)
                    dsdePl1(i+ncomp+1, ncomp+1) = dA3_dD(i)
                    dsdePl1(i+ncomp+1, ncomp+2) = dA3_dP(i)
                    dsdePl1(i+ncomp+1, j+ncomp+2) = dA3_dBeta(i,j)
                
                    dsdePl1(2*ncomp+2, j) = dA4_dSig(j)
                    dsdePl1(2*ncomp+2, ncomp+1) = dA4_dD
                    dsdePl1(2*ncomp+2, ncomp+2) = dA4_dP
                    dsdePl1(2*ncomp+2, j+ncomp+2) = 0
                    end do
                end do
c ****** solve the linear system
          CALL obrat(2*ncomp+2, dsdePl1, dsdePl1)
c          invCe = matmul(dsdePl1, eyeBig)             
c ****** update internals               
         ijac = matmul(dsdePl1, eyeBig)                               
         do i = 1, ncomp
          do j = 1, ncomp
              dsdePl(i,j) =ijac(i,j+1) 
          end do
         end do 
          dsdePl = matmul(dsdePl,dsdeEl1)*Intt
          deps_pl = plasticMult * N2
          if (rrrr.EQ.1000) then
              do i=1,ncomp
                  ustatev(i+20) = Strain(i)
                  ustatev(i+24) = dstrain(i)
                  ustatev(i+28) = epsPl(i)
              end do
          end if    
          epsPl = epsPl + deps_pl   
          Str = THIRD * (epsPl(1) + epsPl(2) + epsPl(3)) 
          StrDev(1) = epsPl(1) - Str
          StrDev(2) = epsPl(2) - Str
          StrDev(3) = epsPl(3) - Str
          StrDev(4) = epsPl(4) * 0.5 
          If (ncomp .EQ. 4) then
                StrDev(5) = 0
                StrDev(6) = 0
          end if
          If (ncomp .EQ. 6) then
                StrDev(5) = epsPl(5) * 0.5
                StrDev(6) = epsPl(6) * 0.5
          end if      
          pleq = 
     &  StrDev(1) * StrDev(1) + StrDev(2) * StrDev(2) +
     &  StrDev(3) * StrDev(3) +
     &  TWO*(StrDev(4) * StrDev(4) + StrDev(5) * StrDev(5) + 
     &  StrDev(6) * StrDev(6))
          pleq = sqrt(pleq/ONEHALF)
          epseq  = pleq         
          DO i = 1 , ncomp
               stress(i) = sigg(i)
          END DO 
          sig_hyd = THIRD * (stress(1) + stress(2) + stress(3))        
c *** initial damage
          ustatev(8) = w
          J22 = sum(beta**2*PhysToEn) / TWO 
          qEl = sqrt(3.0d0*J22)
          ustatev(7) = beta(2)
c *** isotropic hardening internal variable
          ustatev(9) = R0
          ustatev(10) = q_tr/Int0
c *** back stress
          do i=1,ncomp
            ustatev(i+10) = beta(i)
          end do 
c *** true vM stress
c *** effective yield stress
          ustatev(1)  = R0 + abs(plasticMult)
          do i=1,ncomp
            ustatev(i+1) = epsPl(i)
          end do           
      else
          IF (young .EQ. 1000) then
              do i=1, ncomp
                  stress(i) = sigElp(i)
                  strain(i)=strainEl(i)
              end do    
          else          
              do i=1,ncomp
                 do j=1,ncomp        
                     dsdePl(i,j) = Int0 * dsdeEl(i,j)
                 end do
              end do
              do i=1, ncomp
                  stress(i) = sigElp_tr(i)
              end do
          end if
          !IF (young .EQ. 1) then
          !    !pleq = zero
          !    !epsPl = zero
          !else    
          !    pleq = R0
          !end if    
           pleq = R0
c *** true stress           
c *** initial damage
          ustatev(8) = w0
          J22 = sum(beta0**2*PhysToEn) / TWO 
          qEl = sqrt(3.0d0*J22)
          ustatev(7) = beta(2)
c *** isotropic hardening internal variable
          ustatev(9) = qEl
          ustatev(10) = fsig
c *** back stress
          do i=1,ncomp
            ustatev(i+10) = beta0(i)
          end do 
c *** true vM stress
c *** effective yield stress
          ustatev(1) = pleq
          do i=1,ncomp
            ustatev(i+1) = epsPl(i)
          end do           
      end if
c      ustatev(15) = zero
 
  800 continue          
      ustatev(7) = LimSurfVal
      return    
      end
*deck,usermatbm    USERDISTRIB  parallel                                gal
      subroutine usermatbm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMATBM"::usermatbm 

#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,  cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, ONE, TWO, SMALL, SQTWOTHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 3
     &                 )

c
      EXTERNAL         vmove, vzero, vapb1, vamb1,get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 cpleq, dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep,
     &                 c1, c2, c3, c4, c5
      DOUBLE PRECISION pv(3)
      data pv/TWOTHIRD, TWO, TWO/

      return
      end
*deck,usermatps    USERDISTRIB  parallel                                gal
      subroutine usermatps(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5, 
     &                   var6, var7)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERMATPS"::usermatps 

#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp), sigi(ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          NEWTON, mcomp
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, 
     &                 SQTWOTHIRD, SQTWO1,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 SQTWO1     = 0.707106769084930420d0,
     &                 NEWTON     = 20,
     &                 mcomp      = 6
     &                 )

      EXTERNAL         vmove, vzero, vapb1, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      INTEGER          i, j, k
      DOUBLE PRECISION pleq_t,  sigy_t , sigy,
     &                 dpleq,   pleq,    twoG,    et,
     &                 young, posn,    sigy0,   dsigdep, tEo1pm,
     &                 gamma, dgamma,  dfdga,   dplga, 
     &                 funcFb,funcFb2, funcf,   dFdep,   fratio,
     &                 con1,  con2,    con3,  con4,
     &                 con2p1, ocon2p1,
     &                 ocon2p2, con4p1, ocon4p1, ocon4p2,
     &                 c1, c2, c3,c4, c5,dperr(3)

      return
      end

*deck,usermat_harm    USERDISTRIB  parallel                    jmgerken
      subroutine usermat_harm(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nProp,
     &                   freq,dfreq,Temp,stress,jacobi,tsstif,
     &                   strain,prop,coords)

#include "impcom.inc"

      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nProp
      DOUBLE PRECISION freq,dfreq, Temp
      DOUBLE PRECISION 
     &                 stress  (ncomp,2),
     &                 jacobi  (ncomp,ncomp,2),
     &                 strain  (ncomp,2),
     &                 prop    (nProp),
     &                 coords  (3),
     &                 tsstif  (2)
c
c***************** Local *************************************
c

      INTEGER i, j, iprony ! loop counter
      INTEGER nprnsh,nprnvl ! number of prony terms
      INTEGER nvar  ! index
      DOUBLE PRECISION c1, c2, c3 ! temporary
      DOUBLE PRECISION alpha, tau  ! prony parameters
      DOUBLE PRECISION Omega  ! frequency
      DOUBLE PRECISION smG,lmG ! storage and loss shear modulus factors
      DOUBLE PRECISION smK,lmK ! storage and loss bulk modulus factors
      DOUBLE PRECISION young, posn, G, bulk  ! elastic constants

      DOUBLE PRECISION pi
      PARAMETER (pi=3.14159265358979324D0)


        

      return
      end
