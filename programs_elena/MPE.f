      Program MPE!  E. Y. Vedmedenko

      implicit real*8(a-h,o-x)
      implicit complex*16(y,z)
      character*20 f1,f2,f3,f4,f6      !!file with j,m;resut file;file with coordinates
      character*20 s3,s4
      EXTERNAL On_Intrp
      logical Fl_Abort, RFlag
      data Fl_Abort/.FALSE./, RFlag/.FALSE./
      dimension Lj(10),jm(10,3),jmn(400,3),Zjmn(400),Zln(400), s3(20)
      dimension Kln(400,2),Mln(400,2),ZMln(400),xyz(10000,3),se(10000)
      dimension Euler(10000,3), ncalcspin(10000), f2(20), f3(20)
      common/pHmp/pi,z1,Mln,ZMln,LnMax    !! common for main programm and 2 subprograms
      common/pEmp2/ZE,se,xyz,Euler,T,Nj,Njm,Ns  !
      common/file/f2,f3             !! common blocks for Hmp and main program
      common/file1/f1  
      common/sign/Fl_Abort
      common/cnt/iCnt

      open(7,file='MPE.d',status='old')

      read(7,*)f6                ! history
      read(7,*)To,dT,nT          ! To=init.temp.dT=increm.nT=nmbr.of temp.:dT
      read(7,*) (s3(i),i=1,nT)
      read(7,*) s4               !temporary file for MP configurations

      read(7,*)f1                ! Qlmn results
      read(7,*)Nj                ! Number of all possible l
      read(7,*) (f2(i),i=1,Nj)   ! j-s
      read(7,*)Njm               ! Number of jm par composing a MP
      read(7,*) (f3(i),i=1,Njm)  ! j-s and m-s

      read(7,*)f4                ! koordinates of MP

      read(7,*) Ns               ! Number of MP

      read(7,*)nmcs              ! numb.of MC steps
      read(7,*)nmcr              ! numb.of MC steps for reading (writing?)
      read(7,*)init              ! initial random number
      read(7,*) iconfig          ! iconfig=0, random distr.| iconfig=1, given configuration.
      close(7)

      pi=3.14159265350D0
      z1=(0.00D+00,-1.00D+00)

c Initial configuration of the random generator...
c      do i=1,init
         x=drandm(init)
c      enddo

      open(4,file=f4,status='old')
      if(iconfig.eq.0) then         !! New random conf. 0,
                                    !! homogeneous distribution of angles
         do i=1,Ns
            read(4,*) (xyz(i,j),j=1,3)!coord,Eu Winkels & energies
            Euler(i,1) = 0.d0                            !! Euler angles - Alpha
            Euler(i,2) = acos(drandm(0) * 2.d0 - 1.d0)  !!       - Beta
            Euler(i,3) = drandm(0) * 2.d0 * pi          !!       - Gamma
            se(i) = 0.
            open(3,file=f6,status='unknown') !! History.will be from begin
         enddo
      else                          !! Continuation from a given step,
                                    !! angles, energies and coordinated are read from the file
         do i=1,Ns
            read(4,*) (xyz(i,j),j=1,3),(Euler(i,j),j=1,3) !coord,Eu Winkels & energies
            open(3,file=f6,status='unknown') !! History.will be continued..
         enddo
      endif
      close(4)

      write(3,*)'MultiPole ->Vector Spins on Triang.lattice.Sample=Disc'
c      write(3,*) 'Here exchange+dipolar+Aniso'
c      write(3,*)'No cutoff : Interactions all over the sample'
      write(3,*) 'ns=',ns,' To=',To, 'random=',RFlag
c      write(3,*) 'Aniso1=',A1,'Aniso2=',A2,' Field Hz_|_plane =',Hz
c      write(3,*) 'AnisoIN=',Ax,' Field Hin||_plane =',Hin,'alpha=',alpha
c      write(3,*) 'AnisoIN=',Ax,' Field Hin||_plane =',Hin
      write(3,*) 'init=',init, ', nmcs=',nmcs, ', tempfile=',s4
      write(3,*) nT,' Results_Files :'
      write(3,*) (s3(i),i=1,nT)
c      write(3,*) ' Config.current File : ',s4
      close(3)

c      open(1,file=f1,status='unknown')
   
      nsig = signal( 2, On_Intrp, -1)           !! Look for signal SIGINT..
      nsig = signal(15, On_Intrp, -1)           !! ...and signal SIGTERM as well

      call Emp_ini ()                           !! Initialization- 
                                                !! fill on the arrays...

C ---------- Cycle for T (BEGIN) ------------
      do iT=1,Nt
         T = To + dT * (iT - 1)
         if (T.lt.0.05)  T = 1.d-2
        open(3,file=f6,status='old',access='append')
        write(3,*) ' '
        write(3,*) 'Temperature=',T
        write(3,*) 'Iter. Energy/spin  <sz2>/spin'
        close(3)

C ---------- Cycle for MC-steps (BEGIN) -----------
         do iMC=1, NMCs
            Ei  = 0.d0
c			Ei2 = 0.d0
            sz2 = 0.d0
			sx = 0.d0
			sy = 0.d0
			sx2 = 0.d0
			sy2 = 0.d0
C ---------- Cycle for i-spins (BEGIN) ------------
            if(RFlag) then             !! If we had flag "-r"...
C Initialization of the array of calculated spins...            
               do ispin=1,ns
                  ncalcspin(ispin) = ispin
               enddo
            endif
C  =====> MAIN CYCLE(ns) <=====
C  Here we define whether the MP will be addressed randomly or in the list order..
            do ii = 1, Ns
               if(RFlag) then          !! If flag is "-r", then i-th MP
                                       !! is taken randomly..
                  nspin = ns-ii+1
C       Now we take a MP, which has not been calculated yet...
                  ispin = drandm(0) * nspin + 1.d0
                  i = ncalcspin(ispin)
C       Reordering of the MP, which were not used yet:                 
                  konspin = nspin - 1  !! For optimization of DO ...
                  do j=ispin,konspin
                     ncalcspin(j) = ncalcspin(j+1)
                  enddo
c                 print *, 'i =',i     !! Just control - should be deleted!
               else                    !!...if not, 
                  i = ii               !! going over the MP in the list order..
               endif

C   Generation of random angles and Call Emp..
               Alph = 1.d0
               Bet  = acos(drandm(0) * 2.d0 - 1.d0)
               Gamm = drandm(0) * 2.d0 * pi
               call Emp(i, Alph, Bet, Gamm)
c               se(i) = ZE                       !! just for any case
               
c		   if((xyz(i,1).gt.147).and.(xyz(i,1).lt.153).and.
c     1 (xyz(i,2).gt.147).and.(xyz(i,2).lt.153)) then
			   
			   Ei = Ei + se(i)        !! Energy and vertcal component calc.
c			   Ei2 = Ei2 +se(i)*se(i)
               sz2 = sz2 + (dcos(Euler(i,2)) * dcos(Euler(i,2)))
			   sx = sx +(dsin(Euler(i,2))*dcos(Euler(i,3)))
			   sy = sy +(dsin(Euler(i,2))*dsin(Euler(i,3)))
			   sx2 = sx2 +(dsin(Euler(i,2)) * dcos(Euler(i,3)))
     1 *(dsin(Euler(i,2)) * dcos(Euler(i,3)))
			   sy2 = sy2 +(dsin(Euler(i,2)) * dsin(Euler(i,3)))
     1 *(dsin(Euler(i,2)) * dsin(Euler(i,3)))
	 
c	        endif
            enddo
C ---------- Cycle for i-spins (END) ------------

            if((mod(imc,nmcr).eq.0) .OR. Fl_Abort) then  !! nmcr= numb.of MC steps for reading
			   EiAV = Ei/Ns/2.d0
c			   Ei2AV = Ei2/Ns/4.d0
			   EiAV2=EiAV**2
			   sxAV=sx/Ns
			   syAV=sy/Ns
			   sxAV2=sxAV*sxAV
			   syAV2=syAV*syAV
			   sx2AV=sx2/Ns
			   sy2AV=sy2/Ns
c			   Cv=(Ei2AV-EiAV2)/(T**2)
c			   SusX=(sx2AV-sxAV2)/T
c			   SusY=(sy2AV-syAV2)/T
               IMC1=(it-1)*nmcs+imc
               open(3,file=f6,status='old',access='append')
               write(3,5)IMC1, EiAV, sx2AV, sy2AV, sz2AV
               close(3) 

               open(4,file=s4,status='unknown')
c              write(4,*)'Here current config.Final config.will be in :',
c     +                                                          s3(it)
c              write(4,*)'Iter. Energy/spin  <sz2>/spin'
c              write(4,3)imc,ei/ns/2.d0,sz2/ns
c              write(4,*)' '
               do i=1,ns        !! Temp copy of data (on any key)
                  write(4,3) (xyz(i,ii),ii=1,3),(Euler(i,j),j=1,3),se(i)
               enddo
               close(4)
            endif  ! end of loop on nmcr

            if(Fl_Abort) exit    !! Stop-signal was obtained, exit
                                 !! cycle on i and go to write the obtained configutration..
          enddo
C ---------- Cycle for MC-steps (END) ------------

C ++++++  Results  ++++++
         open(30,file=s3(it),status='unknown')
         do i=1,Ns
             write(30,3) (xyz(i,ii),ii=1,3), (Euler(i,j),j=1,3), se(i)
         enddo
         close(30)
         if(Fl_Abort) goto 666   !! .. i ends it job
      enddo
C ---------- Cycle for T (END) ------------

C ++++++  Results  ++++++
      open(3,file=f6,status='old',access='append')
      write(3,*)'Energy/spin of final configuration: ',Ei/Ns/2.d0
      close(3)
      close(1)
      print*,iCnt

1     format(F10.6,F10.6'*i')
2     format(3I5,F10.6,F10.6'*i')
3     format(7E15.7)
5     format(I5,6E15.7)
666   end


C++++++ Function nfact - factorial calculation (quick, without recursion!) +++++

      integer Function nfact(n)
      integer m,n

      m=1
      do i=1,n
         m=m*i
      enddo

      nfact=m

      end


C+++++++++++++++++ Function for  Wigner coefficients ++++++++++

      real*8 Function Set(j,m,n,beta)
      implicit real*8(a-h,o-z)
      pi=3.1415926535
      Set=0.0  

      do k=0,j+n

         et=0.0

         nsig=j+m+k
         n1=nfact(j+n)
         n2=nfact(j-n)
         n3=nfact(j+m)
         n4=nfact(j-m)
         nd1=j+m-k
         nd2=k
         nd3=j+n-k
         nd4=k-n-m

         if((nd1.lt.0).or.(nd2.lt.0)
     1  .or.(nd3.lt.0).or.(nd4.lt.0)) then
            et=0.0
         else

            nd1=nfact(nd1)
            nd2=nfact(nd2)
            nd3=nfact(nd3)
            nd4=nfact(nd4)

            a=n1*n2*n3*n4
            anum=dsqrt(a)
            nden=nd1*nd2*nd3*nd4
            nspow=2*j-2*k+n+m ! sin power
            ncpow=2*k-n-m

            if((beta.eq.0.).and.(nspow.eq.0)) then
               if(n.eq.m) then
                  et=1.0
               else
                  et=0.0
               endif
               goto 25 

            else

               if((beta.eq.pi).and.(ncpow.eq.0)) then
                  if((n+m).eq.0) then
                     et=(-1.)**(j+m)
                  else
                     et=0.0
                  endif
                  goto 25

               else

                  et=(-1)**nsig*(anum/nden)*dcos(beta/2)**ncpow
     1                                     *dsin(beta/2)**nspow
               endif
            endif
         endif
25       Set=Set+et
C        print*,et,Set
      enddo

      end


C++++++++++++++++++++++++++Functions for transformation of+++++++++++++++++++++++++
C ++++++++++++++++++++++++XX,YY,ZZ   into   Theta and Phi++++++++++++++++++++++++

      real*8 Function Theta(zz,r)
      implicit real*8(a-h,o-z)
      pi=3.1415926535

      theta=acos(zz/r)

      end


      real*8 Function Phi(xx,yy,zz,r)
      implicit real*8(a-h,o-z)
      pi=3.1415926535

      sintheta=dsqrt(1.d0-(zz/r)**2)
      if(sintheta.eq.0.0) then
         if(xx.gt.0.0) then
            cosphi=1.d0
         else
            cosphi=-1.d0
         endif
      else
         cosphi=xx/r/sintheta
      endif

      if(cosphi.gt.1.0)  cosphi=1.d0
      if(cosphi.lt.-1.0) cosphi=-1.d0
      yphi=dacos(cosphi)
c      if(yy.eq.0.d0) then
c         phi=0.0
      if((xx.eq.0.d0).and.(yy.gt.0.d0)) then
         phi=pi/2.
      elseif((xx.eq.0.d0).and.(yy.lt.0.d0)) then
         phi=3.*pi/2.
c      elseif((yy.gt.0.d0).and.(xx.gt.0.d0)) then
c         phi=yphi
c      elseif((yy.gt.0.d0).and.(xx.lt.0.d0)) then
c         phi=pi-yphi
      elseif((yy.lt.0.d0).and.(xx.gt.0.d0)) then
         phi=2*pi-yphi
      elseif((yy.lt.0.d0).and.(xx.lt.0.d0)) then
         phi=2*pi-yphi
      else
         phi=yphi
      endif
c     print*,phi
      end


C++++++++++++++++++++++++++Spherical harm1++++++++++++++++++++++++++++++++++++++

      complex*16 Function Ylm(theta,phi,lalb,mamb)
      implicit real*8(a-h,o-z)
      complex*16 z1,Y1
      pi=3.1415926535
      z1=(0.00D+00,1.00D+00)

      if((lalb.eq.1).and.(mamb.eq.1))then
         Ylm=-dsqrt(3./8./pi)*sin(theta)*
     1                      cdexp(z1*phi)

      elseif((lalb.eq.1).and.(mamb.eq.0))then
         Ylm=dsqrt(3./4./pi)*cos(theta)

      elseif((lalb.eq.2).and.(mamb.eq.2))then
         Ylm=0.25d0*dsqrt(15./2./pi)*sin(theta)*
     1      sin(theta)*cdexp(2.d0*z1*phi)

      elseif((lalb.eq.2).and.(mamb.eq.1))then
         Ylm=-dsqrt(15./8./pi)*sin(theta)*
     1      cos(theta)*cdexp(z1*phi)

      elseif((lalb.eq.2).and.(mamb.eq.0))then
         Ylm=dsqrt(5./4./pi)*((3./2.)*cos(theta)*
     1      cos(theta)-0.5d0)

      elseif((lalb.eq.2).and.(mamb.eq.-1))then
         Y1=-dsqrt(15./8./pi)*sin(theta)*
     1      cos(theta)*cdexp(z1*phi)
         Ylm=conjg(Y1)*(-1.d0)

      elseif((lalb.eq.2).and.(mamb.eq.-2))then
         Y1=0.25d0*dsqrt(15./2./pi)*sin(theta)*
     1      sin(theta)*cdexp(2.d0*z1*phi)
         Ylm=conjg(Y1)
      
      elseif((lalb.eq.4).and.(mamb.eq.-4))then
         Ylm=(3.d0*dsqrt(35.d0/2.d0/pi)*dsin(theta)**4)
     1   /16.d0/cdexp(4.d0*z1*phi)

      elseif((lalb.eq.4).and.(mamb.eq.-3))then
         Ylm= (3.d0*dsqrt(35.d0/pi)*dcos(theta)*dsin(theta)**3)   
     1   /8.d0/cdexp(3.d0*z1*phi)

      elseif((lalb.eq.4).and.(mamb.eq.-2))then
         Ylm= (3.d0*dsqrt(5.d0/2.d0/pi)*(-1.d0+7.d0*dcos(theta)**2)
     1   *dsin(theta)**2)/8.d0/cdexp(2.d0*z1*phi) 
     
      elseif((lalb.eq.4).and.(mamb.eq.-1))then
         Ylm= (3.d0*dsqrt(5.d0/pi)*dcos(theta)*
     1    (-3.d0+7.d0*dcos(theta)**2)*dsin(theta))/8.d0/cdexp(z1*phi)

      elseif((lalb.eq.4).and.(mamb.eq.0))then
         Ylm= 3.d0*(3.d0-30.d0*dcos(theta)**2+
     1    35.d0*dcos(theta)**4)/16.d0/dsqrt(pi)

      elseif((lalb.eq.4).and.(mamb.eq.1))then
         Ylm= (-3.d0*cdexp(z1*phi)*dsqrt(5.d0/pi)*dcos(theta)*
     1    (-3.d0+7.d0*dcos(theta)**2)*dsin(theta))/8.d0

      elseif((lalb.eq.4).and.(mamb.eq.2))then
         Ylm= (3.d0*cdexp(2.d0*z1*phi)*dsqrt(5.d0/2.d0/pi)
     1    *(-1.d0+7.d0*dcos(theta)**2)*dsin(theta)**2)/8.d0

      elseif((lalb.eq.4).and.(mamb.eq.3))then
         Ylm= (-3.d0*cdexp(3.d0*z1*phi)*dsqrt(35.d0/pi)*dcos(theta)
     1   *dsin(theta)**3)/8.d0
     
      elseif((lalb.eq.4).and.(mamb.eq.4))then
         Ylm=(3.d0*cdexp(4.d0*z1*phi)*dsqrt(35.d0/2.d0/pi)
     1    *dsin(theta)**4)/16.d0
      elseif((lalb.eq.3).and.(mamb.eq.3))then
         Ylm=-(cdexp(3.d0*z1*phi)*dsqrt(35.d0/pi)
     1    *dsin(theta)**3)/8.d0
      elseif((lalb.eq.3).and.(mamb.eq.2))then
         Ylm=(cdexp(2.d0*z1*phi)*dsqrt(105.d0/2./pi)*dcos(theta)
     1    *dsin(theta)**2)/4.d0
      elseif((lalb.eq.3).and.(mamb.eq.1))then
         Ylm=-(cdexp(z1*phi)*dsqrt(21.d0/pi)*(-1.+5.*dcos(theta)**2)
     1    *dsin(theta))/8.d0
      elseif((lalb.eq.3).and.(mamb.eq.0))then
         Ylm=dsqrt(7.d0/pi)*(-3.d0*dcos(theta)
     1    +5.d0*dcos(theta)**3)/4.d0
      elseif((lalb.eq.3).and.(mamb.eq.-1))then
         Ylm=(dsqrt(21.d0/pi)*(-1.+5.*dcos(theta)**2)
     1    *dsin(theta))/8.d0/cdexp(z1*phi)
      elseif((lalb.eq.3).and.(mamb.eq.-2))then
         Ylm=(dsqrt(105.d0/2./pi)*dcos(theta)
     1    *dsin(theta)**2)/4.d0/cdexp(2.d0*z1*phi)
      elseif((lalb.eq.3).and.(mamb.eq.-3))then
         Ylm=(dsqrt(35.d0/pi)
     1    *dsin(theta)**3)/8.d0/cdexp(3.d0*z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.6))then
         Ylm=(cdexp(6.d0*z1*phi)*dsqrt(3003.d0/pi)
     1    *dsin(theta)**6)/64.d0
      elseif((lalb.eq.6).and.(mamb.eq.5))then
         Ylm=-(3.d0*cdexp(5.d0*z1*phi)*dsqrt(1001.d0/pi)*dcos(theta)
     1    *dsin(theta)**5)/32.d0
      elseif((lalb.eq.6).and.(mamb.eq.4))then
         Ylm=(3.d0*cdexp(4.d0*z1*phi)*dsqrt(91.d0/2./pi)*(-1.+11.
     1    *dcos(theta)**2)*dsin(theta)**4)/32.d0
      elseif((lalb.eq.6).and.(mamb.eq.3))then
         Ylm=-(cdexp(3.d0*z1*phi)*dsqrt(1365.d0/pi)*dcos(theta)
     1    *(-3.+11.*dcos(theta)**2)*dsin(theta)**3)/32.d0
      elseif((lalb.eq.6).and.(mamb.eq.2))then
         Ylm=(cdexp(2.d0*z1*phi)*dsqrt(1365.d0/pi)*(1.-18.
     1    *dcos(theta)**2+33.*dcos(theta)**4)*dsin(theta)**2)/64.d0
      elseif((lalb.eq.6).and.(mamb.eq.1))then
         Ylm=-(cdexp(z1*phi)*dsqrt(273.d0/2./pi)*dcos(theta)*(5.-30.
     1    *dcos(theta)**2+33.*dcos(theta)**4)*dsin(theta))/16.d0
      elseif((lalb.eq.6).and.(mamb.eq.0))then
         Ylm=(dsqrt(13.d0/pi)*(-5.+105.*dcos(theta)**2-315.*
     1    *dcos(theta)**4+231.*dcos(theta)**6))/32.d0
      elseif((lalb.eq.6).and.(mamb.eq.-1))then
         Ylm=(dsqrt(273.d0/2./pi)*dcos(theta)*(5.-30.
     1    *dcos(theta)**2+33.*dcos(theta)**4)*dsin(theta))
     2    /16.d0/cdexp(z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.-2))then
         Ylm=(dsqrt(1365.d0/pi)*(1.-18.*dcos(theta)**2+
     1   33.*dcos(theta)**4)*dsin(theta)**2)/64.d0/cdexp(2.d0*z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.-3))then
         Ylm=(dsqrt(1365.d0/pi)*dcos(theta)*(-3.+11.
     1    *dcos(theta)**2)*dsin(theta)**3)/32.d0/cdexp(3.d0*z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.-4))then
         Ylm=(3.d0*dsqrt(91.d0/2./pi)*(-1.+11.
     1    *dcos(theta)**2)*dsin(theta)**4)/32.d0/cdexp(4.d0*z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.-5))then
         Ylm=(3.d0*dsqrt(1001.d0/pi)*dcos(theta)
     1    *dsin(theta)**5)/32.d0/cdexp(5.d0*z1*phi)
      elseif((lalb.eq.6).and.(mamb.eq.-6))then
         Ylm=(dsqrt(3003.d0/pi)
     1    *dsin(theta)**6)/64.d0/cdexp(6.d0*z1*phi)
         
      else
         Ylm=0.d0
      endif
      end

C+++++++++++ Calculation of Qln++++++++++++++++++++++++++++++++++++
C   Returns array Mln with all possible combinations of l,n
C   and complex array ZMln with Qln corresponding to each l,n pair

      subroutine Hmp(Nj,Njm)
      implicit real*8(a-h,o-x)
      implicit complex*16(y,z)
      character*20 f2,f3
      dimension Lj(10),jm(10,3),jmn(400,3),Zjmn(400),Zln(400)
      dimension Kln(400,2),Mln0(400,2),ZMln0(400),f2(20),f3(20)
      dimension Mln(10000,400,2),ZMln(10000,400),LnMax(10000)
      common/pHmp/pi,z1,Mln,ZMln,LnMax
      common/file/f2,f3 ! obshii dlya Hmp i osnovnoi programmi
      common/pEmp1/alpha,beta,gamma,ZMln0,lnmax0,Mln0

C     open(2,file=f2,status='old') 

      IS=0

      do ij=1,Nj                          !! cycle on all l
         read(f2(ij),*) Lj(ij)
         j=Lj(ij)

C        open(3,file=f3,status='old') 
         do ijm=1,Njm                     !! cycle po m
            read(f3(ijm),*) (jm(ijm,ii),ii=1,3) !! j,m,Qjm

            if(j.eq.jm(ijm,1))then
               m=jm(ijm,2)
               Qjm=jm(ijm,3)

               do n=-j,j            !! cycle to calculate Dnm(l) of 1 MP
                  i=j-n+1           !! serial number inside of the cycle
                  IS=IS+1           !! general serial number
  
                  Zwigner=Set(j,m,n,beta)*cdexp(z1*(m*alpha+n*gamma)) !Wigner koeff

                  jmn(IS,1)=j
                  jmn(IS,2)=m
                  jmn(IS,3)=n
                  Zjmn(IS)=Zwigner*Qjm

               enddo
            endif
         enddo                      !! end of cycle on m

C        close(3)
      enddo                         !! end of cycle by j (l)
      ISmax=IS

C     close(2)

!!!!!! Creation of an array, in which Q with identical ln1 are summed

      ln=0           ! initial length of the array ln
      Kln(1,1)=0     !
      Kln(1,2)=0     ! initialization of the first list element
      Zln(1)=0.      ! initialization of the complex array	
      Mln0(1,1)=0     !
      Mln0(1,2)=0     ! initialization of the first element of Mln0
      ZMln0(1)=0.     ! initialization of the complex array

!!! lookup all possible lmn combinations
!!! ISmax - number of all possible lmn combinations
c     print*,ISmax
      do is=1,ISmax !cycle
         Zqln=Zjmn(is) !Qln for the first ln-kombinacion
         do js=1,ISmax

            if(is.ne.js)then
               if((jmn(is,1).eq.jmn(js,1)).and.
     1            (jmn(is,3).eq.jmn(js,3))) then

                     Zqln=Zqln+Zjmn(js)! summed Qln for identical ln
                   
               endif
            endif  
         enddo

         Kln(is,1)=jmn(is,1)     ! l 
         Kln(is,2)=jmn(is,3)     ! n
         Zln(is)=Zqln            ! Total Qln for identical ln

c        print*,Kln(is,1),Kln(is,2),Zln(is)

!!!here identical lines can appear
!!!!!!  this part delete them

         do ii=1,ln
           if((Kln(is,1).eq.Mln0(ii,1)).and.
     *        (Kln(is,2).eq.Mln0(ii,2))) then
               goto 33
            endif
         enddo
         Mln0(ii,1)=Kln(is,1)   ! l 
         Mln0(ii,2)=Kln(is,2)    ! n
         ZMln0(ii)=Zln(is)       ! Total Qln for identical ln
c        print*,Mln0(ii,1),Mln0(ii,2),ZMln0(ii)
         ln=ii

33    enddo
      lnmax0=ln

      return
2     format(3I5,F10.6,F10.6'*i')
      end

C  ++++++++++++++++++++++++++    EMP     +++++++++++++++++++++++++++
C   +++++   To calculate the interaction energy
C    +++    Parameters - number of MP, and 3 new random angles
C     +     trial (old angles are in Euler() )..

      subroutine Emp(i,Alph,Bet,Gamm)

      implicit real*8(a-h,o-x)
      implicit complex*16(y,z)
      character*20 f1
      dimension Mln0(400,2),xyz(10000,3),se(10000),ForRecalc(10000,4)
      dimension Euler(10000,3)
      dimension Mln1(400,2),ZMln0(400),ZMln1(400),Mln2(400,2),ZMln2(400)
      dimension Mln(10000,400,2),ZMln(10000,400),LnMax(10000)
      common/pHmp/pi,z1,Mln,ZMln,LnMax       ! common part with Hmp and with main prog.
      common/pEmp1/alpha,beta,gamma,ZMln0,lnmax0,Mln0    ! common with Hmp only
      common/pEmp2/ZE,se,xyz,Euler,T,Nj,Njm,Ns ! common with main only
      common/file1/f1
      common/recalc/ForRecalc    ! common with "EmpRecalc"

      ZE_new = 0.

C  New values for i-th MP using Hmp..
      alpha = Alph
      beta  = Bet
      gamma = Gamm
      call Hmp(Nj,Njm)
      do l=1,lnmax0
         Mln1(l,1)=Mln0(l,1)
         Mln1(l,2)=Mln0(l,2)
         ZMln1(l)=ZMln0(l)
      enddo

C  ...and calculate with all j-th MP..
      do j=1,Ns

         if(i.eq.j) cycle        !! avoid self-interaction
         lnmax0=LnMax(j)
         do l=1,lnmax0
            Mln2(l,1)=Mln(j,l,1)
            Mln2(l,2)=Mln(j,l,2)
            ZMln2(l)=ZMln(j,l)
         enddo
C  All distances are calculated only 1(!) time
C  for the actual angles as well as for new angles. Ekonomiya!! ;-)
         rxx=xyz(j,1)-xyz(i,1)
         ryy=xyz(j,2)-xyz(i,2)
         rzz=xyz(j,3)-xyz(i,3)
         r2=rxx*rxx+ryy*ryy+rzz*rzz
         r=dsqrt(r2)

         th=Theta(rzz,r)
         ph=Phi(rxx,ryy,rzz,r)
C  Memorize "r", "th" i "ph"..
         ForRecalc(j,1) = r
         ForRecalc(j,2) = th
         ForRecalc(j,3) = ph

         ZEln_new=0.

         do l=1,lnmax0
            do n=1,lnmax0

C        Calculate all factorials under the sum

               Lab_new=Mln1(l,1)+Mln2(n,1)         !! La+Lb
               Nab_new=Mln1(l,2)+Mln2(n,2)         !! Na+Nb
               LNab1_new=nfact(Lab_new+Nab_new)
               LNab2_new=nfact(Lab_new-Nab_new)
               LNa1_new=nfact(Mln1(l,1)+Mln1(l,2)) !! La+Na
               LNa2_new=nfact(Mln1(l,1)-Mln1(l,2)) !! La-Na
               LNb1=nfact(Mln2(n,1)+Mln2(n,2))     !! Lb+Nb
               LNb2=nfact(Mln2(n,1)-Mln2(n,2))     !! Lb-Nb
               LNab_new=LNa1_new*LNa2_new*LNb1*LNb2

               c1_new=LNab1_new*LNab2_new/LNab_new
               c2_new=dsqrt(c1_new)
               c3=(-1.)**Mln2(n,1)
               Zc4_new=c2_new*c3
               ZY_new=Ylm(th,ph,Lab_new,Nab_new)

               AL_new=dsqrt(4.*pi/(2.*Mln1(l,1)+1.))
               BL=dsqrt(4.*pi/(2.*Mln2(n,1)+1.))
               ABL_new=dsqrt(4.*pi/(2.*Mln1(l,1)+2.*Mln2(n,1)+1.))
               ZAB_new=AL_new*BL*ABL_new

               rr_new=r**(Mln1(l,1)+Mln2(n,1)+1.)

               Zqq_new=ZMln1(l)*ZMln2(n)*Zc4_new*ZAB_new*ZY_new

               ZEln_new=ZEln_new+Zqq_new/rr_new
               ForRecalc(j,4) = ZEln_new        !! memorize
                                                !! for possible re-calculation.
            enddo
         enddo

         ZE_new=ZE_new+ZEln_new
      enddo
c        write(1,3) ZE 

C  Now we have 2 energy values (ZE_cur i ZE_ new) for current anf new
C   Alpha, Beta i Gamma, compare them and sccept lower one..

c      se(i) = ZE_cur                  !! for any case
c      se_old = ZE_cur
      se_new = ZE_new
c      ZE = ZE_cur
      if (se_new .le. se(i)) then      !! New energy is lower, memorise it
                                       !! and correponding angles
         call EmpRecalc(i)             !! rewrite and recalculate all components of
         se(i)      = ZE_new           !!  se() for all  j-th Q..
         Euler(i,1) = Alph
         Euler(i,2) = Bet
         Euler(i,3) = Gamm
         do l=1,lnmax0
            Mln(i,l,1)=Mln1(l,1)
            Mln(i,l,2)=Mln1(l,2)
            ZMln(i,l) =ZMln1(l)
         enddo
         LnMax(i)=lnmax0
c         ZE=ZE_new
      else                            !! check Boltzmann probability
         x = dlog(drandm(0))
         p = (se(i) - se_new) / T
         if (p.gt.x) then             !! accept...
            call EmpRecalc(i)         !! re-address komponents of i-th Q...
            se(i)      = ZE_new
            Euler(i,1) = Alph
            Euler(i,2) = Bet
            Euler(i,3) = Gamm
            do l=1,lnmax0
               Mln(i,l,1)=Mln1(l,1)
               Mln(i,l,2)=Mln1(l,2)
               ZMln(i,l) =ZMln1(l)
            enddo
            LnMax(i)=lnmax0
c            ZE=ZE_new
         endif
      endif

      return
2     format(2I5,F10.6,F10.6'*I')
3     format(6E15.7,F10.6,F10.6'*I')        
      end

C +++++++++++++RANDOM GENERATOR by E.V. ++++++++++++++++++++++
C ------------------------------------------------------------
      real*8 Function drand_v(xxx)

      REAL brray(1)
      INTEGER IJ, KL, Input
      SAVE Input

      IJ = 1802
      KL = 9373

      If (Input .eq. 0) then
         call rmarin(IJ,KL)
         Input = 1
      EndIf

      call RANMAR(brray, 1)
      drand_v = brray(1)

      end

      subroutine RMARIN(IJ,KL)
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient
C length to complete an entire calculation with. For example, if sveral
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random
C number generator can create 900 million different subsequences -- with
C each subsequence having a length of approximately 10^30.
C
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0


      real U(97), C, CD, CM
      integer I97, J97
      logical TEST
      data TEST /.FALSE./
      common /raset1/ U, C, CD, CM, I97, J97, TEST

      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.
     *    KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value
     *between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and
     *30081'
            stop
      endif

      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169)

      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
3        continue
         U(ii) = s
2     continue

      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0

      I97 = 97
      J97 = 33

      TEST = .TRUE.
      return
      end

      subroutine ranmar(RVEC, LEN)
C This is the random number generator proposed by George Marsaglia in
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.
      REAL RVEC(*)
      real U(97), C, CD, CM
      integer I97, J97
      logical TEST
      common /raset1/ U, C, CD, CM, I97, J97, TEST

      integer ivec

      if( .NOT. TEST ) then
         print '(A)',' Call the init routine (RMARIN) before calling RAN
     *MAR'
         stop
      endif

      do 100 ivec = 1, LEN
         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         uni = uni - C
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RVEC(ivec) = uni
100   continue
      return
      end

!!
!!  Signal handler routine
!!
      INTEGER*4 FUNCTION On_Intrp (sig_num)
      INTEGER*4 sig_num
      logical Fl_Abort
      common/sign/Fl_Abort

C      nsig = signal(2,,1)            !! Alle weitere Sygnale werden ignorieren sein - M.V.
      Fl_Abort = .TRUE.
      On_Intrp = 1
      END

C
C   +++++  EMP_ini initial calculation of energy       ++++++
C +++++++++ of starting configuration, and initial values of all arrays  ++++++++++
C
      subroutine Emp_ini ()

      implicit real*8(a-h,o-x)
      implicit complex*16(y,z)
      dimension Mln0(400,2),ZMln0(400),xyz(10000,3),se(10000)
      dimension Euler(10000,3), Mln2(400,2),Zmln2(400)
      dimension Mln(10000,400,2),ZMln(10000,400),LnMax(10000)
      common/pHmp/pi,z1,Mln,ZMln,LnMax       ! common with Hmp, as well as with main
      common/pEmp1/alpha,beta,gamma,ZMln0,lnmax0,Mln0 ! common with Hmp
      common/pEmp2/ZE,se,xyz,Euler,T,Nj,Njm,Ns  ! common with main only

C    "Hmp" and its coefficients...
      do i=1,Ns
         alpha=Euler(i,1)
         beta =Euler(i,2)
         gamma=Euler(i,3)
         call Hmp(Nj,Njm)
         do l=1,lnmax0
            Mln(i,l,1)=Mln0(l,1)
            Mln(i,l,2)=Mln0(l,2)
            ZMln(i,l) =ZMln0(l)
         enddo
         LnMax(i)=lnmax0
      enddo

C  Cycle on all sites to calculate their energies se(i)...

      do i=1,Ns

      se(i) = 0.

C  Values of Mln i ZMln for all Q we take from corresponding arrays (without calling Hmp!)..
      lnmax0=LnMax(i)
      do l=1,lnmax0
         Mln0(l,1)=Mln(i,l,1)
         Mln0(l,2)=Mln(i,l,2)
         ZMln0(l)=ZMln(i,l)
      enddo

C  ...calculate with all j-th Q..
      do j=1,Ns

         if(i.eq.j) cycle        !! no self-energy
         lnmax0=LnMax(j)
         do l=1,lnmax0
            Mln2(l,1)=Mln(j,l,1)
            Mln2(l,2)=Mln(j,l,2)
            ZMln2(l)=ZMln(j,l)
         enddo
         rxx=xyz(j,1)-xyz(i,1)
         ryy=xyz(j,2)-xyz(i,2)
         rzz=xyz(j,3)-xyz(i,3)
         r2=rxx*rxx+ryy*ryy+rzz*rzz
         r=dsqrt(r2)

         th=Theta(rzz,r)
         ph=Phi(rxx,ryy,rzz,r)

         ZEln_cur=0.

         do l=1,lnmax0
            do n=1,lnmax0

C        All factorials under the summation

               Lab_cur=Mln0(l,1)+Mln2(n,1)          !! La+Lb
               Nab_cur=Mln0(l,2)+Mln2(n,2)          !! Na+Nb
               LNab1_cur=nfact(Lab_cur+Nab_cur)             !!
               LNab2_cur=nfact(Lab_cur-Nab_cur)
               LNa1_cur=nfact(Mln0(l,1)+Mln0(l,2))  !! La+Na
               LNa2_cur=nfact(Mln0(l,1)-Mln0(l,2))  !! La-Na
               LNb1=nfact(Mln2(n,1)+Mln2(n,2))      !! Lb+Nb
               LNb2=nfact(Mln2(n,1)-Mln2(n,2))      !! Lb-Nb
               LNab_cur=LNa1_cur*LNa2_cur*LNb1*LNb2

               c1_cur=LNab1_cur*LNab2_cur/LNab_cur
               c2_cur=dsqrt(c1_cur)
               c3=(-1.)**Mln2(n,1)
               Zc4_cur=c2_cur*c3
               ZY_cur=Ylm(th,ph,Lab_cur,Nab_cur)

               AL_cur=dsqrt(4.*pi/(2.*Mln0(l,1)+1.))
               BL=dsqrt(4.*pi/(2.*Mln2(n,1)+1.))
               ABL_cur=dsqrt(4.*pi/(2.*Mln0(l,1)+2.*Mln2(n,1)+1.))
               ZAB_cur=AL_cur*BL*ABL_cur

               rr_cur=r**(Mln0(l,1)+Mln2(n,1)+1.)

               Zqq_cur=ZMln0(l)*ZMln2(n)*Zc4_cur*ZAB_cur*ZY_cur

               ZEln_cur=ZEln_cur+Zqq_cur/rr_cur

            enddo
         enddo
C  Summing se(i)...
         se(i) = se(i) + ZEln_cur
      enddo    !! => "j"

      enddo    !! => "i"
      end

C  ++++++++++++++++++++++++    EmpRecalc     +++++++++++++++++++++++++
C   +++++ Re-calculation of all se(j) for all Q if the angles of i-th Q were changed
C     +
      subroutine EmpRecalc(i)

      implicit real*8(a-h,o-x)
      implicit complex*16(y,z)
      dimension Mln0(400,2),xyz(10000,3),se(10000),ForRecalc(10000,4)
      dimension Euler(10000,3)
      dimension Mln1(400,2),ZMln0(400),ZMln1(400),Mln2(400,2),ZMln2(400)
      dimension Mln(10000,400,2),ZMln(10000,400),LnMax(10000)
      common/pHmp/pi,z1,Mln,ZMln,LnMax       ! common Hmp as well as with main
      common/pEmp1/alpha,beta,gamma,ZMln0,lnmax0,Mln0    ! common  Hmp only
      common/pEmp2/ZE,se,xyz,Euler,T,Nj,Njm,Ns ! common with main
      common/recalc/ForRecalc    !common with "EmpRecalc"
      common/cnt/iCnt

      iCnt = iCnt+1
C  Current values for i-th MP berem from arrays..
      lnmax0=LnMax(i)
      do l=1,lnmax0
         Mln0(l,1)=Mln(i,l,1)
         Mln0(l,2)=Mln(i,l,2)
         ZMln0(l)=ZMln(i,l)
      enddo

C  ...calculate with all j-th MP..
      do j=1,Ns

         if(i.eq.j) cycle        !! except of itself
         lnmax0=LnMax(j)
         do l=1,lnmax0
            Mln2(l,1)=Mln(j,l,1)
            Mln2(l,2)=Mln(j,l,2)
            ZMln2(l)=ZMln(j,l)
         enddo
C  All distances and angles between i-th and j-th MP were already calculated
C  and saved in the array ForRecalc, therefore
C  we take "r", "th" and "ph" from there..  Ekonomiya!! ;-)
c         rxx=xyz(j,1)-xyz(i,1)
c         ryy=xyz(j,2)-xyz(i,2)
c         rzz=xyz(j,3)-xyz(i,3)
c         r2=rxx*rxx+ryy*ryy+rzz*rzz
c         r=dsqrt(r2)
         r=ForRecalc(j,1)
c         th=Theta(rzz,r)
         th=ForRecalc(j,2)
c         ph=Phi(rxx,ryy,rzz,r)
         ph=ForRecalc(j,3)

         ZEln_cur=0.

         do l=1,lnmax0
            do n=1,lnmax0

C        Factorials under the summation

               Lab_cur=Mln0(l,1)+Mln2(n,1)          !! La+Lb
               Nab_cur=Mln0(l,2)+Mln2(n,2)          !! Na+Nb
               LNab1_cur=nfact(Lab_cur+Nab_cur)             !!
               LNab2_cur=nfact(Lab_cur-Nab_cur)
               LNa1_cur=nfact(Mln0(l,1)+Mln0(l,2))  !! La+Na
               LNa2_cur=nfact(Mln0(l,1)-Mln0(l,2))  !! La-Na
               LNb1=nfact(Mln2(n,1)+Mln2(n,2))      !! Lb+Nb
               LNb2=nfact(Mln2(n,1)-Mln2(n,2))      !! Lb-Nb
               LNab_cur=LNa1_cur*LNa2_cur*LNb1*LNb2

               c1_cur=LNab1_cur*LNab2_cur/LNab_cur
               c2_cur=dsqrt(c1_cur)
               c3=(-1.)**Mln2(n,1)
               Zc4_cur=c2_cur*c3
               ZY_cur=Ylm(th,ph,Lab_cur,Nab_cur)

               AL_cur=dsqrt(4.*pi/(2.*Mln0(l,1)+1.))
               BL=dsqrt(4.*pi/(2.*Mln2(n,1)+1.))
               ABL_cur=dsqrt(4.*pi/(2.*Mln0(l,1)+2.*Mln2(n,1)+1.))
               ZAB_cur=AL_cur*BL*ABL_cur

               rr_cur=r**(Mln0(l,1)+Mln2(n,1)+1.)

               Zqq_cur=ZMln0(l)*ZMln2(n)*Zc4_cur*ZAB_cur*ZY_cur

               ZEln_cur=ZEln_cur+Zqq_cur/rr_cur

            enddo
         enddo

C  Edit se(j) - minus old value of ZEln_cur and plus
C  calculated in Emp new value (saved in ForRecalc(j,4) )...
         se(j) = se(j) - ZEln_cur + ForRecalc(j,4)
      enddo

      end
