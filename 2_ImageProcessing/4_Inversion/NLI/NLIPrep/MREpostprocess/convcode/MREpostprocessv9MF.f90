! Run this program in the inv directory to generate convergence information 
! and delete old iteration files. 

program MREpostprocessv8

  implicit none
  
  
  character*500 fname,mtrf,stm,meshindf,dispsetf
  character*2 mstr
  character*1 vstr
  integer ii,jj,kk,ll,numval,nset,maxitr,ind(3),ijunk,N,imtr,numprop,nummf
  logical fexist
  real*4, dimension(:,:), allocatable :: val,sval
  real*8 rjunk
  integer meshind(2,3)
  integer, dimension(:), allocatable :: valind,sense,reg
  real*8, dimension(:,:,:), allocatable :: meanv,stdv,medv,maxv,minv,lowpcnt,highpcnt
  real*8, dimension(:,:), allocatable :: avmtr
  integer nav
  logical avfinalitr,rcnind(6),mfind(100)
  
  
  nav=8
  avfinalitr=.false.
  mfind(:)=.false.
  mfind(1)=.true.
  
  nset=iargc()
  print *,'iargc() = ',iargc()
  
  if(nset.gt.0) then
    do ii=1,nset
      call getarg(ii,mtrf) ! Currently only processes one set
    enddo
  else 
    call system('ls *.RE.*0001.prop.01.mf.mtr')
    print *,'Enter name of file set to be processed >> '  
    nset=1 
    read(*,*) mtrf
  endif
 
  print *,'trim(mtrf) =',trim(mtrf)

  ind(1)=index(mtrf, '.RE..',.true.)  
  print *,'ind =',ind(1)
  stm=mtrf(1:ind(1)-1)
  print *,'File stem = ',trim(stm)

  !find out how many files there are  
  fexist=.true.
  maxitr=1
  do while (fexist)
    write(mtrf,'(A,A,i4.4,A)') trim(stm),'.RE..',maxitr,'.prop.01.mf.mtr'
    inquire(file=mtrf, exist=fexist)
    if (.not. fexist) then
      print *,'File ', trim(mtrf), "' not found, exiting loop"
    else
      !print *,'File: ', trim(mtrf),' found.'
      maxitr = maxitr + 1    
    endif
  enddo
  maxitr=maxitr-1
  
  
  !find number of properties
  print *,'finding number of properties'
  fexist=.true.
  numprop=2
  nummf=1
  do while (fexist)
    write(mtrf,'(A,A,i2.2,A)') trim(stm),'.RE..0001.prop.',numprop,'.mtr'
    inquire(file=mtrf, exist=fexist)
    if (.not. fexist) then
     write(mtrf,'(A,A,i2.2,A)') trim(stm),'.RE..0001.prop.',numprop,'.mf.mtr'
     inquire(file=mtrf, exist=fexist)
	 if (fexist) then
	  mfind(numprop)=.true.
	  nummf=nummf+1
	 endif
	endif	
	if (.not. fexist) then
      print *,'File ', trim(mtrf), "' not found, exiting loop"
    else
      !print *,'File: ', trim(mtrf),' found.'
      numprop = numprop + 1    
    endif
  enddo
  numprop=numprop-1
  print *,'numprop = ',numprop
  print *,'nummf = ',nummf

  if(maxitr.ge.nav) avfinalitr=.true.

  print *,'maxitr = ',maxitr
  allocate(valind(numprop),meanv(numprop+2*nummf,2,maxitr),stdv(numprop+2*nummf,2,maxitr),medv(numprop+2*nummf,2,maxitr))
  allocate(maxv(numprop+2*nummf,2,maxitr),minv(numprop+2*nummf,2,maxitr),lowpcnt(numprop+2*nummf,2,maxitr),highpcnt(numprop+2*nummf,2,maxitr))
  
  
  valind(1)=1
  do ii=1,numprop
	if (ii.lt.numprop) then
	  if(mfind(ii)) then
	    valind(ii+1)=valind(ii)+3
	  else
	    valind(ii+1)=valind(ii)+1
	  endif
	endif
  enddo 
  print *,'valind = ',(valind(ii),ii=1,numprop)
  meanv(:,:,:)=0.d0
  medv(:,:,:)=0.d0
  stdv(:,:,:)=0.d0
  maxv(:,:,:)=0.d0
  minv(:,:,:)=0.d0
  highpcnt(:,:,:)=0.d0
  lowpcnt(:,:,:)=0.d0
    
  ! Read mesh index
  write(fname,'(A,A)') trim(stm),'.meshind'
  open(unit=9,file=fname,status='old',position='append')
  backspace(9)
  backspace(9)
    read(9,*) (meshind(1,jj),jj=1,numprop)
    read(9,*) (meshind(2,jj),jj=1,numprop)
  close(9)
  print *,'meshind = ',meshind
  
  ! Read recon index
  write(fname,'(A,A)') trim(stm),'.reconind'
  open(unit=9,file=fname,status='old')
    read(9,*) (rcnind(ii),ii=1,2*numprop)
  close(9)
  print *,'rcnind = ',rcnind
  
  
  !Loop over all files
  do imtr=1,numprop
    do jj=1,2
      ! Read in node sensitivity
      write(fname,'(A,A,i2.2,A)') trim(stm),'.mtrmesh.',meshind(jj,imtr),'.nod'
      print *,'nodf = ',trim(fname)
      open(unit=9,file=fname,position='append')
        backspace(9)
	    read(9,*) N,rjunk,rjunk,rjunk,ijunk,ijunk
	    print *,'N=',N
        allocate(sense(N),reg(N),val(N,3),sval(N,3),avmtr(N,3))
	    rewind(9)
	    do kk=1,N
	     read(9,*) ijunk,rjunk,rjunk,rjunk,sense(kk),reg(kk)
	    enddo
      close(9)
      
      
      !Loop over all material files
      avmtr(:,:)=0.d0
      do ii=1,maxitr
	 if (mfind(imtr)) then
	  numval=3
          if(jj.eq.1) then
           write(mtrf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.RE..',ii,'.prop.',imtr,'.mf.mtr'
          else
           write(mtrf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.IM..',ii,'.prop.',imtr,'.mf.mtr'  
          endif
	 else
	  numval=1
          if(jj.eq.1) then
           write(mtrf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.RE..',ii,'.prop.',imtr,'.mtr'
          else
           write(mtrf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.IM..',ii,'.prop.',imtr,'.mtr'  
          endif
	 endif
        !print *,'mtrf = ',trim(mtrf)
      
        open(unit=9,file=mtrf,status='old')
        do kk=1,N
	     read(9,*) ijunk,(val(kk,ll),ll=1,numval)
	     if(ii.gt.maxitr-nav) then
	      do ll=1,numval
		   avmtr(kk,ll)=avmtr(kk,ll)+val(kk,ll)/dble(nav)
	      enddo
	     endif	
	     do ll=1,numval
		  val(kk,ll)=abs(val(kk,ll))*real(sense(kk))  ! It is easier to just take the absolute value.
	     enddo
	enddo
        close(9)
	if(rcnind(2*(imtr-1)+jj)) then
         !print *,'Calling sort: ','mtr ',imtr,' part ',jj
         do ll=1,numval
	  call sort(val(:,ll),N,sval(:,ll))
	 enddo         
	else
	 !print *,'mtr ',imtr,' part ',jj,' not reconstructed'
	endif

	do ll=1,numval
	 minv(valind(imtr)-1+ll,jj,ii)=0.d0
	 meanv(valind(imtr)-1+ll,jj,ii)=0.d0
     stdv(valind(imtr)-1+ll,jj,ii)=0.d0
     ind(ll)=0
	 do kk=1,N
	  if(sval(kk,ll).gt.0.d0) then
	    ind(ll)=ind(ll)+1
	    meanv(valind(imtr)-1+ll,jj,ii)=meanv(valind(imtr)-1+ll,jj,ii)+dble(sval(kk,ll))
	    if(ind(ll).eq.1) minv(valind(imtr)-1+ll,jj,ii)=dble(sval(kk,ll))
	    if(sval(kk,ll).lt.minv(valind(imtr)-1+ll,jj,ii)) minv(valind(imtr)-1+ll,jj,ii)=dble(sval(kk,ll))
	    stdv(valind(imtr)-1+ll,jj,ii)=stdv(valind(imtr)-1+ll,jj,ii)+dble(sval(kk,ll))*dble(sval(kk,ll))
	  endif
	 enddo
	
	 if(ind(ll).gt.0) then
	  meanv(valind(imtr)-1+ll,jj,ii)=meanv(valind(imtr)-1+ll,jj,ii)/dble(ind(ll))
	 else
	  meanv(valind(imtr)-1+ll,jj,ii)=0.d0
	 endif
	 maxv(valind(imtr)-1+ll,jj,ii)=sval(N,ll)
	 highpcnt((valind(imtr)-1+ll),jj,ii)=sval((N-nint(0.05*real(ind(ll)))),ll) ! 95th percentile
	 lowpcnt((valind(imtr)-1+ll),jj,ii)=sval((N-nint(0.95*real(ind(ll)))),ll) ! 5th percentile
	 medv((valind(imtr)-1+ll),jj,ii)=sval((N-nint(0.5*real(ind(ll)))),ll) ! 50th percentile
	 if(ind(ll).gt.0) then
	  stdv(valind(imtr)-1+ll,jj,ii)=dsqrt(abs(stdv(valind(imtr)-1+ll,jj,ii)/dble(ind(ll))-meanv(valind(imtr)-1+ll,jj,ii)*meanv(valind(imtr)-1+ll,jj,ii)))
	 else
	  stdv(valind(imtr)-1+ll,jj,ii)=0.d0
	 endif
	 !if(isnan(stdv(imtr,jj,ii))) stdv(imtr,jj,ii)=0.d0
	enddo !end ll loop over MF numvals
	
	!delete iteration files
	if(ii.lt.maxitr) then
		call system('rm '//trim(mtrf))
		if(jj.eq.1) then
		  write(meshindf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.RE..',ii,'.prop.',imtr,'.meshind'
		else
		  write(meshindf,'(A,A,i4.4,A,i2.2,A)') trim(stm),'.IM..',ii,'.prop.',imtr,'.meshind'		  
		endif
		
		call system('rm '//trim(meshindf))
		
		if((imtr.eq.1).and.(jj.eq.1)) then
		  write(dispsetf,'(A,A,i4.4,A)') trim(stm),'.',ii,'.dispset.1.dsp'
		  inquire(file=trim(dispsetf), exist=fexist)
           if(fexist) then
		    call system('rm '//trim(dispsetf))
		   endif
		endif
	endif
		
      enddo !end ii maxiter looop
      
      if(avfinalitr) then
        if (mfind(imtr)) then
	 if(jj.eq.1) then
          write(mtrf,'(A,A,i2.2,A,i4.4,A,i2.2,A)') trim(stm),'.avlast',nav,'.RE..',maxitr,'.prop.',imtr,'.mf.mtr'
         else
          write(mtrf,'(A,A,i2.2,A,i4.4,A,i2.2,A)') trim(stm),'.avlast',nav,'.IM..',maxitr,'.prop.',imtr,'.mf.mtr'  
         endif
	else
	 if(jj.eq.1) then
          write(mtrf,'(A,A,i2.2,A,i4.4,A,i2.2,A)') trim(stm),'.avlast',nav,'.RE..',maxitr,'.prop.',imtr,'.mtr'
         else
          write(mtrf,'(A,A,i2.2,A,i4.4,A,i2.2,A)') trim(stm),'.avlast',nav,'.IM..',maxitr,'.prop.',imtr,'.mtr'  
         endif
	endif
        open(unit=9,file=mtrf,status='unknown')
          do kk=1,N
	    write(9,*) kk,(avmtr(kk,ll),ll=1,numval)
	  enddo
        close(9)
      endif   
      deallocate(sense,reg,val,sval,avmtr)
      
    enddo
  enddo
  
  print *,'minv(1,1,:) ',(minv(jj,1,1),jj=1,numprop+2*nummf)
  
  
  ! Write output files
  do ii=1,9
    print *,meanv(1,1,ii),'+/-',stdv(1,1,ii),' max/min ',maxv(1,1,ii),'/',minv(3,1,ii)
  enddo
  
  do ii=1,9
    print *,'lowpcnt, medv, highpcnt ',lowpcnt(1,1,ii),medv(1,1,ii),highpcnt(1,1,ii)
  enddo

   
  do imtr=1,numprop
   if (mfind(imtr)) then
    write(mstr,'(i2.2)') imtr
    do ll=1,3
    write(vstr,'(i1.1)') ll
    open(unit=9,file=trim(stm) // '.prop.' // mstr // '.mf.' // vstr // '.RE.convinfo',status='unknown')
      print *,maxitr,imtr
      do ii=1,maxitr
        write(9,'(i7,7(1x,es11.4))') ii,meanv(valind(imtr)-1+ll,1,ii),stdv(valind(imtr)-1+ll,1,ii),minv(valind(imtr)-1+ll,1,ii),lowpcnt(valind(imtr)-1+ll,1,ii)&
	&,medv(valind(imtr)-1+ll,1,ii),highpcnt(valind(imtr)-1+ll,1,ii),maxv(valind(imtr)-1+ll,1,ii)
      enddo      
    close(9)
    open(unit=9,file=trim(stm) // '.prop.' // mstr // '.mf.' // vstr // '.IM.convinfo',status='unknown')
      print *,maxitr,imtr
      do ii=1,maxitr
        write(9,'(i7,7(1x,es11.4))') ii,meanv(valind(imtr)-1+ll,2,ii),stdv(valind(imtr)-1+ll,2,ii),minv(valind(imtr)-1+ll,2,ii),lowpcnt(valind(imtr)-1+ll,2,ii)&
	&,medv(valind(imtr)-1+ll,2,ii),highpcnt(valind(imtr)-1+ll,2,ii),maxv(valind(imtr)-1+ll,2,ii)
      enddo      
    close(9)
    enddo
   else
    ll=1
    write(mstr,'(i2.2)') imtr
    open(unit=9,file=trim(stm) // '.prop.' // mstr // '.RE.convinfo',status='unknown')
      print *,maxitr,imtr
      do ii=1,maxitr
        write(9,'(i7,7(1x,es11.4))') ii,meanv(valind(imtr)-1+ll,1,ii),stdv(valind(imtr)-1+ll,1,ii),minv(valind(imtr)-1+ll,1,ii),lowpcnt(valind(imtr)-1+ll,1,ii)&
	&,medv(valind(imtr)-1+ll,1,ii),highpcnt(valind(imtr)-1+ll,1,ii),maxv(valind(imtr)-1+ll,1,ii)
      enddo      
    close(9)
    open(unit=9,file=trim(stm) // '.prop.' // mstr // '.IM.convinfo',status='unknown')
      print *,maxitr,imtr
      do ii=1,maxitr
        write(9,'(i7,7(1x,es11.4))') ii,meanv(valind(imtr)-1+ll,2,ii),stdv(valind(imtr)-1+ll,2,ii),minv(valind(imtr)-1+ll,2,ii),lowpcnt(valind(imtr)-1+ll,2,ii)&
	&,medv(valind(imtr)-1+ll,2,ii),highpcnt(valind(imtr)-1+ll,2,ii),maxv(valind(imtr)-1+ll,2,ii)
      enddo      
    close(9)
   endif
  enddo
  
  deallocate(valind,meanv,stdv,medv,maxv,minv,lowpcnt,highpcnt)
  

end program MREpostprocessv8



      SUBROUTINE SORT(X,N,Y)
!C
!C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
!C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
!C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
!C              SINGLE PRECISION VECTOR Y.
!C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
!C                                OBSERVATIONS TO BE SORTED. 
!C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
!C                                IN THE VECTOR X. 
!C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
!C                                INTO WHICH THE SORTED DATA VALUES
!C                                FROM X WILL BE PLACED.
!C     OUTPUT--THE SINGLE PRECISION VECTOR Y
!C             CONTAINING THE SORTED
!C             (IN ASCENDING ORDER) VALUES
!C             OF THE SINGLE PRECISION VECTOR X.
!C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
!C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU 
!C                   (DEFINED AND USED INTERNALLY WITHIN
!C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
!C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
!C                   IF IL AND IU EACH HAVE DIMENSION K,
!C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
!C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
!C                   OF IL AND IU HAVE BEEN SET TO 36,
!C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
!C                   APPROXIMATELY 137 BILLION.
!C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
!C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
!C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
!C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
!C                   THEN THERE IS NO PRACTICAL RESTRICTION
!C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
!C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE 
!C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
!C                   INTO THIS SUBROUTINE.)
!C     OTHER DATAPA!C   SUBROUTINES NEEDED--NONE.
!C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
!C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
!C     LANGUAGE--ANSI FORTRAN. 
!C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
!C              WILL BE PLACED IN THE FIRST POSITION
!C              OF THE VECTOR Y,
!C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
!C              WILL BE PLACED IN THE SECOND POSITION
!C              OF THE VECTOR Y, ETC.
!C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
!C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
!C              THIS IS DONE BY HAVING THE SAME
!C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE. 
!C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
!C              CALL SORT(X,N,X)
!C              IS ALLOWABLE AND WILL RESULT IN
!C              THE DESIRED 'IN-PLACE' SORT.
!C     COMMENT--THE SORTING ALGORTHM USED HEREIN
!C              IS THE BINARY SORT.
!C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
!C              FOLLOWING TIME TRIALS INDICATE.
!C              THESE TIME TRIALS WERE CARRIED OUT ON THE
!C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
!C              IN AUGUST OF 1974.
!C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
!C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
!C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
!C              ALSO BEEN INCLUDED--
!C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
!C               NUMBERS SORTED
!C                N = 10                 .002 SE!C          .002 SEC
!C                N = 100                .011 SE!C          .045 SEC
!C                N = 1000               .141 SE!C         4.332 SEC
!C                N = 3000               .476 SE!C        37.683 SEC
!C                N = 10000             1.887 SE!C      NOT COMPUTED
!C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
!C                 BY RICHARD C. SINGLETON).
!C               --CACM JANUARY 1970, PAGE 54.
!C               --CACM OCTOBER 1970, PAGE 624.
!C               --JACM JANUARY 1961, PAGE 41.
!C     WRITTEN BY--JAMES J. FILLIBEN
!C                 STATISTICAL ENGINEERING LABORATORY (205.03)
!C                 NATIONAL BUREAU OF STANDARDS
!C                 WASHINGTON, D. C. 20234
!C                 PHONE--301-921-2315
!C     ORIGINAL VERSION--JUNE      1972. 
!C     UPDATED         --NOVEMBER  1975. 
!C
!C---------------------------------------------------------------------
!C
      DIMENSION X(1),Y(1)
      DIMENSION IU(36),IL(36) 
!C
      IPR=6
!C
!C     CHECK THE INPUT ARGUMENTS FOR ERRORS
!C
      IF(N.LT.1) GOTO 50
      IF(N.EQ.1) GOTO 55
      HOLD=X(1)
      DO 60 I=2,N
      IF(X(I).NE.HOLD) GOTO 90
   60 CONTINUE
      !print *,'***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT (A VECTOR) TO THE SORT SUBROUTINE HAS ALL ELEMENTS = ',HOLD
      DO 61 I=1,N
      Y(I)=X(I)
   61 CONTINUE
      RETURN
   50 print *,'***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE SORT SUBROUTINE IS NON-POSITIVE *****'
      print *, '***** THE VALUE OF THE ARGUMENT IS ',N
      RETURN
   55 print *, '***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****'
      Y(1)=X(1)
      RETURN
   90 CONTINUE
!C
!C-----START POINT-----------------------------------------------------
!C
!C     COPY THE VECTOR X INTO THE VECTOR Y
      DO 100 I=1,N
      Y(I)=X(I)
  100 CONTINUE
!C
!C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
!C
      NM1=N-1
      DO 200 I=1,NM1
      IP1=I+1
      IF(Y(I).LE.Y(IP1)) GOTO 200
      GOTO 250
  200 CONTINUE
      RETURN
  250 M=1 
      I=1 
      J=N 
  305 IF(I.GE.J) GOTO 370
  310 K=I 
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I).LE.AMED) GOTO 320 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
  320 L=J 
      IF(Y(J).GE.AMED) GOTO 340 
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I).LE.AMED) GOTO 340 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO 340
  330 Y(L)=Y(K)
      Y(K)=TT
  340 L=L-1
      IF(Y(L).GT.AMED) GOTO 340 
      TT=Y(L)
  350 K=K+1
      IF(Y(K).LT.AMED) GOTO 350 
      IF(K.LE.L) GOTO 330
      LMI=L-I
      JMK=J-K
      IF(LMI.LE.JMK) GOTO 360
      IL(M)=I
      IU(M)=L
      I=K 
      M=M+1
      GOTO 380
  360 IL(M)=K
      IU(M)=J
      J=L 
      M=M+1
      GOTO 380
  370 M=M-1
      IF(M.EQ.0)RETURN
      I=IL(M)
      J=IU(M)
  380 JMI=J-I
      IF(JMI.GE.11) GOTO 310
      IF(I.EQ.1) GOTO 305
      I=I-1
  390 I=I+1
      IF(I.EQ.J) GOTO 370
      AMED=Y(I+1)
      IF(Y(I).LE.AMED) GOTO 390 
      K=I 
  395 Y(K+1)=Y(K)
      K=K-1
      IF(AMED.LT.Y(K)) GOTO 395 
      Y(K+1)=AMED
      GOTO 390
      END 

