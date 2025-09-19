c23456789012345678901234567890123456789012345678901234567890123456789012

      program anmcentroid
   
      integer npar1,npar2,npar3
      PARAMETER (npar1=15000)
      PARAMETER (npar2=10000)
      PARAMETER (npar3=20)

    
      character*3 restyp(npar1),restypp(npar1),resnam
      character*6 dummy
      character*3 atname(npar1),resatom(npar1,100),at
      character*1 chain,chainold,reschain(npar1),ch
      character*1 let,letold
      character*4 dumm
      character*2 aa(npar1)
      character(len=1000) :: buf
       
      integer resnum,res,rno,nresidue,ires
      integer resn(npar1),natoms(npar1),resnumo(npar1)
      integer ee,dd,atnum

      real xxnew(npar1),yynew(npar1),zznew(npar1)
      real vxa(npar1,100),vya(npar1,100),vza(npar1,100)
      real vx(npar1,100),vy(npar1,100),vz(npar1,100)      
      real xx(npar1,100),yy(npar1,100),zz(npar1,100)
      real x,y,z,bb
      real mass(npar1),fc(npar1),eig(npar2)
      
      real e(npar2),ecum(npar2),f(npar2),r2
      real bfactorr(npar1),rcut,r,b,bfac(npar1)
       real ddx,ddy,ddz,bx,by,bz
       real sum(100),contact
       real cx(npar1),cy(npar1),cz(npar1)
        real cxx(npar1),cyx(npar1),czx(npar1)
        
	open(16,file='7aku-apo.pdb')
	
        open(60,file='coordmin1.pdb')

	open(18,file='connectivity.out')



C***********************************************************************
c   rcutt= the atomic cutoff distance, rno=number of residues
c   nresidue=number of atoms

	rcutt=4.5
         rno=610
         nresidue=4702
       	jres=3*rno
		
C***********************************************************************

 50     FORMAT(I6,1x,20(f8.5,1x))
 51     FORMAT(I6,1x,10(f10.8,1x))
 52   Format (A6)
 53   Format(A6,I5,2X,a2,2x,A3,1X,A1,I4,4X,3F8.3,6X,f6.2)
 54    Format(I5,2x,3F8.3)

 57     Format(A6,I5,2X,a2,2x,A3,1X,A1,I4,4X,3F8.3,2x,f4.2,1x,f5.1)
 58     Format(A6,I5,2X,a2,2x,A3,2X,I4,4X,3F8.3,1x,f4.2,2x,f6.2)
 59     Format(A6,I5,2X,a2,2x,A3,1X,A1,I4,4X,3F8.3,1x,f4.2,2x,f6.2)
 77   format (A6,I5,2X,a2,2x,A3,1X,A1,I4,A1,3X,3F8.3,6X,f6.2)
C***********************************************************************

c---------------------------------------------------------------

c           read all the atomic coordinates and sort


	i=0.
	ires=0.
	

	do  k=1,nresidue

 31      read(16,52)dummy
      if (dummy.eq.'ATOM  ') then
      backspace (16)

      
      read(16,77)dummy,atnum,at,resnam,ch,
     : resnum,let,x,y,z,bb

       else
       goto 31
       endif
       
	i=i+1
	
	
        if(atnum.eq.1) then
         resold=resnum
         chainold=ch
         letold=let
         ires=1
         ichain=1
      endif

      if(resnum.ne.resold.or.letold.ne.let)then
         natoms(ires)=i-1
         ires=ires+1
         i=1
         resold=resnum
         letold=let
      endif

      if(ch.ne.chainold)then
         ichain=ichain+1
         chainold=ch
      endif


      xx(ires,i)=x
      yy(ires,i)=y
      zz(ires,i)=z



        restypp(ires)=resnam
        resnumo(ires)=resnum
	resatom(ires,i)=at
        bfac(ires)=bb
        reschain(ires)=ch


	enddo

c 33   maxat=atnum

       natoms(ires)=i
       write(6,*)ires,rno

	if(ires.ne.rno) write(6,*)'check the residue number!'


	do i=1,rno
        cx(i)=0.
        cy(i)=0.
        cz(i)=0.
        cxx(i)=0.
        cyx(i)=0.
        czx(i)=0.
        
	do j=1,rno

        enddo
	enddo

         do i=1,rno
	    do k=1,natoms(i)
	    cx(i)=cx(i)+xx(i,k)
	    cy(i)=cy(i)+yy(i,k)
	    cz(i)=cz(i)+zz(i,k)


	    enddo

            cxx(i)=cx(i)/natoms(i)
            cyx(i)=cy(i)/natoms(i)
            czx(i)=cz(i)/natoms(i)

       write (60,53)'ATOM  ',i,'CA',restypp(i),reschain(i),resnumo(i),
     :  cxx(i),cyx(i),czx(i),bfac(i)

            enddo


     
c       calculate the local interaction strength between residues


	do i=1,rno
	do j=i+1,rno
	     contact=0.
	     
	do k=1,natoms(i)
	do l=1,natoms(j)


	dx=xx(i,k)-xx(j,l)
	dy=yy(i,k)-yy(j,l)
	dz=zz(i,k)-zz(j,l)

	r2=sqrt(dx*dx+dy*dy+dz*dz)
	
	
	if(i.ne.j.and.r2.le.rcutt) then
	contact=contact+1

	endif

	enddo
	enddo
           if(contact.ne.0) then
           write(18,*)i,j,contact/(natoms(i)**0.5)/(natoms(j)**0.5)
           endif
           
	enddo
	enddo
	




	
	stop
	end
	
