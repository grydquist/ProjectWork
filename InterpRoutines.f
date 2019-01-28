! This is the sb Type you'll need     
	  TYPE sb
      ! Size of searchbox
         REAL(KIND=8) :: step(nsd)
      ! Searchbox dimensions
         REAL(KIND=8) :: dim(nsd*2)
      ! Elements contained in searchbox
         INTEGER, ALLOCATABLE :: els(:)
      END TYPE sb

! These are just variables I used in the examples below.
	  INTEGER  split(nsd)
	  REAL(KIND=8) :: pint,uint(nsd),xpt(nsd)
      INTEGER      :: ID
	  ! Domain split into sb's
      TYPE(sb), ALLOCATABLE :: sbdom(:)
	  
	  ! split searchboxes this many times (can change to whatever)
      split = (/10,10,10/)
	  
	  ! Split domain into searchboxes once
	  CALL SBDomain(split,sbdom)
	  
	  ! Call these functions as much as needed after calling SBDomain
	  CALL SBx(xpt,sbdom,split,ID)
      CALL INTERP(xpt,sbdom(ID),pint,uint)
	  
!#################################################################### SBDOMAIN
      SUBROUTINE SBDomain(split,sbdom)
      INTEGER, INTENT(IN) :: split(nsd)
      TYPE(sb), INTENT(OUT), ALLOCATABLE :: sbdom(:)
      INTEGER :: ii,jj,cnt2,kk
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd), elbox(2*nsd,nEl)
      INTEGER, ALLOCATABLE :: sbel(:)

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sbdom(split(1)*split(2)*split(3)))

      ALLOCATE(sbel(nEl))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(split(3)*split(2)),seq2(split(3)*split(1)),seq3(split(2)*split(1)))

      ! Domain ranges
      diff(1) = MAXVAL(x(1,:)) - MINVAL(x(1,:))
      diff(2) = MAXVAL(x(2,:)) - MINVAL(x(2,:))
      diff(3) = MAXVAL(x(3,:)) - MINVAL(x(3,:))
      ! Size of sb
      do ii = 1,1!(split(1)*split(2)*split(3))
         sbdom(ii)%step = diff/split
      end do

      seq1=(/(ii, ii=0, split(2)*split(3)-1, 1)/)*split(1)+1
      cnt2=0
      do ii=1,split(1)*split(3)
            seq2(ii) = ii + cnt2*(split(2) - 1)*split(1)
            if (MOD(ii,split(1)).eq.0) cnt2 = cnt2 + 1
      end do
      seq3=(/(ii, ii=0, split(1)*split(2) - 1, 1)/) + 1

      ! Allocating sbdim with min and max dimensions
      do ii=1,split(1)
         sbdom(seq1 + ii - 1)%dim(1) = MINVAL(x(1,:)) + sbdom(1)%step(1)*(ii)
      end do

      do ii=1,split(2)
         sbdom(seq2 + ( ii - 1)*split(1))%dim(3) = MINVAL(x(2,:)) + sbdom(1)%step(2)*(ii)
      end do

      do ii=1,split(3)
         sbdom(seq3 + (ii - 1)*split(1)*split(2))%dim(5) = &
     &    MINVAL(x(3,:)) + sbdom(1)%step(3)*(ii)
      end do

      sbdom%dim(2) = sbdom%dim(1) + sbdom(1)%step(1)
      sbdom%dim(4) = sbdom%dim(3) + sbdom(1)%step(2)
      sbdom%dim(6) = sbdom%dim(5) + sbdom(1)%step(3)

      ! Making boxes surrounding elements
      do ii=1,Nel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(x(jj,IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(x(jj,IEN(:,ii)))
         end do
      end do

      do ii=1,split(1)*split(2)*split(3)
         cnt2=1
         sbel=0
         outer: do jj=1,Nel
            ! Check if elements are completely outside searchbox
            inner: do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).gt.sbdom(ii)%dim(2*kk  )) cycle outer
               if (elbox(2*kk  ,jj).lt.sbdom(ii)%dim(2*kk-1)) cycle outer
            enddo inner

            sbel(cnt2) = jj
            cnt2=cnt2+1
         enddo outer
         ALLOCATE(sbdom(ii)%els(cnt2-1))
         sbdom(ii)%els=sbel(1:cnt2-1)
      end do
      
      END SUBROUTINE SBDomain
	  
!#################################################################### SBx
      ! Find all Searchboxes x is in
      SUBROUTINE SBx(x,sbdom,split,sbid)
      IMPLICIT NONE
      TYPE(sb), INTENT(IN), ALLOCATABLE :: sbdom(:)
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER, INTENT(IN) :: split(nsd)
      INTEGER, INTENT(OUT) :: sbid
      REAL(KIND=8) :: xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Set domain back to zero
      xzero(1) = x(1) - minval(sbdom%dim(1))
      xzero(2) = x(2) - minval(sbdom%dim(3))
      xzero(3) = x(3) - minval(sbdom%dim(5))

      ! Number of searchbox steps in x,y,and z
      xsteps=FLOOR(xzero/(sbdom(1)%step))
      ! ID of sb x is in
      sbid = xsteps(1) + split(1)*xsteps(2) + split(1)*split(2)*xsteps(3) + 1

      RETURN
      END SUBROUTINE SBx
	  
!#################################################################### INTERP
      ! Interpolates pressure/velocity
      SUBROUTINE INTERP(xpt,sbdom,pint,uint)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: xpt(nsd)
      REAL(KIND=8), INTENT(OUT):: pint, uint(nsd)
      TYPE(sb), INTENT(IN)     :: sbdom
      INTEGER :: ii,jj,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,eNoN), &
       &               shps(eNoN), prntx(nsd)
    
      xXi = 0D0
      pint = 0D0
      uint = 0D0

      do ii=1,size(sbdom%els)

      IF (nsd .EQ. 2) THEN
      !
      ! 2D not done
      !
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         
         DO a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         END DO


      ! 3D case
      ELSE

         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,IEN(a,sbdom%els(ii)))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,IEN(a,sbdom%els(ii)))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,IEN(a,sbdom%els(ii)))*Nxi(3,a)
         END DO
         
      ! Inverting matrix
         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)&
     &       + xXi(1,2)*xXi(2,3)*xXi(3,1)&
     &       + xXi(1,3)*xXi(2,1)*xXi(3,2)&
     &       - xXi(1,1)*xXi(2,3)*xXi(3,2)&
     &       - xXi(1,2)*xXi(2,1)*xXi(3,3)&
     &       - xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac
         
         ! Finding coordinates in parent domain
         DO a=1, nsd
       prntx(a) =      xiX(a,1)*( xpt(1) - x(1,IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,2)*( xpt(2) - x(2,IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,3)*( xpt(3) - x(3,IEN(4,sbdom%els(ii))))
         END DO
      END IF

      ! Finding shape functions
       shps(1) =  prntx(1)
       shps(2) =  prntx(2)
       shps(3) =  prntx(3)
       shps(4) =  1 - prntx(1) - prntx(2) - prntx(3)
      ! If shape functions positive, we've found the correct element
      IF (ALL(shps.gt.0D0)) then
         EXIT
      END IF

   end do
   
   ! Using shape functions to interpolate value at coordinate
   do jj=1,eNoN
      pint = pint + pres(IEN(jj,sbdom%els(ii)))*shps(jj)
      uint = uint + vel(:,IEN(jj,sbdom%els(ii)))*shps(jj)
      print *, vel(1,IEN(jj,sbdom%els(ii)))
   end do
         
   END SUBROUTINE INTERP

	  
	  