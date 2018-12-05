      PROGRAM READVTK
      IMPLICIT NONE
!     This is for writing binary files
      CHARACTER, PARAMETER :: eol = CHAR(10)
!     This is the standard length for all the strings used in this code
      INTEGER, PARAMETER :: stdL = 256
!     Assuming element_number_of_node = 4, degrees_of_freedom = 4, and
!     number_of_spatial_dimensions = 3
      INTEGER, PARAMETER :: eNoN=4, dof=4, nsd=3, nG = 4, Np=2
      LOGICAL isBinary, foundVar(2)
      INTEGER i, a, e, fid, m, iPos, nNo, nEl, tmpI(eNoN+1), g
      REAL(KIND=8) N(eNoN,nG), xi(nsd,eNoN), Nxi(nsd,eNoN), s, w(nG), t
      CHARACTER c
      CHARACTER(LEN=stdL) rLine, tmp, fName
      INTEGER, ALLOCATABLE :: IEN(:,:)
      REAL, ALLOCATABLE :: tmpS(:), tmpV(:,:)
!Grant Test
      REAL(KIND=8) :: xg(nsd,eNoN), Nx(nsd,eNoN), Jac, ks(nsd,nsd),time
      INTEGER cnt, split(nsd),j, k

      !prtcollide tests
      !REAL(KIND=8) :: x1(nsd), x2(nsd), v1(nsd), v2(nsd), test1(nsd),tester(nsd), test2(nsd)

      TYPE prt
      ! Properties
         REAL(KIND=8) :: mp, dp = 1D0, rhop = 0.1D0
      ! Flow characteristics
         REAL(KIND=8) :: x(nsd), vel(nsd), prntx(nsd), shps(eNoN), remdtp
      ! Searechbox/ element location
         INTEGER :: sbid(2**nsd), elid, near(Np)
      ! Collisions with other particles
         LOGICAL :: collided=.false.
      ! Restitution Coefficient
         REAL(KIND=8) :: k
      END TYPE prt

      ! Collection of particles
      type(prt) :: prts(Np)

      TYPE sb
      ! Searchbox dimensions
         REAL(KIND=8) :: dim(nsd*2)
      ! Elements contained in searchbox
         INTEGER, ALLOCATABLE :: els(:)
      END TYPE sb

      ! Domain spliited into sb's
      type(sb), ALLOCATABLE :: sbdom(:)

      REAL, PARAMETER :: pi=3.1415926535897932384626433

      REAL(KIND=8), ALLOCATABLE :: vel(:,:), pres(:), X(:,:)

      prts%mp = pi*prts%rhop/6D0*prts%dp**3D0
!End Grant Test

      w = 1D0/24D0
      s = (5D0 + 3D0*SQRT(5D0))/2D1
      t = (5D0 -     SQRT(5D0))/2D1
      xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
      xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
      xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
      xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
      DO g=1, nG
         N(1,g) = xi(1,g)
         N(2,g) = xi(2,g)
         N(3,g) = xi(3,g)
         N(4,g) = 1D0 - xi(1,g) - xi(2,g) - xi(3,g)
      END DO
      Nxi(1,1) =  1D0
      Nxi(2,1) =  0D0
      Nxi(3,1) =  0D0
      Nxi(1,2) =  0D0
      Nxi(2,2) =  1D0
      Nxi(3,2) =  0D0
      Nxi(1,3) =  0D0
      Nxi(2,3) =  0D0
      Nxi(3,3) =  1D0
      Nxi(1,4) = -1D0
      Nxi(2,4) = -1D0
      Nxi(3,4) = -1D0

!     Reading the name of the vtk file
      fName="pipe_example_100.vtk"
      !CALL GETARG(1,fName)
      !fName = ADJUSTL(fName)
      fid = 1
      print * , fName
      OPEN(fid, FILE=fName, STATUS='OLD')
      
      isBinary = .FALSE.
      DO i=1, 20
         READ (fid,'(a)') rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:6) .EQ. 'BINARY') THEN
            isBinary = .TRUE.
            CLOSE(fid)
            OPEN (fid, FILE=fName, STATUS='UNKNOWN', ACCESS='STREAM',&
     &         FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
            iPos = 1
            EXIT
         END IF
      END DO
      REWIND(fid)
!####################################################################
!     Position vector      
!####################################################################
!     Skipping stuff until I get to a point which contains position
!     vector
      DO
         IF (isBinary) THEN
!     Since reading from STREAM, we need to find the end of a line
!     manually
            rLine = ''
            DO i=1, stdL
               READ (fid,END=001,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=001) rLine
         END IF
         rLine = ADJUSTL(rLine)
         IF (rLine(1:6) .EQ. 'POINTS') EXIT
      END DO
      rLine = rLine(7:)
      READ (rLine,*) nNo
      PRINT *, "Number of nodes:", nNo
!     Allocating space for position vector, vel, and pres
      ALLOCATE(x(nsd,nNo), vel(nsd,nNo), pres(nNo), tmpS(nNo), &
     &   tmpV(nsd,nNo))
      IF (isBinary) THEN
         DO a=1, nNo
            READ(fid,END=001) tmpV(:,a)
            iPos = iPos + nsd*KIND(tmpV)
         END DO
         x = tmpV
      ELSE
         DO a=1, nNo
            READ(fid,*,END=001) x(:,a)
         END DO
      END IF
!####################################################################
!     Connectivity
!####################################################################
!     Skipping stuff until I get to a point which contains connectivity
      DO
         IF (isBinary) THEN
            rLine = ''
            DO i=1, stdL
               READ (fid,END=001,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=001) rLine
         END IF
         rLine = ADJUSTL(rLine)
         IF (rLine(1:5) .EQ. 'CELLS') EXIT
      END DO
      rLine = rLine(6:)
      READ (rLine,*) nEl
      PRINT *, "Number of elements:", nEl
!     Allocating space for IEN array
      ALLOCATE(IEN(eNoN,nEl))
      IF (isBinary) THEN
         DO e=1, nEl
            READ(fid,END=001) tmpI
            iPos = iPos + (eNoN+1)*KIND(IEN)
            IEN(:,e) = tmpI(2:eNoN+1)
            IF (tmpI(1).NE.4) STOP "This code is designed for TET only"
         END DO
      ELSE
         DO e=1, nEl
            READ(fid,*,END=001) i, IEN(:,e)
         END DO
      END IF
!####################################################################
!     Velocity and pressure  
!####################################################################
!     Skipping stuff until I get to a point which contains data
      DO
         IF (isBinary) THEN
            rLine = ''
            DO i=1, stdL
               READ (fid,END=001,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=001) rLine
         END IF
         rLine = ADJUSTL(rLine)
         IF (rLine(1:10) .EQ. 'POINT_DATA') EXIT
      END DO
!     Neither of variables are found by default
      foundVar = .FALSE.
      DO
         IF (isBinary) THEN
            rLine = ''
            DO i=1, stdL
               READ (fid,END=003,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=003) rLine
         END IF
         rLine = ADJUSTL(rLine)
!     Continuing until to get to a block of nodal data
         IF (rLine(1:7) .EQ. 'SCALARS') THEN
            m = 1
!     To skip lookup table statment            
            IF (isBinary) THEN
               DO i=1, stdL
                  READ (fid,END=003,POS=iPos) c
                  iPos = iPos + 1
                  IF (c .EQ. eol) EXIT
               END DO
            ELSE
               READ (fid,'(a)',END=001) tmp
            END IF
         ELSE IF (rLine(1:7) .EQ. 'VECTORS') THEN
            m = nsd
         ELSE 
            CYCLE
         END IF
!     Skipping VECTOR/SCALAR
         rLine = rLine(9:)
         i = LEN(TRIM(rLine))
!     Skipping 'float'
         rLine = rLine(1:i-6)
!     Expecting to find Navier-Stokes
         IF (rLine(1:2) .NE. 'NS') THEN
            PRINT *, "Skipping equation <"//rLine(1:2)//"> in " &
     &         //TRIM(fName)
            CYCLE
         END IF
         rLine = rLine(4:)
!     Finding output based on the name of variable
         IF (rLine .EQ. 'Velocity') THEN
            PRINT *, "Reading velocity"
            IF (m .NE. 3) STOP "m.NE.3 for velocity"
            IF (isBinary) THEN
               DO a=1, nNo
                  READ(fid,END=001) tmpV(:,a)
                  iPos = iPos + m*KIND(tmpV)
               END DO
               vel = tmpV
            ELSE
               DO a=1, nNo
                  READ(fid,*,END=001) vel(:,a)
               END DO
            END IF
            foundVar(1) = .TRUE.
         ELSE IF (rLine .EQ. 'Pressure') THEN
            PRINT *, "Reading pressure"
            IF (m .NE. 1) STOP "m.NE.1 for pressure"
            IF (isBinary) THEN
               DO a=1, nNo
                  READ(fid,END=001) tmpS(a)
                  iPos = iPos + m*KIND(tmpS)
               END DO
               pres = tmpS
            ELSE
               DO a=1, nNo
                  READ(fid,*,END=001) pres(a)
               END DO
            END IF
            foundVar(2) = .TRUE.
         END IF
      END DO
 003  IF (ANY(.NOT.foundVar)) STOP "Did not find some varibales" 
      CLOSE(fid)
      PRINT *, TRIM(fName)//" was read successfully"
!####################################################################
!     Your own implementation: nNo, nEl, x, vel, pres, IEN are given
!####################################################################
      pres = 1D0
      !ALLOCATE(pGrad(nsd,nNo))
      IEN=IEN+1
      cnt=1

      !Test particle velocity
      prts(1)%vel(1)=0
      prts(1)%vel(2)=0
      prts(1)%vel(3)=1

      prts(2)%vel(1)=0
      prts(2)%vel(2)=0
      prts(2)%vel(3)=-1D0

      !Test particle position

      prts(1)%x(1) = 0D0
      prts(1)%x(2) = 0D0
      prts(1)%x(3) = 10D0

      prts(2)%x(1) = 0D0
      prts(2)%x(2) = 0.7D0
      prts(2)%x(3) = 11.5D0

      ! split searchboxes this many times
      split = (/10,10,10/)
      
      ! Resititution Coefficient
      prts%k = 1D0

      vel = 0D0
      vel(3,:) = 0.1D0
      
      CALL SBDomain(x,split,sbdom)
      CALL xSB(prts(1),sbdom,split)
      CALL xEl(sbdom(prts(1)%sbid(1)),prts(1),x)

      time=0d0
      open(88,file='pos.txt')

      do a=1,10000
         do i=1,Np
            CALL xSB(prts(i),sbdom,split)
            CALL xEl(sbdom(prts(i)%sbid(1)),prts(i),x)
            if (i.eq.2) then
               k=-1
            else
               k=1
            end if
            CALL prtAdvance(prts(i),k*vel,x)   
         end do

         do i=1,Np
               ! Collisions
            do j=1,Np
               ! Check if the particle collides with any other particles and hasn't collided. Advance if so.
               if ((i.ne.j).and.(.not.(prts(i)%collided))) CALL prtCollide(prts(i),prts(j),0.01D0)
            end do

            ! If particles haven't collided, advance by vel*dtp
            if (.not.(prts(i)%collided)) then
               prts(i)%x = 0.01D0*prts(i)%vel +prts(i)%x
            else
               prts(i)%collided = .false.
            end if

            print *, prts(i)%x,time
         end do

         !!! Need searchbox implementsation (shouldn't be too hard)

         !!! Need to fix timestep so I can bring it into prtCollide. Should just do it in the overall program

         !!! Idea: multiple collisions: get minimum tcr, do that one first
         !!! Only advance to this collision, then make new dt = dt - tcr
         !!! When checking for particle collisions, advance other particles by tcr
         !!! And Boom! Should work. Will need some testing
         !!! Honestly, probably don't need mult in one time step though


         write(88,*) prts(1)%x,prts(2)%x
      end do

      close(88)

      do a=1,nEl
      xg(:,1) = x(:,IEN(1,a))
      xg(:,2) = x(:,IEN(2,a))
      xg(:,3) = x(:,IEN(3,a))
      xg(:,4) = x(:,IEN(4,a))
      CALL GNN(xg,Nx,Jac,ks,prts(1)%prntx,prts(1)%x,prts(1)%shps)
      if ((ALL(prts(1)%shps.gt.0))) EXIT
      cnt=cnt+1
      end do
      print *, cnt,prts(1)%shps
      
      STOP 
 001  STOP "A block of data is missing"

      CONTAINS
!#################################################################### GRAD      
      
      !Be careful here with IEN here
      SUBROUTINE grad(s, v)
      REAL(KIND=8), INTENT(IN) :: s(nNo)
      REAL(KIND=8), INTENT(OUT) :: v(nsd,nNo)

      INTEGER e, g, a, Ac
      REAL(KIND=8) Nx(nsd,eNoN), xl(nsd,eNoN), Jac, &
     &   vl(nsd)
      REAL(KIND=8), ALLOCATABLE :: sA(:), sF(:,:)

      ALLOCATE(sA(nNo), sF(nsd,nNo))

      sA = 0D0
      sF = 0D0
      DO e=1, nEl
         DO a=1, eNoN
            Ac = IEN(a,e)
            xl(:,a) = x(:,Ac)
         END DO
         !CALL GNN(xl, Nx, Jac, ksix)
         DO g=1, nG
            vl = 0D0
            DO a=1, eNoN
               Ac = IEN(a,e)
               vl = vl + Nx(:,a)*s(Ac)
            END DO

      !Mapping Tau into the nodes by assembling it into a local vector
            DO a=1, eNoN
               Ac       = IEN(a,e)
               sA(Ac)   = sA(Ac)   + w(g)*Jac*N(a,g)
               sF(:,Ac) = sF(:,Ac) + w(g)*Jac*N(a,g)*vl
            END DO
         END DO
      END DO
      stop
      
      DO Ac=1, nNo
         v(:,Ac) = sF(:,Ac)/sA(Ac)
      END DO

      DEALLOCATE(sA, sF)

      RETURN
      END SUBROUTINE grad

!#################################################################### GNN
      PURE SUBROUTINE GNN(x, Nx, Jac, ks,prntx,prtx,shps)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: x(nsd,eNoN), prtx(nsd)
      REAL(KIND=8), INTENT(OUT) :: Nx(nsd,eNoN), Jac, ks(nsd,nsd) 
      REAL(KIND=8), INTENT(OUT) :: prntx(nsd),shps(eNoN)
      INTEGER :: a
      REAL(KIND=8) xXi(nsd,nsd), xiX(nsd,nsd)

      prntx = 0D0
      Nx  = 0D0
      xXi = 0D0
      IF (nsd .EQ. 2) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
         ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
         ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
         ks(2,1) = ks(1,2)
         
         DO a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         END DO
      ELSE
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         END DO
         
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

         ks(1,1) = xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1)+xiX(3,1)*xiX(3,1)
         ks(1,2) = xiX(1,2)*xiX(1,1)+xiX(2,2)*xiX(2,1)+xiX(3,2)*xiX(3,1)
         ks(1,3) = xiX(1,3)*xiX(1,1)+xiX(2,3)*xiX(2,1)+xiX(3,3)*xiX(3,1)
         ks(2,2) = xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2)+xiX(3,2)*xiX(3,2)
         ks(2,3) = xiX(1,2)*xiX(1,3)+xiX(2,2)*xiX(2,3)+xiX(3,2)*xiX(3,3)
         ks(3,3) = xiX(1,3)*xiX(1,3)+xiX(2,3)*xiX(2,3)+xiX(3,3)*xiX(3,3)
         ks(2,1) = ks(1,2)
         ks(3,1) = ks(1,3)
         ks(3,2) = ks(2,3)
         
         DO a=1, eNoN
            Nx(1,a) = Nx(1,a) + Nxi(1,a)*xiX(1,1)  &
     &                        + Nxi(2,a)*xiX(2,1)  &
     &                        + Nxi(3,a)*xiX(3,1)
            
            Nx(2,a) = Nx(2,a) + Nxi(1,a)*xiX(1,2) &
     &                        + Nxi(2,a)*xiX(2,2) &
     &                        + Nxi(3,a)*xiX(3,2)
            
            Nx(3,a) = Nx(3,a) + Nxi(1,a)*xiX(1,3) &
     &                        + Nxi(2,a)*xiX(2,3) &
     &                        + Nxi(3,a)*xiX(3,3)

         END DO
         DO a=1, nsd
            prntx(a) = xiX(a,1)*(prtx(1) - x(1,4)) + &
     &                 xiX(a,2)*(prtx(2) - x(2,4)) + &
     &                 xiX(a,3)*(prtx(3) - x(3,4))
         END DO

         shps(1)=prntx(1)
         shps(2)=prntx(2)
         shps(3)=prntx(3)
         shps(4)=1-prntx(1)-prntx(2)-prntx(3)

      END IF

      RETURN
      END SUBROUTINE GNN

!#################################################################### SBDOMAIN
      SUBROUTINE SBDomain(x,split,sbdom)
      REAL(KIND=8), INTENT(IN) :: x(nsd,nNo)
      INTEGER, INTENT(IN) :: split(nsd)
      TYPE(sb), INTENT(OUT), ALLOCATABLE :: sbdom(:)
      INTEGER :: ii,jj,cnt2,kk
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd),step(nsd), elbox(2*nsd,nEl)
      INTEGER, ALLOCATABLE :: sbel(:)

      ! dim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sbdom(split(1)*split(2)*split(3)))

      ALLOCATE(sbel(nEl))

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(split(3)*split(2)),seq2(split(3)*split(1)),seq3(split(2)*split(1)))

      ! Domain ranges
      diff(1)=MAXVAL(x(1,:))-MINVAL(x(1,:))
      diff(2)=MAXVAL(x(2,:))-MINVAL(x(2,:))
      diff(3)=MAXVAL(x(3,:))-MINVAL(x(3,:))
      ! Size of sb
      step=diff/((split+1)/2)

      seq1=(/(ii, ii=0, split(2)*split(3)-1, 1)/)*split(1)+1
      cnt2=0
      do ii=1,split(1)*split(3)
            seq2(ii)=ii+cnt2*(split(2)-1)*split(1)
            if (MOD(ii,split(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, split(1)*split(2)-1, 1)/)+1

      ! Allocating sbdim, such that they overlap by 50%
      do ii=1,split(1)
         sbdom(seq1+ii-1)%dim(1)=MINVAL(x(1,:))+step(1)*(ii-1)/2
      end do

      do ii=1,split(2)
         sbdom(seq2+(ii-1)*split(1))%dim(3)=MINVAL(x(2,:))+step(2)*(ii-1)/2
      end do

      do ii=1,split(3)
         sbdom(seq3+(ii-1)*split(1)*split(2))%dim(5)=MINVAL(x(3,:))+step(3)*(ii-1)/2
      end do

      sbdom%dim(2)=sbdom%dim(1)+step(1)
      sbdom%dim(4)=sbdom%dim(3)+step(2)
      sbdom%dim(6)=sbdom%dim(5)+step(3)

      ! Making boxes surrounding elements
      do ii=1,Nel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(x(jj,IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(x(jj,IEN(:,ii)))
         end do
      end do

      !! I think I might want to change the order of one of the "if (elbox(" statements below

      do ii=1,split(1)*split(2)*split(3)
         cnt2=1
         sbel=0
         do jj=1,Nel
            ! Check if elements are completely outside searchbox
            do kk=1,nsd
               ! Cycle if min value elbox .gt. max value searchbox & vice-verse
               if (elbox(2*kk-1,jj).lt.sbdom(ii)%dim(2*kk  )) cycle
               if (elbox(2*kk  ,jj).gt.sbdom(ii)%dim(2*kk-1)) cycle
            end do

            sbel(cnt2) = jj
            cnt2=cnt2+1
            !end if
         end do
         sbdom(ii)%els=sbel
      end do
      
      END SUBROUTINE SBDomain

!#################################################################### XSB
      ! Find all Searchboxes x is in
      SUBROUTINE xSB(myprt,sbdom,split)
      IMPLICIT NONE
      TYPE(sb), INTENT(IN), ALLOCATABLE :: sbdom(:)
      TYPE(prt), INTENT(INOUT) :: myprt
      INTEGER, INTENT(IN) :: split(nsd)
      REAL(KIND=8) :: step(nsd),xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Searchbox dimensions

      !! add in steps to sb type
      step(1) = sbdom(1)%dim(2) - sbdom(1)%dim(1)
      step(2) = sbdom(1)%dim(4) - sbdom(1)%dim(3)
      step(3) = sbdom(1)%dim(6) - sbdom(1)%dim(5)

      ! Set domain back to zero
      xzero(1) = myprt%x(1) - minval(sbdom%dim(1))
      xzero(2) = myprt%x(2) - minval(sbdom%dim(3))
      xzero(3) = myprt%x(3) - minval(sbdom%dim(5))

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps=FLOOR(xzero/step)
      ! furthest searchbox in front
      myprt%sbid(1) = xsteps(1)+split(1)*xsteps(2)+split(1)*split(2)*xsteps(3)+1
      ! previous sb in x
      myprt%sbid(2) = myprt%sbid(1)-1
      ! previous sb's in y
      myprt%sbid(3) = myprt%sbid(1)-split(1)
      myprt%sbid(4) = myprt%sbid(3)-1
      ! Next sb's in z (if available)
      if (nsd.eq.3) then
         myprt%sbid(5) = myprt%sbid(1) - split(1)*split(2)
         myprt%sbid(6) = myprt%sbid(5) - 1
         myprt%sbid(7) = myprt%sbid(5) - split(1)
         myprt%sbid(8) = myprt%sbid(7) - 1
      end if
      

      END SUBROUTINE xSB

!#################################################################### XEL
      ! Finds element particle of position x is in is in
      SUBROUTINE xEl(sbdom,myprt,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: x(nsd,nNo)
      TYPE(sb), INTENT(IN) :: sbdom
      TYPE(prt), INTENT(INOUT) :: myprt
      INTEGER :: ii,cnt,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,eNoN)
      cnt=1
      myprt%elid=0

      do ii=1,nEl

      xXi = 0D0

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


      ELSE

         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,IEN(a,sbdom%els(ii)))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,IEN(a,sbdom%els(ii)))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,IEN(a,sbdom%els(ii)))*Nxi(3,a)
         END DO
         
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
         
         DO a=1, nsd
      myprt%prntx(a) = xiX(a,1)*(myprt%x(1) - x(1,IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,2)*(myprt%x(2) - x(2,IEN(4,sbdom%els(ii)))) + &
     &                 xiX(a,3)*(myprt%x(3) - x(3,IEN(4,sbdom%els(ii))))
         END DO
      END IF

      myprt%shps(1) = myprt%prntx(1)
      myprt%shps(2) = myprt%prntx(2)
      myprt%shps(3) = myprt%prntx(3)
      myprt%shps(4) = 1 - myprt%prntx(1) - myprt%prntx(2) - myprt%prntx(3)

      IF (ALL(myprt%shps.gt.0D0)) then
         myprt%elid=sbdom%els(ii)
         EXIT
      END IF
         
      end do

      if (myprt%elid.eq.0) print *, 'outside domain'

      END SUBROUTINE xEl

!#################################################################### PRTDRAG
      ! Find acceleration on one particle from drag
      SUBROUTINE prtDrag(myprt,vel,apd,taupo)
      TYPE(prt), INTENT(IN) :: myprt
      REAL(KIND=8), INTENT(IN)  :: vel(nsd,nNo)
      REAL(KIND=8), INTENT(OUT) :: apd(nsd)
      REAL(KIND=8) :: fvel(nsd),taup
      INTEGER ii,jj
      REAL(KIND=8), INTENT(OUT), OPTIONAL :: taupo

      !fluid Parameters
      REAL(KIND=8) :: rho, mu
      ! Derived from flow
      REAL(KIND=8) :: fSN, magud, Rep, relvel(nsd)


      ! Fluid density
      rho=1D0
      ! Fluid viscosity
      mu=0.01D0

      ! Particle relaxation time
      taup=myprt%rhop * myprt%dp**2D0/mu/18D0
      if (present(taupo)) taupo=taup

      ! Interpolate velocity at particle point
      fvel=0D0
      do ii=1,nsd
         do jj=1,eNoN
            fvel(ii) = fvel(ii) + vel(ii,IEN(jj,myprt%elid))*myprt%shps(jj)
         end do
      end do

      ! Relative velocity
      relvel = fvel-myprt%vel
      ! Relative velocity magnitude
      magud = SUM(relvel**2D0)**0.5D0
      ! Reynolds Number
      Rep = myprt%dp*magud*myprt%rhop/mu
      ! Schiller-Neumann (finite Re) correction
      fSN = 1D0 + 0.15D0*Rep**0.687D0
      ! Stokes corrected drag force
      apD = fSN/taup*relvel
      
      END SUBROUTINE prtDrag

!#################################################################### PRTADVANCE
      ! Advance 1 Particle through flow
      SUBROUTINE prtAdvance(myprt,vel,x)
      IMPLICIT NONE
      TYPE(prt), INTENT(INOUT) :: myprt
      TYPE(prt) :: tmpprt
      REAL(KIND=8), INTENT(IN) :: vel(nsd,nNo), x(nsd,nNo)
      INTEGER ii

      !Particle/fluid Parameters
      REAL(KIND=8) :: g(nsd), rho, dtp,maxdtp,sbdt(nsd)
      ! Derived from flow
      REAL(KIND=8) :: apT(nsd), apd(nsd), apdpred(nsd), apTpred(nsd), taup
      ! RK2 stuff
      REAL(KIND=8) :: prtxpred(nsd), pvelpred(nsd)

      tmpprt = myprt

      ! Gravity
      g=0D0
      !g(3)=1D0
      ! Fluid density
      rho=1D0

      !! Time step (still need to match to overall flow solver)
      ! Maxdt for overall solver
      maxdtp = 0.01D0

      ! Get drag acceleration/particle relaxation time
      CALL prtDrag(myprt,vel,apd,taup)

      ! Separate into sb sizes (could be optimized here)
      do ii=1,nsd
         sbdt(ii)=(sbdom(1)%dim(2*ii)-sbdom(1)%dim(2*ii-1))/abs(myprt%vel(ii))
      end do
      
      ! dtp is minimum between time to travel half searchbox, 1/10 relaxation time, and maxdtp
      dtp = min(maxdtp,0.5*minval(sbdt),taup/10)

      ! Total acceleration (just drag and buoyancy now)

      apT = apd + g*(1D0 - rho/myprt%rhop)

      ! 2nd order advance (Heun's Method)
      ! Predictor
      pvelpred = myprt%vel + dtp*apT
      prtxpred = myprt%x + dtp*myprt%vel

      tmpprt%vel = pvelpred
      tmpprt%x   = prtxpred

      CALL xSB(tmpprt,sbdom,split)
      CALL xEl(sbdom(tmpprt%sbid(1)),tmpprt,x)
      CALL prtDrag(tmpprt,vel,apdpred)

      apTpred = apdpred + g*(1D0 - rho/myprt%rhop)

      ! Corrector
      myprt%vel = myprt%vel + 0.5D0*dtp*(apT+apTpred)

      ! Collisions work. Just need to figure out how to change velocity of both particles
      !! Maybe get only velocities for particles, bring those out, then advance all particles through flow in
      !!    collision solver?
      
      !prtx = prtx + 0.5D0*dtp*(pvel+pvelpred)
      time=time+dtp

      END SUBROUTINE prtAdvance

!#################################################################### PRTCOLLIDE
      ! Detects and enacts collisions
      !! Only between particles right now
      SUBROUTINE prtCollide(prt1,prt2,dtp)
      IMPLICIT NONE
      TYPE(prt), INTENT(INOUT) :: prt1, prt2
      REAL(KIND=8), INTENT(IN) :: dtp

      ! Calculating distance coefficient
      REAL(KIND=8) :: a, b, c, d, e, f, qa, qb, qc, zeros(2), tcr
      REAL(KIND=8) :: n1(nsd), n2(nsd), t1(nsd), t2(nsd)
      REAL(KIND=8) :: vpar1, vpar2, vperp1, vperp2
      ! Coefficients to make calculating parallel/perp vel easier
      REAL(KIND=8) :: pa, pb

      prt1%collided = .false.
      prt1%collided = .false.

      ! First, check if particles will collide at current trajectory
      a = prt1%x(1)   - prt2%x(1)
      b = prt1%vel(1) - prt2%vel(1)
      c = prt1%x(2)   - prt2%x(2)
      d = prt1%vel(2) - prt2%vel(2)
      if(nsd.eq.3) then
         e = prt1%x(3)   - prt2%x(3)
         f = prt1%vel(3) - prt2%vel(3)
      else
         e=0D0
         f=0D0
      end if
      
      qa = b**2D0 + d**2D0 + f**2D0
      qb = 2D0*(a*b + c*d +e*f)
      qc = a**2D0 + c**2D0 + e**2D0 - ((prt1%dp + prt2%dp)/2D0)**2D0

      ! Imaginary zeros means particles won't collide
      if ((qb**2D0-4D0*qa*qc).lt.0) RETURN

      ! Zeros are when the particle either enters or leaves vicinity of other particle
      zeros(1) = (-qb + sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)
      zeros(2) = (-qb - sqrt(qb**2D0-4D0*qa*qc))/(2D0*qa)

      ! Negative zeros mean the particle would collide previously in time
      if (ANY(zeros.le.0D0)) RETURN

      tcr = minval(zeros)

      ! Exit function if collision won't occur during timestep
      if (tcr.gt.dtp) RETURN

      ! particle locations at point of collision
      prt1%x = prt1%vel*tcr + prt1%x
      prt2%x = prt2%vel*tcr + prt2%x

      ! Vector parallel and pependicular to collision tangent line
      n1 = (prt1%x - prt2%x)/((prt1%dp + prt2%dp)/2)
      n2 = -n1
      t1 = cross(cross(n1,prt1%vel),n1)
      t2 = cross(cross(n2,prt2%vel),n2)

      ! Rare case with no perpendicular velocity
      if (ANY(ISNAN(t1))) t1 = 0D0
      if (ANY(ISNAN(t2))) t2 = 0D0
      
      ! Get precollision parallel and perpendicular velocities
      vperp1 = sum(t1*prt1%vel)
      vpar1  = sum(n1*prt1%vel)
      vperp2 = sum(t2*prt2%vel)
      vpar2  = sum(n2*prt2%vel)

      ! Note that perpendicular velocities don't change, so we only need to calculate parallel
      pa = prt1%mp*vpar1 - prt2%mp*vpar2
      pb = (-vpar1 - vpar2)*prt1%k

      vpar2 = (pa - prt1%mp*pb)/(prt1%mp + prt2%mp)
      vpar1 = pb + vpar2
      vpar2 = -vpar2

      ! V here is split into just two velocities, so just add them as vector

      prt1%vel = vpar1*n1 + vperp1*t1
      prt2%vel = vpar2*n2 + vperp2*t2

      !!! Needs to be extended for multiple collisions per time step (will probably be here)
      !! Doesn't work if particle is still because of cross products
      ! Advance particle the rest of the time step at this velocity.
      prt1%x = prt1%x + prt1%vel*(dtp - tcr)
      prt2%x = prt2%x + prt2%vel*(dtp - tcr)

      prt1%collided = .true.
      prt2%collided = .true.

      prt1%remdtp = dtp-tcr
      prt2%remdtp = dtp-tcr


      END SUBROUTINE prtCollide

!#################################################################### CROSS
      ! I use cross products a couple times above. Also normalizes to unit vector
      FUNCTION cross(v1,v2)
      IMPLICIT NONE
      REAL(KIND=8) :: cross(nsd)
      REAL(KIND=8), INTENT(IN) :: v1(nsd), v2(nsd)

      cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
      cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
      cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

      cross = cross/sqrt(cross(1)**2D0 + cross(2)**2D0 + cross(3)**2D0)
      
      END FUNCTION cross


      END PROGRAM READVTK
