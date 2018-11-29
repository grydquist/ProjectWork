      PROGRAM READVTK
      IMPLICIT NONE
!     This is for writing binary files
      CHARACTER, PARAMETER :: eol = CHAR(10)
!     This is the standard length for all the strings used in this code
      INTEGER, PARAMETER :: stdL = 256
!     Assuming element_number_of_node = 4, degrees_of_freedom = 4, and
!     number_of_spatial_dimensions = 3
      INTEGER, PARAMETER :: eNoN=4, dof=4, nsd=3, nG = 4
      LOGICAL isBinary, foundVar(2)
      INTEGER i, a, e, fid, m, iPos, nNo, nEl, tmpI(eNoN+1), g
      REAL(KIND=8) N(eNoN,nG), xi(nsd,eNoN), Nxi(nsd,eNoN), s, w(nG), t
      CHARACTER c
      CHARACTER(LEN=stdL) rLine, tmp, fName
      INTEGER, ALLOCATABLE :: IEN(:,:)
      REAL, ALLOCATABLE :: tmpS(:), tmpV(:,:)
!Grant Test
      REAL(KIND=8) :: xg(nsd,eNoN), Nx(nsd,eNoN), Jac, ks(nsd,nsd)
      REAL(KIND=8) :: prntx(nsd), prtx(nsd),shps(eNoN)
      INTEGER cnt
      INTEGER, ALLOCATABLE :: test (:,:)
      REAL(KIND=8),ALLOCATABLE ::  test2(:,:)
      INTEGER :: sbid,elid

      REAL(KIND=8), ALLOCATABLE :: vel(:,:), pres(:), X(:,:)

      cnt=1
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

      prtx(1) = 1D0
      prtx(2) = 1.5D0
      prtx(3) = 20D0
      
      CALL SBDomain(x,(/10,10,10/),test,test2)
      CALL xSB(prtx,test2,(/10,10,10/),sbid)
      CALL xEl(test(sbid,:),prtx,x,elid)

      print *, elid, shape(IEN),minval(x(3,:))

      do a=1,nEl
      xg(:,1) = x(:,IEN(1,a))
      xg(:,2) = x(:,IEN(2,a))
      xg(:,3) = x(:,IEN(3,a))
      xg(:,4) = x(:,IEN(4,a))
      CALL GNN(xg,Nx,Jac,ks,prntx,prtx,shps)
      if ((ALL(shps.gt.0))) EXIT
      cnt=cnt+1
      end do
      print *, cnt,shps
      
      STOP 
 001  STOP "A block of data is missing"

      CONTAINS
!--------------------------------------------------------------------         
      
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

!     Mapping Tau into the nodes by assembling it into a local vector
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

!####################################################################
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

!####################################################################
      PURE SUBROUTINE SBDomain(x,split,sbel,sbdim)
      REAL(KIND=8), INTENT(IN) :: x(nsd,nNo)
      INTEGER, INTENT(IN) :: split(nsd)
      REAL(KIND=8), INTENT(OUT),ALLOCATABLE :: sbdim(:,:)
      INTEGER, INTENT(OUT), ALLOCATABLE :: sbel(:,:)
      INTEGER :: ii,jj,cnt2,kk
      !LOGICAL :: inbox(eNoN)
      INTEGER, ALLOCATABLE :: seq1(:),seq2(:),seq3(:)
      REAL(KIND=8) :: diff(nsd),step(nsd), elbox(2*nsd,nEl)

      ! sbdim is the dimensions of each of the search boxes, with minx,maxx,miny,maxy,minz,maxz
      ALLOCATE(sbdim(nsd*2,split(1)*split(2)*split(3)))
      sbdim=0D0
      ! sbel is the id of the elements with a searchbox
      ALLOCATE(sbel(split(1)*split(2)*split(3),nEl))
      sbel=0D0

      ! these sequences are just for allocating sbdim
      ALLOCATE(seq1(split(3)*split(2)),seq2(split(3)*split(1)),seq3(split(2)*split(1)))

      ! Domain ranges
      diff(1)=MAXVAL(x(1,:))-MINVAL(x(1,:))
      diff(2)=MAXVAL(x(2,:))-MINVAL(x(2,:))
      diff(3)=MAXVAL(x(3,:))-MINVAL(x(3,:))
      ! Absolute searchbox dimensions
      step=diff/split

      seq1=(/(ii, ii=0, split(2)*split(3)-1, 1)/)*split(1)+1
      cnt2=0
      do ii=1,split(1)*split(3)
            seq2(ii)=ii+cnt2*(split(2)-1)*split(1)
            if (MOD(ii,split(1)).eq.0) cnt2=cnt2+1
      end do
      seq3=(/(ii, ii=0, split(1)*split(2)-1, 1)/)+1

      ! Allocating sbdim
      do ii=1,split(1)
         sbdim(1,seq1+ii-1)=MINVAL(x(1,:))+step(1)*(ii-1)
      end do

      do ii=1,split(2)
         sbdim(3,seq2+(ii-1)*split(1))=MINVAL(x(2,:))+step(2)*(ii-1)
      end do

      do ii=1,split(3)
         sbdim(5,seq3+(ii-1)*split(1)*split(2))=MINVAL(x(3,:))+step(3)*(ii-1)
      end do

      sbdim(2,:)=sbdim(1,:)+step(1)
      sbdim(4,:)=sbdim(3,:)+step(2)
      sbdim(6,:)=sbdim(5,:)+step(3)

      ! Making boxes surrounding elements
      do ii=1,Nel
         do jj=1,nsd
            elbox(2*jj-1,ii) = MINVAL(x(jj,IEN(:,ii)))
            elbox(2*jj  ,ii) = MAXVAL(x(jj,IEN(:,ii)))
         end do
      end do

      do ii=1,split(1)*split(2)*split(3)
         cnt2=1
         do jj=1,Nel
            ! Check if elements are completely outside searchbox
            do kk=1,nsd
               ! Cycle if min value elbox .gt. max value seachbox & vice-verse
               if (elbox(2*kk-1,jj).gt.sbdim(2*kk-1,ii)) cycle
               if (elbox(2*kk  ,jj).gt.sbdim(2*kk  ,ii)) cycle
            end do

            ! Old, dated method

            ! Checks if node values of element jj is in searchbox ii
            !inbox=((x(1,IEN(:,jj)).gt.sbdim(1,ii)).and.(x(1,IEN(:,jj)).lt.sbdim(2,ii)).and. &
            !          &   (x(2,IEN(:,jj)).gt.sbdim(3,ii)).and.(x(2,IEN(:,jj)).lt.sbdim(4,ii)).and. &
            !          &   (x(3,IEN(:,jj)).gt.sbdim(5,ii)).and.(x(3,IEN(:,jj)).lt.sbdim(6,ii))) 
            ! Puts element into searchbox 
            !if (any(inbox)) then
            sbel(ii,cnt2) = jj
            cnt2=cnt2+1
            !end if
         end do
      end do
      !
      ! Still need to work on corners
      !
      END SUBROUTINE SBDomain

!####################################################################
      ! Find which Searchbox x is in
      SUBROUTINE xSB(x,sbdim,split,sbid)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN),ALLOCATABLE:: sbdim(:,:)
      REAL(KIND=8), INTENT(IN) :: x(nsd)
      INTEGER, INTENT(IN) :: split(nsd)
      INTEGER, INTENT(OUT) :: sbid
      REAL(KIND=8) :: step(nsd),xzero(nsd)
      INTEGER :: xsteps(nsd)

      ! Searchbox dimensions
      step(1)=sbdim(2,1)-sbdim(1,1)
      step(2)=sbdim(4,1)-sbdim(3,1)
      step(3)=sbdim(6,1)-sbdim(5,1)

      ! Set domain back to zero
      xzero(1)=x(1)-minval(sbdim(1,:))
      xzero(2)=x(2)-minval(sbdim(3,:))
      xzero(3)=x(3)-minval(sbdim(5,:))

      ! Find which searchbox the particle is in
      ! Number of searchbox steps in x,y,and z
      xsteps=FLOOR(xzero/step)
      sbid=xsteps(1)+split(1)*xsteps(2)+split(1)*split(2)*xsteps(3)+1
      ! Could likely move xEl up here
      !
      !
      END SUBROUTINE xSB

!####################################################################
      ! Finds element x is in
      PURE SUBROUTINE xEl(sbel,prtx,x,elid)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: prtx(nsd), x(nsd,nNo)
      INTEGER, INTENT(IN) :: sbel(nEl)
      INTEGER, INTENT(OUT) :: elid
      INTEGER :: ii,cnt,a
      REAL(KIND=8) :: Jac,prntx(nsd),xXi(nsd,nsd), xiX(nsd,nsd),shps(eNoN),Nx(nsd,eNoN)
      cnt=1

      do ii=1,nEl

      prntx = 0D0
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
            xXi(:,1) = xXi(:,1) + x(:,IEN(a,sbel(ii)))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,IEN(a,sbel(ii)))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,IEN(a,sbel(ii)))*Nxi(3,a)
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
            prntx(a) = xiX(a,1)*(prtx(1) - x(1,IEN(4,sbel(ii)))) + &
     &                 xiX(a,2)*(prtx(2) - x(2,IEN(4,sbel(ii)))) + &
     &                 xiX(a,3)*(prtx(3) - x(3,IEN(4,sbel(ii))))
         END DO
      END IF

      shps(1)=prntx(1)
      shps(2)=prntx(2)
      shps(3)=prntx(3)
      shps(4)=1-prntx(1)-prntx(2)-prntx(3)

      IF (ALL(shps.gt.0D0)) then
         elid=sbel(ii)
         EXIT
      END IF
         
      end do

      END SUBROUTINE xEl

      END PROGRAM READVTK
