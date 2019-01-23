      ! Finds element particle of position x is in is in
      SUBROUTINE INTERP(xpt,pint,uint)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: xpt(nsd)
      REAL(KIND=8), INTENT(OUT):: pint, uint(nsd)
      INTEGER :: ii,jj,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,eNoN), &
       &               Nxi(nsd,eNoN), shps(eNoN), prntx(nsd)
    
       ! Setting up matrix for inversion
      Nxi(1,1) =  1D0
      Nxi(2,1) =  0D0
      Nxi(3,1) =  0D0
      Nxi(1,2) =  0D0
      Nxi(2,2) =  1D0
      Nxi(3,2) =  0D0
      Nxi(1,3) =  0D0
      Nxi(2,3) =  0D0
      Nxi(3,3) =  1D02
      Nxi(1,4) = -1D0
      Nxi(2,4) = -1D0
      Nxi(3,4) = -1D0

      xXi = 0D0
      pint = 0D0
      uint = 0D0

      do ii=1,nEl

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
            xXi(:,1) = xXi(:,1) + x(:,IEN(a,ii))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,IEN(a,ii))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,IEN(a,ii))*Nxi(3,a)
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
       prntx(a) =      xiX(a,1)*( xpt(1) - x(1,IEN(4,ii))) + &
     &                 xiX(a,2)*( xpt(2) - x(2,IEN(4,ii))) + &
     &                 xiX(a,3)*( xpt(3) - x(3,IEN(4,ii)))
         END DO
      END IF

      ! Finding shape functions
       shps(1) =  prntx(1)
       shps(2) =  prntx(2)
       shps(3) =  prntx(3)
       shps(4) = 1 -  prntx(1) -  prntx(2) -  prntx(3)
      ! If shape functions positive, we've found the correct element
      IF (ALL(shps.gt.0D0)) then
         EXIT
      END IF

   end do
   
   ! Using shape functions to interpolate value at coordinate
   do jj=1,eNoN
      pint = pint + pres(IEN(jj,ii))*shps(jj)
      uint = uint +vel(:,IEN(jj,ii))*shps(jj)
   end do
         
   END SUBROUTINE INTERP