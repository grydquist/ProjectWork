      PROGRAM GAUSSER
      IMPLICIT NONE
      INTEGER :: nsd=3,eNoN=4
      REAL(KIND=8) ::Nxi(3,4)
      


      REAL(KIND=8) :: x(3,4)
      REAL(KIND=8) :: Nx(3,4), Jac, ks(3,4), prntx(3), prtx(3)

      INTEGER :: a
      REAL(KIND=8) :: xXi(3,3), xiX(3,3)

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

      x(1,1) =  1D0
      x(2,1) = -3D0
      x(3,1) =  2D0
      x(1,2) =  5D0
      x(2,2) =  2D0
      x(3,2) =  1D0
      x(1,3) =  1D0
      x(2,3) =  0D0
      x(3,3) =  5D0
      x(1,4) = -1D0
      x(2,4) = -2D0
      x(3,4) = -3D0

      prtx(1) = 0.2D0
      prtx(2) = -1D0
      prtx(3) = -1D0

      nsd=3
      eNoN=4

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
         
         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3) &
     &       + xXi(1,2)*xXi(2,3)*xXi(3,1) &
     &       + xXi(1,3)*xXi(2,1)*xXi(3,2) &
     &       - xXi(1,1)*xXi(2,3)*xXi(3,2) &
     &       - xXi(1,2)*xXi(2,1)*xXi(3,3) &
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

         prntx = 0D0

         DO a=1, nsd
            prntx(a) = xiX(a,1)*(prtx(1) - x(1,4)) + &
            &          xiX(a,2)*(prtx(2) - x(2,4)) + &
            &          xiX(a,3)*(prtx(3) - x(3,4))
         END DO
         
print *, prntx
      END IF

      END PROGRAM GAUSSER

