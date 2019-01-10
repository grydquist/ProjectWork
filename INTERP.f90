      ! Finds element particle of position x is in is in
      SUBROUTINE INTERP(IEN,el,x,xpt,uin,u)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: IEN(eNoN,nEl)
      INTEGER,      INTENT(IN) :: el
      REAL(KIND=8), INTENT(IN) :: x(nsd,eNoN), xpt(nsd), uin(nsd,eNoN)
      REAL(KIND=8), INTENT(OUT):: u(nsd)
      INTEGER :: ii,jj,a
      REAL(KIND=8) :: Jac,xXi(nsd,nsd), xiX(nsd,nsd),Nx(nsd,eNoN),&
                      Nxi(nsd,eNoN), shps(eNoN), prntx(nsd)
    
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
            xXi(:,1) = xXi(:,1) + x(:,IEN(a,el))*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,IEN(a,el))*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,IEN(a,el))*Nxi(3,a)
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
       prntx(a) =      xiX(a,1)*( xpt(1) - x(1,IEN(4,el))) + &
     &                 xiX(a,2)*( xpt(2) - x(2,IEN(4,el))) + &
     &                 xiX(a,3)*( xpt(3) - x(3,IEN(4,el)))
         END DO
      END IF

       shps(1) =  prntx(1)
       shps(2) =  prntx(2)
       shps(3) =  prntx(3)
       shps(4) = 1 -  prntx(1) -  prntx(2) -  prntx(3)

     do ii=1,nsd
        do jj=1,eNoN
           u(ii) = u(ii) + uin(ii,jj)*shps(jj)
        end do
     end do
         
      END SUBROUTINE INTERP