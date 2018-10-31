!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!-----------------------------------------------------------------------
!
!     Here, all functions for treating immersed boundaries are defined.
!
!--------------------------------------------------------------------

!     This routine reads IB mesh data
      SUBROUTINE IB_READMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL :: flag
      INTEGER :: i, j, iM, iFa, a, b, Ac, e
      REAL(KIND=8) :: scaleF, fibN(nsd), rtmp
      CHARACTER(LEN=stdL) :: ctmp
      TYPE(listType), POINTER :: lPtr, lPM
      TYPE(fileType) :: fTmp

      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:), gX(:,:)

      ib%nMsh  = list%srch("Add IB",ll=1)
      std = " Number of immersed boundaries: "//ib%nMsh
      ALLOCATE (ib%msh(ib%nMsh), gX(0,0))

      ib%tnNo = 0
      DO iM=1, ib%nMsh
         lPM => list%get(ib%msh(iM)%name, "Add IB", iM)
         lPtr => lPM%get(ib%msh(iM)%lShl, "Set mesh as shell")

         std  = " Reading IB mesh <"//CLR(TRIM(ib%msh(iM)%name))//">"
         CALL READSV(lPM, ib%msh(iM))
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READCCNE(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READNRB(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READGAMBIT(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            err = " Failed to identify format of the IB mesh"
         END IF

         std = " Number of IB nodes: "//ib%msh(iM)%gnNo
         std = " Number of IB elements: "//ib%msh(iM)%gnEl

!     Making sure face names are unique
         DO iFa=1, ib%msh(iM)%nFa
            ib%msh(iM)%fa(iFa)%iM = iM
            ctmp = ib%msh(iM)%fa(iFa)%name
            DO i=1, iM
               DO j=1, ib%msh(i)%nFa
                  IF (ctmp .EQ. ib%msh(i)%fa(j)%name .AND.
     2               (i.NE.iM .OR. j.NE.iFa)) THEN
                     err = " Repeating face names is not allowed"
                  END IF
               END DO
            END DO
         END DO

!     To scale the mesh, while attaching x to gX
         scaleF = 1D0
         lPtr => lPM%get(scaleF, "Mesh scale factor", lb=0D0)
         a = ib%tnNo + ib%msh(iM)%gnNo
         IF (iM .GT. 1) THEN
            ALLOCATE(tmpX(nsd,ib%tnNo))
            tmpX = gX
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
            gX(:,1:ib%tnNo) = tmpX
            DEALLOCATE(tmpX)
         ELSE
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
         END IF
         gX(:,ib%tnNo+1:a) = ib%msh(iM)%x * scaleF
         ib%tnNo           = a
         DEALLOCATE(ib%msh(iM)%x)

         ib%msh(iM)%dx = 0.1D0
         lPtr => lPM%get(rtmp,"Mesh edge size")
         IF (ASSOCIATED(lPtr)) ib%msh(iM)%dx = rtmp
      END DO
      ALLOCATE(ib%x(nsd,ib%tnNo))
      ib%x = gX
      DEALLOCATE(gX)

!     Checks for shell elements
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
            IF (ib%msh(iM)%eType .NE. eType_TRI) THEN
               err = "Immersed shell elements can only be triangles"
            END IF
         END IF
      END DO

!     Setting msh%gN, msh%lN parameter
      b = 0
      DO iM=1, ib%nMsh
         ib%msh(iM)%nNo = ib%msh(iM)%gnNo
         ALLOCATE(ib%msh(iM)%gN(ib%msh(iM)%nNo), ib%msh(iM)%lN(ib%tnNo))
         ib%msh(iM)%gN = 0
         ib%msh(iM)%lN = 0
         DO a=1, ib%msh(iM)%nNo
            b = b + 1
            ib%msh(iM)%gN(a) = b
            ib%msh(iM)%lN(b) = a
         END DO
      END DO
      IF (b .NE. ib%tnNo) err =
     2   " Mismatch in ib%tnNo. Correction needed"

!     Remap msh%gIEN array
      DO iM=1, ib%nMsh
         ib%msh(iM)%nEl = ib%msh(iM)%gnEl
         ALLOCATE(ib%msh(iM)%IEN(ib%msh(iM)%eNoN,ib%msh(iM)%nEl))
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%gIEN(a,e)
               Ac = ib%msh(iM)%gN(Ac)
               ib%msh(iM)%IEN(a,e) = Ac
            END DO
         END DO
         DEALLOCATE(ib%msh(iM)%gIEN)
      END DO

!     Re-arranging fa structure
      b = 0
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            ALLOCATE(ib%msh(iM)%fa(iFa)%lN(ib%tnNo))
            ib%msh(iM)%fa(iFa)%lN = 0
            DO a=1, ib%msh(iM)%fa(iFa)%nNo
               Ac = ib%msh(iM)%fa(iFa)%gN(a)
               Ac = ib%msh(iM)%gN(Ac)
               ib%msh(iM)%fa(iFa)%gN(a)  = Ac
               ib%msh(iM)%fa(iFa)%lN(Ac) = a
            END DO
            DO e=1, ib%msh(iM)%fa(iFa)%nEl
               DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                  Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                  Ac = ib%msh(iM)%gN(Ac)
                  ib%msh(iM)%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         END DO
      END DO

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId.
      DO iM=1, nMsh
         lPM => list%get(ib%msh(iM)%name,"Add IB",iM)

         lPtr => lPM%get(fTmp,"Domain (IB) file path")
         IF (ASSOCIATED(lPtr)) CALL SETDMNIDFF(ib%msh(iM), fTmp%open())

         lPtr => lPM%get(i,"Domain (IB)",ll=0,ul=BIT_SIZE(ib%dmnId)-1)
         IF (ASSOCIATED(lPtr)) CALL SETDMNID(ib%msh(iM), i)
         IF (.NOT.ALLOCATED(ib%msh(iM)%eId)) err =
     2      " Immersed bodies require domain ID parameter to be set"
      END DO
      ALLOCATE(ib%dmnId(ib%tnNo))
      ib%dmnId = 0
      DO iM=1, ib%nMsh
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               ib%dmnId(Ac) = IOR(ib%dmnId(Ac),ib%msh(iM)%eId(e))
            END DO
         END DO
      END DO

!     Read fiber orientation
      ib%nFn = 0
      flag   = .FALSE.
      DO iM=1, ib%nMsh
         lPM => list%get(ib%msh(iM)%name,"Add IB",iM)
         j = lPM%srch("Fiber direction file path")
         IF (j .EQ. 0) j = lPM%srch("Fiber direction")
         IF (ib%nFn .LT. j) ib%nFn = j
      END DO
      IF (ib%nFn .GT. 0) flag = .TRUE.

      IF (flag) THEN
         DO iM=1, ib%nMsh
            lPM => list%get(ib%msh(iM)%name,"Add IB",iM)

            j = lPM%srch("Fiber direction file path")
            IF (j .NE. 0) THEN
               ALLOCATE(ib%msh(iM)%x(ib%nFn*nsd,ib%msh(iM)%gnNo))
               ib%msh(iM)%x = 0D0
               DO i=1, j
                  lPtr => lPM%get(cTmp, "Fiber direction file path", i)
                  CALL READVTUPDATA(ib%msh(iM), cTmp, "FIB_DIR", nsd, i)
               END DO
            ELSE
               j = lPM%srch("Fiber direction")
               IF (j .NE. 0) THEN
                  ALLOCATE(ib%msh(iM)%x(ib%nFn*nsd,msh(iM)%gnNo))
                  ib%msh(iM)%x = 0D0
                  DO i=1, j
                     lPtr => lPM%get(fibN, "Fiber direction", i)
                     DO a=1, ib%msh(iM)%gnNo
                        ib%msh(iM)%x((i-1)*nsd+1:i*nsd,a) = fibN(1:nsd)
                     END DO
                  END DO
               END IF
            END IF
         END DO

         ALLOCATE(ib%fN(ib%nFn*nsd,ib%tnNo))
         ib%fN = 0D0
         DO iM=1, ib%nMsh
            IF (.NOT.ALLOCATED(ib%msh(iM)%x)) CYCLE
            DO a=1, ib%msh(iM)%gnNo
               Ac = ib%msh(iM)%gN(a)
               DO i=1, ib%nFn
                  fibN(:) = ib%msh(iM)%x((i-1)*nsd+1:i*nsd,a)
                  rtmp = SQRT(NORM(fibN))
                  IF (.NOT.ISZERO(rtmp)) fibN(:) = fibN(:)/rtmp
                  ib%fN((i-1)*nsd+1:i*nsd,Ac) = fibN
               END DO
            END DO
            DEALLOCATE(ib%msh(iM)%x)
         END DO
      ELSE
         ib%nFn = 1
      END IF

      IF (ib%nMsh .GT. 1) THEN
         std = " Total number of IB nodes: "//ib%tnNo
         std = " Total number of IB elements: "//SUM(ib%msh%nEl)
      END IF

      std = CLR(" IB mesh data imported successfully",3)

      RETURN
      END SUBROUTINE IB_READMSH
!####################################################################
!     This routine reads IB options
      SUBROUTINE IB_READOPTS(list)
      USE COMMOD
      USE LISTMOD
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL lShl
      INTEGER iM
      CHARACTER(LEN=stdL) :: ctmp
      TYPE(listType), POINTER :: lIBs, lPtr

      lIBs => list%get(ctmp, "IB settings")
      IF (.NOT.ASSOCIATED(lIBs)) err = "IB settings not provided"

      lShl = .FALSE.
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
            lShl = .TRUE.
            EXIT
         END IF
      END DO

      lPtr => lIBs%get(ctmp, "Method")
      IF (ASSOCIATED(lPtr)) THEN
         SELECT CASE (ctmp)
         CASE ("SSM")
            ib%mthd = ibMthd_SSM
         CASE ("Penalty")
            ib%mthd = ibMthd_Penalty
         CASE ("Nitsche")
            ib%mthd = ibMthd_Nitsche
         CASE ("IFEM")
            ib%mthd = ibMthd_IFEM
         CASE DEFAULT
            err = " Invalid IB method"
         END SELECT
      ELSE
         std = " Choosing default IB method "//CLR("SSM")
         ib%mthd = ibMthd_SSM
      END IF

      lPtr => lIBs%get(ib%fbFlag, "Use feedback forcing")
      IF (.NOT.ASSOCIATED(lPtr)) THEN
         IF (ib%mthd .EQ. ibMthd_Penalty) ib%fbFlag = .TRUE.
      END IF

      lPtr => lIBs%get(ib%fcFlag, "Use finite cell integration")
      IF (.NOT.ASSOCIATED(lPtr)) THEN
         IF (ib%mthd .EQ. ibMthd_Penalty .OR.
     2       ib%mthd .EQ. ibMthd_Nitsche) ib%fcFlag = .TRUE.
      END IF

      IF (lShl) THEN
         IF (ib%mthd .NE. ibMthd_Penalty) err = "Immersed shells "//
     2      "are treated using penalty formulation only"
         IF (.NOT.ib%fbFlag) err = "Feedback forcing should be "//
     2      "applied for immersed shells"
         IF (ib%fcFlag) err = "Finite cell method is not applicable "//
     2      "for immersed shells"
      END IF

      IF (ib%mthd .EQ. ibMthd_IFEM) THEN
         ib%fbFlag = .FALSE.
         ib%fcFlag = .FALSE.
      END IF

      RETURN
      END SUBROUTINE IB_READOPTS
!####################################################################
!     This routine reads IB domain properties and BCs in a given Eq
      SUBROUTINE IB_READEQ(lEq, list, eqName)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE

      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=stdL), INTENT(IN) :: eqName

      INTEGER, PARAMETER :: maxOutput = 5

      LOGICAL flag
      INTEGER i, fid, iBc, propL(maxNProp), iDmn, iProp, prop, nDOP(4),
     2   outputs(maxOutput)
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPD, lSFp, lPBC

      IF (.NOT.ALLOCATED(ib%msh)) err = " No IB mesh is read yet"
      IF (.NOT.ALLOCATED(ib%dmnId)) err = " No IB domain is read yet"

      lEq%nDmnIB = list%srch("Domain (IB)")
      flag = .FALSE.
      IF (lEq%nDmnIB .EQ. 0) THEN
         IF (ib%mthd .EQ. ibMthd_IFEM) err =
     2      "No IB domain provided for IFEM method"
         flag = .TRUE.
         lEq%nDmnIB = 1
         ALLOCATE(lEq%dmnIB(lEq%nDmnIB))
         lEq%dmnIB(1)%Id = -1
         lEq%dmnIB(1)%phys = phys_NA
      ELSE
         ALLOCATE(lEq%dmnIB(lEq%nDmnIB))
      END IF

!--------------------------------------------------------------------
!     Searching through each IB domain properties
      DO iDmn=1, lEq%nDmnIB
         IF (flag) EXIT
         lPD => list%get(lEq%dmnIB(iDmn)%Id,"Domain (IB)",iDmn,
     2      ll=0,ul=(BIT_SIZE(ib%dmnId)-1))
         DO i=1, iDmn-1
            IF (lEq%dmnIB(iDmn)%Id .EQ. lEq%dmnIB(i)%Id) THEN
               err = TRIM(list%ping("Domain (IB)",lPD))//
     2            " Repeated IB domain ID"
            END IF
         END DO

!        IB Equation being solved: shell/struct/lElas
         lPtr => lPD%get(ctmp, "Equation", 1)
         propL = prop_NA
         SELECT CASE(TRIM(ctmp))
         CASE("shell")
            lEq%dmnIB(iDmn)%phys = phys_shell
            IF (nsd .NE. 3) err = "Shell mechanics can be solved only"//
     2         " in 3-dimensions"
            propL(1) = solid_density
            propL(2) = damping
            propL(3) = elasticity_modulus
            propL(4) = poisson_ratio
            propL(5) = shell_thickness
            propL(6) = f_x
            propL(7) = f_y
            propL(8) = f_z

            nDOP = (/2,1,0,0/)
            outPuts(1) = out_displacement
            outPuts(2) = out_velocity

         CASE("struct")
            lEq%dmnIB(iDmn)%phys = phys_struct
            propL(1) = solid_density
            propL(2) = elasticity_modulus
            propL(3) = poisson_ratio
            propL(4) = f_x
            propL(5) = f_y
            IF (nsd .EQ. 3) propL(6) = f_z

            nDOP = (/5,1,0,0/)
            outPuts(1) = out_displacement
            outPuts(2) = out_velocity
            outPuts(3) = out_pressure
            outPuts(4) = out_strainInv
            outPuts(5) = out_acceleration

         CASE DEFAULT
            err = TRIM(lPD%ping("Equation (IB)",lPtr))//
     2         "IB must be solved using shell/struct equation"
         END SELECT

!        Domain properties
         DO iProp=1, maxNProp
            rtmp = 0D0
            prop = propL(iProp)
            SELECT CASE (prop)
            CASE (prop_NA)
               EXIT
            CASE (solid_density)
               lPtr => lPD%get(rtmp,"Density",1,ll=0D0)
            CASE (elasticity_modulus)
               lPtr => lPD%get(rtmp,"Elasticity modulus",1,lb=0D0)
            CASE (poisson_ratio)
               lPtr => lPD%get(rtmp,"Poisson ratio",1,ll=0D0,ub=5D-1)
            CASE (f_x)
               lPtr => lPD%get(rtmp,"Force_X")
            CASE (f_y)
               lPtr => lPD%get(rtmp,"Force_Y")
            CASE (f_z)
               lPtr => lPD%get(rtmp,"Force_Z")
            CASE (source_term)
               lPtr => lPD%get(rtmp,"Source term")
            CASE (damping)
               lPtr => lPD%get(rtmp,"Mass damping")
            CASE (shell_thickness)
               lPtr => lPD%get(rtmp,"Shell thickness",1,lb=0D0)
            CASE DEFAULT
               err = "Undefined properties (IB)"
            END SELECT
            lEq%dmnIB(iDmn)%prop(prop) = rtmp
         END DO

         IF (lEq%dmnIB(iDmn)%phys .EQ. phys_shell) THEN
            lSFp => lPD%get(ctmp,"Follower pressure load")
            IF (ASSOCIATED(lSFp)) THEN
               ALLOCATE(lEq%dmnIB(iDmn)%shlFp)
               CALL READShlFp(lSFp, lEq%dmnIB(iDmn)%shlFp, ctmp)
            END IF
         END IF

         IF (lEq%dmnIB(iDmn)%phys .EQ. phys_struct) THEN
            CALL READMATMODEL(lEq%dmnIB(iDmn), lPD)
         END IF
      END DO

      IF (.NOT.flag) THEN
!        Read IB outputs
         CALL IB_READOUTPUTS(lEq, nDOP, outPuts, list)

!        IB linear solver properties:
!        Default: CG, relTol=1E-2; absTol=1E-10; maxItr=1000)
         CALL FSILS_LS_CREATE(lEq%lsIB, LS_TYPE_CG)
!         SUBROUTINE FSILS_LS_CREATE(ls, LS_type, relTol, absTol, &
!     &      maxItr, dimKry, relTolIn, absTolIn, maxItrIn)
      END IF

!--------------------------------------------------------------------
!     Searching for BCs on immersed bodies
      lEq%nBcIB = list%srch("Add BC (IB)")
      ALLOCATE(lEq%bcIB(lEq%nBcIB))
      std = " Number of imposed BC for equation <"//TRIM(eqName)//
     2   ">: "//lEq%nBcIB
      DO iBc=1, lEq%nBcIB
         lPBC => list%get(ctmp,"Add BC (IB)",iBc)
         CALL IB_FINDFACE(ctmp, lEq%bcIB(iBc)%iM, lEq%bcIB(iBc)%iFa)
         CALL IB_READBC(lEq%bcIB(iBc), lPBC)
      END DO

      RETURN
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE READShlFp(list, lFp, ctmp)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list
      TYPE(shlFpType), INTENT(INOUT) :: lFp
      CHARACTER(LEN=stdL), INTENT(IN) :: ctmp

      INTEGER i, j
      TYPE(fileType) fTmp

      SELECT CASE (TRIM(ctmp))
      CASE ("steady")
         lFp%bType = IBSET(lFp%bType, bType_std)
         lPtr => list%get(lFp%p,"Value",1)
      CASE ("unsteady")
         lFp%bType = IBSET(lFp%bType, bType_ustd)
         ALLOCATE(lFp%pt)
         lPtr => list%get(ftmp, "Temporal values file path")
         IF (ASSOCIATED(lPtr)) THEN
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lFp%pt%n = j
            ALLOCATE(lFp%pt%r(j))
            ALLOCATE(lFp%pt%i(j))
            CALL FFT(fid, i, lFp%pt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
            fid = fTmp%open()
            READ (fid,*) lFp%pt%ti
            READ (fid,*) lFp%pt%T
            READ (fid,*) lFp%pt%qi
            READ (fid,*) lFp%pt%qs
            READ (fid,*) j
            lFp%pt%n = j
            ALLOCATE(lFp%pt%r(j))
            ALLOCATE(lFp%pt%i(j))
            DO i=1, j
               READ (fid,*) lFp%pt%r(i), lFp%pt%i(i)
            END DO
            CLOSE(fid)
         END IF
      CASE DEFAULT
         err = "Undefined follower load type"
      END SELECT

      RETURN
      END SUBROUTINE READShlFp
!--------------------------------------------------------------------
      END SUBROUTINE IB_READEQ
!####################################################################
!     This subroutine is to read from input file on how to process the
!     output quantities: for VTK files or for surface/volume integrated
!     quantities. nDOP(1) is the total number of outputs, nDOP(2) is the
!     default number of outputs for VTK files, nDOP(3) is for boundaries
!     nDOP(4) is for volume
      SUBROUTINE IB_READOUTPUTS(lEq, nDOP, outputs, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER, INTENT(IN) :: nDOP(4), outputs(nDOP(1))
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER nOut, iOut, i, j
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPO

      lEq%nOutIB = nDOP(1)
      ALLOCATE(lEq%outIB(nDOP(1)))

      DO iOut=1, nDOP(1)
         SELECT CASE (outputs(iOut))
         CASE (out_velocity)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Velocity"
         CASE (out_pressure)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = nsd
            lEq%outIB(iOut)%l    = 1
            lEq%outIB(iOut)%name = "Pressure"
         CASE (out_acceleration)
            lEq%outIB(iOut)%grp  = outGrp_A
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Acceleration"
         CASE (out_displacement)
            lEq%outIB(iOut)%grp  = outGrp_D
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Displacement"
         CASE (out_strainInv)
            lEq%outIB(iOut)%grp  = outGrp_stInv
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Strain_invariants"
         CASE DEFAULT
            err = "Internal output undefined"
         END SELECT
      END DO

!     These are the default values, we use the first nDef/nBDef outputs
      DO j=1, 3
         lEq%outIB(1:nDOP(j+1))%wtn(j) = .TRUE.
      END DO

!     First reading the outputs for VTK files and then for boundaries
!     and last for the volume
      nOut = list%srch("Output (IB)")
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output (IB)",iOut)
         SELECT CASE(TRIM(ctmp))
         CASE("Spatial")
            j = 1
         CASE("B_INT")
            j = 2
         CASE("V_INT")
            j = 3
         CASE DEFAULT
            j = -1
            err = TRIM(list%ping("Output (IB)",lPO))//
     2         " Undefined keyword"
         END SELECT
         DO i=1, lEq%nOutIB
            lPtr => lPO%get(lEq%outIB(i)%wtn(j),lEq%outIB(i)%name)
            IF (TRIM(lEq%outIB(i)%name) .EQ. "Vortex") THEN
               IF (nsd .ne. maxNSD) THEN
                  lEq%outIB(i)%wtn(j) = .FALSE.
               END IF
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE IB_READOUTPUTS
!####################################################################
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE IB_FINDFACE(faName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE

      CHARACTER(LEN=stdL) :: faName
      INTEGER, INTENT(OUT) :: iM, iFa

      iFa = 0
      MY_LOOP : DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
            IF (ib%msh(iM)%name .EQ. faName) THEN
               iFa = 0
               EXIT
            END IF
         END IF
         DO iFa=1, ib%msh(iM)%nFa
            IF (ib%msh(iM)%fa(iFa)%name .EQ. faName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP

      IF (iM .GT. ib%nMsh) err =
     2   " Unable to find face (IB) <"//TRIM(faName)//">"

      RETURN
      END SUBROUTINE IB_FINDFACE
!####################################################################
!     This routine reads BC for immersed bodies
      SUBROUTINE IB_READBC(lBc, list)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE

      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL ltmp
      INTEGER a, b, i, j, Ac, iM, iFa, fid, nNo
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) fTmp

      INTEGER, ALLOCATABLE :: ptr(:)

      iM  = lBc%iM
      iFa = lBc%iFa
      IF (iFa .NE. 0) THEN
         nNo = msh(iM)%fa(iFa)%nNo
      ELSE
         nNo = msh(iM)%gnNo
      END IF

!     Reading the type: Dir/Neu
      lPtr => list%get(ctmp,"Type")
      SELECT CASE (ctmp)
      CASE ("Dirichlet","Dir")
         lBc%bType = IBSET(lBc%bType,bType_Dir)
      CASE ("Neumann","Neu")
         lBc%bType = IBSET(lBc%bType,bType_Neu)
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected BC type"
      END SELECT

!     Weak Dirichlet BC for fluid/FSI equations
      lBc%weakDir = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (ib%mthd .EQ. ibMthd_Nitsche) lBc%weakDir = .TRUE.
         lPtr => list%get(ltmp, "Weakly applied")
         IF (ASSOCIATED(lPtr)) lBc%weakDir = ltmp

         IF (lBc%weakDir .AND. ib%mthd.NE.ibMthd_Nitsche) err =
     2      "Dir BC can be applied weakly with only Nitsche method"

!        Read penalty values
         lPtr => list%get(rtmp, "Penalty constant")
         IF (ASSOCIATED(lPtr)) lBc%tauB(:) = rtmp
         lPtr => list%get(rtmp, "Penalty constant (tangential)")
         IF (ASSOCIATED(lPtr)) lBc%tauB(1) = rtmp
         lPtr => list%get(rtmp, "Penalty constant (normal)")
         IF (ASSOCIATED(lPtr)) lBc%tauB(2) = rtmp

!        Read feedback parameter
         lPtr => list%get(rtmp, "Feedback force constant")
         IF (ASSOCIATED(lPtr)) THEN
            IF (.NOT.ib%fbFlag) wrn = "Feedback forcing option is "//
     2         "turned off"
            lBc%tauF = rtmp
            lPtr => list%get(lBc%fbN, "Feedback along normal direction")
         END IF

         IF (ib%msh(iM)%lShl .AND. iFa.EQ.0) THEN
            IF (lBc%weakDir) err = "Dir cannot be weakly applied "//
     2         "on immersed shells"
            IF (.NOT.ib%fbFlag) err = "Feedback forcing should be "//
     2         "applied for immersed shells"
            IF (.NOT.lBc%fbN) err = "Feedback forcing should be "//
     2         "applied along normal direction for immersed shells"
         END IF
      END IF

!     Read BCs for shells with triangular elements. Not necessary for
!     NURBS elements
      lPtr => list%get(ctmp,"Shell BC type")
      IF (ASSOCIATED(lPtr)) THEN
         SELECT CASE (ctmp)
         CASE ("Fixed", "fixed", "Clamped", "clamped")
            lBc%bType = IBSET(lBc%bType,bType_fix)
            IF (.NOT.BTEST(lBc%bType,bType_Dir)) err = "Fixed BC "//
     2         "can be applied for Dirichlet boundaries only"
         CASE ("Hinged", "hinged")
            lBc%bType = IBSET(lBc%bType,bType_hing)
            IF (.NOT.BTEST(lBc%bType,bType_Dir)) err = "Hinged BC "//
     2         "can be applied for Dirichlet boundaries only"
         CASE ("Free", "free")
            lBc%bType = IBSET(lBc%bType,bType_free)
            IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "Free BC "//
     2         "can be applied for Neumann boundaries only"
!       Symm BC needs to be verified
c         CASE ("Symm", "symm", "Symmetric", "symmetric")
c            lBc%bType = IBSET(lBc%bType,bType_symm)
c            IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "Symm BC "//
c     2         "can be applied for Neumann boundaries only"
         CASE DEFAULT
            err = TRIM(list%ping("Shell BC type",lPtr))//
     2         " Unexpected Shell BC type"
         END SELECT
      END IF

      ALLOCATE(lBc%eDrn(nsd))
      lBc%eDrn = 0
      lPtr => list%get(lBc%eDrn,"Effective direction")
      ctmp = "Steady"
      lPtr => list%get(ctmp,"Time dependence")
      SELECT CASE (ctmp)
      CASE ('Steady')
         lBc%bType = IBSET(lBc%bType,bType_std)
         lPtr => list%get(lBc%g,"Value",1)
      CASE ('Unsteady')
         lBc%bType = IBSET(lBc%bType,bType_ustd)
         ALLOCATE(lBc%gt)
         lPtr => list%get(fTmp,"Temporal values file path")
         IF (ASSOCIATED(lPtr)) THEN
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lBc%gt%n = j
            ALLOCATE(lBc%gt%r(j))
            ALLOCATE(lBc%gt%i(j))
            CALL FFT(fid, i, lBc%gt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
            fid = fTmp%open()
            READ (fid,*) lBc%gt%ti
            READ (fid,*) lBc%gt%T
            READ (fid,*) lBc%gt%qi
            READ (fid,*) lBc%gt%qs
            READ (fid,*) j
            lBc%gt%n = j
            ALLOCATE(lBc%gt%r(j))
            ALLOCATE(lBc%gt%i(j))
            DO i=1, j
               READ (fid,*) lBc%gt%r(i), lBc%gt%i(i)
            END DO
            CLOSE(fid)
         END IF
      CASE ('Coupled')
         err = " Cannot apply Coupled BCs for immersed bodies"
      CASE ('Resistance')
         err = " Cannot apply Resistance BCs for immersed bodies"
      CASE ('General')
         lBc%bType = IBSET(lBc%bType,bType_gen)
         lPtr =>list%get(fTmp,"Temporal and spatial values file path",1)
         fid = fTmp%open()
         READ (fid,*) i, j, a
         IF (a .NE. nNo) THEN
            IF (iFa .NE. 0) THEN
               err = "Number of nodes does not match between "//
     2            TRIM(ib%msh(iM)%fa(iFa)%name)//" and "//fTmp%fname
            ELSE
               err = "Number of nodes does not match between "//
     2            TRIM(ib%msh(iM)%name)//" and "//fTmp%fname
            END IF
         END IF
         IF (i.LT.1 .OR. i.GT.nsd) err = "0 < dof <= "//nsd//
     2      " is violated in "//fTmp%fname

         ALLOCATE(lBc%gm)
         ALLOCATE(lBc%gm%t(j), lBc%gm%d(i,a,j), ptr(ib%msh(iM)%gnNo))
!     I am seting all the nodes to zero just in case a node is not set
         lBc%gm%d   = 0D0
         lBc%gm%dof = i
         lBc%gm%nTP = j
         ptr        = 0
!     Preparing the pointer array
         IF (iFa .NE. 0) THEN
            DO a=1, ib%msh(iM)%fa(iFa)%nNo
               Ac = ib%msh(iM)%fa(iFa)%gN(a)
               Ac = ib%msh(iM)%lN(Ac)
               IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2            "detected for BC. Mesh: "//TRIM(ib%msh(iM)%name)//
     3            ", Face: "//TRIM(ib%msh(iM)%fa(iFa)%name)//
     4            ", Node: "//STR(a)//" gN: "//
     5            STR(ib%msh(iM)%fa(iFa)%gN(a))
               ptr(Ac) = a
            END DO
         ELSE
            DO a=1, ib%msh(iM)%gnNo
               ptr(a) = a
            END DO
         END IF

         DO i=1, j
            READ (fid,*) rtmp
            lBc%gm%t(i) = rtmp
            IF (i .EQ. 1) THEN
               IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2            " should be zero in <"//TRIM(ftmp%fname)//">"
            ELSE
               rtmp = rtmp - lBc%gm%t(i-1)
               IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err = "Non-in"//
     2            "creasing time trend is found in <"//
     3            TRIM(ftmp%fname)//">"
            END IF
         END DO

         lBc%gm%period = lBc%gm%t(j)
         DO b=1, nNo
            READ (fid,*) Ac
            IF (Ac.GT.ib%msh(iM)%gnNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//ftmp%fname
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry "//b//" not found in "//ftmp%fname
            END IF
            DO i=1, j
               READ (fid,*) lBc%gm%d(:,a,i)
            END DO
         END DO
         CLOSE(fid)
         DEALLOCATE(ptr)
      CASE DEFAULT
         err = TRIM(list%ping("Time dependence",lPtr))//
     2      " Unexpected type"
      END SELECT

!     To zero-out perimeter or not. Default is .true. for Dir
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) ltmp = .TRUE.
      lPtr => list%get(ltmp,"Zero out perimeter")
      lBc%bType = IBCLR(lBc%bType,bType_zp)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_zp)

!     Impose BC on the state variable or its integral
      lPtr => list%get(ltmp,"Impose on state variable integral")
      lBc%bType = IBCLR(lBc%bType,bType_impD)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_impD)

!     Reading the spatial profile: flat/para/ud
      ctmp = "Flat"
      lPtr => list%get(ctmp,"Profile")
      SELECT CASE (ctmp)
      CASE ('Flat')
         lBc%bType = IBSET(lBc%bType,bType_flat)
      CASE ('Parabolic')
         IF (iFa .EQ. 0) err = "Parabolic profile not yet set up for "//
     2      "immersed shells"
         lBc%bType = IBSET(lBc%bType,bType_para)
      CASE ('D_dependent')
         lBc%bType = IBSET(lBc%bType,bType_ddep)
      CASE ('User_defined')
         lBc%bType = IBSET(lBc%bType,bType_ud)
         lPtr => list%get(fTmp,"Spatial profile file path",1)
         fid = fTmp%open()
         ALLOCATE(lBc%gx(ib%msh(iM)%fa(iFa)%nNo), ptr(ib%msh(iM)%gnNo))
!     I am seting all the nodes to zero just in case a node is not set
         ptr = 0
!     Preparing the pointer array
         IF (iFa .NE. 0) THEN
            DO a=1, ib%msh(iM)%fa(iFa)%nNo
               lBc%gx(a) = 0D0
               Ac = ib%msh(iM)%fa(iFa)%gN(a)
               Ac = ib%msh(iM)%lN(Ac)
               IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2            "detected for BC. Mesh: "//TRIM(ib%msh(iM)%name)//
     3            ", Face: "//TRIM(ib%msh(iM)%fa(iFa)%name)//
     4            ", Node: "//STR(a)//" gN: "//
     5            STR(ib%msh(iM)%fa(iFa)%gN(a))
               ptr(Ac) = a
            END DO
         ELSE
            DO a=1, ib%msh(iM)%gnNo
               lBc%gx(a) = 0D0
               ptr(a) = a
            END DO
         END IF

         DO b=1, nNo
            READ (fid,*) Ac, rtmp
            IF (Ac.GT.ib%msh(iM)%gnNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//fTmp%fname
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry <"//b//" not found in face "//fTmp%fname
            END IF
            lBc%gx(a) = rtmp
         END DO
         CLOSE(fid)
         DEALLOCATE(ptr)
      CASE DEFAULT
         err = TRIM(list%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

      RETURN
      END SUBROUTINE IB_READBC
!####################################################################
!     This routine initializes IB solution and FSILS data structures
      SUBROUTINE IB_INIT(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER a, iM, iFa, iEq, iDmn, iBc, nnz
      REAL(KIND=8) lD(nsd,tnNo)
      CHARACTER(LEN=stdL) sOut
      TYPE(FSILS_commuType) lsCom

      ib%callD = 0D0
      ib%callD(1) = CPUT()

      std = " ================================="
      std = " Initializing IB data structures.."

      lD = 0D0
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     Initialize IB face normals
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            ib%msh(iM)%fa(iFa)%iM = iM
            CALL IB_FACEINI(ib%msh(iM)%fa(iFa))
         END DO
      END DO

!     Initialize IB face BC profile
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBcIB
            iFa = eq(iEq)%bcIB(iBc)%iFa
            iM  = eq(iEq)%bcIB(iBc)%iM
            IF (iFa .NE. 0) THEN
               CALL IB_BCINI(eq(iEq)%bcIB(iBc), ib%msh(iM)%fa(iFa))
               IF (ib%msh(iM)%lShl) THEN
                  CALL SHLBCINI(eq(iEq)%bcIB(iBc), ib%msh(iM)%fa(iFa),
     2               ib%msh(iM))
               END IF
            ELSE
               CALL IB_SBCINI(eq(iEq)%bcIB(iBc), ib%msh(iM))
            END IF
         END DO
      END DO

!     Initialize shell structures (extended IEN and shell normal). This
!     call has to be made before calling LHSA.
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) CALL IB_SHLINI(ib%msh(iM))
      END DO

!     Create sparse matrix data structures
      CALL IB_LHSA(nnz)
      std = "    Non-zeros in LHS matrix (IB): "//nnz

!     To compute IB traces, calculate node/element adjacency for
!     background lumen mesh
      DO iM=1, nMsh
         msh(iM)%iGC = 0
         CALL GETEADJCNCY(msh(iM))
         CALL GETNADJCNCY(msh(iM))
      END DO

!     Calculate node/element adjacency for immersed bodies
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            CALL GETEADJCNCY(ib%msh(iM))
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               CALL GETEADJCNCY(ib%msh(iM)%fa(iFa))
            END DO
         END IF
      END DO

!     Set iblank field for immersed boundaries
      CALL IB_SETIBLANK(lD)

!     Compute the traces of IB nodes on the background mesh
      DO iM=1, ib%nMsh
         CALL IB_INITTRACES(ib%msh(iM), lD)
      END DO

!     Initialize IB communication structure
      CALL IB_SETCOMMU()

!     Identify ghost cells for immersed boundaries
      CALL IB_SETIGHOST()

      ib%callD(1) = CPUT() - ib%callD(1)
      WRITE(sOut,'(F6.2)') ib%callD(1)
      WRITE(sOut,'(A)') "    Immersed boundary setting time: "//
     2   TRIM(sOut)//" sec"
      std = TRIM(sOut)

!     Set Dirichlet BCs
      CALL IB_SETBCDIR(ib%Ao, ib%Yo, ib%Uo)

      std = " IB initialization complete.."
      std = " ================================="

      RETURN
      END SUBROUTINE IB_INIT
!####################################################################
!     Initializing immersed boundary faces
      SUBROUTINE IB_FACEINI(lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER e, a, Ac, g, iM
      REAL(KIND=8) tmp, area, n(nsd), Jac

      REAL(KIND=8), ALLOCATABLE :: sV(:,:)

!     Calculating the center of the face, diameter and its area
      iM = lFa%iM
      IF (ALLOCATED(lFa%nV)) DEALLOCATE(lFa%nV)
      ALLOCATE(lFa%nV(nsd,lFa%nNo))

      ALLOCATE(sV(nsd,ib%tnNo))
      sV   = 0D0
      area = 0D0
      DO e=1, lFa%nEl
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(ib%msh(iM), lFa, e)
         DO g=1, lFa%nG
            CALL GNNIB(lFa, e, g, n)
            Jac  = SQRT(NORM(n))
            area = area + Jac*lFa%w(g)
            DO a=1, lFa%eNoN
               Ac       = lFa%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + n*lFa%N(a,g)*lFa%w(g)
            END DO
         END DO
      END DO
      lFa%area = area
      std = "    Area of face <"//TRIM(lFa%name)//"> is "//STR(area)
      IF (ISZERO(area)) THEN
         IF (cm%mas()) wrn = " <"//TRIM(lFa%name)//"> area is zero"
      END IF

      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         tmp = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(tmp)) THEN
            wrn = " Skipping normal calculation of node "//a//
     2         " in face <"//TRIM(lFa%name)//">"
            lFa%nV(:,a) = 0D0
            lFa%nV(1,a) = 1D0
            CYCLE
         END IF
         lFa%nV(:,a) = sV(:,Ac)/tmp
      END DO

      RETURN
      END SUBROUTINE IB_FACEINI
!--------------------------------------------------------------------
!     Set BC spatial profile on the face
      SUBROUTINE IB_BCINI(lBc, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER iM, iFa, jFa, i, a, b, Ac, j
      REAL(KIND=8) tmp, nV(nsd), center(nsd), maxN

      INTEGER, ALLOCATABLE :: gNodes(:)
      REAL(KIND=8), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = " Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      iM  = lFa%iM
      iFa = lBc%iFa
      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lFa%nNo))

      ALLOCATE(s(ib%tnNo))
      s = 0D0
      IF (BTEST(lBc%bType,bType_flat)) THEN
!     Just a constant value for Flat profile
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = 1D0
         END DO
      ELSE IF (BTEST(lBc%bType,bType_para)) THEN
!     Here is the method that is used for imposing parabolic profile:
!     1- Find the coordinate of the points on the boundary 2- find unit
!     vector from center to each of points on the boundary: ew
!     3- maximize ew(i).e where e is the unit vector from current
!     point to the center 4- Use the point i as the diam here
         center = 0D0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            center(:) = center(:) + ib%x(:,Ac)
         END DO
         center(:) = center(:)/REAL(lFa%nNo,KIND=8)
         ALLOCATE(gNodes(ib%tnNo), sV(nsd,ib%tnNo))
!     gNodes is one if a node located on the boundary (beside iFa)
         gNodes = 0
         DO jFa=1, ib%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ib%msh(iM)%fa(jFa)%nNo
               Ac         = ib%msh(iM)%fa(jFa)%gN(a)
               gNodes(Ac) = 1
            END DO
         END DO
!     "j" is a counter for the number of nodes that are located on the
!     boundary of lFa and sVl contains the list of their coordinates
         j  = 0
         sV = 0D0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (gNodes(Ac) .EQ. 1) THEN
               j = j + 1
               sV(:,j) = ib%x(:,Ac)
            END IF
         END DO
         IF (cm%mas() .AND. j.EQ.0) err = " No perimeter"//
     2      " found for face "//lFa%name
!     sVl will keep the normal unit vector from center to perimeter
         ALLOCATE(sVl(nsd,j))
         DO a=1, j
            sV(:,a)  = sV(:,a) - center
            sVl(:,a) = sV(:,a)/SQRT(NORM(sV(:,a)))
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            nV = ib%x(:,Ac) - center
            maxN = NORM(nV, sVl(:,1))
            i = 1
            DO b=2, j
               tmp = NORM(nV, sVl(:,b))
               IF (tmp .GT. maxN) THEN
                  maxN = tmp
                  i = b
               END IF
            END DO
            s(Ac) = 1D0 - NORM(nV)/NORM(sV(:,i))
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ud)) THEN
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = lBc%gx(a)
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ddep)) THEN
         s = 1D0
      END IF

!     Now correcting the inlet BC for the inlet ring
      IF (BTEST(lBc%bType,bType_zp)) THEN
         DO jFa=1, ib%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ib%msh(iM)%fa(jFa)%nNo
               Ac    = ib%msh(iM)%fa(jFa)%gN(a)
               s(Ac) = 0D0
            END DO
         END DO
      END IF

      DO a=1, lFa%nNo
         Ac        = lFa%gN(a)
         lBc%gx(a) = s(Ac)
      END DO

      RETURN
      END SUBROUTINE IB_BCINI
!--------------------------------------------------------------------
!     Set BC spatial profile on the shell
      SUBROUTINE IB_SBCINI(lBc, lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(mshType), INTENT(IN) :: lM

      INTEGER a, Ac, iFa

      REAL(KIND=8), ALLOCATABLE :: s(:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = " Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lM%nNo))

      ALLOCATE(s(ib%tnNo))
      s = 0D0
      IF (BTEST(lBc%bType,bType_flat)) THEN
!     Just a constant value for Flat profile
         DO a=1, lM%nNo
            Ac    = lM%gN(a)
            s(Ac) = 1D0
         END DO
      ELSE IF (BTEST(lBc%bType,bType_para)) THEN
         err = "Parabolic BC not set up yet for immersed shells"

      ELSE IF (BTEST(lBc%bType,bType_ud)) THEN
         DO a=1, lM%nNo
            Ac    = lM%gN(a)
            s(Ac) = lBc%gx(a)
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ddep)) THEN
         s = 1D0
      END IF

!     Zero perimeter - all faces of the shell are zeroed
      IF (BTEST(lBc%bType,bType_zp)) THEN
         DO iFa=1, lM%nFa
            DO a=1, lM%fa(iFa)%nNo
               Ac    = lM%fa(iFa)%gN(a)
               s(Ac) = 0D0
            END DO
         END DO
      END IF

      DO a=1, lM%nNo
         Ac        = lM%gN(a)
         lBc%gx(a) = s(Ac)
      END DO

      RETURN
      END SUBROUTINE IB_SBCINI
!--------------------------------------------------------------------
!     Initializing shells
      SUBROUTINE IB_SHLINI(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER :: a, b, e, f, g, Ac, Bc, nNo, nEl, eNoN, ep(2,3)
      LOGICAL :: flag
      REAL(KIND=8) :: Jac, area, nV(nsd), tmpR(nsd,nsd-1)

      INTEGER, ALLOCATABLE :: incN(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), sV(:,:)

      IF (lM%eType .EQ. eType_NRB) THEN
         ALLOCATE(lM%eIEN(0,0), lM%sbc(lM%eNoN,lM%nEl))
         lM%sbc = 0
         RETURN
      END IF

      ep   = RESHAPE((/2,3,3,1,1,2/), SHAPE(ep))
      nNo  = lM%nNo
      nEl  = lM%nEl
      eNoN = lM%eNoN
      ALLOCATE(incN(eNoN), lM%eIEN(eNoN,nEl), lM%sbc(eNoN,nEl))

      lM%eIEN = 0
      lM%sbc  = 0
      DO e=1, nEl
         DO a=1, eNoN
            Ac = lM%IEN(ep(1,a),e)
            Bc = lM%IEN(ep(2,a),e)
            DO f=1, nEl
               IF (e .EQ. f) CYCLE
               incN = 0
               DO b=1, eNoN
                  IF (lM%IEN(b,f).EQ.Ac .OR. lM%IEN(b,f).EQ.Bc)
     2               incN(b) = incN(b) + 1
               END DO
               IF (SUM(incN) .EQ. 2) THEN
                  DO b=1, eNoN
                     IF (incN(b) .EQ. 0) THEN
                        lM%eIEN(a,e) = lM%IEN(b,f)
                        EXIT
                     END IF
                  END DO
                  EXIT
               END IF
            END DO

            IF (lM%eIEN(a,e) .EQ. 0) THEN
               lM%sbc(a,e) = IBSET(lM%sbc(a,e), bType_free)
            END IF
         END DO
      END DO
      DEALLOCATE(incN)

!     Compute shell director (normal)
      ALLOCATE(xl(nsd,eNoN), sV(nsd,ib%tnNo), lM%nV(nsd,nNo))
      sV   = 0D0
      area = 0D0
      DO e=1, nEl
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            xl(:,a) = ib%x(:,Ac)
         END DO

         DO g=1, lM%nG
            CALL GNNS(eNoN, lM%Nx(:,:,g), xl, nV, tmpR, tmpR)
            Jac = SQRT(NORM(nV))
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + lM%w(g)*lM%N(a,g)*nV
               area     = area     + lM%w(g)*lM%N(a,g)*Jac
            END DO
         END DO
      END DO

      std = "    Area of the shell surface <"//TRIM(lM%name)//"> is "//
     2   STR(area)

      flag = .TRUE.
      DO a=1, nNo
         Ac  = lM%gN(a)
         Jac = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(Jac)) THEN
            IF (flag) THEN
               wrn = "Skipping normal calculation of node "//a//
     2            " in mesh <"//TRIM(lM%name)//">"
               flag = .FALSE.
            END IF
            lM%nV(:,a) = 0D0
            lM%nV(1,a) = 1D0
            CYCLE
         END IF
         lM%nV(:,a) = sV(:,Ac)/Jac
      END DO

      DEALLOCATE(xl, sV)

      RETURN
      END SUBROUTINE IB_SHLINI
!####################################################################
!     Form the LHS sparse data structures for immersed bodies
      SUBROUTINE IB_LHSA(nnz)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: nnz

      INTEGER i, j, a, b, e, mnnzeic, rowN, colN, iM

      INTEGER, ALLOCATABLE :: uInd(:,:)

      mnnzeic = 10*MAXVAL(ib%msh(:)%eNoN)
 003  mnnzeic = mnnzeic + MAX(5,mnnzeic/5)

      IF (ALLOCATED(uInd)) DEALLOCATE(uInd)
      ALLOCATE (uInd(mnnzeic,ib%tnNo))
      uInd = 0
      DO iM=1, ib%nMsh
!        Treat shell with triangular elements separately
         IF (ib%msh(iM)%lShl .AND. ib%msh(iM)%eType.EQ.eType_TRI) CYCLE
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               rowN = ib%msh(iM)%IEN(a,e)
               DO b=1, ib%msh(iM)%eNoN
                  colN = ib%msh(iM)%IEN(b,e)
                  DO i=1, mnnzeic
                     IF (uInd(i,rowN) .EQ. 0) THEN
                        uInd(i,rowN) = colN
                        EXIT
                     END IF
                     IF (colN .GT. uInd(i,rowN)) CYCLE
                     IF (colN .EQ. uInd(i,rowN)) EXIT
                     IF (uInd(mnnzeic,rowN) .NE. 0) GOTO 003
                     DO j=mnnzeic, i+1, -1
                        uInd(j,rowN) = uInd(j-1,rowN)
                     END DO
                     uInd(i,rowN) = colN
                     EXIT
                  END DO
                  IF (i .GT. mnnzeic) GOTO 003
               END DO
            END DO
         END DO
      END DO

!     Treat shells with triangular elements here
      DO iM=1, ib%nMsh
         IF (.NOT.ib%msh(iM)%lShl) CYCLE
         IF (ib%msh(iM)%eType .EQ. eType_NRB) CYCLE
         DO e=1, ib%msh(iM)%nEl
            DO a=1, 2*ib%msh(iM)%eNoN
               IF (a .LE. ib%msh(iM)%eNoN) THEN
                  rowN = ib%msh(iM)%IEN(a,e)
               ELSE
                  rowN = ib%msh(iM)%eIEN(a-ib%msh(iM)%eNoN,e)
               END IF
               IF (rowN .EQ. 0) CYCLE
               DO b=1, 2*ib%msh(iM)%eNoN
                  IF (b .LE. ib%msh(iM)%eNoN) THEN
                     colN = ib%msh(iM)%IEN(b,e)
                  ELSE
                     colN = ib%msh(iM)%eIEN(b-ib%msh(iM)%eNoN,e)
                  END IF
                  IF (colN .EQ. 0) CYCLE
                  DO i=1, mnnzeic
                     IF (uInd(i,rowN) .EQ. 0) THEN
                        uInd(i,rowN) = colN
                        EXIT
                     END IF
                     IF (colN .GT. uInd(i,rowN)) CYCLE
                     IF (colN .EQ. uInd(i,rowN)) EXIT
                     IF (uInd(mnnzeic,rowN) .NE. 0) GOTO 003
                     DO j=mnnzeic, i+1, -1
                        uInd(j,rowN) = uInd(j-1,rowN)
                     END DO
                     uInd(i,rowN) = colN
                     EXIT
                  END DO
                  IF (i .GT. mnnzeic) GOTO 003
               END DO
            END DO
         END DO
      END DO

!--------------------------------------------------------------------
!     Finding number of non-zeros in colPtr vector
      nnz = 0
      DO rowN=1, ib%tnNo
         IF (uInd(1,rowN) .EQ. 0) THEN
            err = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               nnz = nnz + 1
            END IF
         END DO
      END DO

!--------------------------------------------------------------------
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE (ib%colPtr(nnz), ib%rowPtr(ib%tnNo+1))
      j  = 1
      ib%rowPtr(1) = 1
      DO rowN=1, ib%tnNo
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               ib%colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         ib%rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      END SUBROUTINE IB_LHSA
!####################################################################
!     Set iblank field
!        iblank is set only immersed solids and not set for thin shells
!        iblank(A) = 1   =>   node A is inside the immersed solid
!        iblank(A) = 0   =>   node A is outside the immersed solid
      SUBROUTINE IB_SETIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER :: a, Ac, i, iM, itmp
      LOGICAL :: flag, incL
      REAL(KIND=8) :: minb(nsd), maxb(nsd), xp(nsd), dx

      LOGICAL, ALLOCATABLE :: chck(:)

!     Initialize iblank field
      iblank(:) = 0

!     Return if all the immersed meshes are shells
      flag = .FALSE.
      DO iM=1, ib%nMsh
         IF (.NOT.ib%msh(iM)%lShl) THEN
            flag = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT.flag) RETURN

!     Fit a bounding box around IB and probe only those nodes lying
!     inside the box
      dx = TINY(dx)
      DO iM=1, ib%nMsh
         IF (dx .LT. ib%msh(iM)%dx) THEN
            dx = ib%msh(iM)%dx
         END IF
      END DO

      DO i=1, nsd
         minb(i) = MINVAL(ib%x(i,:) + ib%Uo(i,:)) - dx
         maxb(i) = MAXVAL(ib%x(i,:) + ib%Uo(i,:)) + dx
      END DO

      ALLOCATE(chck(tnNo))
      chck = .FALSE.
      DO Ac=1, tnNo
         xp(:) = x(:,Ac) + lD(:,Ac)
         itmp = 0
         DO i=1, nsd
            IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i)) THEN
               itmp = itmp + 1
            END IF
         END DO
         IF (itmp .EQ. nsd) chck(Ac) = .TRUE.
      END DO

      DO iM=1, nMsh
!        Probe each node if it is inside or outside the immersed body
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            IF (.NOT.chck(Ac)) CYCLE

            xp(:) = x(:,Ac) + lD(:,Ac)
            incL  = .FALSE.
            CALL IB_CHECKINOUT(xp, incL)
            IF (incL) iblank(Ac) = 1
            chck(Ac) = .FALSE.
         END DO
      END DO

      RETURN
      END SUBROUTINE IB_SETIBLANK
!--------------------------------------------------------------------
!     Update iblank field around the neighborhood of previous iblank. In
!     case iblank is 0 entirely on the current process, redo iblank
!     initialization to look for any fresh iblank node
      SUBROUTINE IB_UPDATEIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL flag
      INTEGER a, b, i, Ac, Bc, iswp, iM, nNo
      REAL(KIND=8) xp(nsd)

      INTEGER, ALLOCATABLE :: incNd(:), ptr(:)

!     Just in case a fresh iblank appears in the current process
      IF (SUM(iblank) .EQ. 0) CALL IB_SETIBLANK(lD)

      ALLOCATE(incNd(tnNo))
      incNd(:) = iblank(:)

!     Do multiple sweeps to get the neighborhood of iblank field. Probe
!     is performed only on these neighborhood nodes.
      iswp = 0
      DO WHILE (iswp .LE. 2)
         DO iM=1, nMsh
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               IF (incNd(Ac) .EQ. 1) THEN
                  DO i=msh(iM)%nAdj%prow(a), msh(iM)%nAdj%prow(a+1)-1
                     b = msh(iM)%nAdj%pcol(i)
                     Bc = msh(iM)%gN(b)
                     incNd(Bc) = 1
                  END DO
               END IF
            END DO
         END DO
         iswp = iswp + 1
      END DO

      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) a = a + 1
      END DO
      nNo = a

!     Transfer all the nodes to be probed to a local pointer array
      ALLOCATE(ptr(nNo))
      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) THEN
            a = a + 1
            ptr(a) = Ac
         END IF
      END DO

!     Probe each node if it is inside or outside the immersed body
      iblank = 0
      DO a=1, nNo
         Ac = ptr(a)
         xp = x(:,Ac) + lD(:,Ac)
         flag  = .FALSE.
         CALL IB_CHECKINOUT(xp, flag)
         IF (flag) iblank(Ac) = 1
      END DO
      DEALLOCATE(incNd, ptr)

      RETURN
      END SUBROUTINE IB_UPDATEIBLANK
!--------------------------------------------------------------------
!     Checks if a probe lies inside or outside an immersed boundary
      SUBROUTINE IB_CHECKINOUT(xp, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: xp(nsd)
      LOGICAL, INTENT(INOUT) :: flag

      INTEGER :: a, e, Ac, Ec, iM, iFa, jM, jFa
      REAL(KIND=8) :: dS, minS, Jac, nV(nsd), xb(nsd), dotP

!     Find the closest immersed face centroid from the probe
      minS = HUGE(minS)
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            DO e=1, ib%msh(iM)%fa(iFa)%nEl
               xb = 0D0
               DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                  Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                  xb = xb + ib%x(:,Ac) + ib%Uo(:,Ac)
               END DO
               xb = xb / REAL(ib%msh(iM)%fa(iFa)%eNoN, KIND=8)
               dS = SQRT( SUM( (xp(:)-xb(:))**2 ) )
               IF (dS .LT. minS) THEN
                  minS = dS
                  Ec   = e
                  jM   = iM
                  jFa  = iFa
               END IF
            END DO
         END DO
      END DO

!     Compute the normal of the face element
      CALL GNNIB(ib%msh(jM)%fa(jFa), Ec, 1, nV)
      Jac = SQRT(NORM(nV))
      nV  = nV(:)/Jac

!     Check the sign of op.n
      xb = 0D0
      DO a=1, ib%msh(jM)%fa(jFa)%eNoN
         Ac = ib%msh(jM)%fa(jFa)%IEN(a,Ec)
         xb = xb + ib%x(:,Ac) + ib%Uo(:,Ac)
      END DO
      xb   = xb / REAL(ib%msh(jM)%fa(jFa)%eNoN, KIND=8)
      dotP = NORM(xp-xb, nV)

      IF (dotP .LT. -1D-9) THEN
!        probe lies inside IB
         flag = .TRUE.
      ELSE IF (dotP .GT. 1D-9) THEN
!        probe lies outside IB
         flag = .FALSE.
      ELSE
!        if the probe is along the tangent, perform sign check with one
!        of the vertices of the closest node instead of the centroid
         Ac   = ib%msh(jM)%fa(jFa)%IEN(1,Ec)
         xb   = ib%x(:,Ac) + ib%Uo(:,Ac)
         dotP = NORM(xp-xb,nV)
         IF (ABS(dotP) .LT. 1D-9) THEN
            flag = .TRUE.
         ELSE
            flag = .FALSE.
         END IF
      END IF

      RETURN
      END SUBROUTINE IB_CHECKINOUT
!####################################################################
!     Initialize traces of IB mesh integration points on the
!     background mesh elements and tag those traces as ghost elements.
      SUBROUTINE IB_INITTRACES(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER e, g, i, j, Ec, iM, iFa

      INTEGER, ALLOCATABLE :: incEl(:), ePtr(:,:,:)

!     For shells and IFEM method, all the mesh elements are included for
!     trace search. For solids, only boundary/face nodes are included
!     for trace search.
      IF (lM%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
         ALLOCATE(incEl(lM%nEl), ePtr(2,lM%nG,lM%nEl))
         incEl = 1
         ePtr  = 0
         CALL IB_FINDMSHTRACES(lM, lD, incEl, ePtr)
         DEALLOCATE(incEl, ePtr)
      ELSE
         DO iFa=1, lM%nFa
            ALLOCATE(incEl(lM%fa(iFa)%nEl),
     2         ePtr(2,lM%fa(iFa)%nG,lM%fa(iFa)%nEl))
            incEl = 1
            ePtr  = 0
            CALL IB_FINDFACETRACES(lM%fa(iFa), lD, incEl, ePtr)
            DEALLOCATE(incEl, ePtr)
         END DO
      END IF

      RETURN
      END SUBROUTINE IB_INITTRACES
!--------------------------------------------------------------------
!     Update IB tracers
      SUBROUTINE IB_UPDATE(Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER :: a, iM, iFa
      REAL(KIND=8), ALLOCATABLE :: lD(:,:)

      ALLOCATE(lD(nsd,tnNo))
      lD = 0D0
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     As the IB may have moved update iblank field
      CALL IB_UPDATEIBLANK(lD)

!     Update IB traces on background fluid mesh and update ghost cells
      DO iM=1, nMsh
         msh(iM)%iGC(:) = 0
      END DO
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            CALL IB_UPDATETRACESMSH(ib%msh(iM), lD)
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               CALL IB_UPDATETRACESFACE(ib%msh(iM)%fa(iFa), lD)
            END DO
         END IF
      END DO

!     Reset IB communicator
      CALL IB_SETCOMMU()

!     Now that we have updated iblank field and IB traces, time to
!     update ghost cells
      CALL IB_SETIGHOST()

      RETURN
      END SUBROUTINE IB_UPDATE
!--------------------------------------------------------------------
!     Update traces of IB mesh integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_UPDATETRACESMSH(lM, lD)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: ipass
      INTEGER :: a, e, g, i, j, Ac, Ec, iM, iswp, ne
      REAL(KIND=8) :: xp(nsd), xi(nsd)

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER, ALLOCATABLE :: incEl(:), eList(:), ePtr(:,:,:),
     2   gPtr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:)

!     It is not necessary to find traces for all nodes. Instead, use
!     previously determined traces and update them in their neighborhood
!     elements. Likewise, neighborhood nodes of IB are also added for
!     trace search.
      ALLOCATE(ePtr(2,lM%nG,lM%nEl), incEl(lM%nEl))
      ePtr  = 0
      incEl = 0
      DO i=1, lM%trc%n
         e = lM%trc%gE(1,i)
         g = lM%trc%gE(2,i)
         incEl(e)    = 1
         ePtr(:,g,e) = lM%trc%ptr(:,i)
      END DO

!     Ignore if there are no tracers earlier. However, some elements
!     could move to the neighboring process and create fresh tracer
!     pointers. These will be tracked later.
      ipass = .TRUE.
      IF (SUM(incEl) .EQ. 0) ipass = .FALSE.

!     For elements that already have tracer pointers, search in the
!     neighborhood of previous element tracers
      ALLOCATE(gPtr(2,lM%nG,lM%nEl), ichk(lM%nEl), xl(nsd,lM%eNoN))
      gPtr = 0
      iswp = 1
      ichk = .FALSE.
      DO WHILE (iswp.LE.3 .AND. ipass)
         DO e=1, lM%nEl
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO

            DO g=1, lM%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. ichk(e)) CYCLE

               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               xp = 0D0
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
               END IF

               gPtr(1,g,e) = Ec
               gPtr(2,g,e) = iM
               msh(iM)%iGC(Ec) = 1

!              Initialize the neighbors using the current trace
               DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
                  j = lM%eAdj%pcol(i)
                  IF (j.EQ.e .OR. incEl(j).EQ.0 .OR. ichk(j)) CYCLE
                  ePtr(1,g,j) = Ec
                  ePtr(2,g,j) = iM
               END DO
            END DO
            ichk(e) = .TRUE.
         END DO
         iswp = iswp + 1
      END DO
      DEALLOCATE(ePtr)

!     Reset incEl for the elements whose traces are newly determined
      incEl = 1
      DO e=1, lM%nEl
         i = 0
         DO g=1, lM%nG
            IF (gPtr(1,g,e) .NE. 0) i = i + 1
         END DO
         IF (i .EQ. lM%nG) incEl(e) = 0
      END DO

!     Perform a general search for all the remaining nodes including
!     the ones that may have migrated into other process domains
      CALL IB_FINDMSHTRACES(lM, lD, incEl, gPtr)
      DEALLOCATE(incEl, gPtr)

      RETURN
      END SUBROUTINE IB_UPDATETRACESMSH
!--------------------------------------------------------------------
!     Update traces of IB face integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_UPDATETRACESFACE(lFa, lD)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE

      TYPE(faceType), INTENT(INOUT) :: lFa
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: ipass
      INTEGER :: a, e, g, i, j, Ac, Ec, iM, iswp, ne
      REAL(KIND=8) :: xp(nsd), xi(nsd)

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER, ALLOCATABLE :: incEl(:), eList(:), ePtr(:,:,:),
     2   gPtr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:)

!     It is not necessary to find traces for all nodes. Instead, use
!     previously determined traces and update them in their neighborhood
!     elements. Likewise, neighborhood nodes of IB are also added for
!     trace search.
      ALLOCATE(ePtr(2,lFa%nG,lFa%nEl), incEl(lFa%nEl))
      ePtr  = 0
      incEl = 0
      DO i=1, lFa%trc%n
         e = lFa%trc%gE(1,i)
         g = lFa%trc%gE(2,i)
         incEl(e)    = 1
         ePtr(:,g,e) = lFa%trc%ptr(:,i)
      END DO

!     Ignore if there are no tracers earlier. However, some elements
!     could move to the neighboring process and create fresh tracer
!     pointers. These will be tracked later.
      ipass = .TRUE.
      IF (SUM(incEl) .EQ. 0) ipass = .FALSE.

!     For elements that already have tracer pointers, search in the
!     neighborhood of previous element tracers
      ALLOCATE(gPtr(2,lFa%nG,lFa%nEl), ichk(lFa%nEl), xl(nsd,lFa%eNoN))
      gPtr = 0
      iswp = 0
      ichk = .FALSE.
      DO WHILE (iswp.LE.1 .AND. ipass)
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,Ac) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO

            DO g=1, lFa%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. ichk(e)) CYCLE

               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               xp = 0D0
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
               END IF

               gPtr(1,g,e) = Ec
               gPtr(2,g,e) = iM
               msh(iM)%iGC(Ec) = 1

!              Initialize the neighbors using the current trace
               DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
                  j = lFa%eAdj%pcol(i)
                  IF (j.EQ.e .OR. incEl(j).EQ.0 .OR. ichk(j)) CYCLE
                  ePtr(1,g,j) = Ec
                  ePtr(2,g,j) = iM
               END DO
            END DO
            ichk(e) = .TRUE.
         END DO
         iswp = iswp + 1
      END DO
      DEALLOCATE(ePtr)

!     Reset incEl for the elements whose traces are newly determined
      incEl = 1
      DO e=1, lFa%nEl
         i = 0
         DO g=1, lFa%nG
            IF (gPtr(1,g,e) .NE. 0) i = i + 1
         END DO
         IF (i .EQ. lFa%nG) incEl(e) = 0
      END DO

!     Perform a general search for all the remaining nodes including
!     the ones that may have migrated into other process domains
      CALL IB_FINDFACETRACES(lFa, lD, incEl, gPtr)
      DEALLOCATE(incEl, gPtr)

      RETURN
      END SUBROUTINE IB_UPDATETRACESFACE
!--------------------------------------------------------------------
!     Find traces of IB integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_FINDMSHTRACES(lM, lD, srchEl, ePtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: srchEl(lM%nEl)
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER, INTENT(INOUT) :: ePtr(2,lM%nG,lM%nEl)

      INTEGER :: a, b, e, g, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=8) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd)

      TYPE(queueType) :: probeElQ
      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER, ALLOCATABLE :: incEl(:), rootEl(:), eList(:), sCount(:),
     2   disps(:), masEList(:), ptr(:), gptr(:), tmpI(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), xpL(:,:)

!     Create a list of all the mesh nodes
      ALLOCATE(xpL(nsd,lM%nNo))
      xpL = 0D0
      DO a=1, lM%nNo
         Ac = lM%gN(a)
         xpL(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - lM%dx
      maxb(:) = maxb(:) + lM%dx

!     Loop over all possible background mesh and find traces of each IB
!     integration point onto the elements of the background mesh.
      DO iM=1, nMsh
!        Identify the IB elements which are owned by the current
!        process. Search for traces is performed only these elements
         ALLOCATE(incEl(lM%nEl))
         incEl(:) = srchEl(:)
         DO e=1, lM%nEl
            b = 0
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               Ac = lM%lN(Ac)
               xp = xpL(:,Ac)
               j = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               j = j + 1
               END DO
               IF (j .EQ. nsd) b = b + 1
            END DO
            IF (b .LT. lM%eNoN) incEl(e) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           Reset eList if the there is no overlap with iblank field
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b  = b + iblank(Ac)
            END DO
            IF (b .EQ. 0) eList(e) = 0
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            CALL IB_PARTMSH(lM, msh(iM), lD, eList, incEl)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incEl) .EQ. 0) THEN
            DEALLOCATE(incEl, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IB nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lM%nEl), rootEl(lM%nEl))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed element using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed element and form a queue for
!        subsequent search. The root element for the seed element is
!        assigned as a seed search element for the queued neighboring
!        elements
         ALLOCATE(xl(nsd,lM%eNoN))
         Ec = 0
         E_LOOP: DO e=1, lM%nEl
            ichk(e) = .TRUE.
            IF (incEl(e) .EQ. 0) CYCLE
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lM%nG
               xp = 0D0
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               ichk(e) = .FALSE.
               EXIT E_LOOP
            END DO
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
            CYCLE
         END IF

         ePtr(1,g,e) = Ec
         ePtr(2,g,e) = iM
         rootEl(e)   = Ec
         DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
            j = lM%eAdj%pcol(i)
            IF (incEl(j) .EQ. 1) THEN
               CALL ENQUEUE(probeElQ, j)
               rootEl(j) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeElQ, e))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(e)) CYCLE
            ichk(e) = .TRUE.

            DO g=1, lM%nG
               DO a=1, lM%eNoN
                  Ac = lM%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
               END DO

               xp = 0D0
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO

               Ec = rootEl(e)
               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!              If a trace is not found, then include the neighborhood
!              elements for search. Repeat this twice. If not found yet,
!              search using master element list. If master list search
!              also fails, continue
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) THEN
                     CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2                  Ec, xi)
                     IF (Ec .EQ. 0) CYCLE
                  END IF
               END IF

!              At this point, the trace is definitely found. Save it.
               ePtr(1,g,e) = Ec
               ePtr(2,g,e) = iM

!              Use the trace as a seed search element for the neigboring
!              IB elements in the following trace search
               rootEl(e)  = Ec
               DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
                  j = lM%eAdj%pcol(i)
                  IF (incEl(j).EQ.1 .AND. .NOT.ichk(j)) THEN
                     CALL ENQUEUE(probeElQ, j)
                     rootEl(j) = Ec
                  END IF
               END DO
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeElQ)

!        Perform a brute search on any missed elements
         DO e=1, lM%nEl
            IF (.NOT.ichk(e) .AND. incEl(e).EQ.1) THEN
               DO g=1, lM%nG
                  DO a=1, lM%eNoN
                     Ac = lM%IEN(a,e)
                     xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
                  END DO

                  xp = 0D0
                  DO a=1, lM%eNoN
                     xp = xp + xl(:,a)*lM%N(a,g)
                  END DO

                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END DO
            END IF
         END DO

         DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
      END DO

      DEALLOCATE(xpL)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lM%trc)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%ptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)  = e
               lM%trc%gE(2,i)  = g
               lM%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Check if all traces are found and return if yes
      i = cm%reduce(lM%trc%n)
      j = lM%nEl*lM%nG
      IF (i .GE. j) RETURN

!     If all traces are not found, add a fool-proof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lM%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      j = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(ptr(4*lM%trc%n), gptr(4*j))
      DO i=1, lM%trc%n
         ptr(4*i-3) = lM%trc%gE(1,i)
         ptr(4*i-2) = lM%trc%gE(2,i)
         ptr(4*i-1) = lM%trc%ptr(1,i)
         ptr(4*i)   = lM%trc%ptr(2,i)
      END DO

      i = 4*lM%trc%n
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lM%nG,lM%nEl), ptr(lM%nEl))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            e  = gptr(4*i-3)
            g  = gptr(4*i-2)
            Ec = gptr(4*i-1)
            iM = gptr(4*i)
            tmpI(1,g,e) = Ec
            tmpI(2,g,e) = iM
         END DO

         DO e=1, lM%nEl
            DO g=1, lM%nG
               Ec = tmpI(1,g,e)
               iM = tmpI(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(e) = 1
            END DO
         END DO

         ne = SUM(ptr)
         ALLOCATE(incEl(ne))
         i = 0
         DO e=1, lM%nEl
            IF (ptr(e) .EQ. 1) THEN
               i = i + 1
               incEl(i) = e
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF

!     Share the element list to all processes
      CALL cm%bcast(ne)
      IF (cm%slv()) ALLOCATE(incEl(ne))
      CALL cm%bcast(incEl)

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(xl(nsd,lM%eNoN))
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        IB bounding box position
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, ne
            e = incEl(i)
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lM%nG
               xp = 0D0
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .NE. 0) THEN
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END IF
            END DO
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(xl)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lM%trc)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%ptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)  = e
               lM%trc%gE(2,i)  = g
               lM%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Abort simulation if all traces are still not found
      i = cm%reduce(lM%trc%n)
      j = lM%nEl * lM%nG
      IF (i .LT. j) CALL DEBUGIBMSHTRC(lM)

      RETURN
      END SUBROUTINE IB_FINDMSHTRACES
!--------------------------------------------------------------------
!     Find traces of IB integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_FINDFACETRACES(lFa, lD, srchEl, ePtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(INOUT) :: lFa
      INTEGER, INTENT(IN) :: srchEl(lFa%nEl)
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER, INTENT(INOUT) :: ePtr(2,lFa%nG,lFa%nEl)

      LOGICAL lShl
      INTEGER :: a, b, e, g, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=8) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd)

      TYPE(queueType) :: probeElQ
      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER, ALLOCATABLE :: incEl(:), rootEl(:), eList(:), sCount(:),
     2   disps(:), masEList(:), ptr(:), gptr(:), tmpI(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), xpL(:,:)

!     Create a list of all the face nodes
      ALLOCATE(xpL(nsd,lFa%nNo))
      xpL = 0D0
      DO a=1, lFa%nNo
         Ac = lFa%gN(a)
         xpL(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - ib%msh(lFa%iM)%dx
      maxb(:) = maxb(:) + ib%msh(lFa%iM)%dx

!     Determine if there are any shell surfaces
      lShl = .FALSE.
      DO i=1, ib%nMsh
         IF (ib%msh(i)%lShl) THEN
            lShl = .TRUE.
            EXIT
         END IF
      END DO

!     Loop over all possible background mesh and find traces of each IB
!     integration point onto the elements of the background mesh.
      DO iM=1, nMsh
!        Identify the IB elements which are owned by the current
!        process. Search for traces is performed only these elements
         ALLOCATE(incEl(lFa%nEl))
         incEl(:) = srchEl(:)
         DO e=1, lFa%nEl
            b = 0
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               Ac = lFa%lN(Ac)
               xp = xpL(:,Ac)
               j = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               j = j + 1
               END DO
               IF (j .EQ. nsd) b = b + 1
            END DO
            IF (b .LT. lFa%eNoN) incEl(e) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           For solids, reset eList if the element has no overlap with
!           IB based on iblank
            IF (.NOT.lShl) THEN
               b = 0
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  b  = b + iblank(Ac)
               END DO
               IF (b .EQ. 0) eList(e) = 0
            END IF
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            CALL IB_PARTFACE(lFa, msh(iM), lD, eList, incEl)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incEl) .EQ. 0) THEN
            DEALLOCATE(incEl, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IB nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lFa%nEl), rootEl(lFa%nEl))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed element using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed element and form a queue for
!        subsequent search. The root element for the seed element is
!        assigned as a seed search element for the queued neighboring
!        elements
         ALLOCATE(xl(nsd,lFa%eNoN))
         Ec = 0
         E_LOOP: DO e=1, lFa%nEl
            ichk(e) = .TRUE.
            IF (incEl(e) .EQ. 0) CYCLE
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lFa%nG
               xp = 0D0
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               ichk(e) = .FALSE.
               EXIT E_LOOP
            END DO
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
            CYCLE
         END IF

         ePtr(1,g,e) = Ec
         ePtr(2,g,e) = iM
         rootEl(e)   = Ec
         DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
            j = lFa%eAdj%pcol(i)
            IF (incEl(j) .EQ. 1) THEN
               CALL ENQUEUE(probeElQ, j)
               rootEl(j) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeElQ, e))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(e)) CYCLE
            ichk(e) = .TRUE.

            DO g=1, lFa%nG
               DO a=1, lFa%eNoN
                  Ac = lFa%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
               END DO

               xp = 0D0
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO

               Ec = rootEl(e)
               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!              If a trace is not found, then include the neighborhood
!              elements for search. Repeat this twice. If not found yet,
!              search using master element list. If master list search
!              also fails, continue
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) THEN
                     CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2                  Ec, xi)
                     IF (Ec .EQ. 0) CYCLE
                  END IF
               END IF

!              At this point, the trace is definitely found. Save it.
               ePtr(1,g,e) = Ec
               ePtr(2,g,e) = iM

!              Use the trace as a seed search element for the neigboring
!              IB elements in the following trace search
               rootEl(e)  = Ec
               DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
                  j = lFa%eAdj%pcol(i)
                  IF (incEl(j) .EQ. 1) THEN
                     CALL ENQUEUE(probeElQ, j)
                     rootEl(j) = Ec
                  END IF
               END DO
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeElQ)

!        Perform a brute search on any missed elements
         DO e=1, lFa%nEl
            IF (.NOT.ichk(e) .AND. incEl(e).EQ.1) THEN
               DO g=1, lFa%nG
                  DO a=1, lFa%eNoN
                     Ac = lFa%IEN(a,e)
                     xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
                  END DO

                  xp = 0D0
                  DO a=1, lFa%eNoN
                     xp = xp + xl(:,a)*lFa%N(a,g)
                  END DO

                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END DO
            END IF
         END DO

         DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
      END DO

      DEALLOCATE(xpL)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lFa%trc)
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lFa%trc%n = i
      ALLOCATE(lFa%trc%gE(2,i), lFa%trc%ptr(2,i))
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lFa%trc%gE(1,i)  = e
               lFa%trc%gE(2,i)  = g
               lFa%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Check if all traces are found and return if yes
      i = cm%reduce(lFa%trc%n)
      j = lFa%nEl*lFa%nG
      IF (i .GE. j) RETURN

!     If all traces are not found, add a fool-proof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lFa%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      j = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(ptr(4*lFa%trc%n), gptr(4*j))
      DO i=1, lFa%trc%n
         ptr(4*i-3) = lFa%trc%gE(1,i)
         ptr(4*i-2) = lFa%trc%gE(2,i)
         ptr(4*i-1) = lFa%trc%ptr(1,i)
         ptr(4*i)   = lFa%trc%ptr(2,i)
      END DO

      i = 4*lFa%trc%n
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lFa%nG,lFa%nEl), ptr(lFa%nEl))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            e  = gptr(4*i-3)
            g  = gptr(4*i-2)
            Ec = gptr(4*i-1)
            iM = gptr(4*i)
            tmpI(1,g,e) = Ec
            tmpI(2,g,e) = iM
         END DO

         DO e=1, lFa%nEl
            DO g=1, lFa%nG
               Ec = tmpI(1,g,e)
               iM = tmpI(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(e) = 1
            END DO
         END DO

         ne = SUM(ptr)
         ALLOCATE(incEl(ne))
         i = 0
         DO e=1, lFa%nEl
            IF (ptr(e) .EQ. 1) THEN
               i = i + 1
               incEl(i) = e
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF

!     Share the element list to all processes
      CALL cm%bcast(ne)
      IF (cm%slv()) ALLOCATE(incEl(ne))
      CALL cm%bcast(incEl)

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(xl(nsd,lFa%eNoN))
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        IB bounding box position
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, ne
            e = incEl(i)
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lFa%nG
               xp = 0D0
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .NE. 0) THEN
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END IF
            END DO
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(xl)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lFa%trc)
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lFa%trc%n = i
      ALLOCATE(lFa%trc%gE(2,i), lFa%trc%ptr(2,i))
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lFa%trc%gE(1,i)  = e
               lFa%trc%gE(2,i)  = g
               lFa%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Abort simulation if all traces are still not found
      i = cm%reduce(lFa%trc%n)
      j = lFa%nEl * lFa%nG
      IF (i .LT. j) CALL DEBUGIBFATRC(lFa)

      RETURN
      END SUBROUTINE IB_FINDFACETRACES
!--------------------------------------------------------------------
      SUBROUTINE IB_PARTMSH(lM, gM, lD, eList, incEl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM, gM
      INTEGER, INTENT(IN) :: eList(gM%nEl)
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER, INTENT(INOUT) :: incEl(lM%nEl)

      INTEGER :: itr, a, b, e, Ac, Bc, ierr
      REAL(KIND=8) :: f, tol, dS, minS, xp(nsd), xb(nsd)

      INTEGER, ALLOCATABLE :: part(:), tmpI(:), gN(:)

      ALLOCATE(gN(tnNo))
      gN = 0
      DO e=1, gM%nEl
         IF (eList(e) .EQ. 0) CYCLE
         DO a=1, gM%eNoN
            Ac = gM%IEN(a,e)
            gN(Ac) = 1
         END DO
      END DO

      itr = 0
      f   = 0.05D0
      ALLOCATE(part(lM%nNo), tmpI(lM%nNo))
 001  part = 0
      tmpI = 0
      itr  = itr + 1
      f    = 2.0D0*f
      tol  = (1.0D0 + f)*lM%dx
      DO a=1, lM%nNo
         IF (part(a) .NE. 0) CYCLE
         Ac   = lM%gN(a)
         xp   = ib%x(:,Ac) + ib%Uo(:,Ac)
         minS = HUGE(minS)
         DO b=1, gM%nNo
            Bc = gM%gN(b)
            IF (gN(Bc) .EQ. 0) CYCLE
            xb = x(:,Bc) + lD(:,Bc)
            dS = SQRT( SUM( (xp(:)-xb(:))**2 ))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            part(a) = cm%tF()
         END IF
      END DO

      CALL MPI_ALLREDUCE(part, tmpI, lM%nNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      b = 0
      DO a=1, lM%nNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. lM%nNo) THEN
         wrn = "Found only "//STR(b)//" nodes in pass "//STR(itr)//
     2      " out of "//STR(lM%nNo)//" nodes"
         IF (itr .GT. 5) err = "Could not distribute all nodes in "//
     2      STR(itr)//" passes. Try changing mesh edge size."
         GOTO 001
      END IF

!     While all the IB nodes are partitioned, the ones which do not lie
!     within the current process are ignored
      DO e=1, lM%nEl
         b = 0
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            IF (part(Ac) .NE. cm%tF()) b = b + 1
         END DO
         IF (b .EQ. lM%eNoN) incEl(e) = 0
      END DO

      DEALLOCATE(gN, part, tmpI)

      RETURN
      END SUBROUTINE IB_PARTMSH
!--------------------------------------------------------------------
      SUBROUTINE IB_PARTFACE(lFa, gM, lD, eList, incEl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      TYPE(mshType), INTENT(IN) :: gM
      INTEGER, INTENT(IN) :: eList(gM%nEl)
      REAL(KIND=8), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER, INTENT(INOUT) :: incEl(lFa%nEl)

      INTEGER :: itr, a, b, e, Ac, Bc, ierr
      REAL(KIND=8) :: f, tol, dS, minS, xp(nsd), xb(nsd)

      INTEGER, ALLOCATABLE :: part(:), tmpI(:), gN(:)

      ALLOCATE(gN(tnNo))
      gN = 0
      DO e=1, gM%nEl
         IF (eList(e) .EQ. 0) CYCLE
         DO a=1, gM%eNoN
            Ac = gM%IEN(a,e)
            gN(Ac) = 1
         END DO
      END DO

      itr = 0
      f   = 0.05D0
      ALLOCATE(part(lFa%nNo), tmpI(lFa%nNo))
 001  part = 0
      tmpI = 0
      itr  = itr + 1
      f    = 2.0D0*f
      tol  = (1.0D0 + f)*ib%msh(lFa%iM)%dx
      DO a=1, lFa%nNo
         IF (part(a) .NE. 0) CYCLE
         Ac   = lFa%gN(a)
         xp   = ib%x(:,Ac) + ib%Uo(:,Ac)
         minS = HUGE(minS)
         DO b=1, gM%nNo
            Bc = gM%gN(b)
            IF (gN(Bc) .EQ. 0) CYCLE
            xb = x(:,Bc) + lD(:,Bc)
            dS = SQRT( SUM( (xp(:)-xb(:))**2 ))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            part(a) = cm%tF()
         END IF
      END DO

      CALL MPI_ALLREDUCE(part, tmpI, lFa%nNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      b = 0
      DO a=1, lFa%nNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. lFa%nNo) THEN
         wrn = "Found only "//STR(b)//" nodes in pass "//STR(itr)//
     2      " out of "//STR(lFa%nNo)//" nodes"
         IF (itr .GT. 5) err = "Could not distribute all nodes in "//
     2      STR(itr)//" passes. Try changing mesh edge size."
         GOTO 001
      END IF

!     While all the IB nodes are partitioned, the ones which do not lie
!     within the current process are ignored
      DO e=1, lFa%nEl
         b = 0
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            IF (part(Ac) .NE. cm%tF()) b = b + 1
         END DO
         IF (b .EQ. lFa%eNoN) incEl(e) = 0
      END DO

      DEALLOCATE(gN, part, tmpI)

      RETURN
      END SUBROUTINE IB_PARTFACE
!--------------------------------------------------------------------
      SUBROUTINE IB_FPSRCH(xp, lM, lD, ne, eSrch, itMax, Ec, xi)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: xp(nsd), lD(nsd,tnNo)
      INTEGER, INTENT(IN)  :: ne, eSrch(ne), itMax
      INTEGER, INTENT(OUT) :: Ec
      REAL(KIND=8), INTENT(OUT) :: xi(nsd)
      TYPE(mshType), INTENT(IN) :: lM

      LOGICAL flag
      INTEGER :: a, e, i, El, En, nEl, nEn, iter, is, ie

      LOGICAL, ALLOCATABLE :: eChk(:)
      INTEGER, ALLOCATABLE :: eList(:), tmpI(:)

      nEn = ne
      ALLOCATE(eList(nEn), eChk(lM%nEl))
      eList = eSrch
      eChk  = .FALSE.
      iter  = 0

      DO WHILE (iter .LE. itMax)
         ALLOCATE(tmpI(lM%nEl))
         tmpI = 0
         nEl  = nEn
         nEn  = 0
         DO e=1, nEl
            El = eList(e)
            eChk(El) = .TRUE.
         END DO
         DO e=1, nEl
            El = eList(e)
            is = lM%eAdj%prow(El)
            ie = lM%eAdj%prow(El+1) - 1
            DO i=is, ie
               En = lM%eAdj%pcol(i)
               IF (Ec.EQ.En .OR. eChk(En)) CYCLE
               flag = .TRUE.
               DO a=1, nEn
                  IF (En .EQ. tmpI(a)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  nEn = nEn + 1
                  tmpI(nEn) = En
               END IF
            END DO
         END DO
         DEALLOCATE(eList)
         ALLOCATE(eList(nEn))
         eList(:) = tmpI(1:nEn)
         CALL FINDE(xp, lM, x, lD, tnNo, nEn, eList, Ec, xi)
         DEALLOCATE(tmpI)
         IF (Ec .NE. 0) THEN
            DEALLOCATE(eChk, eList)
            RETURN
         ELSE
            iter = iter + 1
         END IF
      END DO

      DEALLOCATE(eChk, eList)

      RETURN
      END SUBROUTINE IB_FPSRCH
!####################################################################
!     Communication structure for IB is initialized here. Here we
!     create a list of nodal traces on master that are local to other
!     processes. The master then gathers all the data, projects flow
!     variables (acceleration, velocity and pressure) and broadcasts to
!     all the processes.
      SUBROUTINE IB_SETCOMMU()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER i, a, e, n, Ac, iM, iFa, ierr, tag, sReq

      INTEGER, ALLOCATABLE :: incNd(:), ptr(:), rA(:,:), rReq(:)

!     Free memory if already allocated
      CALL DESTROY(ib%cm)

!     Map trace pointers local to a process into a nodal vector. Note,
!     however, that the traces point to integration point of an element.
!     Therefore, we set all the nodes of an element that get contribution
!     from a valid trace.
      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            DO i=1, ib%msh(iM)%trc%n
               e = ib%msh(iM)%trc%gE(1,i)
               DO a=1, ib%msh(iM)%eNoN
                  Ac = ib%msh(iM)%IEN(a,e)
                  incNd(Ac) = 1
               END DO
            END DO
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               DO i=1, ib%msh(iM)%fa(iFa)%trc%n
                  e = ib%msh(iM)%fa(iFa)%trc%gE(1,i)
                  DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                     Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                     incNd(Ac) = 1
                  END DO
               END DO
            END DO
         END IF
      END DO

!     All the included IB nodes are now mapped to a local vector
      n = SUM(incNd)
      ALLOCATE(ptr(n))
      i = 0
      DO a=1, ib%tnNo
         IF (incNd(a) .NE. 0) THEN
            i = i + 1
            ptr(i) = a
         END IF
      END DO

!     Set IB comm data structures for sequential run
      IF (cm%seq()) THEN
         ALLOCATE(ib%cm%n(1), ib%cm%gE(n))
         ib%cm%n(1)  = n
         ib%cm%gE(:) = ptr
         DEALLOCATE(incNd, ptr)
         RETURN
      END IF

!     Gather no of included nodes from each process on master. Data at
!     these nodes will later be gathered by master from other processes
      ALLOCATE(ib%cm%n(cm%np()))
      ib%cm%n = 0
      CALL MPI_GATHER(n, 1, mpint, ib%cm%n, 1, mpint, master, cm%com(),
     2   ierr)

!     The processes that do not have any nodal traces return
      IF (.NOT.cm%mas() .AND. n.EQ.0) THEN
         DEALLOCATE(incNd, ptr)
         ALLOCATE(ib%cm%gE(0))
         RETURN
      END IF

!     Master receives list of all nodal traces from other processes
      IF (cm%mas()) THEN
         n = SUM(ib%cm%n)
         a = MAXVAL(ib%cm%n)
         ALLOCATE(ib%cm%gE(n), rA(a,cm%np()), rReq(cm%np()))
         ib%cm%gE  = 0
         rA(:,:)   = 0
         DO i=1, cm%np()
            n = ib%cm%n(i)
            IF (n .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               rA(1:n,i) = ptr(:)
            ELSE
               tag = i*100
               CALL MPI_IRECV(rA(1:n,i), n, mpint, i-1, tag, cm%com(),
     2            rReq(i), ierr)
            END IF
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ib%cm%n(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO
         DEALLOCATE(rReq)

         a = 0
         DO i=1, cm%np()
            n = ib%cm%n(i)
            ib%cm%gE(a+1:a+n) = rA(1:n,i)
            a = a + n
         END DO
      ELSE
         ALLOCATE(ib%cm%gE(0), rA(0,0))
         tag = cm%tF() * 100
         CALL MPI_SEND(ptr, n, mpint, master, tag, cm%com(), ierr)
      END IF

      DEALLOCATE(incNd, ptr, rA)

      RETURN
      END SUBROUTINE IB_SETCOMMU
!####################################################################
!     Set ighost field
!        ighost is set for both solids and shells
!        For solids, ighost(A) = 1   =>  iblank(A)=1 and is connected to
!           atleast one iblank=0 node
!        For shells, ighost(A) = 1   =>  node belongs to a ghost cell
      SUBROUTINE IB_SETIGHOST()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER :: a, b, e, Ac, iM, jM

      REAL(KIND=8), ALLOCATABLE :: rG(:)

      ALLOCATE(rG(tnNo))
      rG = 0D0

!     Update mesh and nodal ghost cell pointers based on iblank field
      DO iM=1, nMsh
         DO e=1, msh(iM)%nEl
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b = b + iblank(Ac)
            END DO
            IF (b.GT.0 .AND. b.LT.msh(iM)%eNoN) THEN
               msh(iM)%iGC(e) = 1
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  rG(Ac) = REAL(iblank(Ac), KIND=8)
               END DO
            END IF
         END DO
      END DO

!     For Shells, use IB traces to identify ghost cells
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
            DO b=1, ib%msh(iM)%trc%n
               e  = ib%msh(iM)%trc%ptr(1,b)
               jM = ib%msh(iM)%trc%ptr(2,b)
               DO a=1, msh(jM)%eNoN
                  Ac = msh(jM)%IEN(a,e)
                  rG(Ac) = 1.0D0
               END DO
            END DO
         END IF
      END DO
      CALL COMMU(rG)

      ighost = 0
      DO a=1, tnNo
         IF (rG(a) .GT. 1D-6) ighost(a) = 1
      END DO

      RETURN
      END SUBROUTINE IB_SETIGHOST
!####################################################################
!     Write IB solution to a vtu file
      SUBROUTINE IB_WRITEVTUS(lA, lY, lU)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: lA(nsd,ib%tnNo), lY(nsd+1,ib%tnNo),
     2   lU(nsd,ib%tnNo)

      TYPE(dataType) :: d(ib%nMsh)
      TYPE(vtkXMLType) :: vtu

      CHARACTER(LEN=stdL) :: fName
      INTEGER :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo, nEl, s, l, ie,
     2   is, nSh, oGrp, outDof, nOut, cOut

      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)
      INTEGER, ALLOCATABLE :: outS(:), tmpI(:,:)
      REAL(KIND=8), ALLOCATABLE :: tmpV(:,:)

      IF (ib%mthd .NE. ibMthd_IFEM) RETURN

      IF (cm%slv()) THEN
         ib%savedOnce = .TRUE.
         RETURN
      END IF

      nOut   = 1
      outDof = nsd
      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutIB
            IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
            nOut   = nOut + 1
            outDof = outDof + eq(iEq)%outIB(iOut)%l
         END DO
      END DO

      ALLOCATE(outNames(nOut), outS(nOut+1))

!     Prepare all solultions in to dataType d
      nNo = 0
      nEl = 0
      DO iM=1, nMsh
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""

         IF (ib%msh(iM)%eType .EQ. eType_NRB) err =
     2      " Outputs for NURBS data is under development"

         d(iM)%nNo     = ib%msh(iM)%nNo
         d(iM)%nEl     = ib%msh(iM)%nEl
         d(iM)%eNoN    = ib%msh(iM)%eNoN
         d(iM)%vtkType = ib%msh(iM)%vtkType

         ALLOCATE(d(iM)%x(outDof,d(iM)%nNo),
     2      d(iM)%IEN(d(iM)%eNoN,d(iM)%nEl))
         DO a=1, ib%msh(iM)%nNo
            Ac = ib%msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = ib%x(:,Ac)
         END DO

         DO e=1, ib%msh(iM)%nEl
            d(iM)%IEN(:,e) = ib%msh(iM)%IEN(:,e)
         END DO

         DO iEq=1, nEq
            DO iOut=1, eq(iEq)%nOutIB
               IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
               l  = eq(iEq)%outIB(iOut)%l
               s  = eq(iEq)%s + eq(iEq)%outIB(iOut)%o
               e  = s + l - 1

               cOut = cOut + 1
               is   = outS(cOut)
               ie   = is + l - 1
               outS(cOut+1)   = ie + 1
               outNames(cOut) = "IB_"//TRIM(eq(iEq)%outIB(iOut)%name)

               oGrp = eq(iEq)%outIB(iOut)%grp
               SELECT CASE (oGrp)
                  CASE (outGrp_NA)
                  err = "Undefined output grp in VTK"
               CASE (outGrp_A)
                  DO a=1, ib%msh(iM)%nNo
                     Ac = ib%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lA(s:e,Ac)
                  END DO
               CASE (outGrp_Y)
                  DO a=1, ib%msh(iM)%nNo
                     Ac = ib%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lY(s:e,Ac)
                  END DO
               CASE (outGrp_D)
                  DO a=1, ib%msh(iM)%nNo
                     Ac = ib%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lU(s:e,Ac)
                  END DO
               CASE DEFAULT
                  err = "Undefined output "//
     2               TRIM(eq(iEq)%outIB(iOut)%name)
               END SELECT
            END DO
         END DO

         ALLOCATE(d(iM)%xe(d(iM)%nEl,1))
         IF (.NOT.savedOnce) THEN
            IF (ALLOCATED(ib%dmnID)) THEN
               d(iM)%xe(:,1) = ib%msh(iM)%eId(:)
            ELSE
               d(iM)%xe(:,1) = 1
            END IF
         END  IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
      END DO

      ALLOCATE(tmpV(maxnsd,nNo))

!     Writing to vtu file (master only)
      IF (cTS .GE. 1000) THEN
         fName = STR(cTS)
      ELSE
         WRITE(fName,'(I3.3)') cTS
      END IF

      fName = TRIM(saveName)//"_ib_"//TRIM(ADJUSTL(fName))//".vtu"
      dbg = "Writing VTU"

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

!     Writing the position data
      iOut = 1
      s    = outS(iOut)
      e    = outS(iOut+1)-1
      nSh  = 0
      tmpV = 0D0
      DO iM=1, ib%nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%x(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO
      CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

!     Writing the connectivity data
      nSh = -1
      DO iM=1, ib%nMsh
         ALLOCATE(tmpI(d(iM)%eNoN,d(iM)%nEl))
         DO e=1, d(iM)%nEl
            tmpI(:,e) = d(iM)%IEN(:,e) + nSh
         END DO
         CALL putVTK_elemIEN(vtu, tmpI, d(iM)%vtkType, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (ien)"
         DEALLOCATE(tmpI)
         nSh = nSh + d(iM)%nNo
      END DO

!     Writing all solutions
      DO iOut=2, nOut
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         s = outS(iOut)
         e = outS(iOut+1) - 1
         l = e - s + 1
         ALLOCATE(tmpV(l, nNo))
         nSh = 0
         DO iM=1, ib%nMsh
            DO a=1, d(iM)%nNo
               tmpV(:,a+nSh) = d(iM)%x(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (point data)"
      END DO

!     Write element-based variables
      IF (.NOT.savedOnce .OR. mvMsh) THEN
         ib%savedOnce = .TRUE.
         ALLOCATE(tmpI(1,nEl))
!     Write the domain ID
         IF (ALLOCATED(ib%dmnID)) THEN
            Ec = 0
            DO iM=1, ib%nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = d(iM)%xe(e,1)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (dom id)"
         END IF

!     Write the mesh ID
         IF (ib%nMsh .GT. 1) THEN
            Ec = 0
            DO iM=1, ib%nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = iM
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Mesh_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (mesh id)"
         END IF
         DEALLOCATE(tmpI)
      END IF

      DO iM=1, nMsh
         CALL DESTROY(d(iM))
      END DO

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE IB_WRITEVTUS
!####################################################################
!     Apply Dirichlet boundary conditions on the immersed faces
      SUBROUTINE IB_SETBCDIR(Ab, Yb, Ub)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: Ab(nsd,ib%tnNo), Yb(nsd+1,ib%tnNo),
     2   Ub(nsd,ib%tnNo)

      INTEGER :: iFa, iM, iEq, iBc, a, Ac, nNo, lDof, i
      LOGICAL :: eDir(maxnsd)

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: nV(:,:), tmpA(:,:), tmpY(:,:)

      DO iEq=1, nEq
         IF (eq(iEq)%phys .NE. phys_fluid .AND.
     2       eq(iEq)%phys .NE. phys_FSI) CYCLE
         DO iBc=1, eq(iEq)%nBcIB
            IF (.NOT.BTEST(eq(iEq)%bcIB(iBc)%bType,bType_Dir)) CYCLE
            eDir = .FALSE.
            lDof = 0
            DO i=1, nsd
               IF (eq(iEq)%bcIB(iBc)%eDrn(i) .NE. 0) THEN
                  eDir(i) = .TRUE.
                  lDof = lDof + 1
               END IF
            END DO
            IF (lDof .EQ. 0) lDof = nsd

            iFa = eq(iEq)%bcIB(iBc)%iFa
            iM  = eq(iEq)%bcIB(iBc)%iM

!        Prepare a pointer list and normals for a face or a shell
            IF (iFa .EQ. 0) THEN
               IF (.NOT.ib%msh(iM)%lShl) err = "Trying to apply Dir "//
     2            "BC on a non-shell mesh surface"
               nNo = ib%msh(iM)%nNo
               ALLOCATE(ptr(nNo), nV(nsd,nNo))
               ptr = 0
               DO a=1, nNo
                  ptr(a)  = ib%msh(iM)%gN(a)
                  nV(:,a) = ib%msh(iM)%nV(:,a)
               END DO
            ELSE
               nNo = ib%msh(iM)%fa(iFa)%nNo
               ALLOCATE(ptr(nNo), nV(nsd,nNo))
               ptr = 0
               DO a=1, nNo
                  ptr(a)  = ib%msh(iM)%fa(iFa)%gN(a)
                  nV(:,a) = ib%msh(iM)%fa(iFa)%nV(:,a)
               END DO
            END IF

            ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
            CALL IB_SETBCDIRL(eq(iEq)%bcIB(iBc), nNo, lDof, nV, tmpA,
     2         tmpY)

            IF (ANY(eDir)) THEN
               DO a=1, nNo
                  Ac = ptr(a)
                  IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                     DO i=1, nsd
                        lDof = 0
                        IF (eDir(i)) THEN
                           lDof = lDof + 1
                           Yb(i,Ac) = tmpA(lDof,a)
                           Ub(i,Ac) = tmpY(lDof,a)
                        END IF
                     END DO
                  ELSE
                     DO i=1, nsd
                        lDof = 0
                        IF (eDir(i)) THEN
                           lDof = lDof + 1
                           Ab(i,Ac) = tmpA(lDof,a)
                           Yb(i,Ac) = tmpY(lDof,a)
                        END IF
                     END DO
                  END IF
               END DO
            ELSE
               DO a=1, nNo
                  Ac = ptr(a)
                  IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                     Yb(1:nsd,Ac) = tmpA(:,a)
                     Ub(1:nsd,Ac) = tmpY(:,a)
                  ELSE
                     Ab(1:nsd,Ac) = tmpA(:,a)
                     Yb(1:nsd,Ac) = tmpY(:,a)
                  END IF
               END DO
            END IF

            DEALLOCATE(ptr, nV, tmpA, tmpY)
         END DO
      END DO

      RETURN
      END SUBROUTINE IB_SETBCDIR
!--------------------------------------------------------------------
      SUBROUTINE IB_SETBCDIRL(lBc, nNo, lDof, nvL, lA, lY)
      USE COMMOD
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      INTEGER, INTENT(IN) :: lDof, nNo
      REAL(KIND=8), INTENT(IN) :: nvL(nsd,nNo)
      REAL(KIND=8), INTENT(INOUT) :: lA(lDof,nNo), lY(lDof,nNo)

      INTEGER :: a, i
      REAL(KIND=8) :: dirY, dirA, nV(nsd)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (lDof .NE. lBc%gm%dof) err = "Inconsistent DOF to apply "//
     2      "Gen BC"
         IF (nNo .NE. SIZE(lBc%gm%d,2)) err = "Inconsistent nodes "//
     2      "to apply Gen BC"
         CALL IGBC(lBc%gm, lY, lA)
         RETURN
      ELSE IF (BTEST(lBc%bType,bType_ustd)) THEN
         CALL IFFT(lBc%gt, dirY, dirA)
      ELSE ! std / cpl
         dirA = 0D0
         dirY = lBc%g
      END IF

      IF (lDof .EQ. nsd) THEN
         DO a=1, nNo
            nV      = nvL(:,a)
            lA(:,a) = dirA*lBc%gx(a)*nV
            lY(:,a) = dirY*lBc%gx(a)*nV
         END DO
      ELSE
         DO a=1, nNo
            DO i=1, lDof
               lA(i,a) = dirA*lBc%gx(a)
               lY(i,a) = dirY*lBc%gx(a)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE IB_SETBCDIRL
!####################################################################
!     Project velocity from IB to ghost/finite cells using distance-
!     weighted interpolation. This type of projection is used for SSM
!     method only.
      SUBROUTINE IB_PRJCTU(Yg, Dg, Yb)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: Yg(tDof,tnNo)
      REAL(KIND=8), INTENT(IN) :: Dg(tDof,tnNo), Yb(nsd+1,ib%tnNo)

      INTEGER a, b, i, g, e, s, Ac, Bc, Ec, iM, jM, iFa, iEq, eNoN
      REAL(KIND=8) :: w, xe(nsd), xp(nsd), yp(nsd), ctime

      REAL(KIND=8), ALLOCATABLE :: sA(:), sV(:,:), xl(:,:), yl(:,:),
     2   N(:)

      ctime = CPUT()

!     SSM (stair-step method): interpolate velocities on background mesh
!     node based on the traces of IB. We apply an inverse-distance
!     weighted interpolation operator on the ghost nodes.
      ALLOCATE(sA(tnNo), sV(nsd,tnNo))
      sV = 0D0
      sA = 0D0
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
            eNoN = ib%msh(iM)%eNoN
            ALLOCATE(xl(nsd,eNoN), yl(nsd,eNoN), N(eNoN))
            DO i=1, ib%msh(iM)%trc%n
               e = ib%msh(iM)%trc%gE(1,i)
               g = ib%msh(iM)%trc%gE(2,i)
               DO a=1, eNoN
                  Ac = ib%msh(iM)%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
                  yl(:,a) = Yb(1:nsd,Ac)
                  N(a)    = ib%msh(iM)%N(a,g)
               END DO

               xp = 0D0
               yp = 0D0
               DO a=1, eNoN
                  xp = xp + N(a)*xl(:,a)
                  yp = yp + N(a)*yl(:,a)
               END DO

               Ec = ib%msh(iM)%trc%ptr(1,i)
               jM = ib%msh(iM)%trc%ptr(2,i)
               DO b=1, msh(jM)%eNoN
                  Bc = msh(jM)%IEN(b,Ec)
                  IF (ighost(Bc) .EQ. 1) THEN
                     xe = x(:,Bc)
                     IF (mvMsh) xe = xe + Dg(nsd+2:2*nsd+1,Bc)
                     w = SQRT(SUM( (xe - xp)**2 ))
                     IF (ISZERO(w)) w = eps
                     w = 1D0/w
                     sA(Bc) = sA(Bc) + w
                     sV(:,Bc) = sV(:,Bc) + w*yp(:)
                  END IF
               END DO
            END DO
            DEALLOCATE(xl, yl, N)
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               eNoN = ib%msh(iM)%fa(iFa)%eNoN
               ALLOCATE(xl(nsd,eNoN), yl(nsd,eNoN), N(eNoN))
               DO i=1, ib%msh(iM)%fa(iFa)%trc%n
                  e = ib%msh(iM)%fa(iFa)%trc%gE(1,i)
                  g = ib%msh(iM)%fa(iFa)%trc%gE(2,i)
                  DO a=1, eNoN
                     Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                     xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
                     yl(:,a) = Yb(1:nsd,Ac)
                     N(a)    = ib%msh(iM)%fa(iFa)%N(a,g)
                  END DO

                  xp = 0D0
                  yp = 0D0
                  DO a=1, eNoN
                     xp = xp + N(a)*xl(:,a)
                     yp = yp + N(a)*yl(:,a)
                  END DO

                  Ec = ib%msh(iM)%fa(iFa)%trc%ptr(1,i)
                  jM = ib%msh(iM)%fa(iFa)%trc%ptr(2,i)
                  DO b=1, msh(jM)%eNoN
                     Bc = msh(jM)%IEN(b,Ec)
                     IF (ighost(Bc) .EQ. 1) THEN
                        xe = x(:,Bc)
                        IF (mvMsh) xe = xe + Dg(nsd+2:2*nsd+1,Bc)
                        w = SQRT(SUM( (xe - xp)**2 ))
                        IF (ISZERO(w)) w = eps
                        w = 1D0/w
                        sA(Bc) = sA(Bc) + w
                        sV(:,Bc) = sV(:,Bc) + w*yp(:)
                     END IF
                  END DO
               END DO
               DEALLOCATE(xl, yl, N)
            END DO
         END IF
      END DO

      CALL COMMU(sA)
      CALL COMMU(sV)

      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_fluid .OR.
     2       eq(iEq)%phys .EQ. phys_FSI) EXIT
      END DO

      s = eq(iEq)%s
      e = eq(iEq)%e
      IF (eq(iEq)%dof .EQ. nsd+1) e = e - 1
      DO a=1, tnNo
         IF (.NOT.ISZERO(sA(a)) .AND. ighost(a).EQ.1) THEN
            Yg(s:e,a) = sV(:,a)/sA(a)
         END IF
      END DO
      DEALLOCATE(sA, sV)

      ib%callD(1) = CPUT() - ctime

      RETURN
      END SUBROUTINE IB_PRJCTU
!####################################################################
!     Add contribution from IB to the residue (RHS)
      SUBROUTINE IB_CONSTRUCT(Yg, Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER a

      IF (eq(cEq)%phys .NE. phys_fluid .AND.
     2    eq(cEq)%phys .NE. phys_FSI) RETURN

      IF (ib%mthd .EQ. ibMthd_Penalty .OR.
     2    ib%mthd .EQ. ibMthd_Nitsche) THEN
         CALL IB_SETBCDIRA(Yg, Dg)

      ELSE IF (ib%mthd .EQ. ibMthd_IFEM) THEN
         IF (.NOT.ALLOCATED(Rib)) err =
     2      "Correction needed IFEM construct"
         DO a=1, tnNo
            R(1:nsd,a) = R(1:nsd,a) - Rib(:,a)
         END DO
      END IF

      RETURN
      END SUBROUTINE IB_CONSTRUCT
!####################################################################
!     Treat IB with Dir BCs using augmented Lagrange multiplier method
!     (extended Nitsche), or using pure Penalty method w-/wo- feedback
      SUBROUTINE IB_SETBCDIRA(Yg, Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER :: iBc, iFa, iM
      REAL(KIND=8) ctime

      ctime = CPUT()
      DO iBc=1, eq(cEq)%nBcIB
         IF (.NOT.BTEST(eq(cEq)%bcIB(iBc)%bType, bType_Dir)) CYCLE
         iM  = eq(cEq)%bcIB(iBc)%iM
         iFa = eq(cEq)%bcIB(iBc)%iFa
         IF (iFa .NE. 0) THEN
            CALL IB_SETBCDIRAL(eq(cEq)%bcIB(iBc), ib%msh(iM)%fa(iFa),
     2         Yg, Dg)
         ELSE
            CALL IB_SETBCDIRASL(eq(cEq)%bcIB(iBc), ib%msh(iM), Yg, Dg)
         END IF
      END DO
      ib%callD(1) = CPUT() - ctime

      RETURN
      END SUBROUTINE IB_SETBCDIRA
!--------------------------------------------------------------------
      SUBROUTINE IB_SETBCDIRAL(lBc, lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      LOGICAL :: flag
      INTEGER :: a, e, g, i, Ac, Ec, nNo, nEl, nG, eNoN, eNoNb, cPhys,
     2   iM, nGP
      REAL(KIND=8) :: w, Jac, fb(nsd), xp(nsd), xi(nsd), nV(nsd),
     2   ub(nsd), tauB(2), Ks(nsd,nsd)

      INTEGER, ALLOCATABLE :: ePtr(:,:,:), ptr(:)
      REAL(KIND=8), ALLOCATABLE :: lD(:,:), N(:), Nb(:), Nxi(:,:),
     2   Nx(:,:), xl(:,:), xbl(:,:), yl(:,:), ubl(:,:), fbl(:,:),
     3   lR(:,:), lK(:,:,:)

      nNo   = lFa%nNo
      nEl   = lFa%nEl
      nG    = lFa%nG
      eNoNb = lFa%eNoN

!     Get penalty constants
      tauB  = lBc%tauB

!     Project IB face nodal values with Dirichlet BCs to global array
      ALLOCATE(lD(nsd,tnNo))
      lD = 0D0
      DO a=1, tnNo
         IF (mvMsh) lD(:,a) = Dg(nsd+2:2*nsd+1,a)
      END DO

!     Load all the ghost cell pointers into a global array
      ALLOCATE(ePtr(2,nG,nEl))
      ePtr = 0
      DO i=1, lFa%trc%n
         e = lFa%trc%gE(1,i)
         g = lFa%trc%gE(2,i)
         ePtr(:,g,e) = lFa%trc%ptr(:,i)
      END DO

      ALLOCATE(Nb(eNoNb), xbl(nsd,eNoNb), ubl(nsd,eNoNb),fbl(nsd,eNoNb))

!     Loop over all the elements and assemble to the residue
      nGP = 0
      DO e=1, nEl
         DO a=1, eNoNb
            Ac = lFa%IEN(a,e)
            xbl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            ubl(:,a) = ib%Yo(1:nsd,Ac)
            fbl(:,a) = ib%R(:,Ac)
         END DO

         DO g=1, nG
            Ec = ePtr(1,g,e)
            iM = ePtr(2,g,e)
            IF (Ec .EQ. 0) CYCLE

!           Keep track of the total trace count
            nGP = nGP + 1

            CALL GNNIB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
!           Normal pointing outwards from solid to fluid
            nV = nV / Jac
            w  = lFa%w(g)*Jac
            Nb = lFa%N(:,g)

!           Coordinates, IB velocity and feedback forcing at the
!           integration point
            xp  = 0D0
            ub  = 0D0
            fb  = 0D0
            DO a=1, eNoNb
               xp = xp + Nb(a)*xbl(:,a)
               ub = ub + Nb(a)*ubl(:,a)
               fb = fb + Nb(a)*fbl(:,a)
            END DO

            cDmn  = DOMAIN(msh(iM), cEq, Ec)
            cPhys = eq(cEq)%dmn(cDmn)%phys
!           Domain physics has to be fluid
            IF (cPhys .NE. phys_fluid) err = "IB formulation is for "//
     2         "fluid phys only"

            eNoN = msh(iM)%eNoN
            ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN), xl(nsd,eNoN),
     2         yl(tDof,eNoN), ptr(eNoN), lR(dof,eNoN),
     3         lK(dof*dof,eNoN,eNoN))

!           Transfer to local arrays
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,Ec)
               ptr(a) = Ac
               xl(:,a) = x(:,Ac) + lD(:,Ac)
               yl(:,a) = Yg(:,Ac)
            END DO

!           Get the parametric coordinates and shape functions of the
!           trace w.r.t. background mesh
            xi = 0D0
            DO i=1, msh(iM)%nG
               xi = xi + msh(iM)%xi(:,i)
            END DO
            xi(:) = xi(:) / REAL(msh(iM)%nG,KIND=8)

            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge. "//
     2         "Invalid trace pointer"

!           Get the shape functions and their derivatives at the
!           integration point, as a function of the basis defined on the
!           background mesh
            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)
!           Check if the shape functions are correct
            i = 0
            DO a=1, eNoN
               IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
            END DO
            IF (i .NE. eNoN) err =
     2         "Error in computing IB trace shape functions"

!           Compute spatial derivatives of shape functions
            CALL GNN(eNoN, Nxi, xl, Nx, Jac, Ks)

!           Compute the residue and stiffness matrices
            lR = 0D0
            lK = 0D0
            IF (nsd .EQ. 3) THEN
               IF (ib%mthd .EQ. ibMthd_Penalty) THEN
                  CALL BPFLUID3D(eNoN, w, N, yl, ub, fb, nV, tauB,
     2               lBc%fbN, lR, lK)
               ELSE IF (ib%mthd .EQ. ibMthd_Nitsche) THEN
!                 Flip the normal to point away from fluid
                  nV = -nV
                  CALL BWFLUID3D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR,
     2               lK)
               END IF

            ELSE
               IF (ib%mthd .EQ. ibMthd_Penalty) THEN
                  CALL BPFLUID2D(eNoN, w, N, yl, ub, fb, nV, tauB,
     2               lBc%fbN, lR, lK)
               ELSE IF (ib%mthd .EQ. ibMthd_Nitsche) THEN
!                 Flip the normal to point away from fluid domain
                  nV = -nV
                  CALL BWFLUID2D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR,
     2               lK)
               END IF
            END IF

!           Assemble
#ifdef WITH_TRILINOS
            IF (useTrilinosAssemAndLS) THEN
               CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
            ELSE
#endif
               CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
      END IF
#endif
            DEALLOCATE(N, Nxi, Nx, xl, yl, ptr, lR, lK)
         END DO
      END DO

      nGP = cm%reduce(nGP)
      i   = nEl*nG
      IF (nGP .NE. i) THEN
         err = " Error: found only "//STR(nGP)//" integration point "//
     2      "traces out of "//STR(i)
      END IF
      DEALLOCATE(lD, ePtr, Nb, xbl, ubl, fbl)

      RETURN
      END SUBROUTINE IB_SETBCDIRAL
!--------------------------------------------------------------------
!     Treat shells here
      SUBROUTINE IB_SETBCDIRASL(lBc, lM, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      LOGICAL :: flag
      INTEGER :: a, e, g, i, Ac, Ec, nNo, nEl, nG, eNoN, eNoNb, cPhys,
     2   iM, nGP
      REAL(KIND=8) :: w, Jac, fb(nsd), xp(nsd), xi(nsd), nV(nsd),
     2   ub(nsd), tauB(2), Ks(nsd,nsd), gC(nsd,nsd-1)

      INTEGER, ALLOCATABLE :: ePtr(:,:,:), ptr(:)
      REAL(KIND=8), ALLOCATABLE :: lD(:,:), N(:), Nb(:), Nxi(:,:),
     2   Nx(:,:), xl(:,:), xbl(:,:), yl(:,:), ubl(:,:), fbl(:,:),
     3   lR(:,:), lK(:,:,:)

      IF (ib%mthd .NE. ibMthd_Penalty) err = "Only penalty method "//
     2   "is allowed for immersed shells"

      nNo   = lM%nNo
      nEl   = lM%nEl
      nG    = lM%nG
      eNoNb = lM%eNoN

!     Get penalty and feedback constants
      tauB  = lBc%tauB

!     Project IB face nodal values with Dirichlet BCs to global array
      ALLOCATE(lD(nsd,tnNo))
      lD = 0D0
      DO a=1, tnNo
         IF (mvMsh) lD(:,a) = Dg(nsd+2:2*nsd+1,a)
      END DO

!     Load all the ghost cell pointers into a global array
      ALLOCATE(ePtr(2,nG,nEl))
      ePtr = 0
      DO i=1, lM%trc%n
         e = lM%trc%gE(1,i)
         g = lM%trc%gE(2,i)
         ePtr(:,g,e) = lM%trc%ptr(:,i)
      END DO

      ALLOCATE(Nb(eNoNb), xbl(nsd,eNoNb), ubl(nsd,eNoNb),fbl(nsd,eNoNb))

!     Loop over all the elements and assemble to the residue
      nGP = 0
      DO e=1, nEl
         DO a=1, eNoNb
            Ac = lM%IEN(a,e)
            xbl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            ubl(:,a) = ib%Yo(1:nsd,Ac)
            fbl(:,a) = ib%R(:,Ac)
         END DO

         DO g=1, nG
            Ec = ePtr(1,g,e)
            iM = ePtr(2,g,e)
            IF (Ec .EQ. 0) CYCLE

!           Keep track of the total trace count
            nGP = nGP + 1

!           IB shell normal
            CALL GNNS(eNoNb, lM%Nx(:,:,g), xbl, nV, gC, gC)
            Jac = SQRT(NORM(nV))
!           Normal pointing outwards from solid to fluid
            nV  = nV/Jac
            w   = lM%w(g) * Jac
            Nb  = lM%N(:,g)

!           Coordinates, IB velocity and feedback forcing at the
!           integration point
            xp  = 0D0
            ub  = 0D0
            fb  = 0D0
            DO a=1, eNoNb
               xp = xp + Nb(a)*xbl(:,a)
               ub = ub + Nb(a)*ubl(:,a)
               fb = fb + Nb(a)*fbl(:,a)
            END DO

            cDmn  = DOMAIN(msh(iM), cEq, Ec)
            cPhys = eq(cEq)%dmn(cDmn)%phys
!           Domain physics has to be fluid
            IF (cPhys .NE. phys_fluid) err = "IB formulation is for "//
     2         "fluid phys only"

            eNoN = msh(iM)%eNoN
            ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN), xl(nsd,eNoN),
     2         yl(tDof,eNoN), ptr(eNoN), lR(dof,eNoN),
     3         lK(dof*dof,eNoN,eNoN))

!           Transfer to local arrays
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,Ec)
               ptr(a) = Ac
               xl(:,a) = x(:,Ac) + lD(:,Ac)
               yl(:,a) = Yg(:,Ac)
            END DO

!           Get the parametric coordinates and shape functions of the
!           trace w.r.t. background mesh
            xi = 0D0
            DO i=1, msh(iM)%nG
               xi = xi + msh(iM)%xi(:,i)
            END DO
            xi(:) = xi(:) / REAL(msh(iM)%nG,KIND=8)

            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge. "//
     2         "Invalid trace pointer"

!           Get the shape functions and their derivatives at the
!           integration point, as a function of the basis defined on the
!           background mesh
            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)
!           Check if the shape functions are correct
            i = 0
            DO a=1, eNoN
               IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
            END DO
            IF (i .NE. eNoN) err =
     2         "Error in computing IB trace shape functions"

!           Compute spatial derivatives of shape functions
            CALL GNN(eNoN, Nxi, xl, Nx, Jac, Ks)

!           Compute the residue and stiffness matrices
            lR = 0D0
            lK = 0D0
            IF (nsd .EQ. 3) THEN
               IF (ib%mthd .EQ. ibMthd_Penalty) THEN
                  CALL BPFLUID3D(eNoN, w, N, yl, ub, fb, nV, tauB,
     2               lBc%fbN, lR, lK)
               ELSE IF (ib%mthd .EQ. ibMthd_Nitsche) THEN
!                 Flip the normal to point away from fluid
                  nV = -nV
                  CALL BWFLUID3D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR,
     2               lK)
               END IF

            ELSE
               IF (ib%mthd .EQ. ibMthd_Penalty) THEN
                  CALL BPFLUID2D(eNoN, w, N, yl, ub, fb, nV, tauB,
     2               lBc%fbN, lR, lK)
               ELSE IF (ib%mthd .EQ. ibMthd_Nitsche) THEN
!                 Flip the normal to point away from fluid domain
                  nV = -nV
                  CALL BWFLUID2D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR,
     2               lK)
               END IF
            END IF

!           Assemble
#ifdef WITH_TRILINOS
            IF (useTrilinosAssemAndLS) THEN
               CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
            ELSE
#endif
               CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
      END IF
#endif
            DEALLOCATE(N, Nxi, Nx, xl, yl, ptr, lR, lK)
         END DO
      END DO

      nGP = cm%reduce(nGP)
      i   = nEl*nG
      IF (nGP .NE. i) THEN
         err = " Error: found only "//STR(nGP)//" integration point "//
     2      "traces out of "//STR(i)
      END IF
      DEALLOCATE(lD, ePtr, Nb, xbl, ubl, fbl)

      RETURN
      END SUBROUTINE IB_SETBCDIRASL
!####################################################################
!     Computes the feedback forcing term for IBs with Dirichlet BCs
!     treated using Penalty method
      SUBROUTINE IB_SETFBF()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER iEq, iBc, iM, iFa, a, i
      REAL(KIND=8) ctime

      REAL(KIND=8), ALLOCATABLE :: sA(:), sF(:,:)

      IF (.NOT.ib%fbFlag) RETURN

      ctime = CPUT()
      ALLOCATE(sA(ib%tnNo), sF(nsd,ib%tnNo))
      sA = 0D0
      sF = 0D0
      DO iEq=1, nEq
         IF (eq(iEq)%phys .NE. phys_fluid .AND.
     2       eq(iEq)%phys .NE. phys_FSI) CYCLE
         DO iBc=1, eq(iEq)%nBcIB
            IF (.NOT.BTEST(eq(iEq)%bcIB(iBc)%bType, bType_Dir)) CYCLE
            iM  = eq(iEq)%bcIB(iBc)%iM
            iFa = eq(iEq)%bcIB(iBc)%iFa
            IF (iFa .NE. 0) THEN
               CALL IB_SETFBFL(eq(iEq)%bcIB(iBc), ib%msh(iM)%fa(iFa),
     2            Yn, Dn, sF, sA)
            ELSE
               CALL IB_SETFBFSL(eq(iEq)%bcIB(iBc), ib%msh(iM), Yn, Dn,
     2            sF, sA)
            END IF
         END DO
      END DO

!     Synchronize ib%R across all the processes
c      CALL IB_SYNC(sF)
c      CALL IB_SYNC(sA)

      DO a=1, ib%tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            DO i=1, nsd
               ib%R(i,a) = ib%R(i,a) + (sF(i,a)/sA(a))
            END DO
         END IF
      END DO

      DEALLOCATE(sA, sF)

      ib%callD(1) = ib%callD(1) + CPUT() - ctime

      RETURN
      END SUBROUTINE IB_SETFBF
!--------------------------------------------------------------------
      SUBROUTINE IB_SETFBFL(lBc, lFa, Yg, Dg, sF, sA)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      REAL(KIND=8), INTENT(INOUT) :: sA(tnNo), sF(nsd,tnNo)

      LOGICAL :: flag
      INTEGER :: a, e, i, g, Ac, Ec, nNo, nEl, nG, eNoN, eNoNb, iM, nGP
      REAL(KIND=8) :: w, Jac, fb(nsd), xp(nsd), xi(nsd), nV(nsd),u(nsd),
     2   ub(nsd), tauF

      INTEGER, ALLOCATABLE :: ePtr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: lD(:,:), N(:), Nb(:), Nxi(:,:),
     2   xl(:,:), xbl(:,:), ul(:,:), ubl(:,:)

      nNo   = lFa%nNo
      nEl   = lFa%nEl
      nG    = lFa%nG
      eNoNb = lFa%eNoN

!     Get penalty constants
      tauF  = lBc%tauF

      ALLOCATE(lD(nsd,tnNo))
      lD = 0D0
      DO a=1, tnNo
         IF (mvMsh) lD(:,a) = Dg(nsd+2:2*nsd+1,a)
      END DO

!     Load all the ghost cell pointers into a global array
      ALLOCATE(ePtr(2,nG,nEl))
      ePtr = 0
      DO i=1, lFa%trc%n
         e = lFa%trc%gE(1,i)
         g = lFa%trc%gE(2,i)
         ePtr(:,g,e) = lFa%trc%ptr(:,i)
      END DO

      ALLOCATE(Nb(eNoNb), xbl(nsd,eNoNb), ubl(nsd,eNoNb))

!     Loop over all the elements and assemble to the residue
      nGP = 0
      DO e=1, nEl
         DO a=1, eNoNb
            Ac = lFa%IEN(a,e)
            xbl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            ubl(:,a) = ib%Yo(1:nsd,Ac)
         END DO

         DO g=1, nG
            Ec = ePtr(1,g,e)
            iM = ePtr(2,g,e)
            IF (Ec .EQ. 0) CYCLE

!           Keep track of the total trace count
            nGP = nGP + 1

            CALL GNNIB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
!           Normal pointing outwards from solid to fluid
            nV = nV / Jac
            w  = lFa%w(g)*Jac
            Nb = lFa%N(:,g)

!           Coordinates, IB velocity and feedback forcing at the
!           integration point
            xp  = 0D0
            ub  = 0D0
            DO a=1, eNoNb
               xp = xp + Nb(a)*xbl(:,a)
               ub = ub + Nb(a)*ubl(:,a)
            END DO

            cDmn = DOMAIN(msh(iM), cEq, Ec)
            eNoN = msh(iM)%eNoN
            ALLOCATE(N(eNoN), Nxi(nsd,eNoN), xl(nsd,eNoN), ul(nsd,eNoN))

!           Transfer to local arrays
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,Ec)
               xl(:,a) = x(:,Ac) + lD(:,Ac)
               ul(:,a) = Yg(1:nsd,Ac)
            END DO

!           Get the parametric coordinates and shape functions of the
!           trace w.r.t. background mesh
            xi = 0D0
            DO i=1, msh(iM)%nG
               xi = xi + msh(iM)%xi(:,i)
            END DO
            xi(:) = xi(:) / REAL(msh(iM)%nG,KIND=8)

            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge. "//
     2         "Invalid trace pointer"

!           Get the shape functions and their derivatives at the
!           integration point, as a function of the basis defined on the
!           background mesh
            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)
!           Check if the shape functions are correct
            i = 0
            DO a=1, eNoN
               IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
            END DO
            IF (i .NE. eNoN) err =
     2         "Error in computing IB trace shape functions"

!           Get the fluid velocity at the integration point
            u = 0D0
            DO a=1, eNoN
               u(:) = u(:) + N(a)*ul(:,a)
            END DO

!           Compute the feedback forcing with penalty forces along
!           normal direction
            IF (lBc%fbN) THEN
               fb(:) = tauF * NORM(u-ub, nV)
            ELSE
               fb(:) = tauF * (u(:)-ub(:))
            END IF

!           Transfer back to global vector
            DO a=1, eNoNb
               Ac = lFa%IEN(a,e)
               sF(:,Ac) = sF(:,Ac) + w*Nb(a)*fb(:)
               sA(Ac) = sA(Ac) + w*Nb(a)
            END DO

            DEALLOCATE(N, Nxi, xl, ul)
         END DO
      END DO

      nGP = cm%reduce(nGP)
      i   = nEl*nG
      IF (nGP .NE. i) THEN
         err = " Error: found only "//STR(nGP)//" integration point "//
     2      "traces out of "//STR(i)
      END IF
      DEALLOCATE(ePtr, lD, Nb, xbl, ubl)

      RETURN
      END SUBROUTINE IB_SETFBFL
!--------------------------------------------------------------------
      SUBROUTINE IB_SETFBFSL(lBc, lM, Yg, Dg, sF, sA)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)
      REAL(KIND=8), INTENT(INOUT) :: sA(tnNo), sF(nsd,tnNo)

      LOGICAL :: flag
      INTEGER :: a, e, i, g, Ac, Ec, nNo, nEl, nG, eNoN, eNoNb, iM, nGP
      REAL(KIND=8) :: w, Jac, fb(nsd), xp(nsd), xi(nsd), nV(nsd),u(nsd),
     2   ub(nsd), tauF, gC(nsd,nsd-1)

      INTEGER, ALLOCATABLE :: ePtr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: lD(:,:), N(:), Nb(:), Nxi(:,:),
     2   xl(:,:), xbl(:,:), ul(:,:), ubl(:,:)

      nNo   = lM%nNo
      nEl   = lM%nEl
      nG    = lM%nG
      eNoNb = lM%eNoN

!     Get penalty constants
      tauF  = lBc%tauF

      ALLOCATE(lD(nsd,tnNo))
      lD = 0D0
      DO a=1, tnNo
         IF (mvMsh) lD(:,a) = Dg(nsd+2:2*nsd+1,a)
      END DO

!     Load all the ghost cell pointers into a global array
      ALLOCATE(ePtr(2,nG,nEl))
      ePtr = 0
      DO i=1, lM%trc%n
         e = lM%trc%gE(1,i)
         g = lM%trc%gE(2,i)
         ePtr(:,g,e) = lM%trc%ptr(:,i)
      END DO

      ALLOCATE(Nb(eNoNb), xbl(nsd,eNoNb), ubl(nsd,eNoNb))

!     Loop over all the elements and assemble to the residue
      nGP = 0
      DO e=1, nEl
         DO a=1, eNoNb
            Ac = lM%IEN(a,e)
            xbl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            ubl(:,a) = ib%Yo(1:nsd,Ac)
         END DO

         DO g=1, nG
            Ec = ePtr(1,g,e)
            iM = ePtr(2,g,e)
            IF (Ec .EQ. 0) CYCLE

!           Keep track of the total trace count
            nGP = nGP + 1

!           IB shell normal
            CALL GNNS(eNoNb, lM%Nx(:,:,g), xbl, nV, gC, gC)
            Jac = SQRT(NORM(nV))
!           Normal direction outwards from solid to fluid
            nV  = nV/Jac
            w   = lM%w(g) * Jac
            Nb  = lM%N(:,g)

!           Coordinates and IB velocity at the integration point
            xp = 0D0
            ub = 0D0
            DO a=1, eNoNb
               xp(:) = xp(:) + xbl(:,a)*Nb(a)
               ub(:) = ub(:) + ubl(:,a)*Nb(a)
            END DO

            cDmn = DOMAIN(msh(iM), cEq, Ec)
            eNoN = msh(iM)%eNoN
            ALLOCATE(N(eNoN), Nxi(nsd,eNoN), xl(nsd,eNoN), ul(nsd,eNoN))

!           Transfer to local arrays
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,Ec)
               xl(:,a) = x(:,Ac) + lD(:,Ac)
               ul(:,a) = Yg(1:nsd,Ac)
            END DO

!           Get the parametric coordinates and shape functions of the
!           trace w.r.t. background mesh
            xi = 0D0
            DO i=1, msh(iM)%nG
               xi = xi + msh(iM)%xi(:,i)
            END DO
            xi(:) = xi(:) / REAL(msh(iM)%nG,KIND=8)

            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge. "//
     2         "Invalid trace pointer"

!           Get the shape functions and their derivatives at the
!           integration point, as a function of the basis defined on the
!           background mesh
            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)
!           Check if the shape functions are correct
            i = 0
            DO a=1, eNoN
               IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
            END DO
            IF (i .NE. eNoN) err =
     2         "Error in computing IB trace shape functions"

!           Get the fluid velocity at the integration point
            u = 0D0
            DO a=1, eNoN
               u(:) = u(:) + N(a)*ul(:,a)
            END DO

!           Compute the feedback forcing with penalty forces along
!           normal direction
            IF (lBc%fbN) THEN
               fb(:) = tauF * NORM(u-ub, nV)
            ELSE
               fb(:) = tauF * (u(:)-ub(:))
            END IF

!           Transfer back to global vector
            DO a=1, eNoNb
               Ac = lM%IEN(a,e)
               sF(:,Ac) = sF(:,Ac) + w*Nb(a)*fb(:)
               sA(Ac) = sA(Ac) + w*Nb(a)
            END DO

            DEALLOCATE(N, Nxi, xl, ul)
         END DO
      END DO

      nGP = cm%reduce(nGP)
      i   = nEl*nG
      IF (nGP .NE. i) THEN
         err = " Error: found only "//STR(nGP)//" integration point "//
     2      "traces out of "//STR(i)
      END IF
      DEALLOCATE(ePtr, lD, Nb, xbl, ubl)

      RETURN
      END SUBROUTINE IB_SETFBFSL
!####################################################################
!     Immersed Finite Element Method (IFEM):
!     Two step projection: (a) Project acceleration and velocity from
!     the fluid mesh to the immersed bodies. (b) Compute FSI force on
!     the immersed bodies and project it to the background fluid mesh
      SUBROUTINE IB_GETFFSI(Ag, Yg, Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      REAL(KIND=8), ALLOCATABLE :: Ub(:,:)

      INTEGER a, i, iEq

      ib%cEq = 0
      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_fluid) THEN
            ib%cEq = iEq
            EXIT
         END IF
      END DO
      IF (ib%cEq.EQ.0 .OR. ib%mthd.NE.ibMthd_IFEM) RETURN

      ib%callD(1) = CPUT()
!     First, project fluid variables (velocity, pressure, acceleration)
!     at previous time step from fluid mesh to IB solid mesh
      CALL IB_PROJFVAR(Ag, Yg, Dg, ib%Ao, ib%Yo, ib%Uo, ib%callD(2))

!     Update IB solid displacement from the projected solid velocity
      ALLOCATE(Ub(nsd,ib%tnNo))
      Ub = 0D0
      DO a=1, ib%tnNo
         DO i=1, nsd
            Ub(i,a) = ib%Uo(i,a) + ib%Yo(i,a)*dt
         END DO
      END DO

!     Compute FSI force on IB and project it to the fluid mesh
      CALL IB_PROJFFSI(Ub, Dg)

!     Use new displacement field data to update tracer pointers
      ib%Uo = Ub
      ib%callD(3) = CPUT()
      CALL IB_UPDATE(Dg)
      ib%callD(3) = CPUT() - ib%callD(3)
      DEALLOCATE(Ub)

      ib%callD(1) = CPUT() - ib%callD(1)

      RETURN
      END SUBROUTINE IB_GETFFSI
!--------------------------------------------------------------------
!     Poject fluid velocity, pressure and acceleration onto IB
      SUBROUTINE IB_PROJFVAR(Ag, Yg, Dg, Ab, Yb, Ub, cmtime)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo), Ub(nsd,ib%tnNo)
      REAL(KIND=8), INTENT(OUT) :: Ab(nsd,ib%tnNo), Yb(nsd+1,ib%tnNo),
     2   cmtime

      LOGICAL :: flag
      INTEGER :: a, b, e, g, i, is, ie, Ac, Bc, Ec, iM, jM, eNoN, eNoNb
      REAL(KIND=8) :: w, Jac, xp(nsd), xi(nsd), ap(nsd), yp(nsd+1),
     2   Ks(nsd,nsd)

      REAL(KIND=8), ALLOCATABLE :: N(:), Nb(:), Nx(:,:), Nbx(:,:),
     2   xl(:,:), xbl(:,:), al(:,:), yl(:,:), sA(:)

!     Note that the equation physics is either Fluid or FSI
      is = eq(ib%cEq)%s
      ie = eq(ib%cEq)%e
      IF (ie .EQ. nsd+1) ie = ie - 1

!     Use L2 projection with mass lumping to project flow variables from
!     background fluid mesh to IB
      ALLOCATE(sA(ib%tnNo))
      sA = 0D0
      Ab = 0D0
      Yb = 0D0
!     Loop over each IB mesh
      DO iM=1, ib%nMsh
         eNoNb = ib%msh(iM)%eNoN
         ALLOCATE(Nb(eNoNb), Nbx(nsd,eNoNb), xbl(nsd,eNoNb))
!     Loop over each trace, as we need to first interpolate flow var at
!     the IB integration points based on its trace
         DO i=1, ib%msh(iM)%trc%n
            e  = ib%msh(iM)%trc%gE(1,i)
            g  = ib%msh(iM)%trc%gE(2,i)
            Ec = ib%msh(iM)%trc%ptr(1,i)
            jM = ib%msh(iM)%trc%ptr(2,i)

!        Transfer to local arrays: IB mesh variables
            Nb  = ib%msh(iM)%N(:,g)
            DO a=1, eNoNb
               Ac = ib%msh(iM)%IEN(a,e)
               xbl(:,a) = ib%x(:,Ac) + Ub(:,Ac)
            END DO
            CALL GNN(eNoNb, ib%msh(iM)%Nx(:,:,g), xbl, Nbx, Jac, Ks)
            IF (ISZERO(Jac)) err = " Jac < 0 @ element "//e
            w = ib%msh(iM)%w(g) * Jac

!        Coordinates of the integration point
            xp = 0D0
            DO a=1, eNoNb
               xp = xp + Nb(a)*xbl(:,a)
            END DO

!        Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN), al(nsd,eNoN),
     2         yl(nsd+1,eNoN))
            al = 0D0
            yl = 0D0
            DO b=1, eNoN
               Bc = msh(jM)%IEN(b,Ec)
               al(:,b) = Ag(is:ie,Bc)
               yl(:,b) = Yg(is:ie+1,Bc)
               xl(:,b) = x(:,Bc)
               IF (mvMsh) xl(:,b) = xl(:,b) + Dg(ie+2:ie+nsd+1,Bc)
            END DO

!        Initialize parameteric coordinate for Newton's iterations
            xi = 0D0
            DO a=1, msh(jM)%nG
               xi = xi + msh(jM)%xi(:,a)
            END DO
            xi = xi / REAL(msh(jM)%nG, KIND=8)

            CALL GETXI(msh(jM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge "//
     2         "Invalid trace pointer"

            CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)
!        Make sure that the node lies within the element
            a = 0
            DO b=1, eNoN
               IF (N(b).GT.-1E-4 .AND. N(b).LT.1.0001D0) a = a + 1
            END DO
            IF (a .NE. eNoN) err =" IB tracer pointing to incorrect "//
     2         "fluid element (PROJFVAR)"

!        Use the computed shape functions to interpolate flow var at the
!        IB integration point
            ap = 0D0
            yp = 0D0
            DO b=1, eNoN
               Bc = msh(jM)%IEN(b,Ec)
               ap = ap + N(b)*al(:,b)
               yp = yp + N(b)*yl(:,b)
            END DO

!        Project flow variables to IB nodes
            DO a=1, eNoNb
               Ac = ib%msh(iM)%IEN(a,e)
               Ab(:,Ac) = Ab(:,Ac) + w*Nb(a)*ap(:)
               Yb(:,Ac) = Yb(:,Ac) + w*Nb(a)*yp(:)
               sA(Ac)   = sA(Ac)   + w*Nb(a)
            END DO
            DEALLOCATE(N, Nx, xl, al, yl)
         END DO
         DEALLOCATE(Nb, Nbx, xbl)
      END DO

!     Synchronize Ab, Yb across all the processes
      cmtime = CPUT()
      CALL IB_SYNC(Ab)
      CALL IB_SYNC(Yb)
      CALL IB_SYNC(sA)
      cmtime = CPUT() - cmtime

      DO a=1, ib%tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            Ab(:,a) = Ab(:,a) / sA(a)
            Yb(:,a) = Yb(:,a) / sA(a)
         END IF
      END DO

      DEALLOCATE(sA)

      RETURN
      END SUBROUTINE IB_PROJFVAR
!--------------------------------------------------------------------
!     Compute FSI force on the IB and project it to the background mesh.
!     Note that here we compute solid stresses using updated displacement
!     field. But for projection, we still use old displacement to
!     identify traces to satisfy adjoint property and conserve KE
      SUBROUTINE IB_PROJFFSI(Ub, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Ub(nsd,ib%tnNo), Dg(tDof,tnNo)

      LOGICAL :: flag
      INTEGER :: a, b, e, g, i, Ac, Bc, Ec, iM, jM, eNoN, eNoNb, cPhys
      REAL(KIND=8) :: w, Jac, ksix(nsd,nsd), xp(nsd), xi(nsd), fp(nsd)

      INTEGER, ALLOCATABLE :: incNd(:)
      REAL(KIND=8), ALLOCATABLE :: sA(:), N(:), Nb(:), Nx(:,:), al(:,:),
     2   yl(:,:), ul(:,:), xl(:,:), fNl(:,:), lR(:,:), gR(:,:)

!     TODO: This is temporary. Later we get domain based on each IB node
!     and communicate across IB process boundaries
      cDmn = DOMAIN(msh(1), ib%cEq, 1)

!     We assemble the fluid and structure contribution on the IB nodes.
!     Later we compute the FSI force at the IB integration point and
!     project it to background mesh
      ib%R = 0D0
      DO iM=1, ib%nMsh
         eNoN = ib%msh(iM)%eNoN
         ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN), al(nsd,eNoN),
     2      yl(nsd+1,eNoN), ul(nsd,eNoN), fNl(ib%nFn*nsd,eNoN),
     3      lR(nsd,eNoN))
         DO e=1, ib%msh(iM)%nEl
            DO a=1, eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               al(:,a)  = ib%Ao(:,Ac)
               yl(:,a)  = ib%Yo(:,Ac)
               ul(:,a)  = Ub(:,Ac)
               xl(:,a)  = ib%x(:,Ac)
               IF (ALLOCATED(ib%fN)) fNl(:,a) = ib%fN(:,Ac)
            END DO

!           Update IB domain and physics
            ib%cDmn = IB_DOMAIN(ib%msh(iM), ib%cEq, e)
            cPhys   = eq(ib%cEq)%dmnIB(ib%cDmn)%phys
            IF (cPhys .NE. phys_struct) err = "Only struct phys is "//
     2         "allowed on IB using IFEM method"

!           Update shape functions if using NURBS
            IF (ib%msh(iM)%eType .EQ. eType_NRB) THEN
               CALL NRBNNX(ib%msh(iM), e)
            END IF

            lR = 0D0
            DO g=1, ib%msh(iM)%nG
               IF (g.EQ.1 .OR. .NOT.ib%msh(iM)%lShpF) THEN
                  CALL GNN(eNoN, ib%msh(iM)%Nx(:,:,g), xl, Nx, Jac,
     2               ksix)
               END IF
               IF (ISZERO(Jac)) err = " Jac < 0 @ element "//e
               N  = ib%msh(iM)%N(:,g)
               w  = ib%msh(iM)%w(g) * Jac

!              Compute FSI force
               IF (nsd .EQ. 3) THEN
                  CALL IB_FFSI3D(eNoN, w, N, Nx, al, yl, ul, fNl, lR)
               ELSE
                  CALL IB_FFSI2D(eNoN, w, N, Nx, al, yl, ul, fNl, lR)
               END IF
            END DO ! g

!           Add the FSI force residue contribution to the global vector
            DO a=1, eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               ib%R(:,Ac) = ib%R(:,Ac) + lR(:,a)
            END DO ! a
         END DO ! e
         DEALLOCATE(N, Nx, xl, al, yl, ul, fNl, lR)
      END DO ! iM

!     Project ib%R to background fluid mesh using the same domain of
!     influence that was used to project flow variables.
      ALLOCATE(sA(tnNo), incNd(tnNo))
      sA  = 0D0
      Rib = 0D0
      incNd = 0
      DO iM=1, ib%nMsh
         eNoNb = ib%msh(iM)%eNoN
         ALLOCATE(Nb(eNoNb))
         DO i=1, ib%msh(iM)%trc%n
            e  = ib%msh(iM)%trc%gE(1,i)
            g  = ib%msh(iM)%trc%gE(2,i)
            Ec = ib%msh(iM)%trc%ptr(1,i)
            jM = ib%msh(iM)%trc%ptr(2,i)

            xp = 0D0
            fp = 0D0
            Nb = ib%msh(iM)%N(:,g)
            DO a=1, eNoNb
               Ac = ib%msh(iM)%IEN(a,e)
               xp = xp + Nb(a)*(ib%x(:,Ac) + ib%Uo(:,Ac))
               fp = fp + Nb(a)*ib%R(:,Ac)
            END DO

            eNoN  = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN))
            DO b=1, eNoN
               Bc = msh(jM)%IEN(b,Ec)
               xl(:,b) = x(:,Bc)
               IF (mvMsh) xl(:,b) = xl(:,b) + Dg(nsd+2:2*nsd+1,Bc)
            END DO

            xi = 0D0
            DO g=1, msh(jM)%nG
               xi = xi + msh(jM)%xi(:,g)
            END DO
            xi = xi / REAL(msh(jM)%nG, KIND=8)

            CALL GETXI(msh(jM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge "//
     2         "Invalid trace pointer"

            CALL GETGNN(nsd, msh(jM)%eType, eNoN, xi, N, Nx)
!        Make sure that the node lies within the element
            a = 0
            DO b=1, eNoN
               IF (N(b).GT.-1E-4 .AND. N(b).LT.1.0001D0) a = a + 1
            END DO
            IF (a .NE. eNoN) err =" IB tracer pointing to incorrect "//
     2         "fluid element (PROJFFSI)"

            DO b=1, eNoN
               Bc = msh(jM)%IEN(b,Ec)
               incNd(Bc) = 1
               sA(Bc) = sA(Bc) + N(b)
               Rib(:,Bc) = Rib(:,Bc) + N(b)*fp(:)
            END DO
            DEALLOCATE(N, Nx, xl)
         END DO
         DEALLOCATE(Nb)
      END DO

!     Synchronize Rib across process interfaces
      CALL COMMU(Rib)
      CALL COMMU(sA)

      DO a=1, tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            Rib(:,a) = Rib(:,a) / sA(a)
            incNd(a) = 1
         END IF
      END DO
      DEALLOCATE(sA)

c      CALL DEBUGIBR(incNd)

      RETURN
      END SUBROUTINE IB_PROJFFSI
!--------------------------------------------------------------------
!     Compute the 3D FSI force due to IB in reference configuration
      SUBROUTINE IB_FFSI3D(eNoN, w, N, Nx, al, yl, ul, fNl, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), al(3,eNoN),
     2   yl(4,eNoN), ul(3,eNoN), fNl(ib%nFn*3,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(3,eNoN)

      INTEGER :: a, b, iEq, iDmn, iFn
      REAL(KIND=8) :: wl, Jac, rho_s, rho_f, eM_s, nu_s, mu_f, bf(3),
     2   vd(3), v(3), rM(3), vx(3,3), F(3,3), Fi(3,3), S(3,3), P(3,3),
     3   NxFi(3,eNoN), VxFi(3,3), vVx(3), VxNx(3), NxP(3), CC(3,3,3,3),
     4   fl(3,ib%nFn)
      TYPE(stModelType) :: stModel

!     Define IB struct parameters
      iEq     = ib%cEq
      iDmn    = ib%cDmn
      stModel = eq(iEq)%dmnIB(iDmn)%stM
      rho_s   = eq(iEq)%dmnIB(iDmn)%prop(solid_density)
      eM_s    = eq(iEq)%dmnIB(iDmn)%prop(elasticity_modulus)
      nu_s    = eq(iEq)%dmnIB(iDmn)%prop(poisson_ratio)

!     Define fluid parameters
      rho_f   = eq(iEq)%dmn(cDmn)%prop(fluid_density)
      mu_f    = eq(iEq)%dmn(cDmn)%prop(viscosity)

!     Body force
      bf(1)   = eq(iEq)%dmnIB(iDmn)%prop(f_x)
      bf(2)   = eq(iEq)%dmnIB(iDmn)%prop(f_y)
      bf(3)   = eq(iEq)%dmnIB(iDmn)%prop(f_z)

!     Inertia, body force and deformation tensor (F)
      vd      = -bf
      v       = 0D0
      vx      = 0D0
      fl      = 0D0
      F       = 0D0
      F(1,1)  = 1D0
      F(2,2)  = 1D0
      F(3,3)  = 1D0
      DO a=1, eNoN
!        Acceleration (dv_i/dt)
         vd(1) = vd(1) + N(a)*al(1,a)
         vd(2) = vd(2) + N(a)*al(2,a)
         vd(3) = vd(3) + N(a)*al(3,a)

!        Velocity, v
         v(1) = v(1) + N(a)*yl(1,a)
         v(2) = v(2) + N(a)*yl(2,a)
         v(3) = v(3) + N(a)*yl(3,a)

!        Grad v (v_i,I) w.r.t. reference coordinates
         vx(1,1) = vx(1,1) + yl(1,a)*Nx(1,a)
         vx(1,2) = vx(1,2) + yl(1,a)*Nx(2,a)
         vx(1,3) = vx(1,3) + yl(1,a)*Nx(3,a)

         vx(2,1) = vx(2,1) + yl(2,a)*Nx(1,a)
         vx(2,2) = vx(2,2) + yl(2,a)*Nx(2,a)
         vx(2,3) = vx(2,3) + yl(2,a)*Nx(3,a)

         vx(3,1) = vx(3,1) + yl(3,a)*Nx(1,a)
         vx(3,2) = vx(3,2) + yl(3,a)*Nx(2,a)
         vx(3,3) = vx(3,3) + yl(3,a)*Nx(3,a)

         DO iFn=1, ib%nFn
            b = 3*(iFn-1)
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
            fl(3,iFn) = fl(3,iFn) + N(a)*fNl(b+3,a)
         END DO

!        Deformation tensor, F_{iI}. Here the shape function derivatives
!        are w.r.t. the reference coordinates
         F(1,1)  = F(1,1)  + ul(1,a)*Nx(1,a)
         F(1,2)  = F(1,2)  + ul(1,a)*Nx(2,a)
         F(1,3)  = F(1,3)  + ul(1,a)*Nx(3,a)

         F(2,1)  = F(2,1)  + ul(2,a)*Nx(1,a)
         F(2,2)  = F(2,2)  + ul(2,a)*Nx(2,a)
         F(2,3)  = F(2,3)  + ul(2,a)*Nx(3,a)

         F(3,1)  = F(3,1)  + ul(3,a)*Nx(1,a)
         F(3,2)  = F(3,2)  + ul(3,a)*Nx(2,a)
         F(3,3)  = F(3,3)  + ul(3,a)*Nx(3,a)
      END DO
      Jac = MAT_DET(F, nsd)
      Fi  = MAT_INV(F, nsd)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
      CALL GETPK2CC(stModel, F, ib%nFn, fl, S, CC)

!     1st Piola-Kirchhoff tensor (P)
      P(1,1) = F(1,1)*S(1,1) + F(1,2)*S(2,1) + F(1,3)*S(3,1)
      P(1,2) = F(1,1)*S(1,2) + F(1,2)*S(2,2) + F(1,3)*S(3,2)
      P(1,3) = F(1,1)*S(1,3) + F(1,2)*S(2,3) + F(1,3)*S(3,3)
      P(2,1) = F(2,1)*S(1,1) + F(2,2)*S(2,1) + F(2,3)*S(3,1)
      P(2,2) = F(2,1)*S(1,2) + F(2,2)*S(2,2) + F(2,3)*S(3,2)
      P(2,3) = F(2,1)*S(1,3) + F(2,2)*S(2,3) + F(2,3)*S(3,3)
      P(3,1) = F(3,1)*S(1,1) + F(3,2)*S(2,1) + F(3,3)*S(3,1)
      P(3,2) = F(3,1)*S(1,2) + F(3,2)*S(2,2) + F(3,3)*S(3,2)
      P(3,3) = F(3,1)*S(1,3) + F(3,2)*S(2,3) + F(3,3)*S(3,3)

!     grad N_A (N_A,j) w.r.t. spatial coordinates
      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1) + Nx(3,a)*Fi(3,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2) + Nx(3,a)*Fi(3,2)
         NxFi(3,a) = Nx(1,a)*Fi(1,3) + Nx(2,a)*Fi(2,3) + Nx(3,a)*Fi(3,3)
      END DO

!     grad v (v_i,j)
      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1) + vx(1,3)*Fi(3,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2) + vx(1,3)*Fi(3,2)
      VxFi(1,3) = vx(1,1)*Fi(1,3) + vx(1,2)*Fi(2,3) + vx(1,3)*Fi(3,3)

      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1) + vx(2,3)*Fi(3,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2) + vx(2,3)*Fi(3,2)
      VxFi(2,3) = vx(2,1)*Fi(1,3) + vx(2,2)*Fi(2,3) + vx(2,3)*Fi(3,3)

      VxFi(3,1) = vx(3,1)*Fi(1,1) + vx(3,2)*Fi(2,1) + vx(3,3)*Fi(3,1)
      VxFi(3,2) = vx(3,1)*Fi(1,2) + vx(3,2)*Fi(2,2) + vx(3,3)*Fi(3,2)
      VxFi(3,3) = vx(3,1)*Fi(1,3) + vx(3,2)*Fi(2,3) + vx(3,3)*Fi(3,3)

!     V. grad V
      vVx(1) = v(1)*VxFi(1,1) + v(2)*VxFi(1,2) + v(3)*VxFi(1,3)
      vVx(2) = v(1)*VxFi(2,1) + v(2)*VxFi(2,2) + v(3)*VxFi(2,3)
      vVx(3) = v(1)*VxFi(3,1) + v(2)*VxFi(3,2) + v(3)*VxFi(3,3)

      wl     = w * Jac
      DO a=1, eNoN
!        Viscous stress due to fluid
!        grad N_A . (grad v + grad v^T)
         VxNx(1) = NxFi(1,a)*(VxFi(1,1) + VxFi(1,1)) +
     2             NxFi(2,a)*(VxFi(1,2) + VxFi(2,1)) +
     3             NxFi(3,a)*(VxFi(1,3) + VxFi(3,1))

         VxNx(2) = NxFi(1,a)*(VxFi(2,1) + VxFi(1,2)) +
     2             NxFi(2,a)*(VxFi(2,2) + VxFi(2,2)) +
     3             NxFi(3,a)*(VxFi(2,3) + VxFi(3,2))

         VxNx(3) = NxFi(1,a)*(VxFi(3,1) + VxFi(1,3)) +
     2             NxFi(2,a)*(VxFi(3,2) + VxFi(2,3)) +
     3             NxFi(3,a)*(VxFi(3,3) + VxFi(3,3))

!        Sress on structure: NxP
         NxP(1) = Nx(1,a)*P(1,1) + Nx(2,a)*P(1,2) + Nx(3,a)*P(1,3)
         NxP(2) = Nx(1,a)*P(2,1) + Nx(2,a)*P(2,2) + Nx(3,a)*P(2,3)
         NxP(3) = Nx(1,a)*P(3,1) + Nx(2,a)*P(3,2) + Nx(3,a)*P(3,3)

         rM(1) = (rho_f-rho_s)*N(a)*(vd(1) + vVx(1)) + mu_f*VxNx(1)
         rM(2) = (rho_f-rho_s)*N(a)*(vd(2) + vVx(2)) + mu_f*VxNx(2)
         rM(3) = (rho_f-rho_s)*N(a)*(vd(3) + vVx(3)) + mu_f*VxNx(3)

         lR(1,a) = lR(1,a) + wl*rM(1) - w*NxP(1)
         lR(2,a) = lR(2,a) + wl*rM(2) - w*NxP(2)
         lR(3,a) = lR(3,a) + wl*rM(3) - w*NxP(3)
      END DO

      RETURN
      END SUBROUTINE IB_FFSI3D
!--------------------------------------------------------------------
!     Compute the 2D FSI force due to IB in reference configuration
      SUBROUTINE IB_FFSI2D(eNoN, w, N, Nx, al, yl, ul, fNl, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), al(2,eNoN),
     2   yl(3,eNoN), ul(2,eNoN), fNl(2*ib%nFn,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(2,eNoN)

      INTEGER :: a, b, iEq, iDmn, iFn
      REAL(KIND=8) :: wl, Jac, rho_s, rho_f, eM_s, nu_s, mu_f, bf(2),
     2   vd(2), v(2), rM(2), vx(2,2), F(2,2), Fi(2,2), S(2,2), P(2,2),
     3   NxFi(2,eNoN), VxFi(2,2), vVx(2), VxNx(2), NxP(2), CC(2,2,2,2),
     4   fl(2,ib%nFn)
      TYPE(stModelType) :: stModel

!     Define IB struct parameters
      iEq     = ib%cEq
      iDmn    = ib%cDmn
      stModel = eq(iEq)%dmnIB(iDmn)%stM
      rho_s   = eq(iEq)%dmnIB(iDmn)%prop(solid_density)
      eM_s    = eq(iEq)%dmnIB(iDmn)%prop(elasticity_modulus)
      nu_s    = eq(iEq)%dmnIB(iDmn)%prop(poisson_ratio)

!     Define fluid parameters
      rho_f   = eq(iEq)%dmn(cDmn)%prop(fluid_density)
      mu_f    = eq(iEq)%dmn(cDmn)%prop(viscosity)

!     Body force
      bf(1)   = eq(iEq)%dmnIB(iDmn)%prop(f_x)
      bf(2)   = eq(iEq)%dmnIB(iDmn)%prop(f_y)

!     Inertia, body force and deformation tensor (F)
      vd      = -bf
      v       = 0D0
      vx      = 0D0
      fl      = 0D0
      F       = 0D0
      F(1,1)  = 1D0
      F(2,2)  = 1D0
      DO a=1, eNoN
!        Acceleration (dv_i/dt)
         vd(1) = vd(1) + N(a)*al(1,a)
         vd(2) = vd(2) + N(a)*al(2,a)

!        Velocity, v
         v(1) = v(1) + N(a)*yl(1,a)
         v(2) = v(2) + N(a)*yl(2,a)

!        Grad v (v_i,I) w.r.t. reference coordinates
         vx(1,1) = vx(1,1) + yl(1,a)*Nx(1,a)
         vx(1,2) = vx(1,2) + yl(1,a)*Nx(2,a)
         vx(2,1) = vx(2,1) + yl(2,a)*Nx(1,a)
         vx(2,2) = vx(2,2) + yl(2,a)*Nx(2,a)

         DO iFn=1, ib%nFn
            b = 2*(iFn-1)
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
         END DO

!        Deformation tensor, F_{iI}. Here the shape function derivatives
!        are w.r.t. the reference coordinates
         F(1,1)  = F(1,1)  + ul(1,a)*Nx(1,a)
         F(1,2)  = F(1,2)  + ul(1,a)*Nx(2,a)
         F(2,1)  = F(2,1)  + ul(2,a)*Nx(1,a)
         F(2,2)  = F(2,2)  + ul(2,a)*Nx(2,a)
      END DO
      Jac = MAT_DET(F, nsd)
      Fi  = MAT_INV(F, nsd)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
      CALL GETPK2CC(stModel, F, ib%nFn, fl, S, CC)

!     1st Piola-Kirchhoff tensor (P)
      P(1,1) = F(1,1)*S(1,1) + F(1,2)*S(2,1)
      P(1,2) = F(1,1)*S(1,2) + F(1,2)*S(2,2)
      P(2,1) = F(2,1)*S(1,1) + F(2,2)*S(2,1)
      P(2,2) = F(2,1)*S(1,2) + F(2,2)*S(2,2)

!     grad N_A (N_A,j) w.r.t. spatial coordinates
      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2)
      END DO

!     grad v (v_i,j)
      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2)
      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2)

!     V. grad V
      vVx(1) = v(1)*VxFi(1,1) + v(2)*VxFi(1,2)
      vVx(2) = v(1)*VxFi(2,1) + v(2)*VxFi(2,2)

      wl     = w * Jac
      DO a=1, eNoN
!        Viscous stress due to fluid
!        grad N_A . (grad v + grad v^T)
         VxNx(1) = NxFi(1,a)*(VxFi(1,1) + VxFi(1,1)) +
     2             NxFi(2,a)*(VxFi(1,2) + VxFi(2,1))

         VxNx(2) = NxFi(1,a)*(VxFi(2,1) + VxFi(1,2)) +
     2             NxFi(2,a)*(VxFi(2,2) + VxFi(2,2))

!        Sress on structure: NxP
         NxP(1) = Nx(1,a)*P(1,1) + Nx(2,a)*P(1,2)
         NxP(2) = Nx(1,a)*P(2,1) + Nx(2,a)*P(2,2)

         rM(1) = (rho_f-rho_s)*N(a)*(vd(1) + vVx(1)) + mu_f*VxNx(1)
         rM(2) = (rho_f-rho_s)*N(a)*(vd(2) + vVx(2)) + mu_f*VxNx(2)

         lR(1,a) = lR(1,a) + wl*rM(1) - w*NxP(1)
         lR(2,a) = lR(2,a) + wl*rM(2) - w*NxP(2)
      END DO

      RETURN
      END SUBROUTINE IB_FFSI2D
!####################################################################
!     Write IB call duration
      SUBROUTINE IB_OUTR()
      USE COMMOD
      IMPLICIT NONE

      CHARACTER(LEN=stdL) sOut

      std = REPEAT("-",55)
      WRITE(sOut,'(F6.2)') ib%callD(1)
      WRITE(sOut,'(A)') " IB call duration: "//TRIM(sOut)//' sec'
      WRITE(sOut,'(A)') TRIM(sOut)//" (comm."//
     2   STR(NINT(1D2*ib%callD(2)/ib%callD(1)),3)//"%)"
      WRITE(sOut,'(A)') TRIM(sOut)//", (updt."//
     2   STR(NINT(1D2*ib%callD(3)/ib%callD(1)),3)//"%)"
      std = sOut
      std = REPEAT("-",55)

      RETURN
      END SUBROUTINE IB_OUTR
!####################################################################
!     Debugs FSI force on the IB
      SUBROUTINE DEBUGIBR(incNd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: incNd(ib%tnNo)

      INTEGER :: a, i, Ac, iM, fid, ierr
      REAL(KIND=8) :: s, lo, hi, av
      CHARACTER(LEN=stdL) :: fName

      INTEGER, ALLOCATABLE :: lI(:), gI(:)
      REAL(KIND=8), ALLOCATABLE :: lR(:,:), gR(:,:)

!     DEBUG ib%R
      fid = 1289
      IF (cm%mas()) THEN
         WRITE(fName,'(A)') TRIM(appPath)//"dbg_ibR_hist.dat"
         IF (cTS .EQ. 1) THEN
            OPEN(fid,FILE=TRIM(fName))
            CLOSE(fid,STATUS='DELETE')
         END IF
         OPEN(fid,FILE=TRIM(fName),POSITION='APPEND')
         WRITE(fid,'(A)',ADVANCE='NO') STR(cTS)
         av = 0D0
         lo = 1000000.0D0
         hi = -999999.0D0

         DO a=1, ib%tnNo
            s = 0D0
            DO i=1, nsd
               s = s + ib%R(i,a)**2
            END DO
            s = SQRT(s)
            av = av + s
            IF (s .LT. lo) lo = s
            IF (s .GT. hi) hi = s
         END DO
         av = av / REAL(ib%tnNo, KIND=8)
         WRITE(fid,'(A)',ADVANCE='NO') " "//STR(av)//" "//
     2      STR(hi-lo)
      END IF

!     DEBUG Rib
      DO iM=1, nMsh
         ALLOCATE(lR(nsd,msh(iM)%nNo), lI(msh(iM)%nNo))
         IF (cm%mas()) THEN
            ALLOCATE(gR(nsd,msh(iM)%gnNo), gI(msh(iM)%gnNo))
         ELSE
            ALLOCATE(gR(0,0), gI(0))
         END IF
         lR = 0D0
         lI = 0
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            lR(:,a) = Rib(:,Ac)
            lI(a)   = incNd(Ac)
         END DO
         gR = GLOBAL(msh(iM), lR)
         gI = GLOBAL(msh(iM), lI)
         DEALLOCATE(lR, lI)

         IF (cm%mas()) THEN
            av = 0D0
            lo = 10000000.0D0
            hi = -9999999.0D0
            DO a=1, msh(iM)%gnNo
               IF (gI(a) .EQ. 0) CYCLE
               s = 0D0
               DO i=1, nsd
                  s = s + gR(i,a)**2
               END DO
               s = SQRT(s)
               av = av + s
               IF (s .LT. lo) lo = s
               IF (s .GT. hi) hi = s
            END DO
            av = av / REAL(SUM(gI(:)), KIND=8)
            WRITE(fid,'(A)') " "//STR(av)//" "//STR(hi-lo)
            CLOSE(fid)
         END IF
         DEALLOCATE(gR, gI)
      END DO

      iM = cm%reduce(iM)

      RETURN
      END SUBROUTINE DEBUGIBR
!--------------------------------------------------------------------
!     Debugs IB mesh traces
      SUBROUTINE DEBUGIBMSHTRC(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM

      INTEGER :: i, a, e, g, n, Ac, Ec, iM, fid, ierr
      CHARACTER(LEN=stdL) :: fName
      REAL(KIND=8) xp(nsd)

      INTEGER, ALLOCATABLE :: sCount(:), disps(:), lE(:), gE(:),
     2   eptr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:)

      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lM%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      n = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(lE(4*lM%trc%n), gE(4*n))
      DO i=1, lM%trc%n
         lE(4*i-3) = lM%trc%gE(1,i)
         lE(4*i-2) = lM%trc%gE(2,i)
         lE(4*i-1) = lM%trc%ptr(1,i)
         lE(4*i)   = lM%trc%ptr(2,i)
      END DO

      i = 4*lM%trc%n
      CALL MPI_GATHERV(lE, i, mpint, gE, sCount, disps, mpint, master,
     2   cm%com(), ierr)

      DEALLOCATE(lE, disps, sCount)

      WRITE(fName,'(A)') TRIM(appPath)//"dbg_ib_trc_"//STR(cTS)//".dat"
      IF (cm%mas()) THEN
         ALLOCATE(eptr(2,lM%nG,lM%nEl))
         eptr = 0
         DO i=1, n
            e  = gE(4*i-3)
            g  = gE(4*i-2)
            Ec = gE(4*i-1)
            iM = gE(4*i)
            ePtr(1,g,e) = Ec
            ePtr(2,g,e) = iM
         END DO

         fid = 1289
         OPEN(fid,FILE=TRIM(fName))
         WRITE(fid,'(A)') "List of failed traces on mesh: "//
     2      TRIM(lM%name)
         ALLOCATE(xl(nsd,lM%eNoN))
         DO e=1, lM%nEl
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lM%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
                  xp = 0D0
                  DO a=1, lM%eNoN
                     xp(:) = xp(:) + lM%N(a,g)*xl(:,a)
                  END DO
                  WRITE(fid,'(2X,A)',ADVANCE='NO') STR(g)//" "//STR(e)
                  DO i=1, nsd
                     WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
                  END DO
                  WRITE(fid,'(A)')
               END IF
            END DO
         END DO
         DEALLOCATE(xl)
         CLOSE(fid)
      ELSE
         ALLOCATE(ePtr(0,0,0))
      END IF
      DEALLOCATE(gE, ePtr)

      i = cm%reduce(i)
      err = "ERROR: Failed to detect all the traces on "//
     2   TRIM(lM%name)//". See "//TRIM(fName)//" for more information."

      RETURN
      END SUBROUTINE DEBUGIBMSHTRC
!--------------------------------------------------------------------
!     Debugs IB face traces
      SUBROUTINE DEBUGIBFATRC(lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER :: i, a, e, g, n, Ac, Ec, iM, fid, ierr
      CHARACTER(LEN=stdL) :: fName
      REAL(KIND=8) xp(nsd)

      INTEGER, ALLOCATABLE :: sCount(:), disps(:), lE(:), gE(:),
     2   eptr(:,:,:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:)

      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lFa%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      n = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(lE(4*lFa%trc%n), gE(4*n))
      DO i=1, lFa%trc%n
         lE(4*i-3) = lFa%trc%gE(1,i)
         lE(4*i-2) = lFa%trc%gE(2,i)
         lE(4*i-1) = lFa%trc%ptr(1,i)
         lE(4*i)   = lFa%trc%ptr(2,i)
      END DO

      i = 4*lFa%trc%n
      CALL MPI_GATHERV(lE, i, mpint, gE, sCount, disps, mpint, master,
     2   cm%com(), ierr)

      DEALLOCATE(lE, disps, sCount)

      WRITE(fName,'(A)') TRIM(appPath)//"dbg_ib_trc_"//STR(cTS)//".dat"
      IF (cm%mas()) THEN
         ALLOCATE(eptr(2,lFa%nG,lFa%nEl))
         eptr = 0
         DO i=1, n
            e  = gE(4*i-3)
            g  = gE(4*i-2)
            Ec = gE(4*i-1)
            iM = gE(4*i)
            ePtr(1,g,e) = Ec
            ePtr(2,g,e) = iM
         END DO

         fid = 1289
         OPEN(fid,FILE=TRIM(fName))
         WRITE(fid,'(A)') "List of failed traces on mesh: "//
     2      TRIM(lFa%name)
         ALLOCATE(xl(nsd,lFa%eNoN))
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
            END DO
            DO g=1, lFa%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
                  xp = 0D0
                  DO a=1, lFa%eNoN
                     xp(:) = xp(:) + lFa%N(a,g)*xl(:,a)
                  END DO
                  WRITE(fid,'(2X,A)',ADVANCE='NO') STR(g)//" "//STR(e)
                  DO i=1, nsd
                     WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
                  END DO
                  WRITE(fid,'(A)')
               END IF
            END DO
         END DO
         DEALLOCATE(xl)
         CLOSE(fid)
      ELSE
         ALLOCATE(ePtr(0,0,0))
      END IF
      DEALLOCATE(gE, ePtr)

      i = cm%reduce(i)
      err = "ERROR: Failed to detect all the traces on "//
     2   TRIM(lFa%name)//". See "//TRIM(fName)//" for more information."

      RETURN
      END SUBROUTINE DEBUGIBFATRC
!####################################################################
