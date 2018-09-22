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
!--------------------------------------------------------------------
!
!     This routine is desinged to read the mesh/es that may come in
!     several formats and setup all parameters related to meshes.
!
!--------------------------------------------------------------------

!     The higher level routine that calls other routines
      SUBROUTINE READMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list

      CHARACTER, PARAMETER :: dSym(3) = (/"X","Y","Z"/)

      LOGICAL :: flag
      INTEGER :: i, j, iM, iFa, a, b, Ac, e, lDof, lnNo
      REAL(KIND=8) :: maxX(nsd), minX(nsd), scaleF, fibN(nsd), rtmp
      CHARACTER(LEN=stdL) :: ctmp
      TYPE(listType), POINTER :: lPtr, lPM
      TYPE(stackType) :: avNds
      TYPE(fileType) :: fTmp

      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:), gX(:,:), tmpA(:,:),
     2   tmpY(:,:), tmpD(:,:)

      IF (.NOT.resetSim) THEN
         nMsh  = list%srch("Add mesh",ll=1)
         std = " Number of meshes: "//nMsh
         ALLOCATE (msh(nMsh), gX(0,0))

         gtnNo = 0
         minX  =  HUGE(minX)
         maxX  = -HUGE(minX)
         DO iM=1, nMsh
            lPM => list%get(msh(iM)%name,"Add mesh",iM)
            lPtr => lPM%get(msh(iM)%lShl,"Set mesh as shell")

            std  = " Reading mesh <"//CLR(TRIM(msh(iM)%name))//">"
            CALL READSV(lPM, msh(iM))
            IF (msh(iM)%eType .EQ. eType_NA) THEN
               CALL READCCNE(lPM, msh(iM))
            END IF
            IF (msh(iM)%eType .EQ. eType_NA) THEN
               CALL READNRB(lPM, msh(iM))
            END IF
            IF (msh(iM)%eType .EQ. eType_NA) THEN
               CALL READGAMBIT(lPM, msh(iM))
            END IF
            IF (msh(iM)%eType .EQ. eType_NA) THEN
               err = "Failed to identify format of the mesh"
            END IF

            IF (nMsh .LE. 3) THEN
               std = " Number of nodes: "//msh(iM)%gnNo
               std = " Number of elements: "//msh(iM)%gnEl
            END IF
            ALLOCATE (msh(iM)%gN(msh(iM)%gnNo))
            msh(iM)%gN = 0

!     Making sure face names are unique
            DO iFa=1, msh(iM)%nFa
               msh(iM)%fa(iFa)%iM = iM
               ctmp = msh(iM)%fa(iFa)%name
               DO i=1, iM
                  DO j=1, msh(i)%nFa
                     IF (ctmp .EQ. msh(i)%fa(j)%name .AND.
     2                  (i.NE.iM .OR. j.NE.iFa)) THEN
                        err = "Repeating face names is not allowed"
                     END IF
                  END DO
               END DO
            END DO

!     To scale the mesh, while attaching x to gX
            scaleF = 1D0
            lPtr => lPM%get(scaleF,"Mesh scale factor",lb=0D0)
            a = gtnNo + msh(iM)%gnNo
            IF (iM .GT. 1) THEN
               ALLOCATE(tmpX(nsd,gtnNo))
               tmpX = gX
               DEALLOCATE(gX)
               ALLOCATE(gX(nsd,a))
               gX(:,1:gtnNo) = tmpX
               DEALLOCATE(tmpX)
            ELSE
               DEALLOCATE(gX)
               ALLOCATE(gX(nsd,a))
            END IF
            gX(:,gtnNo+1:a) = msh(iM)%x * scaleF
            gtnNo           = a
            DEALLOCATE(msh(iM)%x)
         END DO
         ALLOCATE(x(nsd,gtnNo))
         x = gX

!        Checks for shell elements
         DO iM=1, nMsh
            IF (msh(iM)%lShl) THEN
               IF (msh(iM)%eType.NE.eType_NRB .AND.
     2             msh(iM)%eType.NE.eType_TRI) THEN
                  err = "Shell elements can be either triangles "//
     2               "or C1-NURBS"
               END IF
               IF (msh(iM)%eType .EQ. eType_NRB) THEN
                  DO i=1, nsd-1
                     IF (msh(iM)%bs(i)%p .LE. 1) err =
     2                  "NURBS for shell elements should be p > 1"
                  END DO
               END IF
               IF (msh(iM)%eType .EQ. eType_TRI) THEN
                  IF (.NOT.cm%seq()) err = "Triangular shell elements"//
     2               " should be run sequentially"
               END IF
            END IF
         END DO
      ELSE
         ALLOCATE(gX(nsd,gtnNo))
         gX = x
         DO iM=1, nMsh
            CALL SELECTELE(msh(iM))
            ALLOCATE(msh(iM)%gN(msh(iM)%gnNo))
            msh(iM)%gN = 0
         END DO
      END IF ! resetSim

!     Examining the existance of projection faces and setting %gN.
!     Reseting gtnNo and recounting nodes that are not duplicated
      gtnNo = 0
      CALL SETPROJECTOR(list, avNds)
      DO iM=1, nMsh
         DO a=1, msh(iM)%gnNo
            IF (msh(iM)%gN(a) .EQ. 0) THEN
               IF (PULLSTACK(avNds,i)) THEN
                  msh(iM)%gN(a) = i
               ELSE
                  gtnNo         = gtnNo + 1
                  msh(iM)%gN(a) = gtnNo
               END IF
            END IF
            IF (ALLOCATED(msh(iM)%gpN))
     2         msh(iM)%gpN(a) = msh(iM)%gN(a)
         END DO
      END DO
      DEALLOCATE(x)
      ALLOCATE(x(nsd,gtnNo))
      IF (avNds%n .NE. 0) err = "There are still nodes in the stack"
      CALL DESTROYSTACK(avNds)

!     Temporarily allocate msh%lN array. This is necessary for BCs and
!     will later be deallocated in DISTRIBUTE
      DO iM=1, nMsh
         ALLOCATE(msh(iM)%lN(gtnNo))
         msh(iM)%lN = 0
         DO a=1, msh(iM)%gnNo
            Ac = msh(iM)%gN(a)
            msh(iM)%lN(Ac) = a
         END DO
      END DO

!     Re-arranging fa and x and finding the size of the entire domain
      b = 0
      DO iM=1, nMsh
         DO a=1, msh(iM)%gnNo
            b  = b + 1
            Ac = msh(iM)%gN(a)
            DO i=1, nsd
               x(i,Ac) = gX(i,b)
               IF (maxX(i) .LT. x(i,Ac)) maxX(i) = x(i,Ac)
               IF (minX(i) .GT. x(i,Ac)) minX(i) = x(i,Ac)
            END DO
         END DO
         DO iFa=1, msh(iM)%nFa
            DO a=1, msh(iM)%fa(iFa)%nNo
               Ac = msh(iM)%fa(iFa)%gN(a)
               Ac = msh(iM)%gN(Ac)
               msh(iM)%fa(iFa)%gN(a) = Ac
            END DO
            DO e=1, msh(iM)%fa(iFa)%nEl
               DO a=1, msh(iM)%fa(iFa)%eNoN
                  Ac = msh(iM)%fa(iFa)%IEN(a,e)
                  Ac = msh(iM)%gN(Ac)
                  msh(iM)%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         END DO
      END DO

      IF (resetSim) THEN
         lDof = size(rmsh%Y0,1)
         lnNo = size(rmsh%Y0,2)
         IF (lnNo .NE. gtnNo) THEN
            ALLOCATE(tmpA(lDof,lnNo))
            ALLOCATE(tmpY(lDof,lnNo))
            ALLOCATE(tmpD(lDof,lnNo))
            tmpA = rmsh%A0
            tmpY = rmsh%Y0
            tmpD = rmsh%D0
            DEALLOCATE(rmsh%A0,rmsh%Y0,rmsh%D0)
            ALLOCATE(rmsh%A0(lDof,gtnNo))
            ALLOCATE(rmsh%Y0(lDof,gtnNo))
            ALLOCATE(rmsh%D0(lDof,gtnNo))
            b = 0
            DO iM=1, nMsh
               DO a=1, msh(iM)%gnNo
                  b = b + 1
                  Ac = msh(iM)%gpN(a)
                  DO i=1, lDof
                     rmsh%A0(i,Ac) = tmpA(i,b)
                     rmsh%Y0(i,Ac) = tmpY(i,b)
                     rmsh%D0(i,Ac) = tmpD(i,b)
                  END DO
               END DO
            END DO
         END IF ! lnNo < gtnNo
      END IF ! resetSim

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId.
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)

         lPtr => lPM%get(fTmp,"Domain file path")
         IF (ASSOCIATED(lPtr)) CALL SETDMNIDFF(msh(iM), fTmp%open())

         lPtr => lPM%get(i,"Domain",ll=0,ul=BIT_SIZE(dmnId)-1)
         IF (ASSOCIATED(lPtr)) CALL SETDMNID(msh(iM),i)
         IF (ALLOCATED(msh(iM)%eId)) flag = .TRUE.
      END DO
      IF (flag) THEN
         ALLOCATE(dmnId(gtnNo))
         dmnId = 0
         DO iM=1, nMsh
            IF (.NOT.ALLOCATED(msh(iM)%eId)) CYCLE
            DO e=1, msh(iM)%gnEl
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%gIEN(a,e)
                  Ac = msh(iM)%gN(Ac)
                  dmnId(Ac) = IOR(dmnId(Ac),msh(iM)%eId(e))
               END DO
            END DO
         END DO
      END IF

!     Read fiber orientation
!     Here I read fiber orientation temporarily into msh(iM)%x
!     Later I integrate them into global variable, fN
      nFn  = 0
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)
         j = lPM%srch("Fiber direction file path")
         IF (j .EQ. 0) j = lPM%srch("Fiber direction")
         IF (nFn .LT. j) nFn = j
      END DO
      IF (nFn .GT. 0) flag = .TRUE.

      IF (flag) THEN
         DO iM=1, nMsh
            lPM => list%get(msh(iM)%name,"Add mesh",iM)

            j = lPM%srch("Fiber direction file path")
            IF (j .NE. 0) THEN
               ALLOCATE(msh(iM)%x(nFn*nsd,msh(iM)%gnNo))
               msh(iM)%x = 0D0
               DO i=1, j
                  lPtr => lPM%get(cTmp,
     2               "Fiber direction file path (vtu)", i)
                  IF (ASSOCIATED(lPtr))
     2               CALL READVTUPDATA(msh(iM), cTmp, "FIB_DIR", nsd, i)
               END DO
            ELSE
               j = lPM%srch("Fiber direction")
               IF (j .NE. 0) THEN
                  ALLOCATE(msh(iM)%x(nFn*nsd,msh(iM)%gnNo))
                  msh(iM)%x = 0D0
                  DO i=1, j
                     lPtr => lPM%get(fibN, "Fiber direction", i)
                     IF (ASSOCIATED(lPtr)) THEN
                        DO a=1, msh(iM)%gnNo
                           msh(iM)%x((i-1)*nsd+1:i*nsd,a) = fibN(1:nsd)
                        END DO
                     END IF
                  END DO
               END IF
            END IF
         END DO

!     Copy fiber directions and normalize them to unit vectors
         std = " Loading fiber directions"
         ALLOCATE(fN(nFn*nsd,gtnNo))
         fN = 0D0
         DO iM=1, nMsh
            IF (.NOT.ALLOCATED(msh(iM)%x)) CYCLE
            DO a=1, msh(iM)%gnNo
               Ac = msh(iM)%gN(a)
               DO i=1, nFn
                  fibN(:) = msh(iM)%x((i-1)*nsd+1:i*nsd,a)
                  rtmp = SQRT(NORM(fibN))
                  IF (.NOT.ISZERO(rtmp)) fibN(:) = fibN(:)/rtmp
                  fN((i-1)*nsd+1:i*nsd,Ac) = fibN
               END DO
            END DO
            DEALLOCATE(msh(iM)%x)
         END DO
      ELSE
         nFn = 1
      END IF

!     Read prestress data
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)

         lPtr => lPM%get(cTmp, "Prestress file path (vtu)")
         IF (ASSOCIATED(lPtr)) THEN
            flag = .TRUE.
            ALLOCATE(msh(iM)%x(nstd,msh(iM)%gnNo))
            msh(iM)%x = 0D0
            CALL READVTUPDATA(msh(iM), cTmp, "PS_Prestress", nstd, 1)
         END IF
      END DO
      IF (flag) THEN
         ALLOCATE(pS0(nstd,gtnNo))
         pS0 = 0D0
         DO iM=1, nMsh
            IF (.NOT.ALLOCATED(msh(iM)%x)) CYCLE
            DO a=1, msh(iM)%gnNo
               Ac = msh(iM)%gN(a)
               pS0(:,Ac) = msh(iM)%x(:,a)
            END DO
            DEALLOCATE(msh(iM)%x)
         END DO
      END IF

!     Read general body force data
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)
         lPtr => lPM%get(fTmp, "Body force file path")
         IF (ASSOCIATED(lPtr)) CALL READBF(msh(iM), fTmp%open())
      END DO

!     Initialize pressure field from file
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)

         lPtr => lPM%get(cTmp,"Initial pressures file path (vtu)")
         IF(ASSOCIATED(lPtr)) THEN
            flag = .TRUE.
            ALLOCATE(msh(iM)%x(1,msh(iM)%gnNo))
            msh(iM)%x = 0D0
            CALL READVTUPDATA(msh(iM), cTmp, "NS_Pressure", 1, 1)
         END IF
      END DO
      IF (flag) THEN
         ALLOCATE(Pinit(gtnNo))
         Pinit = 0D0
         DO iM=1, nMsh
            IF(.NOT. ALLOCATED(msh(iM)%x)) CYCLE
            DO a=1, msh(iM)%gnNo
               Ac = msh(iM)%gN(a)
               Pinit(Ac) = msh(iM)%x(1,a)
            END DO
            DEALLOCATE(msh(iM)%x)
         END DO
      END IF

!     Initialize velocity field from file
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name,"Add mesh",iM)

         lPtr => lPM%get(cTmp,"Initial velocities file path (vtu)")
         IF(ASSOCIATED(lPtr)) THEN
            flag = .TRUE.
            ALLOCATE(msh(iM)%x(nsd,msh(iM)%gnNo))
            msh(iM)%x = 0D0
            CALL READVTUPDATA(msh(iM), cTmp, "NS_Velocity", nsd, 1)
         END IF
      END DO
      IF (flag) THEN
         ALLOCATE(Vinit(nsd,gtnNo))
         Vinit = 0D0
         DO iM=1, nMsh
            IF(.NOT. ALLOCATED(msh(iM)%x)) CYCLE
            DO a=1, msh(iM)%gnNo
               Ac = msh(iM)%gN(a)
               Vinit(:,Ac) = msh(iM)%x(:,a)
            END DO
            DEALLOCATE(msh(iM)%x)
         END DO
      END IF

!     Initializing displacement field from file. Used for CMM inflation
      flag = .FALSE.
      DO iM=1, nMsh
         lPM => list%get(msh(iM)%name, "Add mesh",iM)

         lPtr => lPM%get(ctmp,"Initial displacements file path (vtu)")
         IF(ASSOCIATED(lPtr)) THEN
            flag = .TRUE.
            ALLOCATE(msh(iM)%x(nsd,msh(iM)%gnNo))
            msh(iM)%x = 0D0
            CALL READVTUPDATA(msh(iM), cTmp, "Displacement", nsd, 1)
         END IF
      END DO
      IF (flag) THEN
         ALLOCATE(Dinit(nsd,gtnNo))
         Dinit = 0D0
         DO iM=1, nMsh
            IF(.NOT. ALLOCATED(msh(iM)%x)) CYCLE
            DO a=1, msh(iM)%gnNo
               Ac = msh(iM)%gN(a)
               Dinit(:,Ac) = msh(iM)%x(:,a)
            END DO
            DEALLOCATE(msh(iM)%x)
         END DO
      END IF

      IF (nMsh .GT. 1) THEN
         IF (.NOT.resetSim) THEN
            std = " Total number of nodes: "//gtnNo
            std = " Total number of elements: "//SUM(msh%gnEl)
         ELSE
            std = " Total number of nodes after remesh: "//gtnNo
            std = " Total number of elements after remesh: "//
     2         SUM(msh%gnEl)
         END IF
      END IF

      IF (.NOT.resetSim) THEN
         std = CLR(" Mesh data imported successfully",3)
      ELSE
         std = CLR(" Mesh data configured successfully after remesh",3)
      END IF

      DO i=1, nsd
         std = " Domain size "//dSym(i)//": "//CLR(STR(minX(i))
     2      //" "//STR(maxX(i)),4)
      END DO

!     Read contact model parameters
      IF (.NOT.resetSim) THEN
         iCntct = .FALSE.
         lPM => list%get(ctmp, "Contact model")
         IF (ASSOCIATED(lPM)) THEN
            iCntct = .TRUE.
            SELECT CASE (TRIM(ctmp))
            CASE ("penalty")
               cntctM%cType = cntctM_penalty
               lPtr => lPM%get(cntctM%k,
     2            "Penalty constant (k)", 1, ll=0D0)
               lPtr => lPM%get(cntctM%h,
     2            "Desired separation (h)", 1, lb=0D0)
               lPtr => lPM%get(cntctM%c,
     2            "Closest gap to activate penalty (c)", 1, lb=0D0)
               IF (cntctM%c .LT. cntctM%h) err =
     2            "Choose c > h for proper contact penalization"
               lPtr => lPM%get(cntctM%al,
     2            "Min norm of face normals (alpha)",1,lb=0D0,ub=1D0)
            CASE DEFAULT
               err = "Undefined contact model"
            END SELECT
         END IF
      END IF

      RETURN
      END SUBROUTINE READMSH
!--------------------------------------------------------------------
      SUBROUTINE READBF(lM, fid)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: fid

      INTEGER a, i, j, nTP, nNo, Ac
      REAL(KIND=8) rtmp

      READ (fid,*) i, nTP, nNo
      IF (nNo .GT. lM%gnNo) err = "No. of nodes out of bounds"//
     2   " (body force for mesh <"//TRIM(lM%name)//">)"
      IF (i .NE. nsd) err = "DOF should equal spatial dimension "//nsd//
     2   " (body force for mesh <"//TRIM(lM%name)//">)"

      ALLOCATE(lM%bf)
      lM%bf%dof = i
      lM%bf%nTP = nTP

      ALLOCATE(lM%bf%t(nTP), lM%bf%d(i,gtnNo,nTP))
      lM%bf%d = 0D0

      DO j=1, nTP
         READ(fid,*) rtmp
         lM%bf%t(j) = rtmp
         IF (j .EQ. 1) THEN
            IF (.NOT.ISZERO(rtmp)) err = "First time step should be 0"//
     2         " (body force for mesh <"//TRIM(lM%name)//">)"
         ELSE
            rtmp = rtmp - lM%bf%t(j-1)
            IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err = "Non-increasing "//
     2         "time trend found (body force for mesh <"//
     3         TRIM(lM%name)//">)"
         END IF
      END DO

      lM%bf%period = lM%bf%t(nTP)
      DO a=1, nNo
         READ(fid,*) Ac
         IF (Ac.GT.lM%gnNo .OR. Ac.LE.0) THEN
            err = "Entry "//Ac//" is out of bounds"//
     2         " (body force for mesh <"//TRIM(lM%name)//">)"
         END IF
         Ac = lM%gN(Ac)
         DO j=1, nTP
            READ(fid,*) (lM%bf%d(i,Ac,j), i=1, lM%bf%dof)
         END DO
      END DO
      CLOSE(fid)

      RETURN
      END SUBROUTINE READBF
!####################################################################
!     This routines associates two faces with each other and sets gN
      SUBROUTINE SETPROJECTOR(list, avNds)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list
      TYPE(stackType), INTENT(OUT) :: avNds

      INTEGER i, j, k, iM, jM, kM, iFa, jFa, ia, ja, nPrj, iPrj, nStk
      CHARACTER(LEN=stdL) ctmpi, ctmpj
      TYPE(stackType) lPrj
      TYPE(listType), POINTER :: lPtr, lPP

      TYPE(stackType), ALLOCATABLE :: stk(:)

!     Here calculating an upper limit for the number required stacks
      nStk = 0
      nPrj = list%srch("Add projection")
      IF (nPrj .GT. 0) THEN
         DO iM=1, nMsh
            ALLOCATE(msh(iM)%gpN(msh(iM)%gnNo))
         END DO
      END IF

      DO iPrj=1, nPrj
         lPP => list%get(ctmpi,"Add projection",iPrj)
         CALL FINDFACE(ctmpi, iM, iFa)
         nStk = nStk + msh(iM)%fa(iFa)%nNo
      END DO
      ALLOCATE(stk(nStk))

      DO iPrj=1, nPrj
         lPP => list%get(ctmpi,"Add projection",iPrj)
         CALL FINDFACE(ctmpi, iM, iFa)
         lPtr => lPP%get(ctmpj,"Project from face",1)
         CALL FINDFACE(ctmpj, jM, jFa)
         CALL MATCHFACES(msh(iM)%fa(iFa), msh(jM)%fa(jFa), lPrj)
         DO
!     First in, last out
            IF (.NOT.PULLSTACK(lPrj,ja)) EXIT
            IF (.NOT.PULLSTACK(lPrj,ia)) EXIT
            i = msh(iM)%gN(ia)
            j = msh(jM)%gN(ja)
            IF (i .EQ. 0) THEN
               IF (j .EQ. 0) THEN
!     Since neither of them have value, I will add a new node and both
!     of them to the stack
                  IF (.NOT.PULLSTACK(avNds, k)) THEN
                     gtnNo = gtnNo + 1
                     k     = gtnNo
                  END IF
                  msh(iM)%gN(ia) = k
                  msh(jM)%gN(ja) = k
                  CALL PUSHSTACK(stk(k), (/iM,ia,jM,ja/))
               ELSE
!     This is the case one of them has already been assigned. So just
!     using that value for the other one
                  msh(iM)%gN(ia) = j
                  CALL PUSHSTACK(stk(j), (/iM,ia/))
               END IF
            ELSE
               IF (j .EQ. 0) THEN
                  msh(jM)%gN(ja) = i
                  CALL PUSHSTACK(stk(i), (/jM,ja/))
               ELSE
!     Since they are both already have assigned values, I will move the
!     nodes from stack with bigger ID, j, to the other stack, i.
                  IF (i .EQ. j) CYCLE
                  IF (i .GT. j) THEN
                     k = i
                     i = j
                     j = k
                  END IF
                  DO
                     IF (.NOT.PULLSTACK(stk(j),ja)) EXIT
                     IF (.NOT.PULLSTACK(stk(j),kM)) EXIT
                     msh(kM)%gN(ja) = i
                     CALL PUSHSTACK(stk(i), (/kM,ja/))
                  END DO
                  CALL PUSHSTACK(avNds, j)
               END IF
            END IF
         END DO
!     Since all nodes are added, I remove all the members
         CALL DESTROY(lPrj)
      END DO

      RETURN
      END SUBROUTINE SETPROJECTOR
!--------------------------------------------------------------------
!     This is match isoparameteric faces to each other. Project nodes
!     from two adjacent meshes to each other based on a L2 norm.
      SUBROUTINE MATCHFACES(lFa, pFa, lPrj)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(INOUT) :: lFa, pFa
      TYPE(stackType), INTENT(OUT) :: lPrj

      TYPE blkType
         INTEGER :: n = 0
         INTEGER, ALLOCATABLE :: gN(:)
      END TYPE

      LOGICAL nFlt(nsd)
      INTEGER nBkd, i, a, b, Ac, Bc, iBk, nBk, iM, jM, iSh, jSh
      REAL(KIND=8) ds, xMin(nsd), xMax(nsd), dx(nsd)

      INTEGER, ALLOCATABLE :: nodeBlk(:)
      TYPE(blkType), ALLOCATABLE :: blk(:)
      INTEGER :: cnt

      IF ((lFa%eType .NE. pFa%eType) .OR. (lFa%nNo .NE. pFa%nNo))
     2   err = "Incompatible faces in MATCHFACES"

      iM  = lFa%iM
      jM  = pFa%iM
      iSh = 0
      jSh = 0
      DO i=1, iM-1
         iSh = iSh + msh(i)%gnNo
      END DO
      DO i=1, jM-1
         jSh = jSh + msh(i)%gnNo
      END DO

!     We want to have approximately 1000 nodes in each block. So we
!     calculate nBkd, which is the number of separate blockes in each
!     direction, based on that.
      a    = pFa%nNo
      nBkd = NINT((REAL(a,8)/1D3)**(3.33D-1))
      IF (nBkd .EQ. 0) nBkd = 1
      nBk  = nBkd**nsd
      ALLOCATE(nodeBlk(a), blk(nBk))

!     Finding the extends of the domain and size of each block
      DO i=1, nsd
         xMin(i) = MIN(MINVAL(x(i,iSh+lFa%gN)), MINVAL(x(i,jSh+pFa%gN)))
         xMax(i) = MAX(MAXVAL(x(i,iSh+lFa%gN)), MAXVAL(x(i,jSh+pFa%gN)))
         IF (xMin(i) .LT. 0D0) THEN
            xMin(i) = xMin(i)*(1D0+eps)
         ELSE
            xMin(i) = xMin(i)*(1D0-eps)
         END IF
         IF (xMax(i) .LT. 0D0) THEN
            xMax(i) = xMax(i)*(1D0-eps)
         ELSE
            xMax(i) = xMax(i)*(1D0+eps)
         END IF
      END DO
      dx = (xMax - xMin)/REAL(nBkd,8)
      nFlt = .TRUE.
      DO i=1, nsd
         IF (ISZERO(dx(i))) nFlt(i) = .FALSE.
      END DO

!     First finding an estimation for size of each block
      blk%n = 0
      DO a=1, pFa%nNo
         Ac  = pFa%gN(a) + jSh
         iBk = FINDBLK(x(:,Ac))
         nodeBlk(a) = iBk
         blk(iBk)%n = blk(iBk)%n + 1
      END DO
      DO iBk=1, nBk
         ALLOCATE(blk(iBk)%gN(blk(iBk)%n))
      END DO
      blk%n = 0
      DO a=1, pFa%nNo
         Ac  = pFa%gN(a)
         iBk = nodeBlk(a)
         blk(iBk)%n = blk(iBk)%n + 1
         blk(iBk)%gN(blk(iBk)%n) = Ac
      END DO

!     Doing the calculation for every single node on this face
      cnt = 0
      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         iBk = FINDBLK(x(:,Ac+iSh))
!     Checking every single node on the other face
         DO b=1, blk(iBk)%n
            Bc = blk(iBk)%gN(b)
            IF (iM.EQ.jM .AND. Ac.EQ.Bc) CYCLE
            ds  = SQRT(SUM( (x(:,BC+jSh) - x(:,Ac+iSh))**2 ))
            IF (ds .LT. 1D3*EPSILON(ds)) THEN
               CALL PUSHSTACK(lPrj, (/Ac,Bc/))
               cnt = cnt + 1
            END If
         END DO
      END DO

      IF ( (lFa%nNo.EQ.pFa%nNo) .AND. (cnt.NE.lFa%nNo) ) THEN
         err = " Failed to project faces between <"//TRIM(lFa%name)//
     2      "> and <"//TRIM(pFa%name)//">. Try changing the tolerance"//
     3      " of the closest-node-finding algorithm"
      END IF

      IF (lPrj%n .EQ. 0) err = "Unable to find any matching node"//
     2   " between <"//TRIM(lFa%name)//"> and <"//TRIM(pFa%name)//">"

      IF (lPrj%n/2 .NE. lFa%nNo) err = "Mismatch in no. of "//
     2   "projected nodes between <"//TRIM(lFa%name)//"> and <"//
     3   TRIM(pFa%name)//">"

      RETURN
      CONTAINS

!--------------------------------------------------------------------
      INTEGER FUNCTION FINDBLK(x)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: x(nsd)

      INTEGER i, j, k

      i = 1
      j = 1
      k = 1
      IF (nFlt(1)) i = INT((x(1) - xMin(1))/dx(1))
      IF (nFlt(2)) j = INT((x(2) - xMin(2))/dx(2))
      IF (i .EQ. nBkd) i = nBkd - 1
      IF (j .EQ. nBkd) j = nBkd - 1
      IF (nsd .EQ. 3) THEN
         IF (nFlt(3)) k = INT((x(3) - xMin(3))/dx(3))
         IF (k .EQ. nBkd) k = nBkd - 1
         FINDBLK = k + (j + i*nBkd)*nBkd + 1
      ELSE ! nsd .EQ. 2
         FINDBLK = j + i*nBkd + 1
      END IF

      RETURN
      END FUNCTION FINDBLK

      END SUBROUTINE MATCHFACES
!####################################################################
!     Read domain from a file
      SUBROUTINE SETDMNIDFF(lM, fid)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: fid

      INTEGER i, e, a, btSiz
      CHARACTER(LEN=stdL) ctmp, stmp

      INTEGER, ALLOCATABLE :: iDmn(:)

      btSiz = BIT_SIZE(dmnId)
      IF (btSiz .NE. BIT_SIZE(lM%eId))
     2   err = "Use same type for dmnId and lM%eId"

!     Check to see if I need to increase the size of dmnId
      IF (.NOT.ALLOCATED(lM%eId)) THEN
         ALLOCATE(lM%eId(lM%gnEl))
         lM%eId = 0
      END IF

      ALLOCATE (iDmn(btSiz))
      ctmp = "Domain ID exceeds BIT_SIZE(dmnId). Reduce the number"//
     2   " of domains or increase BIT_SIZE(dmnId)."
      DO e=1, lM%gnEl
         READ (fid,"(A)") stmp
         i = CheckNoNumbers(stmp)
         IF (i .GT. btSiz) err = ctmp
         READ (stmp,*) iDmn(1:i)
         DO a=1, i
            IF (iDmn(a) .GE. btSiz) THEN
               err = ctmp
            ELSE IF (iDmn(a) .LT. 0) THEN
               err = "Domain ID must greater or equal to 0"
            END IF
            lM%eId(e) = IBSET(lM%eId(e),iDmn(a))
         END DO
      END DO
      CLOSE (fid)

      RETURN
      END SUBROUTINE SETDMNIDFF
!####################################################################
!     Check the mesh. If this is flag=nonL=.true., then we only check first 4 nodes
!     of IEN array
      SUBROUTINE CHECKIEN(lM, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(IN) :: flag

      INTEGER Ac, b, i, sn(4), e, a, teNoN

      INTEGER, ALLOCATABLE :: incNodes(:)
      REAL(KIND=8), ALLOCATABLE :: v(:,:), xl(:,:)

      ALLOCATE(v(nsd,lM%eNoN), xl(nsd,lM%eNoN))

      teNoN = lM%eNoN
      IF (flag .AND. lM%eType.EQ.eType_BIQ) teNoN = 4
      ALLOCATE (incNodes(lM%gnNo))
      std = " Checking IEN array structure for mesh <"//
     2 TRIM(lM%name)//">"
      incNodes = 0
      DO e=1, lM%gnEl
         DO a=1, teNoN
            Ac = lM%gIEN(a,e)
            IF (Ac.GT.lM%gnNo .OR. Ac.LE.0) THEN
               err = "Element "//e//
     2            " contains an out of bound node"
            END IF
            incNodes(Ac) = 1
         END DO
      END DO
      DO Ac=1, lM%gnNo
         IF (incNodes(Ac) .EQ. 0) THEN
            err = "Node "//Ac//
     2         " is isolated from the mesh"
         END IF
      END DO
      IF (lM%eType .EQ. eType_BRK) THEN
         std = " Make sure nodes in elements are arranged 1-4 on"
         std = "      one face and 5-8 on the opposite face"
         RETURN
      END IF

      DO e=1, lM%gnEl
!     By default no change
         a = 1; b = 1
         IF (lM%eType.EQ.eType_BIL .OR. lM%eType.EQ.eType_BIQ) THEN
            xl     = lM%x(:,lM%gIEN(1:4,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,3)
            v(:,4) = xl(:,1) - xl(:,4)
            sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))
            sn(2)  = SGN(v(1,2)*v(2,3) - v(2,2)*v(1,3))
            sn(3)  = SGN(v(1,3)*v(2,4) - v(2,3)*v(1,4))
            sn(4)  = SGN(v(1,4)*v(2,1) - v(2,4)*v(1,1))
            i = sn(1) + sn(2) + sn(3) + sn(4)
            IF (i .EQ. 0) THEN
               IF (sn(1) .EQ. sn(2)) THEN
!     Sign is changing between 2-3 and 1-4
                  IF (sn(1) .EQ. 1) THEN
                     a = 1; b = 4
                  ELSE
                     a = 2; b = 3
                  END IF
               ELSE
!     Sign is changing between 1-2 and 3-4
                  IF (sn(1) .EQ. 1) THEN
                     a = 3; b = 4
                  ELSE
                     a = 1; b = 2
                  END IF
               END IF
            ELSE IF (i .EQ. -4) THEN
!     Sign is not changing, but this is going C.W.
               a = 1; b = 3
            ELSE IF (i.EQ.2 .OR. ANY(sn.EQ.0)) THEN
!     Two or more edges are on a straight line
               err = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_WDG) THEN
            xl     = lM%x(:,lM%gIEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,1)
            v(:,4) = CROSS(v(:,1:2))
            sn(1)  = SGN(SUM(v(:,3)*v(:,4)))
            i = sn(1)
            IF (i .EQ. -1) THEN
!     Two nodes must be switched
               a = 1; b = 2
               Ac = lM%gIEN(a,e)
               lM%gIEN(a,e) = lM%gIEN(b,e)
               lM%gIEN(b,e) = Ac
               a = 4; b = 5
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               err = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_TET) THEN
            xl     = lM%x(:,lM%gIEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            v(:,3) = xl(:,4) - xl(:,3)
            v(:,4) = CROSS(v(:,1:2))
            sn(1)  = SGN(SUM(v(:,3)*v(:,4)))
            i = sn(1)
            IF (i .EQ. 1) THEN
!     Two nodes must be switched
               a = 1; b = 2
               IF (e .EQ. 1) std = " Reordering element connectivity"
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               err = "Element "//e//" is distorted"
            END IF
         ELSE IF (lM%eType .EQ. eType_TRI) THEN
            IF (nsd .NE. 2) CYCLE
            xl(:,1:3) = lM%x(:,lM%gIEN(:,e))
            v(:,1) = xl(:,2) - xl(:,1)
            v(:,2) = xl(:,3) - xl(:,2)
            sn(1)  = SGN(v(1,1)*v(2,2) - v(2,1)*v(1,2))
            i = sn(1)
            IF (i .EQ. -1) THEN
!     Two nodes must be switched
               a = 1; b = 2
            ELSE IF (i .EQ. 0) THEN
!     Two or more edges are on a straight line
               err = "Element "//e//" is distorted"
            END IF
         END IF

         Ac = lM%gIEN(a,e)
         lM%gIEN(a,e) = lM%gIEN(b,e)
         lM%gIEN(b,e) = Ac
      END DO

      RETURN
      END SUBROUTINE CHECKIEN
!####################################################################
!     Making sure all the faces contain in range values
      SUBROUTINE CALCNBC(lM, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER Ac, e, Ec, a

      INTEGER, ALLOCATABLE :: incNd(:)

!     Finding the nodes that belongs to this face and counting the
!     number of the nodes
      ALLOCATE (incNd(lM%gnNo))
      incNd   = 0
      lFa%nNo = 0
      DO e=1, lFa%nEl
         Ec = lFa%gE(e)
         IF (Ec.GT.lM%gnEl .OR. Ec.LE.0) THEN
            err = "Global ID of element "//e//" of face <"
     2         //TRIM(lFa%name)//"> is out of range "//Ec
         END IF
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            IF (Ac.GT.lM%gnNo .OR. Ac.LE.0) THEN
               err = "IEN of element "//e//" of face <"
     2            //TRIM(lFa%name)//"> is out of range"
            END IF
            IF (ALL(lM%gIEN(:,Ec) .NE. Ac)) THEN
               err = "Parent element "//e//" of face <"
     2            //TRIM(lFa%name)//"> does not contain same"
     3            //" nodes as the boundary"
            END IF
            IF (incNd(Ac) .EQ. 0) THEN
               lFa%nNo = lFa%nNo + 1
               incNd(Ac) = 1
            END IF
         END DO
      END DO

      ALLOCATE (lFa%gN(lFa%nNo))
      a = 0
      DO Ac=1, lM%gnNo
         IF (incNd(Ac) .NE. 0) THEN
            a = a + 1
            lFa%gN(a) = Ac
         END IF
      END DO

      RETURN
      END SUBROUTINE CALCNBC
!####################################################################
!     Calculate mesh properties to check its quality
      SUBROUTINE CALCMESHPROPS(nMesh, mesh)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nMesh
      TYPE(mshType), INTENT(INOUT) :: mesh(nMesh)

      INTEGER :: iM, cPhys, e, ierr

      LOGICAL, ALLOCATABLE :: gFlag(:,:)

      std = "----------------------------------------------------------"
      DO iM=1, nMesh
         std = " Mesh properties: <"//CLR(TRIM(mesh(iM)%name))//">"
         CALL CALCELEMJAC(mesh(iM), rmsh%flag(iM))

         CALL CALCELEMSKEW(mesh(iM), rmsh%flag(iM))

         CALL CALCELEMAR(mesh(iM), rmsh%flag(iM))
      END DO

      IF (rmsh%isReqd .AND. .NOT.ALL(rmsh%flag(:)) ) THEN
         IF (cTS .EQ. rmsh%fTS) THEN
            rmsh%fTS = rmsh%fTS + rmsh%freq
            DO iM=1, nMsh
               DO e=1, msh(iM)%nEl
                  cDmn = DOMAIN(msh(iM), cEq, e)
                  cPhys = eq(cEq)%dmn(cDmn)%phys
                  IF (cPhys .EQ. phys_fluid) THEN
                     rmsh%flag(iM) = .TRUE.
                     EXIT
                  END IF
               END DO
            END DO
            rmsh%cntr = rmsh%cntr + 1
         END IF
      END IF

      ALLOCATE(gFlag(nMesh,cm%np()))
      CALL MPI_ALLGATHER(rmsh%flag, nMesh, mplog, gFlag, nMesh, mplog,
     2   cm%com(), ierr)

      DO iM=1, nMesh
         IF (ANY(gFlag(iM,:))) THEN
            rmsh%flag(iM) = .TRUE.
         ELSE
            rmsh%flag(iM) = .FALSE.
         END IF
      END DO
      IF (ANY(rmsh%flag(:))) resetSim = .TRUE.

      DEALLOCATE(gFlag)

      RETURN
      END SUBROUTINE CALCMESHPROPS
!####################################################################
!     Calculate element Jacobian of a given mesh
      SUBROUTINE CALCELEMJAC(lM, rflag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(INOUT) :: rflag

      INTEGER :: e, a, Ac, cnt, iDmn, cPhys
      REAL(KIND=8) :: minJ, maxJ, tmp

      REAL(KIND=8), ALLOCATABLE :: Jac(:), xl(:,:), dol(:,:)

      rflag = .FALSE.
      ALLOCATE(Jac(lM%nEl))
      ALLOCATE(xl(nsd,lM%eNoN))
      ALLOCATE(dol(nsd,lM%eNoN))

      cnt = 0
      DO e=1, lM%nEl
         iDmn = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(iDmn)%phys
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            IF(mvMsh) THEN
               dol(:,a) = Do(nsd+2:2*nsd+1,Ac)
            END IF
         END DO

         IF (mvMsh) xl = xl + dol

         Jac(e) = JACOBIAN(nsd,lM%eNoN,xl,lM%Nx(:,:,1))
         IF (Jac(e) .LT. 0D0) THEN
            cnt = cnt + 1
            IF (cPhys .NE. phys_fluid) err = "Negative Jacobian in "//
     2         "non-fluid domain"
         END IF
      END DO ! e

      maxJ = MAXVAL(ABS(Jac(:)))
      maxJ = cm%reduce(maxJ, MPI_MAX)
      Jac(:) = Jac(:)/ABS(maxJ)

      minJ = MINVAL(Jac(:))
      minJ = cm%reduce(minJ, MPI_MIN)
      cnt  = cm%reduce(cnt)
      tmp  = 1D2 * REAL(cnt,KIND=8) / REAL(lM%gnEl,KIND=8)
      IF (minJ .LT. 0D0) rflag = .TRUE.

      IF (.NOT.rflag) THEN
         std = "    Min normalized Jacobian <"//minJ//">"
      ELSE
         rmsh%cntr = rmsh%cntr + 1
         std = " "
         std = "cccccccccccccccccccccccccccccccccccccccccccccccccc"//
     2      "cccccccc"
         std = ""
         std = " Mesh is DISTORTED (Min. Jacobian < 0) at time "//cTS
         std = "    Min normalized Jacobian <"//minJ//">"
         std = "    No. of Elements with Jac < 0: "//cnt//
     2      " ("//tmp//"%)"
         IF (.NOT.rmsh%isReqd) THEN
            err = "Unexpected behavior! Mesh is DISTORTED. "//
     2         " But Remesher is not initialized. Fatal"
         END IF
      END IF
      DEALLOCATE(xl, dol, Jac)

      RETURN
      END SUBROUTINE CALCELEMJAC
!####################################################################
!     Calculate element Skewness of a given mesh
      SUBROUTINE CALCELEMSKEW(lM, rflag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(INOUT) :: rflag

      INTEGER :: e, a, Ac, i, iDmn, cPhys, bins(5)
      REAL(KIND=8) :: maxSk, p1, p2, tmp(5)

      REAL(KIND=8), ALLOCATABLE :: Skw(:), xl(:,:), dol(:,:)

      IF (lM%eType .NE. eType_TET .AND.
     2    lM%eType .NE. eType_TRI) THEN
c         wrn = " Skewness is computed for TRI and TET elements only"
         RETURN
      END IF

      ALLOCATE(Skw(lM%nEl))
      ALLOCATE(xl(nsd,lM%eNoN))
      ALLOCATE(dol(nsd,lM%eNoN))

      bins(:) = 0
      DO e=1, lM%nEl
         iDmn = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(iDmn)%phys
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            IF(mvMsh) THEN
               dol(:,a) = Do(nsd+2:2*nsd+1,Ac)
            END IF
         END DO
         IF (mvMsh) xl = xl + dol
         Skw(e) = SKEWNESS(nsd, lM%eNoN, xl)
         p1 = 0D0
         p2 = 6D-1
         DO i=1, 5
            IF (Skw(e).GT.p1 .AND. Skw(e).LE.p2) THEN
               bins(i) = bins(i) + 1
               EXIT
            END IF
            p1 = p2
            p2 = p1 + 1D-1
         END DO ! i
      END DO ! e

      maxSk = cm%reduce(maxSk, MPI_MAX)
      bins  = cm%reduce(bins)
      tmp(:) = 1D2 * REAL(bins(:),KIND=8) / REAL(lM%gnEl,KIND=8)
      IF (rflag) THEN
         std = "    Max Skewness <"//maxSk//">"
         std = "    Skew [    < 0.6] <"//bins(1)//">  ("//tmp(1)//"%)"
         std = "         [0.6 - 0.7] <"//bins(2)//">  ("//tmp(2)//"%)"
         std = "         [0.7 - 0.8] <"//bins(3)//">  ("//tmp(3)//"%)"
         std = "         [0.8 - 0.9] <"//bins(4)//">  ("//tmp(4)//"%)"
         std = "         [0.9 - 1.0] <"//bins(5)//">  ("//tmp(5)//"%)"
      ELSE
         std = "    Max Skewness <"//maxSk//">"
      END IF
      DEALLOCATE(xl, dol, Skw)

      RETURN
      END SUBROUTINE CALCELEMSKEW
!####################################################################
!     Calculate element Aspect Ratio of a given mesh
      SUBROUTINE CALCELEMAR(lM, rflag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(INOUT) :: rflag

      INTEGER :: e, a, Ac, i, iDmn, cPhys, bins(5)
      REAL(KIND=8) :: maxAR, p1, p2, tmp(5)

      REAL(KIND=8), ALLOCATABLE :: AsR(:), xl(:,:), dol(:,:)

      IF (lM%eType .NE. eType_TET .AND.
     2    lM%eType .NE. eType_TRI) THEN
c         wrn = "AR is computed for TRI and TET elements only"
         RETURN
      END IF

      ALLOCATE(AsR(lM%nEl))
      ALLOCATE(xl(nsd,lM%eNoN))
      ALLOCATE(dol(nsd,lM%eNoN))

      bins(:) = 0
      DO e=1, lM%nEl
         iDmn = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(iDmn)%phys
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            IF(mvMsh) THEN
               dol(:,a) = Do(nsd+2:2*nsd+1,Ac)
            END IF
         END DO
         IF (mvMsh) xl = xl + dol
         AsR(e) = ASPECTRATIO(nsd, lM%eNoN, xl)
         p1 = 0D0
         p2 = 5D0
         DO i=1, 4
            IF (AsR(e).GT.p1 .AND. AsR(e).LE.p2) THEN
               bins(i) = bins(i) + 1
               EXIT
            END IF
            p1 = p2
            p2 = p1 + 5D0
         END DO ! i
         IF (i .GT. 4) bins(5) = bins(5) + 1
      END DO ! e

      maxAR = MAXVAL(AsR(:))
      maxAR = cm%reduce(maxAR, MPI_MAX)
      bins  = cm%reduce(bins)

      tmp(:) = 1D2 * REAL(bins(:),KIND=8) / REAL(lM%gnEl,KIND=8)
      IF (rflag) THEN
         std = "    Max Asp. Ratio <"//maxAR//">"
         std = "    AR [   <  5] <"//bins(1)//">  ("//tmp(1)//"%)"
         std = "       [ 5 - 10] <"//bins(2)//">  ("//tmp(2)//"%)"
         std = "       [10 - 15] <"//bins(3)//">  ("//tmp(3)//"%)"
         std = "       [15 - 20] <"//bins(4)//">  ("//tmp(4)//"%)"
         std = "       [   > 20] <"//bins(5)//">  ("//tmp(5)//"%)"
      ELSE
         std = "    Max Asp. Ratio <"//maxAR//">"
      END IF
      DEALLOCATE(xl, dol, AsR)

      RETURN
      END SUBROUTINE CALCELEMAR
!####################################################################

