! mainSharedNoRMAOuterLocks.f90
! Vereinfachte Darstellung des Kommunikationsmusters mit OpenMP in Phoenix
!
! Mittels Cmd-Args können einige Laufzeitgroessen angepasst werden, um eine
! Entwicklung verfolgen zu können 
!
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

      PROGRAM main         
         USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
         USE mpi
!         USE omp_lib

         USE printsMod
         USE dataTypesMod
         USE computationMod
         USE benchmarkUtilityMod
         IMPLICIT NONE
         !--------------------------------------------------------------
         !---------------Phoenix-ähnliche-Variablen---------------------
         !--------------------------------------------------------------
         TYPE(pve_data)                 :: pve_grid
         TYPE(varContainer)             :: container

         INTEGER,TARGET                 :: n_chars, reduction
         INTEGER,TARGET                 :: nx,ny,nz,ix,iy,iz, idx
         INTEGER                        :: idxTmp

         REAL(kind(prec)),TARGET,ALLOCATABLE,DIMENSION(:) &
                                        :: speicher1D
         REAL(kind(prec)),TARGET,ALLOCATABLE,DIMENSION(:,:,:) &
                                        :: speicher3D
         REAL(kind(prec)),TARGET        :: stelle_max   
         INTEGER,         TARGET        :: zaehler
         !--------------------------------------------------------------
         !---------------eigene Variablen-------------------------------
         !--------------------------------------------------------------
         INTEGER                        :: i, iter, stelle, fenceIndex
         CHARACTER(len=32)              :: arg,arg2
         REAL(kind(prec))               :: start
         REAL(kind(prec))               :: meantime, sd
         REAL(kind(prec))               :: j
         REAL(kind(prec))               :: tmp
         LOGICAL                        :: testBool
         INTEGER                        :: reductionTmp
         !--------------------------------------------------------------
         !---------------MPI-Variablen----------------------------------
         !--------------------------------------------------------------
         INTEGER                        :: ierror
         INTEGER                        :: myRank               = 0
         INTEGER                        :: numProcs             = 1
         INTEGER                        :: targetProc           = 0
         INTEGER, PARAMETER             :: masterProc           = 0
         INTEGER, PARAMETER             :: ONE                  = 1
         INTEGER, PARAMETER             :: schritt1ArrayLen     = 2

         INTEGER                        :: pveWin, containerWin 
         INTEGER                        :: s1Win, s2Win
         INTEGER                        :: s3Win1D, s3Win3D, s3WinMax
         INTEGER, DIMENSION(:), POINTER :: schritt1Array
         INTEGER                        :: schritt2
         INTEGER                        :: startIndex, endIndex
         INTEGER                        :: sizeOfInteger, sizeOfReal
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfDouble
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfArray
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfSchritt2
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfSchritt3D1
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfSchritt3D3
         INTEGER(KIND=MPI_ADDRESS_KIND) :: sizeOfSchritt3Max
         INTEGER(KIND=MPI_ADDRESS_KIND) :: s3TargetDisp
         INTEGER                        :: memModel
         LOGICAL                        :: flag
         INTEGER                        :: lb, extent
         INTEGER                        :: prec_mpi, sizeOfprec_mpi
         INTEGER                        :: shmcomm, sizeOfComm
         TYPE(C_PTR)                    :: s1Ptr, s2Ptr
         TYPE(C_PTR)                    :: s3D1Ptr, s3D3Ptr, s3MaxPtr
         INTEGER, POINTER               :: zaehlerPtr
         REAL(kind(prec)), POINTER      :: stelle_maxPtr
         REAL(kind(prec)), DIMENSION(:), POINTER &
                                        :: speicher1DPtr
         REAL(kind(prec)), DIMENSION(:,:,:), POINTER &
                                        :: speicher3DPtr
         !--------------------------------------------------------------
         !--------------Hauptprogramm-----------------------------------
         !--------------------------------------------------------------        

         !----MPI-Initialisierung----
         CALL MPI_INIT(ierror)
         CALL testForError(ierror,"MPI: INIT")
         CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank,ierror)
         CALL testForError(ierror,"MPI: Rang bestimmen")
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numProcs, ierror)
         CALL testForError(ierror, "MPI: Proc-Gesamtzahl bestimmen")
         !Eigenen MPI-Datentyp mit prec erstellen
         CALL MPI_TYPE_CREATE_F90_REAL(KIND(prec), MPI_UNDEFINED,&
               prec_mpi,ierror)
         CALL testForError(ierror, "MPI: Erzeugung von prec_mpi")

         !Wie gross ist der Bereich fuer shared Memory
         CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,&
               0, MPI_INFO_NULL, shmcomm, ierror)
         CALL testForError(ierror, "MPI: Bestimme max. SharedMem-Comm")
         !Programm ist ohne Anpassungen nicht bzgl. Node-Zahl skalierbar
         CALL MPI_COMM_SIZE(shmcomm, sizeOfComm, ierror)
         CALL testForError(ierror, "MPI: Bestimme sizeOfComm")
         IF(sizeOfComm .NE. numProcs) THEN
            CALL testForError(-1,"MPI: shared memory is distributed")
         END IF         

         SELECT CASE(myRank)
            CASE(0)
               master = 6
            CASE DEFAULT
               master = 7 
               OPEN(UNIT=master, FILE='/dev/null')
               infosDrucken             = .FALSE.
               zeitDrucken              = .FALSE.    
               verifikationDrucken      = .FALSE.
         END SELECT

         !----Cmd-Args einlesen----
         nz = laenge 
         DO i = 1, iargc()
            CALL getarg(i, arg)
            !WRITE (*,*) arg
            SELECT CASE (arg)
               CASE ('-v', '--version')
                  WRITE(master,*) "Version ", version
                  WRITE(master,*) "Aktueller Stand: "
                  WRITE(master,*) "Schritt 1: cricical fertig,"&
                        ," reduction fertig"
                  WRITE(master,*) "Schritt 2: atomic fertig"
                  WRITE(master,*) "Schritt 3: konfliktfrei Schreiben",&
                      "fertig, geteilte Max-/Min-Bestimmung fertig"
                  WRITE(master,*) "Verifikation: in Arbeit"
                  CALL EXIT(0)
               CASE ('--xyz')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) laenge
                  nx = laenge
                  ny = laenge
                  nz = laenge
               CASE ('--l1')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) lasteinstellung(1)
               CASE ('--l2')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) lasteinstellung(2)
               CASE ('--l3')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) lasteinstellung(3)
               CASE ('-t','--threads')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) nthreads
               CASE ('--nz')
                  CALL getarg(i+1,arg2)
                  READ(arg2,*) nz
               CASE ('-h', '--help', '--usage')
                  CALL helpPrint()
                  CALL EXIT(0)
               CASE ('--s1')
                  schritt(1) = .TRUE.
               CASE ('--s2')
                  schritt(2) = .TRUE.
               CASE ('--s3')
                  schritt(3) = .TRUE.
            END SELECT
         END DO
         WRITE(master,*) "Ausführung gestartet mit:"
         WRITE(master,*) " -- Kubusseitenlaenge: ", 2*laenge+1
         WRITE(master,*) " -- nz: ",nz
         WRITE(master,*) " -- Lasteinstellung: ", lasteinstellung
         WRITE(master,*) " -- Iterationszahl: ", iterationen
         WRITE(master,*) " -- Thread-Zahl: irrelevant, da unbeachtet"
         WRITE(master,*) " -- Prozess-Zahl: ", numProcs
         WRITE(master,*) " -- Entferne Extremwerte aus Zeitmessung: "&
                     ,zeitBereinigung
         
            
         !Initialisiere Variablen und weise sie Pointer im Container zu
         !Container erleichtert die Wertausgabe per Subroutine
         n_chars                = 0
         reduction              = 0
         idx                    = 0
         nx                     = laenge
         ny                     = laenge
         !nz                     = laenge
         verifikationsArray     = .FALSE.
         zaehler                = 0
         reductionTmp           = 0
         
         container%n_chars      => n_chars
         container%reduction    => reduction
         container%nxyz         => laenge
         container%nz           => nz
         !container%ix          => ix
         !container%iy          => iy
         !container%iz          => iz
         container%idx          => idx
         container%stelle_max   => stelle_max
         pve_grid%zaehler       => zaehler
         
      
         !schritt1Array(1) = n_chars
         !schritt1Array(2) = reduction
         !Allokiere Schritt3-Speicher fuer nachfolgende Win-Erzeugung
         !Diese beiden dynamischen Speicherbereiche soll nur MasterProc 
         !allokieren
         IF(myRank .EQ. masterProc) THEN
            ALLOCATE(speicher1D((container%nxyz*2+1)**2*(nz*2+1)))
            ALLOCATE(speicher3D(-nx:nx,-ny:ny,-nz:nz))
         END IF
         pve_grid%speicher1D => speicher1D
         pve_grid%speicher3D => speicher3D
      
         !------------MPI-RMA-Initialisierung--------------------------
         
         !Bestimme Groesse der Speicherbereiche
!         CALL MPI_Type_get_extent(MPI_INTEGER, lb, extent,ierror)
         CALL MPI_TYPE_SIZE(MPI_INTEGER, sizeOfInteger, ierror)
         CALL testForError(ierror, "MPI: get Int-Size")
         CALL MPI_TYPE_SIZE(MPI_REAL, sizeOfReal, ierror)
         CALL testForError(ierror, "MPI: get Real-Size")
!         CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, sizeOfDouble, ierror)
!         CALL testForError(ierror, "MPI: get Double-Size")
         CALL MPI_TYPE_SIZE(prec_mpi, sizeOfprec_mpi,ierror)
         CALL testForError(ierror, "MPI: get prec_mpi-Size")
         IF(myRank .EQ. masterProc) THEN
            sizeOfArray         = sizeOfInteger *schritt1ArrayLen
            sizeOfSchritt2      = sizeOfInteger *ONE
            sizeOfSchritt3D1    = sizeOfprec_mpi*&
               ((container%nxyz*2+1)**2*(container%nz*2+1))
            sizeOfSchritt3D3    = sizeOfprec_mpi*&
               ((container%nxyz*2+1)**2*(container%nz*2+1))
            sizeOfSchritt3Max   = sizeOfprec_mpi*ONE
         ELSE
            sizeOfArray         = 0_MPI_ADDRESS_KIND
            sizeOfSchritt2      = 0_MPI_ADDRESS_KIND
            sizeOfSchritt3D1    = 0_MPI_ADDRESS_KIND
            sizeOfSchritt3D3    = 0_MPI_ADDRESS_KIND
            sizeOfSchritt3Max   = 0_MPI_ADDRESS_KIND
         END IF
         !Einrichtung der MPI-Windows, auf die dann geschrieben/gelesen
         !Schritt1
!         CALL MPI_WIN_CREATE(schritt1Array, sizeOfArray,sizeOfInteger,&
!                   MPI_INFO_NULL,MPI_COMM_WORLD,s1Win,ierror)  
!         WRITE(master,*) "Punkt A"
         CALL MPI_WIN_ALLOCATE_SHARED(sizeOfArray, sizeOfInteger, &
            MPI_INFO_NULL, shmcomm, s1Ptr, s1Win, ierror)
         CALL testForError(ierror, "MPI: Create shared Win1")
!         WRITE(master,*) "Punkt B"
         !Get Pointer from masterProc-Mem
         WRITE(master,*) "Punkt C"
         !Schritt2
!         CALL MPI_WIN_CREATE(zaehler, sizeOfSchritt2,sizeOfInteger,&
!                   MPI_INFO_NULL,MPI_COMM_WORLD,s2Win,ierror)
!         CALL testForError(ierror, "MPI: Create s2Win")
         CALL MPI_WIN_ALLOCATE_SHARED(sizeOfSchritt2, sizeOfInteger, &
            MPI_INFO_NULL, shmcomm, s2Ptr, s2Win, ierror)
         CALL testForError(ierror, "MPI: Create shared Win1")
         !Get Pointer from masterProc-Mem
         !Schritt3
!         CALL MPI_WIN_CREATE(pve_grid%speicher1D, sizeOfSchritt3D1,&
!                  sizeOfReal,MPI_INFO_NULL,MPI_COMM_WORLD,s3Win1D, &
!                  ierror)
!         CALL testForError(ierror, "MPI: Create s3Win1D")
!         CALL MPI_WIN_CREATE(speicher1D, sizeOfSchritt3D1,&
!                  sizeOfprec_mpi,MPI_INFO_NULL,MPI_COMM_WORLD,s3Win1D, &
!                  ierror)
         CALL MPI_WIN_ALLOCATE_SHARED(sizeOfSchritt3D1,sizeOfprec_mpi,&
                  MPI_INFO_NULL, shmcomm, s3D1Ptr,s3Win1D, ierror)
         CALL testForError(ierror, "MPI: Create s3Win1D")
         WRITE(master,*) "Punkt C2"
!         CALL MPI_WIN_CREATE(speicher3D, sizeOfSchritt3D3,&
!                  sizeOfprec_mpi,MPI_INFO_NULL,MPI_COMM_WORLD,s3Win3D, &
!                  ierror)
         CALL MPI_WIN_ALLOCATE_SHARED(sizeOfSchritt3D3,sizeOfprec_mpi,&
                  MPI_INFO_NULL, shmcomm, s3D3Ptr,s3Win3D, ierror)
         CALL testForError(ierror, "MPI: Create s3Win3D")
         WRITE(master,*) "Punkt C3"
!         CALL MPI_WIN_CREATE(stelle_max,sizeOfprec_mpi ,&
!                  sizeOfSchritt3Max,MPI_INFO_NULL,MPI_COMM_WORLD,&
!                  s3WinMax, ierror)
         CALL MPI_WIN_ALLOCATE_SHARED(sizeOfSchritt3Max,sizeOfprec_mpi,&
                  MPI_INFO_NULL, shmcomm, s3MaxPtr,s3WinMax, ierror)
         CALL testForError(ierror, "MPI: Create s3WinMax")

         !Nicht-MasterProcs brauchen Ptr
         IF(myRank .NE. masterProc) THEN
            !Schritt1
            CALL MPI_WIN_SHARED_QUERY(s1Win, masterProc,&
                 sizeOfArray, sizeOfInteger, s1Ptr, ierror)
            !Schritt2
            CALL MPI_WIN_SHARED_QUERY(s2Win, masterProc,&
                 sizeOfSchritt2, sizeOfInteger, s2Ptr, ierror)
            !Schritt3
            CALL MPI_WIN_SHARED_QUERY(s3Win1D, masterProc,&
                 sizeOfSchritt3D1, sizeOfprec_mpi, s3D1Ptr, ierror)
            CALL MPI_WIN_SHARED_QUERY(s3Win3D, masterProc,&
                 sizeOfSchritt3D3, sizeOfprec_mpi, s3D3Ptr, ierror)
            CALL MPI_WIN_SHARED_QUERY(s3WinMax, masterProc,&
                 sizeOfSchritt3Max, sizeOfprec_mpi, s3MaxPtr, ierror)
         END IF
         !Associate c-pointer with fortran pointer
         WRITE(master,*) "Punkt C4"
         !Schritt1
         CALL C_F_POINTER(s1Ptr, schritt1Array, (/ 2 /))
         !Schritt2
         CALL C_F_POINTER(s2Ptr, zaehlerPtr)
         !Schritt3
         CALL C_F_POINTER(s3D1Ptr, speicher1DPtr,&
                (/ ((container%nxyz*2+1)**2*(nz*2+1)) /))
         CALL C_F_POINTER(s3D3Ptr, speicher3DPtr,&
!                (/ -nx:nx, -ny:ny, -nz:nz /))
                  (/ 2*nx+1, 2*ny+1, 2*nz+1 /))
         CALL C_F_POINTER(s3MaxPtr, stelle_maxPtr)
         WRITE(master,*) "Fenstererstellung erfolgreich abgeschlossen"
         !-------------------------------------------------------------

!         CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
         !IF(myRank .EQ. masterPROC) THEN
!            CALL MPI_Win_get_attr(s1Win, MPI_WIN_MODEL, memModel, flag,&
!                  ierror)
!            CALL testForError(ierror, "MPI: Get Win-Attr")
!            WRITE(master,"(A)",advance='no') "  -- Speichermodell:"
!            IF(flag .NE. .TRUE. .OR. flag .NE. .FALSE.) THEN
!               WRITE(master,*) "Speichermodell wurde nicht erkannt"
               !CALL EXIT(-1) 
!            END IF
!            IF(memModel .EQ. MPI_WIN_UNIFIED) &
!               WRITE(master,*) "unified"
!            IF(memModel .EQ. MPI_WIN_SEPARATE) &
!               WRITE(master,*) "separate"
         !END IF
         WRITE(master,*)

         IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
               container,"Flags eingelesen")         
         ! Schritt 1 ist eine critical section und nicht-konkurrierendes
         ! Schreiben auf einem Datentyp (letzteres in Schritt 3)
         ! Ausserdem ist eine reduce-Operation vorhanden. Die ist zwar
         ! eigentlich total unwichtig, zeigt aber eine mögliche Umsetzung
         ! mit MPI3-Algorithmen
         ! n_chars wird auf Rank 0 inkrementiert, der Rest auf jedem Proc
         ! reduction muss am Ende von allen Procs an Rank 0 aufsummiert werden
         !--------------------------------------------------------------
         !--------------Wiedergabe von Schritt 1 -----------------------
         !-------------------(Z.310-Z.378)------------------------------
         IF(schritt(1)) THEN
            DO iter=1, iterationen
               !Bestimme den individuellen Abschnitt, auf dem jeder Proc laufen soll
               startIndex = (nz*2+1)/numProcs*myRank-nz
               endIndex = (nz*2+1)/numProcs*(myRank+1)-nz-1
               !trivialer Lastausgleich: letzter Proc nimmt Rest
               IF(myRank+1 .EQ. numProcs) endIndex = nz
!               WRITE(*,*) "Rang ",myRank," hat startIndex",startIndex," und "&
!                  ,"endIndex ", endIndex
!               WRITE(master,*) "Iteration: ", iter
               !Reset
               n_chars          = 0
               idx              = 0
               idxTmp           = 0
               reduction        = 0
               reductionTmp     = 0
               IF(myRank .EQ. masterProc) THEN
                  schritt1Array    = (/0,0/)
               END IF
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!               CALL testForError(ierror, "MPI: Step 1, Barrier before")
               start = gibZeit()
               !Alle Zugriffe geschehen auf masterProc, ausserdem sind nur accumulates
               !beteiligt. Es ist also nicht notwendig, innerhalb der Schleife zu syncen
               !Wegen atomarer ops nur LOCK_SHARED
               !
               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, masterProc, 0,& 
                     s1Win, ierror)
               !CALL testForError(ierror, "MPI: Step 1, Lock")
               !iz-Schleife wird unter Procs aufgeteilt
               DO iz=startIndex,endIndex
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                        CALL compLast(lasteinstellung(1))
                        !Inkrementiere n_chars und speichere alten Wert nach idx
!                        CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,&
!                              masterProc, 0, s1Win, ierror)
                        schritt1Array(1) = schritt1Array(1) + 1
                        idx = schritt1Array(1)
!                        CALL MPI_WIN_UNLOCK(masterProc, s1Win, ierror)
                        reductionTmp = reductionTmp + 1
                     END DO 
                  END DO
               END DO
               !Durchfuehrung der Reduktion
               !IF(myRank .NE. masterProc) THEN
!               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,&
!                     masterProc, 0, s1Win, ierror)
               schritt1Array(2) = schritt1Array(2) + reductionTmp
               !END IF !reduction
               CALL MPI_WIN_UNLOCK(masterProc, s1Win, ierror)
!               CALL testForError(ierror, "MPI: Step 1, Unlock")

               !Zu Testzwecken, werden reduction und n_chars gesetzt
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!               CALL testForError(ierror, "MPI: Step 1, Barrier after")
               zeitArray(iter) = (gibZeit() - start) * 1.0e6 
               IF(myRank .EQ. masterProc) THEN
                  n_chars = schritt1Array(1)
                  reduction = schritt1Array(2)
               !   n_chars       = s1Ptr(1)
               !   reduction     = s1Ptr(2)
               END IF
               IF(myRank .EQ. masterproc) &
                  verifikationsArray(iter) =  &
                     verifiziere(schrittNr(1), container, pve_grid)
            END DO !Iterationen
            IF(infosDrucken) CALL infoPrint(verbose_mittel,pve_grid,&
                  container,"Schritt 1 beendet")         
            IF(verifikationDrucken) CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 1")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt1")
         END IF!schritt1
         ! Schritt 2 ist eine atomic-section
         !--------------------------------------------------------------
         !--------------Wiedergabe von Schritt 2------------------------
         !--------------------(Z.413-Z.486)-----------------------------
         IF(schritt(2)) THEN
            n_chars = (container%nxyz*2+1)**2*(nz*2+1) !notwendig, falls Schritt 1 nicht ausgefuehrt
            verifikationsArray = .FALSE.
            !ALLOCATE(speicher1D(n_chars))
            !pve_grid%speicher1D => speicher1D
            !pve_grid%speicher1D(:) = 0.0d0
            !Bestimme den individuellen Abschnitt, auf dem jeder Proc laufen soll
            startIndex = n_chars/numProcs*myRank+1
            endIndex = n_chars/numProcs*(myRank+1)
            !trivialer Lastausgleich: letzter Proc nimmt Rest
!            CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
            IF(myRank+1 .EQ. numProcs) endIndex = n_chars
!            WRITE(*,*) "Rang ",myRank," hat startIndex",startIndex,&
!                  " und ","endIndex ", endIndex
!            CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
            IF(infosDrucken) CALL infoPrint(verbose_mittel,pve_grid,&
                  container,"Nach Initialisierung für Schritt 2")

            DO iter=1, iterationen
!               WRITE(master,*) "Iteration: ", iter
               IF(myRank .EQ. masterProc) THEN
!                  pve_grid%zaehler = 0.0d0
                  zaehlerPtr = 0.0d0
               ENDIF
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!               CALL testForError(ierror, "MPI: Step 2, First Barrier")
               start = gibZeit()
               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, masterProc,0,&
                     s2Win, ierror)
!               CALL testForError(ierror, "MPI: Step 2, Lock")
               !Aufteilung der Schleife auf Prozesse
               DO i=startIndex,endIndex
                  CALL compLast(lasteinstellung(2))
                     zaehlerPtr = zaehlerPtr + 1
!                  WRITE(*,*) "Rank: ",myRank,"i:",i," zaehler: ",zaehler
!                  pve_grid%zaehler = pve_grid%zaehler+1.0d0
               END DO !OMP
               CALL MPI_WIN_UNLOCK(masterProc, s2Win, ierror)
!               CALL testForError(ierror, "MPI: Step 2, Unlock")
               CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!               CALL testForError(ierror, "MPI: Step 2, After Barrier")
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
               pve_grid%zaehler = zaehlerPtr
               IF(myRank .EQ. masterproc) &
                  verifikationsArray(iter) =  &
                     verifiziere(schrittNr(2), container, pve_grid)
            END DO !Iterationen
            IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
                  container,"Schritt 2 beendet")
            IF(verifikationDrucken) CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 2")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt2")
         END IF !schritt2
         ! SChritt3 enthaelt scheinbar keine kritischen Bereiche
         !--------------------------------------------------------------
         !------------- Wiedergabe von Schritt 3 -----------------------
         !--------------------(Z.526-Z.1083)----------------------------
         IF(schritt(3)) THEN
            n_chars = (container%nxyz*2+1)**2*(nz*2+1) !notwendig, falls Schritt 1 nicht ausgefuehrt
            verifikationsArray = .FALSE.
            !ALLOCATE(speicher1D(n_chars))
            !ALLOCATE(speicher3D(-nx:nx,-ny:ny,-nz:nz))
            !pve_grid%speicher1D => speicher1D
            !pve_grid%speicher3D => speicher3D
            
            !Bestimme den individuellen Abschnitt, auf dem jeder Proc laufen soll
            startIndex = (nz*2+1)/numProcs*myRank-nz
            endIndex = (nz*2+1)/numProcs*(myRank+1)-nz-1
            !trivialer Lastausgleich: letzter Proc nimmt Rest
            IF(myRank+1 .EQ. numProcs) THEN
               endIndex = nz
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!            WRITE(*,*) "Rang ",myRank," hat startIndex",startIndex,&
!                        " und ","endIndex ", endIndex

            DO iter=1, iterationen
               IF(myRank .EQ. masterProc) THEN
!                  speicher1D          = 0.0d0
!                  speicher3D          = 0.0d0               
!                  stelle_max = 0
                  stelle_maxPtr         = 0.0d0
                  speicher1DPtr         = 0.0d0
                  speicher3DPtr         = 0.0d0
               END IF
!               WRITE(master,*) "Iteration: ", iter            
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!               CALL testForError(ierror, "MPI: Step 3, First Barrier")
               start = gibZeit()
               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,masterProc,0,&
                     s3Win1D, ierror)
!               CALL testForError(ierror, "MPI: Step 3D1, Lock")
               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,masterProc,0,&
                     s3Win3D, ierror)
!               CALL testForError(ierror, "MPI: Step 3D3, Lock")
               CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,masterProc,0,&
                     s3WinMax, ierror)
!               CALL testForError(ierror, "MPI: Step 3Max, Lock")
               DO iz=startIndex, endIndex
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                        stelle = 1                              + &
                              (iz+nz)*(2*nx+1)*(2*ny+1)         + &
                              (iy+ny)*(2*ny+1)                  + &
                              (ix+nx)
!                        WRITE(*,*) "Stelle: ",stelle
                        CALL compLast(lasteinstellung(3))
                        !Experiment: Da die geschriebenen Werten im Schleifendurchlauf nicht
                        !zwischen Prozessen geteilt werden muessen, behaelt jeder Prozess sein
                        !Ergebnis => weniger Kommunikation mit masterProc. 
                        !Aktuell: erstmal konservativ einbauen
                        !Eigentlich gehoert hier ein atomic hin...nun dient dieser Schritt 
                        ! aber nur als Testfall fuer konfliktfreies Lesen und Array-Windows
                        !Disp = stelle-1, weil stelle mit 1 beginnt, der Speicher aber bei 0?
                        ! Hinweis: RMA-Code ist sehr inperformant, es wuerde sich eher anbieten,
                        ! die Ergebnisse zu sammeln und dann wie in einer reduction am Ende mit
                        ! einem einigen accumulate pro Proc den MasterProc zu aktualisieren. Darauf 
                        ! wird hier aber bewusst noch verzichtet, solche Optimierungen koennen dann
                        ! in Phoenix umgesetzt werden
      
                        s3targetDisp = stelle-1
!                        WRITE(*,*) "Rang", myRank,"s3targetDisp: ",&
!                              s3targetDisp, "Stelle: ",stelle
                        !1D
!                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
!                            prec_mpi, masterProc, s3targetDisp,1,&
!                            prec_mpi, MPI_SUM, s3Win1D, ierror)
!                        CALL testForError(ierror,"MPI: Step 3, acc 1D")
                        speicher1DPtr(stelle) = speicher1DPtr(stelle)+&
                                                REAL(stelle,kind(prec))
                        !3D
!                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
!                            prec_mpi, masterProc, s3targetDisp,1,&
!                            prec_mpi, MPI_SUM, s3Win3D, ierror)
!                        CALL testForError(ierror,"MPI: Step 3, acc 3D")
                        speicher3DPtr(ix+nx+1,iy+ny+1,iz+nz+1) =       &
                              speicher3DPtr(ix+nx+1,iy+ny+1,iz+nz+1)+  &
                              REAL(stelle,kind(prec))
                        !Maximum
!                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
!                           prec_mpi, masterProc,0_MPI_ADDRESS_KIND,1,&
!                           prec_mpi,&
!                           MPI_MAX,s3WinMax, ierror)
!                        CALL testForError(ierror,"MPI: Step 3, acc Max")
                        stelle_maxPtr = max(stelle_maxPtr,&
                                                REAL(stelle,kind(prec)))
                        !evtl. zu sehr vereinfachte Max-Bestimmung 
!                        WRITE(master,*) "ix:",ix,",iy:",iy,",iz:",iz,&
!                              "1D-Stelle: ",stelle,"Max:",stelle_max
                        !Abwarten, bis lokale Buffer freigegeben sind
!                        CALL MPI_WIN_FLUSH_LOCAL(masterProc, s3Win1D,&
!                              ierror)
!                        CALL MPI_WIN_FLUSH_LOCAL(masterProc, s3Win3D,&
!                              ierror)
!                        CALL MPI_WIN_FLUSH_LOCAL(masterProc, s3WinMax,&
!                              ierror)
                     END DO 
                  END DO
               END DO !omp 
               CALL MPI_WIN_UNLOCK(masterProc, s3Win1D, ierror)
!               CALL testForError(ierror, "MPI: Step 3 1D, Unlock")
               CALL MPI_WIN_UNLOCK(masterProc, s3Win3D, ierror)
!               CALL testForError(ierror, "MPI: Step 3 3D, Unlock")
               CALL MPI_WIN_UNLOCK(masterProc, s3WinMax, ierror)
!               CALL testForError(ierror, "MPI: Step 3 Max, Unlock")
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!               CALL testForError(ierror, "MPI: Step 3, Second Barrier")
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
               !Werte fuer den Verifikationsschritt setzen
               IF(myRank .EQ. masterProc) THEN
                  stelle_max = stelle_maxPtr
                  speicher1D = speicher1DPtr
                  speicher3D = speicher3DPtr
               END IF
               IF(myRank .EQ. masterProc) THEN
                  verifikationsArray(iter) =  &
                     verifiziere(schrittNr(3), container, pve_grid)
               END IF
            END DO !iterationen
!            IF(myRank .EQ. masterProc) THEN
!               WRITE(master,*) "Ausgabe des geschriebenen Arrays"
!               DO i = 1,container%n_chars
!                  WRITE(master,*) "i: ",i,"is: ",pve_grid%speicher1D(i)
!               END DO
!            END IF!masterproc-Print
            IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
                  container,"Schritt 3 beendet")
            IF(verifikationDrucken) CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 3")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt3")
            
         END IF !schritt3
         !Nur der MasterProc hat allokiert
         IF(myRank .EQ. masterProc) THEN
            DEALLOCATE(speicher1D)
            DEALLOCATE(speicher3D)
         END IF

         CALL MPI_WIN_FREE(s1Win, ierror)
         CALL testForError(ierror, "MPI: Free s1Win")
         CALL MPI_WIN_FREE(s2Win, ierror)
         CALL testForError(ierror, "MPI: Free s2Win")
         CALL MPI_WIN_FREE(s3Win1D, ierror)
         CALL testForError(ierror, "MPI: Free s3Win1D")
         CALL MPI_WIN_FREE(s3Win3D, ierror)
         CALL testForError(ierror, "MPI: Free s3Win3D")
         CALL MPI_WIN_FREE(s3WinMax, ierror)
         CALL testForError(ierror, "MPI: Free s3WinMax")
         CALL MPI_FINALIZE(ierror)
         CALL testForError(ierror, "MPI: Finalize")

      END PROGRAM main
