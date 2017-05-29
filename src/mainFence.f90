! mainFence.f90
!Vereinfachte Darstellung des Kommunikationsmusters mit OpenMP in Phoenix
!
!Mittels Cmd-Args können einige Laufzeitgroessen angepasst werden, um eine
!Entwicklung verfolgen zu können 
!
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

      PROGRAM main         
         USE mpi

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
         INTEGER, DIMENSION(schritt1ArrayLen)                          &
                                        :: schritt1Array
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
         INTEGER                        :: prec_mpi
         INTEGER                        :: sizeOfprec_mpi
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
         verifikationsArray     = .FALSE.
         zaehler                = 0
         
         container%n_chars      => n_chars
         container%reduction    => reduction
         container%nxyz         => laenge
         container%nz           => nz
         container%idx          => idx
         container%stelle_max   => stelle_max
         pve_grid%zaehler       => zaehler
      
         schritt1Array(1) = n_chars
         schritt1Array(2) = reduction
         !Allokiere Schritt3-Speicher fuer nachfolgende Win-Erzeugung
         ALLOCATE(speicher1D((container%nxyz*2+1)**2*&
            (container%nz*2+1)))
         ALLOCATE(speicher3D(-nx:nx,-ny:ny,-nz:nz))
         pve_grid%speicher1D => speicher1D
         pve_grid%speicher3D => speicher3D
      
         !------------MPI-RMA-Initialisierung--------------------------
         !Bestimme Groesse der Speicherbereiche
         CALL MPI_TYPE_SIZE(MPI_INTEGER, sizeOfInteger, ierror)
         CALL testForError(ierror, "MPI: get Int-Size")
         CALL MPI_TYPE_SIZE(MPI_REAL, sizeOfReal, ierror)
         CALL testForError(ierror, "MPI: get Real-Size")
         CALL MPI_TYPE_SIZE(prec_mpi, sizeOfprec_mpi,ierror)
         CALL testForError(ierror, "MPI: get prec_mpi-Size")

         sizeOfArray            = sizeOfInteger *schritt1ArrayLen
         sizeOfSchritt2         = sizeOfInteger *ONE
         sizeOfSchritt3D1       = sizeOfprec_mpi*&
                              ((container%nxyz*2+1)**2+container%nz*2+1)
         sizeOfSchritt3D3       = sizeOfprec_mpi*&
                              ((container%nxyz*2+1)**2+container%nz*2+1)
         sizeOfSchritt3Max      = sizeOfprec_mpi*ONE
         !Einrichtung der MPI-Windows, auf die dann geschrieben/gelesen
         !Schritt1
         CALL MPI_WIN_CREATE(schritt1Array, sizeOfArray,sizeOfInteger,&
                   MPI_INFO_NULL,MPI_COMM_WORLD,s1Win,ierror)
         CALL testForError(ierror, "MPI: Create s1Win")
         !Schritt2
         CALL MPI_WIN_CREATE(zaehler, sizeOfSchritt2,sizeOfInteger,&
                   MPI_INFO_NULL,MPI_COMM_WORLD,s2Win,ierror)
         CALL testForError(ierror, "MPI: Create s2Win")
         !Schritt3
         CALL MPI_WIN_CREATE(speicher1D, sizeOfSchritt3D1,&
                  sizeOfprec_mpi,MPI_INFO_NULL,MPI_COMM_WORLD,s3Win1D, &
                  ierror)
         CALL testForError(ierror, "MPI: Create s3Win1D")
         CALL MPI_WIN_CREATE(speicher3D, sizeOfSchritt3D3,&
                  sizeOfprec_mpi,MPI_INFO_NULL,MPI_COMM_WORLD,s3Win3D, &
                  ierror)
         CALL testForError(ierror, "MPI: Create s3Win3D")
         CALL MPI_WIN_CREATE(zaehler, sizeOfSchritt2,&
                  sizeOfprec_mpi, MPI_INFO_NULL,MPI_COMM_WORLD,&
                  s3WinMax, ierror)
         !CALL MPI_WIN_CREATE(stelle_max, sizeOfSchritt3Max,&
         !         sizeOfprec_mpi, MPI_INFO_NULL,MPI_COMM_WORLD,&
         !         s3WinMax, ierror)
         CALL testForError(ierror, "MPI: Create s3WinMax")
         WRITE(master,*) "Fenstererstellung erfolgreich abgeschlossen"
         !-------------------------------------------------------------

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
               !Reset
               n_chars = 0
               idx = 0
               idxTmp = 0
               reduction = 0
               schritt1Array = (/0,0/)
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               start = gibZeit()
               CALL MPI_WIN_FENCE(0,s1Win,ierror)
               !iz-Schleife wird unter Procs aufgeteilt
               DO iz=startIndex,endIndex
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                        CALL compLast(lasteinstellung(1))
                        !Inkrementiere n_chars und speichere alten Wert nach idx
                        CALL MPI_Fetch_and_op(ONE,idx,MPI_INTEGER,  &
                              masterProc,0_MPI_ADDRESS_KIND,           &
                              MPI_SUM,s1Win,ierror)
                        !Wert in idx ist um 1 zu klein
                        idx = idx + ONE
                        schritt1Array(2) = schritt1Array(2) + ONE
                     END DO 
                  END DO
               END DO
               !Durchfuehrung der Reduktion
               IF(myRank .NE. masterProc) THEN
                  CALL MPI_ACCUMULATE(schritt1Array(2), 1, MPI_INTEGER,&
                        masterProc, 1_MPI_ADDRESS_KIND, 1, MPI_INTEGER,&
                        MPI_SUM, s1Win, ierror)
               END IF !reduction
               CALL MPI_WIN_FENCE(0,s1Win,ierror)

               !Zu Testzwecken, werden reduction und n_chars gesetzt
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               zeitArray(iter) = (gibZeit() - start) * 1.0e6 
               n_chars = schritt1Array(1)
               reduction = schritt1Array(2)
               IF(myRank .EQ. masterproc) &
               verifikationsArray(iter) =  &
                     verifiziere(schrittNr(1), container, pve_grid)
            END DO !Iterationen
            IF(infosDrucken) CALL infoPrint(verbose_mittel,pve_grid,&
                  container,"Schritt 1 beendet")         
            IF(verifikationDrucken .AND. myRank .EQ. masterProc)&
                CALL verifikationsPrint(verifikationsArray, "Schritt 1")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt1")
         END IF!schritt1
         ! Schritt 2 ist eine atomic-section
         !--------------------------------------------------------------
         !--------------Wiedergabe von Schritt 2------------------------
         !--------------------(Z.413-Z.486)-----------------------------
         IF(schritt(2)) THEN
            n_chars = (container%nxyz*2+1)**2*(nz*2+1) !notwendig, falls Schritt 1 nicht ausgefuehrt
            verifikationsArray = .FALSE.
            !Bestimme den individuellen Abschnitt, auf dem jeder Proc laufen soll
            startIndex = n_chars/numProcs*myRank+1
            endIndex = n_chars/numProcs*(myRank+1)
            !trivialer Lastausgleich: letzter Proc nimmt Rest
            IF(myRank+1 .EQ. numProcs) endIndex = n_chars
            IF(infosDrucken) CALL infoPrint(verbose_mittel,pve_grid,&
                  container,"Nach Initialisierung für Schritt 2")

            DO iter=1, iterationen
               pve_grid%zaehler = 0.0d0
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               start = gibZeit()
               CALL MPI_WIN_FENCE(0,s2Win,ierror)
               !Aufteilung der Schleife auf Prozesse
               DO i=startIndex,endIndex
                  CALL compLast(lasteinstellung(2))
                  CALL MPI_ACCUMULATE(ONE, 1, MPI_INTEGER,&
                      masterProc, 0_MPI_ADDRESS_KIND, 1, MPI_INTEGER,&
                      MPI_SUM,s2Win, ierror)
               END DO !OMP
               CALL MPI_WIN_FENCE(0,s2Win,ierror)
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
               CALL MPI_WIN_FENCE(0,s2Win,ierror)
               IF(myRank .EQ. masterproc) &
               verifikationsArray(iter) =  &
                     verifiziere(schrittNr(2), container, pve_grid)
            END DO !Iterationen
            IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
                  container,"Schritt 2 beendet")
            IF(verifikationDrucken .AND. myRank .EQ. masterproc)&
                CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 2")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt2")
         END IF !schritt2
         ! SChritt3 enthaelt scheinbar keine kritischen Bereiche
         !--------------------------------------------------------------
         !------------- Wiedergabe von Schritt 3 -----------------------
         !--------------------(Z.526-Z.1083)----------------------------
         IF(schritt(3)) THEN
            n_chars = (container%nxyz*2+1)**2*(nz*2+1)!notwendig, falls Schritt 1 nicht ausgefuehrt
            verifikationsArray = .FALSE.
            
            !Bestimme den individuellen Abschnitt, auf dem jeder Proc laufen soll
            startIndex = (nz*2+1)/numProcs*myRank-nz
            endIndex = (nz*2+1)/numProcs*(myRank+1)-nz-1
            !trivialer Lastausgleich: letzter Proc nimmt Rest
            IF(myRank+1 .EQ. numProcs) THEN
               endIndex = nz
            END IF
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

            DO iter=1, iterationen
               speicher1D          = 0.0d0
               speicher3D          = 0.0d0
               stelle_max = 0
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               start = gibZeit()
               CALL MPI_WIN_FENCE(0,s3Win1D, ierror)
               CALL MPI_WIN_FENCE(0,s3Win3D, ierror)
               CALL MPI_WIN_FENCE(0,s3WinMax,ierror)
               DO iz=startIndex, endIndex
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                        stelle = 1                              + &
                              (iz+nz)*(2*nx+1)*(2*ny+1)         + &
                              (iy+ny)*(2*ny+1)                  + &
                              (ix+nx)
                        CALL compLast(lasteinstellung(3))
                        s3targetDisp = stelle-1
                        !1D
                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
                            prec_mpi, masterProc, s3targetDisp,1,&
                            prec_mpi, MPI_SUM, s3Win1D, ierror)
                        CALL MPI_WIN_FENCE(0,s3Win1D, ierror)
                        !3D
                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
                            prec_mpi, masterProc, s3targetDisp,1,&
                            prec_mpi, MPI_SUM, s3Win3D, ierror)
                        CALL MPI_WIN_FENCE(0,s3Win3D, ierror)
                        !Maximum
                        CALL MPI_ACCUMULATE(REAL(stelle,kind(prec)), 1,&
                           prec_mpi, masterProc,0_MPI_ADDRESS_KIND,1,&
                           prec_mpi,&
                           MPI_MAX,s3WinMax, ierror)
                        !Abwarten, bis lokale Buffer freigegeben sind
                        CALL MPI_WIN_FENCE(0,s3WinMax, ierror)
                     END DO 
                  END DO
               END DO !omp 

               !Last proc uses more Fences, so everyone else needs to match them
               IF(myRank .LT. numProcs-1) THEN
                  fenceIndex=MODULO((nz*2+1),numProcs)&
                              *(2*ny+1)*(2*nx+1)
                  DO i = 1,fenceIndex
                        CALL MPI_WIN_FENCE(0,s3Win1D, ierror)
                        CALL MPI_WIN_FENCE(0,s3Win3D, ierror)
                        CALL MPI_WIN_FENCE(0,s3WinMax, ierror)
                  END DO
               END IF
                  
               CALL MPI_WIN_FENCE(0,s3Win1D, ierror)
               CALL MPI_WIN_FENCE(0,s3Win3D, ierror)
               CALL MPI_WIN_FENCE(0,s3WinMax,ierror)
               CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
               !Werte fuer den Verifikationsschritt setzen
               IF(myRank .EQ. masterproc) &
               verifikationsArray(iter) =  &
                     verifiziere(schrittNr(3), container, pve_grid)
            END DO !iterationen
            IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
                  container,"Schritt 3 beendet")
            IF(verifikationDrucken .AND. myRank .EQ. masterproc)&
                CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 3")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt3")
            
         END IF !schritt3
         DEALLOCATE(speicher1D)
         DEALLOCATE(speicher3D)

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
