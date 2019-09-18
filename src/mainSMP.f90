! mainSMP.f90
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
!         USE mpi
         USE omp_lib

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

         REAL(kind(prec)),TARGET,ALLOCATABLE,DIMENSION(:) &
                                        :: speicher1D
         REAL(kind(prec)),TARGET,ALLOCATABLE,DIMENSION(:,:,:) &
                                        :: speicher3D
         REAL(kind(prec)),TARGET        :: stelle_max
         INTEGER,         TARGET        :: zaehler
         !--------------------------------------------------------------
         !---------------eigene Variablen-------------------------------
         !--------------------------------------------------------------
         INTEGER                        :: i, iter, stelle
         CHARACTER(len=32)              :: arg,arg2
         REAL(kind(prec))               :: start
         REAL(kind(prec))               :: meantime, sd
         REAL(kind(prec))               :: j
         REAL(kind(prec))               :: tmp
         LOGICAL                        :: testBool
         INTEGER*4                      :: iargc
         CHARACTER                      :: argv*10
         !--------------------------------------------------------------
         !--------------Hauptprogramm-----------------------------------
         !--------------------------------------------------------------     
         nz = laenge
         DO i = 1, iargc()
            CALL getarg(i, arg)
            !WRITE (*,*) arg
            SELECT CASE (arg)
               CASE ('-v', '--version')
                  WRITE(*,*) "Version ", version
                  WRITE(*,*) "Aktueller Stand: "
                  WRITE(*,*) "Schritt 1: cricical fertig,"&
                        ," reduction fertig"
                  WRITE(*,*) "Schritt 2: atomic fertig"
                  WRITE(*,*) "Schritt 3: konfliktfrei Schreiben fertig"&
                        ,", geteilte Max-/Min-Bestimmung fertig"
                  WRITE(*,*) "Verifikation: in fertig"
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
         WRITE(*,*) "Ausführung gestartet mit:"
         WRITE(*,*) " -- Kubusseitenlaenge: ", 2*laenge+1
         WRITE(*,*) " -- nz: ",nz
         WRITE(*,*) " -- Lasteinstellung: ", lasteinstellung
         WRITE(*,*) " -- Iterationszahl: ", iterationen
         WRITE(*,*) " -- Threadparameter: ", nthreads
         WRITE(*,*) " -- Entferne Extremwerte aus Zeitmessung: "&
                     ,zeitBereinigung
         WRITE(*,*)

         CALL OMP_SET_NUM_THREADS(nthreads)
!$omp parallel
!$omp master
         nthreads=omp_get_num_threads()
!$omp end master
!$omp end parallel
         WRITE(*,*) "Running OpenMP benchmark on ",nthreads," thread(s)"
         WRITE(*,*) "+++++++++++++++++++++++++++++++++++++++++++++"
         WRITE(*,*)

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

         IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
               container,"Flags eingelesen")         
         ! Schritt 1 ist eine critical section und nicht-konkurrierendes
         ! Schreiben auf einem Datentyp (letzteres in Schritt 3)
         ! Außerdem ist eine reduce-Operation vorhanden. Die ist zwar
         ! eigentlich total unwichtig, zeigt aber eine mögliche Umsetzung
         ! mit MPI3-Algorithmen
         !--------------------------------------------------------------
         !--------------Wiedergabe von Schritt 1 -----------------------
         !-------------------(Z.310-Z.378)------------------------------
         IF(schritt(1)) THEN
            DO iter=1, iterationen
!               WRITE(*,*) "Iteration: ", iter
               !Reset
               n_chars = 0
               idx = 0
               reduction = 0
               start = gibZeit()
!$omp parallel do default(none) collapse(3)     &
!$omp   shared(nx,ny,nz,n_chars)                &
!$omp   private(ix,iy,iz, idx)                  &
!$omp   firstprivate(lasteinstellung)           &
!$omp   reduction(+:reduction)
               DO iz=-nz,nz
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                       CALL compLast(lasteinstellung(1))
!$omp critical
                        n_chars = n_chars + 1 
                        idx = n_chars
!$omp end critical
                        !WRITE(*,*) "ix: ", ix, ", iy: ",iy,", iz: ",iz
                        !WRITE(*,*) n_chars, idx
                        reduction = reduction + 1 
                     END DO 
                  END DO
               END DO 
               zeitArray(iter) = (gibZeit() - start) * 1.0e6 
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
            n_chars = (container%nxyz*2+1)**2*(container%nz*2+1) !notwendig, falls Schritt 1 nicht ausgefuehrt
            verifikationsArray = .FALSE.
            !ALLOCATE(speicher1D(n_chars))
            !pve_grid%speicher1D => speicher1D
            !pve_grid%speicher1D(:) = 0.0d0
            IF(infosDrucken) CALL infoPrint(verbose_mittel,pve_grid,&
                  container,"Nach Initialisierung für Schritt 2")

            DO iter=1, iterationen
!               WRITE(*,*) "Iteration: ", iter
               pve_grid%zaehler = 0.0d0
               start = gibZeit()
!$omp parallel do default(none)                 &
!$omp   firstprivate(n_chars,lasteinstellung)   &
!$omp   shared(zaehler)
               DO i=1,n_chars
                  CALL compLast(lasteinstellung(2))
!$omp atomic
                  zaehler = zaehler + 1
!                  pve_grid%zaehler = pve_grid%zaehler+1.0d0
               END DO !OMP
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
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
            ALLOCATE(speicher1D(n_chars))
            ALLOCATE(speicher3D(-nx:nx,-ny:ny,-nz:nz))
            pve_grid%speicher1D => speicher1D
            pve_grid%speicher3D => speicher3D
            
            DO iter=1, iterationen
               speicher1D          = 0.0d0
               speicher3D          = 0.0d0
               stelle_max = 0
!               WRITE(*,*) "Iteration: ", iter            
               start = gibZeit()
!$omp parallel do default(none) collapse(3)     &
!$omp   shared(nx,ny,nz)                        &
!$omp   shared(pve_grid,stelle_max)              &
!$omp   private(ix,iy,iz,stelle)                &
!$omp   firstprivate(lasteinstellung)           
               DO iz=-nz,nz
                  DO iy=-ny,ny
                     DO ix=-nx,nx
                        stelle = 1                              + &
                              (iz+nz)*(2*nx+1)*(2*ny+1)         + &
                              (iy+ny)*(2*nx+1)                  + &
                              (ix+nx)
                        CALL compLast(lasteinstellung(3))
                        !Eigentlich gehoert hier ein atomic hin...nun dient dieser Schritt 
                        ! aber nur als Testfall fuer konfliktfreies Lesen und Array-Windows
!$omp atomic                        
                        pve_grid%speicher3D(ix,iy,iz) =                &
                              pve_grid%speicher3D(ix,iy,iz)+dble(stelle)
!$omp atomic                              
                        pve_grid%speicher1D(stelle) =                  &
                              pve_grid%speicher1D(stelle) + dble(stelle)
!$omp atomic
                        stelle_max = max(pve_grid%speicher1D(stelle),  &
                              stelle_max)
!                        WRITE(*,*) "ix: ", ix, ", iy: ",iy,", iz: ",iz,&
!                              "1D-Stelle: ",stelle 
                     END DO 
                  END DO
               END DO !omp    
               zeitArray(iter) = (gibZeit() - start) * 1.0e6
               verifikationsArray(iter) =  &
                     verifiziere(schrittNr(3), container, pve_grid)
            END DO !iterationen
            IF(infosDrucken) CALL infoPrint(verbose_voll,pve_grid,&
                  container,"Schritt 3 beendet")
            IF(verifikationDrucken) CALL verifikationsPrint(&
                  verifikationsArray, "Schritt 3")
            IF(zeitDrucken) CALL zeitStatistik(zeitArray,"Schritt3")
            
            DEALLOCATE(speicher1D)
            DEALLOCATE(speicher3D)
         END IF !schritt3

      END PROGRAM main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++ S U B R O U T I N E N +++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
