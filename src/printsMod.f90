! printsMod.f90
! Modul mit gesammelten Printaufrufen
! infoPrint     : zu Beginn nuetzliche Infos
! helpPrint     : Auflistung verfügbarer Flags
! testForError  : Ueberprueft den Fehlerwert
!
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

         MODULE printsMod
         USE dataTypesMod
         CONTAINS
            !Einfache Debugginganweisungen
            SUBROUTINE infoPrint(verbose,pve_grid,container,pos)
               IMPLICIT NONE
               LOGICAL,DIMENSION(:),INTENT(IN)  :: verbose
               TYPE(pve_data),      INTENT(IN)  :: pve_grid
               TYPE(varContainer),  INTENT(IN)  :: container
               CHARACTER(len=*),    INTENT(IN)  :: pos
               INTEGER                          :: ix,iy,iz, i

               LOGICAL, PARAMETER               :: alles = .FALSE.
               INTEGER                          :: nxyz
               IF(verbose(1)) THEN
                  WRITE(*,*) "####### Größenausgabe #######"
                  WRITE(*,*) "pve_data: ",sizeof(pve_grid)
                  WRITE(*,*) "1D-Array: ",sizeof(pve_grid%speicher1D)
                  WRITE(*,*) "3D-Array: ",sizeof(pve_grid%speicher3D)
                  WRITE(*,*)
               END IF
               IF(verbose(2)) THEN
                  WRITE(*,*) "####### Variablenblegung in pve_grid",&
                     " #######"
                  IF(associated(pve_grid%speicher1D) .AND. alles) THEN
                     WRITE(*,*) "1d-Array: "
                     DO i=1,container%n_chars
                        WRITE(*,*) "Eintrag ",i,": ",&
                              pve_grid%speicher1D(i)
                     END DO
                     WRITE(*,*)
                  END IF
                  IF(associated(pve_grid%speicher3D) .AND. alles) THEN
                     nxyz = container%nxyz
                     WRITE(*,*) "3d-Array: "
                     DO iz=-nxyz,nxyz
                        DO iy=-nxyz,nxyz
                           DO ix=-nxyz,nxyz
                              WRITE(*,*) "Eintrag (",ix,iy,iz,"): ",&
                                    pve_grid%speicher3d(ix,iy,iz)
                           END DO
                        END DO
                     END DO
                     WRITE(*,*)
                  END IF
                  WRITE(*,*) "zaehler: ",pve_grid%zaehler
                  WRITE(*,*) "####### Variablenbelegung im ",&
                     "VarContainer #######"
                  WRITE(*,*) "n_chars: ",container%n_chars
                  WRITE(*,*) "reduction: ",container%reduction
                  WRITE(*,*) "nxyz: ", container%nxyz
                  !WRITE(*,*) "ix: ", container%ix, ", iy: ", &
                  !   container%iy, ", iz: ", container%iz
                  WRITE(*,*) "idx: ", container%idx
                  WRITE(*,*) "stelle_max: ", container%stelle_max
                  WRITE(*,*)
               END IF
               IF(verbose(3)) THEN
                  WRITE(*,*) "Ich bin hier: ",pos
               END IF
               WRITE(*,*) "--------------------------------------------"
               WRITE(*,*)
            END SUBROUTINE infoPrint

            SUBROUTINE zeitenPrint(mean, sd, kommentar, unit)
               IMPLICIT NONE
               REAL(kind(prec)),    INTENT(IN)     :: mean, sd
               CHARACTER(len=*), INTENT(IN)     :: kommentar, unit
               
               WRITE(*,*) "$$$$$$$$$$$$$$$ Zeitmessung $$$$$$$$$$$$$$$"
               WRITE(*,*) kommentar
               WRITE(*,*) kommentar," Ergebnis[",unit,"]: ",mean,"+/-",&
                     sd
               WRITE(*,*)               
            END SUBROUTINE zeitenPrint

            !Auflistung verfügbarer Flags
            SUBROUTINE helpPrint()
               WRITE(*,*) "Erlaubte Flags sind:"
               WRITE(*,*) "--xyz n: Seitenlaenge des halben Kubus"
               WRITE(*,*) "--[l1|l2|l3] n: Lastsimulierung in Schritt",&
                     " [1|2|3] mit n Iterationen"
               WRITE(*,*) "--[s1|s2|s3]: Schritt [1|2|3] wird",&
                     " ausgeführt"
               WRITE(*,*) "-t n: Ausführung mit n Threads"
               WRITE(*,*) "--file: std-Output wird in File gespeichert"
               WRITE(*,*) "-v: Zeig die aktuelle Version an"
               WRITE(*,*) 
            END SUBROUTINE helpPrint

            !Ausgabe der Ergebnisse der Verifikationsschritte
            SUBROUTINE verifikationsPrint(verifikationsArray, kommentar)
               IMPLICIT NONE
               LOGICAL, DIMENSION(:), INTENT(IN)   :: verifikationsArray
               CHARACTER(len=*),      INTENT(IN)   :: kommentar
   
               IF(ALL(verifikationsArray)) THEN
                  WRITE(*,*) "Verifikation von ",kommentar,&
                              " war erfolgreich"
               ELSE
                  WRITE(*,*) "Verifikation von ",kommentar,&
                              " war NICHT erfolgreich"
                  CALL EXIT(0)
               END IF
               WRITE(*,*)             
            END SUBROUTINE verifikationsPrint
      
            !Ueberprueft den Fehlerwert und meckert ggfs.
            SUBROUTINE testForError(ierror, msg)
               IMPLICIT NONE
               INTEGER,            INTENT(IN)              :: ierror
               CHARACTER(len=*),   INTENT(IN)              :: msg
         
               IF(ierror .NE. 0) THEN
                  WRITE(*,*) "Fehler im Bereich ", msg
                  CALL MPI_FINALIZE()
                  CALL EXIT(ierror)
               END IF
            END SUBROUTINE testForError

      !Alle Write-Aurufe sollen nur vom Master gedruckt werden
      ! Der Einfachheit halber wird hierfuer eine neue Funktion verwendet
      !@deprecated wird im Main-Programm durch die master-unit geregelt
      SUBROUTINE masterPrint(msg,myRank,nr)
      USE datatypesMod
      IMPLICIT NONE
      
      CHARACTER(len=*), INTENT(IN)              :: msg
      INTEGER,          INTENT(IN)              :: myRank
      INTEGER,OPTIONAL, INTENT(IN)              :: nr
         IF(myRank .EQ. 0) THEN
            IF(PRESENT(nr)) THEN
               WRITE(*,*)  msg,nr
            ELSE
               WRITE(*,*) msg
            END IF
         END IF
      END SUBROUTINE masterPrint

         END MODULE printsMod
