! computationMod.f90
! Mod mit einer Methode, um Rechenlast zu erzeugen
! Das soll einer realistischeren Arbeitsumgebung entsprechen,
! da sonst nur Kommunikation betrieben wird
! Staerke der Rechenlast kann durch Parameter angegeben werden
!      
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

      MODULE computationMod
         CONTAINS
            SUBROUTINE compLastAlt(einstellung)
               INTEGER, INTENT(IN)      :: einstellung
               INTEGER                  :: i,j
               REAL,DIMENSION(einstellung)&
                                        ::  zwischenspeicher
               IF(einstellung .GT. 0) THEN  
                  zwischenspeicher(1) = 0
                  DO i = 1, einstellung
                     zwischenspeicher(i) = real(i)
                     DO j = 1,i
                        zwischenspeicher(i) = zwischenspeicher(i) + &
                              zwischenspeicher(j)/real(einstellung)
                     END DO            
                  END DO
               END IF
            END SUBROUTINE compLastAlt

            SUBROUTINE compLast(einstellung)
               INTEGER, INTENT(IN)      :: einstellung
               INTEGER                  :: i
               REAL(KIND=kind(1.0d0))   :: j
               
               j = 0.0d0
               DO i=1,einstellung
                  j = j + real(i)
               END DO
            END SUBROUTINE compLast

            FUNCTION verifiziere(aktuellerSchritt, container, pve_grid)
         USE dataTypesMod
         IMPLICIT NONE
         INTEGER,                  INTENT(IN)      :: aktuellerSchritt
         TYPE(varContainer),       INTENT(IN)      :: container
         TYPE(pve_data),           INTENT(IN)      :: pve_grid
      
         LOGICAL                                   :: verifiziere
         CHARACTER(32), PARAMETER                  :: einleitung &
                                                = "----> Verifikation: "
         INTEGER                                   :: i,ix,iy,iz
         INTEGER                                   :: counter
         
         verifiziere = .TRUE.
         SELECT CASE(aktuellerSchritt)
            CASE(schrittNr(1))
               IF(container%n_chars .NE. (container%nxyz*2+1)**2*(&
                  container%nz*2+1)) THEN
                  verifiziere = .FALSE.
                  WRITE(*,*) einleitung, "n_chars falsch"
               END IF
              IF(container%reduction .NE. container%n_chars) THEN
                  verifiziere = .FALSE.
                  WRITE(*,*) einleitung, "Reduktion falsch"
               END IF
            !-------------------------------------------------------------
            CASE(schrittNr(2))
               IF(container%n_chars .NE. pve_grid%zaehler) THEN
                  verifiziere = .FALSE.
                  WRITE(*,*) einleitung, "zaehler falsch"
               END IF
            !-------------------------------------------------------------   
            CASE(schrittNr(3))
               counter = 0
               DO iz = -container%nz, container%nz
                  DO iy = -container%nxyz, container%nxyz
                     DO ix = -container%nxyz, container%nxyz
                        counter = counter +1
                        IF(pve_grid%speicher3D(ix,iy,iz) .NE. counter)&
                           THEN
                           verifiziere = .FALSE.
                           EXIT
                        END IF
                     END DO
                  END DO
               END DO
               IF(verifiziere .NEQV. .TRUE.) &
                           WRITE(*,*) einleitung, "3D-Speicher falsch"
               DO i = 1,container%n_chars
                  IF(pve_grid%speicher1D(i) .NE. i) THEN
                     verifiziere = .FALSE.
                     WRITE(*,*) einleitung, "1D-Speicher falsch"
                     EXIT
                  END IF 
               END DO
               !Thread-Sicherheit von max/min im Tracker klären!
               IF(container%stelle_max .NE. container%n_chars) THEN
                     verifiziere = .FALSE.
                     WRITE(*,*) einleitung, "stelle_max falsch"
               END IF
               !-------------------------------------------------------------
            CASE DEFAULT
               WRITE(*,*)"Zulässige Werte für schritt: {1,2,3}"
               CALL EXIT(0)
         END SELECT
      END FUNCTION verifiziere
      END MODULE computationMod
