! benchmarkUtilityMod.f90 
! Modul mit allen Operationen für die Zeitmessung im Benchmarkbereich
! zeitstatistik : Berechnet mean und sd aus Zeitarray
! gib Zeit      : Vereinfachung der Fortran-Zeitmessung
! 
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.
         MODULE benchmarkUtilityMod
         CONTAINS
            SUBROUTINE zeitStatistik(gemesseneZeiten, kommentar)
            USE dataTypesMod
            USE printsMod
            IMPLICIT NONE
   
            REAL(kind(prec)), DIMENSION(iterationen),  INTENT(IN)&
                                                   :: gemesseneZeiten
            CHARACTER(LEN=*), INTENT(IN)           :: kommentar
            !--------------
            REAL(kind(prec))                       :: mean, sd
            REAL(kind(prec))                       :: zwischensumme
            INTEGER                                :: i,j
            REAL(kind(prec)), DIMENSION(:), ALLOCATABLE &        
                                                   :: zeitenArray
            !-------------
            ! Wenn bereinigung wahr ist, werden bei der Statistik das Minimum und
            ! Maximum ignoriert, um system-bedingte Ausreißer besser abfangen
            ! zu können
            !-------------
            mean = 0.0d0
            sd = 0.0d0
            zwischensumme = 0.0d0
            !Entferne Extremwerte
            IF(ZeitBereinigung .AND. iterationen .GT. 2) THEN
               j = 1
               ALLOCATE(zeitenArray(iterationen-2))
               DO i = 1,iterationen         
                  !Ueberspringe die beiden Extremwerte
                  IF((i .EQ. MAXLOC(gemesseneZeiten,1)) .OR. &
                     (i .EQ. MINLOC(gemesseneZeiten,1))) CYCLE
                  zeitenArray(j) = gemesseneZeiten(i)
                  j = j + 1
               END DO !iterationen
            ELSE 
               ALLOCATE(zeitenArray(iterationen))
               zeitenArray = gemesseneZeiten
            END IF !bereinigung
   
            !-------------
            !Mittelwert
            DO i=1, SIZE(zeitenArray)
               mean = mean + zeitenArray(i)
            END DO !iterationen
            mean = mean / real(SIZE(zeitenArray))
            
            !Standardabweichung
            DO i=1, SIZE(zeitenArray)
               zwischensumme = zwischensumme + (zeitenArray(i)-mean)**2
            END DO !iterationen
   
            !--------------
            sd = sqrt(zwischensumme/real(SIZE(zeitenArray)))/real(1000)
            mean = mean/real(1000)
            CALL zeitenPrint(mean, sd, kommentar, "ms")
            DEALLOCATE(zeitenArray)
         END SUBROUTINE zeitStatistik
         !-----------------------------------------------------------------------
   
         FUNCTION gibZeit()
            USE dataTypesMod
            IMPLICIT NONE
            REAL(kind(prec))                       :: gibZeit
            INTEGER                                :: count, rate
   
            !-------------
            CALL system_clock(count,rate)
            gibZeit = dble(count)/dble(rate)
         END FUNCTION
      END MODULE benchmarkUtilityMod
