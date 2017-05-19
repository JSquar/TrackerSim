! datatypesMod.f90
! Modul mit definierten Datentypen
!
! Copyright (C) 2017 Jannek Squar
!
! This software may be modified and distributed under the terms
! of the MIT license.  See the LICENSE file for details.

      MODULE dataTypesMod
            REAL*8, PARAMETER           :: prec=1d-13

         TYPE pve_data
            REAL(kind(prec)), POINTER   :: speicher3D(:,:,:)    =>NULL()
            REAL(kind(prec)), POINTER   :: speicher1D(:)        =>NULL()
            INTEGER,          POINTER   :: zaehler      
         END TYPE pve_data
      
         TYPE varContainer
            INTEGER,POINTER             :: n_chars, reduction
            INTEGER,POINTER             :: nxyz, idx,nz!,ix,iy,iz
            REAL(kind(prec)), POINTER   :: stelle_max
            !Ergebnisse der Verifikation der Ergebnisse
         END TYPE varContainer

         !------------------Steuer-Konstanten---------------------------
         !--------------------------------------------------------------
         !--------------------------------------------------------------
         INTEGER, PARAMETER             :: verbosity           = 3
         INTEGER, PARAMETER             :: iterationen         = 5
         LOGICAL                        :: infosDrucken        = .FALSE.
         LOGICAL                        :: zeitDrucken         = .TRUE.
         LOGICAL                        :: verifikationDrucken = .TRUE.
         LOGICAL, DIMENSION(iterationen):: verifikationsArray  = .FALSE.
         LOGICAL, PARAMETER             :: zeitBereinigung     = .TRUE.
         CHARACTER(len=32),PARAMETER    :: version             = "1.3.0"


         !------------------Steuer-Parameter----------------------------
         !--------------------------------------------------------------
         !--------------------------------------------------------------
         INTEGER                        :: master               = 6
         REAL(kind(prec)), DIMENSION(iterationen) :: zeitArray
         INTEGER,TARGET                 :: laenge               = 32
         INTEGER, DIMENSION(3)          :: lasteinstellung      = 0   
         LOGICAL, DIMENSION(3)          :: schritt              = .FALSE.

         LOGICAL,DIMENSION(verbosity)   :: verbose_voll         = .TRUE.
         LOGICAL,DIMENSION(verbosity)   :: verbose_mittel       = & 
      (/.FALSE.,.TRUE.,.TRUE./)
         LOGICAL,DIMENSION(verbosity)   :: verbose_wenig        = & 
      (/.FALSE.,.FALSE.,.TRUE./)
         !Verbose(1): Datentypgroesse; Verbose(2): Varbelegung
         !Verbose(3): Positionsmeldung
         !extra Abfrage, da pve_data jedes mal an infoPrint uebergeben
         !werden muss...das könnte Zeiten verfälschen
         INTEGER                        :: nthreads     = 1 
         INTEGER,PARAMETER,DIMENSION(3) :: schrittNr    = (/1,2,3/) 

      END MODULE dataTypesMod        
