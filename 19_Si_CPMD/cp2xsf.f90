!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
! modified by Gerardo Ballabio and Carlo Cavazzoni
! on Thu Jul 22 18:57:26 CEST 2004
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! -------------------------------------------------------------------
!   this routine writes the crystal structure in XSF format
!   from a FPMD output files
! -------------------------------------------------------------------
PROGRAM xsf_struct
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER, PARAMETER       :: maxsp = 20
  REAL(kind=DP), PARAMETER :: convert = 0.529177d0, alat = 1.0d0

  iNTEGER                    :: natoms, nsp, nat(maxsp), species(maxsp)
  INTEGER                    :: ounit, cunit, punit, funit
  INTEGER                    :: nr1, nr2, nr3, ispin
  INTEGER, ALLOCATABLE       :: ityp(:)
  REAL(kind=DP), ALLOCATABLE :: rho(:,:,:), tau(:,:), force(:,:)
  REAL(kind=DP)              :: at(3, 3)
  REAL(kind=DP)              :: rhof, rhomax, rhomin, rhoc(6)

  CHARACTER(len=80) :: prefix, filepp, fileout, output
  CHARACTER(len=80) :: filecel, filepos, filefor
  LOGICAL           :: lcharge, lforces, ldynamics
  INTEGER           :: iprefix(256)
  INTEGER           :: nframes
  INTEGER           :: ios
  INTEGER           :: irep, irex, irey, irez, irot, ipbc

  REAL(kind=DP) :: x, y, z, fx, fy, fz
  INTEGER       :: i, j, n, ix, iy, iz

  NAMELIST /inputpp/ prefix, filepp, fileout, output, lcharge, lforces, ldynamics, nr1, nr2, nr3, ispin, nat, nsp, species, nframes

  ! default values
  prefix = 'pp'
  filepp = 'CHARGE_DENSITY'
  fileout = ' '
  output = 'xsf'
  lcharge = .false.
  lforces = .false.
  ldynamics = .false.
  nsp = 0
  nat(:) = 0
  nframes = 1
  irep = 0
  irex = 1
  irey = 1
  irez = 1
  irot = 1
  ipbc = 0

  ! read namelist
  READ(*, inputpp, iostat=ios)
  print*,nat(1),nframes

  ! set file names
  filecel = TRIM(prefix) // '.cel'
  filepos = TRIM(prefix) // '.pos'
  filefor = TRIM(prefix) // '.for'
  IF (fileout == ' ') THEN
     IF (ldynamics) THEN
        fileout = 'out.axsf'
     ELSE
        fileout = 'out.xsf'
     END IF
  END IF
  iprefix(:) = 0
  DO i = 1, LEN(TRIM(prefix))
     iprefix(i) = ichar(prefix(i:i))
  END DO

  IF (ldynamics .AND. nframes < 2) THEN
     WRITE(*,*) 'Error: dynamics requested, but only one frame'
     STOP
  END IF
  IF (.NOT. ldynamics) nframes = 1

  IF (ldynamics .AND. lcharge) THEN
     WRITE(*,*) 'Error: dynamics with charge density non supported'
     STOP
  END IF

  IF (nsp > maxsp) THEN
     WRITE(*,*) 'Error: too many atomic species'
     STOP
  END IF

  natoms = 0
  DO i = 1, nsp
     natoms = natoms + nat(i)
  END DO
  ALLOCATE(tau(3, natoms))
  ALLOCATE(ityp(natoms))
  IF (lforces) ALLOCATE(force(3, natoms))
  IF (lcharge) ALLOCATE(rho(nr1, nr2, nr3))

  ! assign species to each atom
  n = 0
  DO i = 1, nsp
     DO j = 1, nat(i)
        n = n + 1
        ityp(n) = species(i)
        WRITE(*,*) 'ityp: ', ityp(n)
     END DO
  END DO

  ! open files
  ounit = 12
  cunit = 11
  punit = 10
  funit = 13
  OPEN(ounit, file=fileout, status='unknown')
  OPEN(punit, file=filepos, status='old')
  OPEN(cunit, file=filecel, status='old')
  if (lforces) OPEN(funit, file=filefor, status='old')

  IF (ldynamics) WRITE(ounit,*) 'ANIMSTEPS', nframes
  WRITE(ounit,*) 'CRYSTAL'

  DO n = 1, nframes

     ! read cell vectors
     READ(cunit,*)
     DO i = 1, 3
        READ(cunit,*) (at(j, i), j=1,3)
     END DO
     ! rescale
     at(:, :) = at(:, :) * alat * convert
     WRITE(*,'(3f10.6)') ((at(i, j), i=1,3), j=1,3)

     ! read atomic coordinates
     READ(punit,*)
     IF (lforces) READ(funit,*)
     DO i = 1, natoms
        ! positions are in Angstroms
        READ(punit,*) x, y, z
        tau(1, i) = x * alat * convert
        tau(2, i) = y * alat * convert
        tau(3, i) = z * alat * convert

        IF (lforces) THEN
           READ (funit,*) fx, fy, fz
           force(1, i) = fx
           force(2, i) = fy
           force(3, i) = fz
        END IF
     END DO

     CALL write_xsf( ldynamics, lforces, ounit, n, at, natoms, &
                    ityp, tau, force )
     !CALL write_pdb( at, tau, nsp, nat, ityp, irep, irex, irey, &
     !               irez, irot, ipbc, iprefix )
  END DO

  CLOSE(punit)
  CLOSE(cunit)
  CLOSE(funit)

  IF (lcharge) THEN

     ! -------------------------------------------------------------------
     !   this routine writes the 3D scalar field (i.e. uniform mesh of points)
     !   in XSF format using the FFT mesh (i.e. fast write)
     ! -------------------------------------------------------------------
     ! XSF scalar-field header
     WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
     WRITE(ounit,'(a)') '3D_PWSCF'
     WRITE(ounit,'(a)') 'DATAGRID_3D_UNKNOWN'

     OPEN(punit, file=filepp, status='old')
     WRITE(*,*) nr1, nr2, nr3, ispin

     ! read charge density from file
     ! note: must transpose
     DO ix = 1, nr1
        DO iy = 1, nr2
           DO iz = 1, nr3
              READ(punit,*) rhof
              rho(ix, iy, iz) = rhof
           END DO
        END DO
     END DO
     CLOSE(punit)
     rhomax = MAXVAL(rho(:,:,:))
     rhomin = MINVAL(rho(:,:,:))
     WRITE(*,*) nr1, nr2, nr3, ispin, rhomin, rhomax

     ! write data to output file

     ! mesh dimensions
     WRITE(ounit,*) nr1, nr2, nr3
     ! origin
     WRITE(ounit,'(3f10.6)') 0.0, 0.0, 0.0
     ! lattice vectors
     WRITE(ounit,'(3f10.6)') ((at(i, j), i=1,3), j=1,3)
     ! charge density
     WRITE(ounit,'(6e13.5)') &
          (((rho(ix, iy, iz), ix=1,nr1), iy=1,nr2), iz=1,nr3)

     WRITE(ounit,'(a)') 'END_DATAGRID_3D'
     WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_3D'

  END IF

  CLOSE(ounit)

  DEALLOCATE(tau)
  DEALLOCATE(ityp)
  IF (lforces) DEALLOCATE(force)
  IF (lcharge) DEALLOCATE(rho)

  STOP
END PROGRAM xsf_struct

SUBROUTINE write_xsf( ldynamics, lforces, ounit, n, at, natoms, &
                     ityp, tau, force )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  LOGICAL       :: ldynamics, lforces
  INTEGER       :: ounit, n, natoms, ityp(natoms)
  REAL(kind=DP) :: at(3, 3), tau(3, natoms), force(3, natoms)

  INTEGER :: i, j

  ! write cell
  IF (ldynamics) THEN
     WRITE(ounit,*) 'PRIMVEC', n
  ELSE
     WRITE(ounit,*) 'PRIMVEC'
  END IF
  WRITE(ounit,'(2(3f15.9/),3f15.9)') at
  IF (ldynamics) THEN
     WRITE(ounit,*) 'CONVVEC', n
     WRITE(ounit,'(2(3f15.9/),3f15.9)') at
  END IF

  ! write atomic coordinates
  IF (ldynamics) THEN
     WRITE(ounit,*) 'PRIMCOORD', n
  ELSE
     WRITE(ounit,*) 'PRIMCOORD'
  END IF
  WRITE(ounit,*) natoms, 1
  DO i = 1, natoms
     IF (lforces) THEN
        WRITE (ounit,'(i3,3x,3f15.9,1x,3f12.5)') ityp(i), &
              (tau(j, i), j=1,3), (force(j, i), j=1,3)
     ELSE
        WRITE (ounit,'(i3,3x,3f15.9,1x,3f12.5)') ityp(i), &
              (tau(j, i), j=1,3)
     END IF
  END DO
END SUBROUTINE write_xsf

SUBROUTINE write_grd( rho, nr1, nr2, nr3, ds, ht, prefix )
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(0.0d0)

  INTEGER          :: nr1, nr2, nr3, ds
  REAL(kind=DP)    :: rho( nr1, nr2, nr3 ), ht(3, 3)
  CHARACTER(len=*) :: prefix

  INTEGER       :: nr1d, nr2d, nr3d
  REAL(kind=DP) :: celldm(6), ps(3,3)

  REAL(kind=DP), PARAMETER   :: pi = 3.14159265358979323846d0
  REAL(kind=DP), PARAMETER   :: au_to_amg = 0.529177d0
  REAL(kind=DP), ALLOCATABLE :: rho_ds( :, :, : )
  INTEGER :: i, j, k

  DO j = 1, 3
     DO k = 1, j
        ps(j,k) = 0.0d0
        DO i = 1, 3
           ps(j,k) = ps(j,k) + ht(j,i) * ht(k,i)
        END DO
     END DO
  END DO
  DO i = 1, 3
     celldm(i) = sqrt( ps(i,i) ) * au_to_amg
  END DO
  celldm(4) = ps(3,2) / sqrt( ps(3,3) * ps(2,2) )
  celldm(4) = acos( celldm(4) ) * 180.0d0 / pI
  celldm(5) = ps(3,1) / sqrt( ps(3,3) * ps(1,1) )
  celldm(5) = acos( celldm(5) ) * 180.0d0 / pi
  celldm(6) = ps(2,1) / sqrt( ps(2,2) * ps(1,1) )
  celldm(6) = acos( celldm(6) ) * 180.0d0 / pi

  nr1d = nr1
  nr2d = nr2
  nr3d = nr3

  ALLOCATE( rho_ds( nr1d, nr2d, nr3d ) )

  DO i = 1, nr1d
     DO j = 1, nr2d
        DO k = 1, nr3d
           rho_ds(i,j,k) = rho(i,j,k)
        END DO
     END DO
  END DO

  OPEN(UNIT=3, FILE=TRIM(prefix)//'.grd')

  WRITE(3,*) 'xsf converted'
  WRITE(3,*) '(1p,e12.5)'
  WRITE(3,fmt="(6f9.3)" ) ( celldm(i), i=1,6 )
  WRITE(3,fmt="(3i5)" ) nr1d - 1, nr2d - 1, nr3d - 1
  WRITE(3,fmt="(7i5)" ) 1, 0, 0, 0, nr1d - 1, nr2d - 1, nr3d - 1
  DO k = 1, nr3d
     DO j = 1, nr2d
        WRITE(3,fmt='(1p,e12.5)') rho_ds( 1:nr1d, j, k )
     END DO
  END DO

  CLOSE(3)

  DEALLOCATE( rho_ds )

  RETURN
END SUBROUTINE
