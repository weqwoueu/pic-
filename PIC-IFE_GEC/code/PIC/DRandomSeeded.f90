MODULE ModuleSimulationRandomSeed
    IMPLICIT NONE
    INTEGER, SAVE :: DRandomR1 = 1271199957
    INTEGER, SAVE :: DRandomR2 = 1013501921

CONTAINS

    SUBROUTINE InitializeSimulationRandomSeed(seed)
        INTEGER, INTENT(IN) :: seed
        INTEGER :: i, seed_size
        INTEGER, ALLOCATABLE :: seed_values(:)
        INTEGER(8) :: seed64, intrinsic_modulus

        seed64 = ABS(INT(seed, 8))
        intrinsic_modulus = INT(HUGE(0), 8) - 1_8

        CALL RANDOM_SEED(SIZE=seed_size)
        ALLOCATE(seed_values(seed_size))
        DO i = 1, seed_size
            seed_values(i) = INT(MODULO(seed64 + 104729_8 * INT(i, 8), intrinsic_modulus) + 1_8)
        END DO
        CALL RANDOM_SEED(PUT=seed_values)
        DEALLOCATE(seed_values)

        DRandomR1 = INT(MODULO(seed64 * 104729_8 + 1271199957_8, 2147483646_8) + 1_8)
        DRandomR2 = INT(MODULO(seed64 * 130363_8 + 1013501921_8, 2147483646_8) + 1_8)
    END SUBROUTINE InitializeSimulationRandomSeed

END MODULE ModuleSimulationRandomSeed


SUBROUTINE DRandom(randum)
    USE ModuleSimulationRandomSeed, ONLY: r1 => DRandomR1, r2 => DRandomR2
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT) :: randum
    DOUBLE PRECISION :: h1l, h1u, r0, r3, asc, bsc
    INTEGER :: i1, isc

    h1l = 65533.0D0
    h1u = 32767.0D0
    isc = 65536
    asc = DBLE(isc)
    bsc = asc * asc
    i1 = r1 - (r1 / isc) * isc
    r3 = h1l * DBLE(r1) + asc * h1u * DBLE(i1)
    i1 = INT(r3 / bsc)
    r3 = r3 - DBLE(i1) * bsc
    bsc = 0.5D0 * bsc
    i1 = r2 / isc
    isc = r2 - i1 * isc
    r0 = h1l * DBLE(r2) + asc * h1u * DBLE(isc)
    asc = 1.0D0 / bsc
    isc = INT(r0 * asc)
    r2 = INT(r0 - DBLE(isc) * bsc)
    r3 = r3 + DBLE(isc) + 2.0D0 * h1u * DBLE(i1)
    isc = INT(r3 * asc)
    r1 = INT(r3 - DBLE(isc) * bsc)
    randum = (DBLE(r1) + DBLE(r2) * asc) * asc
END SUBROUTINE DRandom
