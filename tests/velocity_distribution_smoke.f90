Program VelocityDistributionSmoke
    Use, Intrinsic :: ieee_arithmetic, Only: ieee_is_finite
    Use ModuleSimulationRandomSeed, Only: InitializeSimulationRandomSeed
    Use ModuleVelocityDistribution
    Implicit None

    Integer, Parameter :: sample_count = 200000
    Real(8), Parameter :: electron_mass = 9.1095D-31
    Real(8), Parameter :: temperature = 11605.0D0
    Integer :: failures

    failures = 0
    Call InitializeSimulationRandomSeed(90210)

    Call CheckInitial('maxwellian', DIST_MAXWELLIAN, 2.0D0, 2.0D0, 0.04D0, failures)
    Call CheckInitial('kappa2', DIST_KAPPA, 2.0D0, 2.0D0, 0.25D0, failures)
    Call CheckInitial('kappa6', DIST_KAPPA, 6.0D0, 2.0D0, 0.06D0, failures)
    Call CheckInitial('poly2', DIST_POLYTROPIC, 2.0D0, 2.0D0, 0.04D0, failures)
    Call CheckInitial('poly3', DIST_POLYTROPIC, 2.0D0, 3.0D0, 0.04D0, failures)

    Call CheckThermalBoundary('maxwellian', DIST_MAXWELLIAN, 2.0D0, 2.0D0, 0.04D0, failures)
    Call CheckThermalBoundary('kappa2', DIST_KAPPA, 2.0D0, 2.0D0, 0.25D0, failures)
    Call CheckThermalBoundary('kappa6', DIST_KAPPA, 6.0D0, 2.0D0, 0.06D0, failures)
    Call CheckThermalBoundary('poly2', DIST_POLYTROPIC, 2.0D0, 2.0D0, 0.04D0, failures)
    Call CheckThermalBoundary('poly3', DIST_POLYTROPIC, 2.0D0, 3.0D0, 0.04D0, failures)

    If (failures /= 0) Then
        Print*, 'velocity distribution smoke test FAILED: ', failures
        Stop 1
    End If
    Print*, 'velocity distribution smoke test PASSED'

Contains

    Subroutine CheckInitial(label, distribution, kappa, poly_gamma, variance_tolerance, failure_count)
        Character(Len=*), Intent(In) :: label
        Integer(4), Intent(In) :: distribution
        Real(8), Intent(In) :: kappa, poly_gamma, variance_tolerance
        Integer, Intent(InOut) :: failure_count
        Integer :: i
        Real(8) :: vx, vy, vz, sum_value, sum_square, mean_value
        Real(8) :: normalized_variance, expected_variance, vmax, largest_speed
        Real(8) :: shape_sum, shape_mean, shape_expected, kappa_width

        sum_value = 0.0D0
        sum_square = 0.0D0
        shape_sum = 0.0D0
        largest_speed = 0.0D0
        Do i = 1, sample_count
            Call SampleInitialVelocity(distribution, electron_mass, temperature, kappa, &
                                       poly_gamma, vx, vy, vz)
            If (.Not. ieee_is_finite(vx) .Or. .Not. ieee_is_finite(vy) .Or. &
                .Not. ieee_is_finite(vz)) Then
                failure_count = failure_count + 1
                Print*, Trim(label), ': non-finite initial velocity'
                Return
            End If
            sum_value = sum_value + vx
            sum_square = sum_square + vx * vx
            largest_speed = Max(largest_speed, Abs(vx))
            If (distribution == DIST_KAPPA) Then
                kappa_width = (2.0D0 * kappa - 3.0D0) * &
                              KB_VD * temperature / electron_mass
                shape_sum = shape_sum + vx * vx / (kappa_width + vx * vx)
            End If
        End Do

        expected_variance = KB_VD * temperature / electron_mass
        mean_value = sum_value / Dble(sample_count)
        normalized_variance = (sum_square / Dble(sample_count) - mean_value**2) / expected_variance
        Print '(A,2(A,ES12.4))', Trim(label), ' mean/vth=', &
              mean_value / Sqrt(expected_variance), ' variance/expected=', normalized_variance

        If (Abs(mean_value) / Sqrt(expected_variance) > 0.04D0) Then
            failure_count = failure_count + 1
            Print*, Trim(label), ': initial mean check failed'
        End If
        If (Abs(normalized_variance - 1.0D0) > variance_tolerance) Then
            failure_count = failure_count + 1
            Print*, Trim(label), ': initial variance check failed'
        End If

        If (distribution == DIST_KAPPA) Then
            shape_mean = shape_sum / Dble(sample_count)
            shape_expected = 1.0D0 / (2.0D0 * kappa)
            Print '(A,2(A,F10.6))', Trim(label), ' Kappa-shape mean=', shape_mean, &
                  ' expected=', shape_expected
            If (Abs(shape_mean - shape_expected) > 0.005D0) Then
                failure_count = failure_count + 1
                Print*, Trim(label), ': Kappa exponent check failed'
            End If
        End If

        If (distribution == DIST_POLYTROPIC) Then
            vmax = Sqrt(2.0D0 * poly_gamma / (poly_gamma - 1.0D0) * expected_variance)
            If (largest_speed > vmax * (1.0D0 + 1.0D-12)) Then
                failure_count = failure_count + 1
                Print*, Trim(label), ': cutoff check failed'
            End If
        End If
    End Subroutine CheckInitial

    Subroutine CheckThermalBoundary(label, distribution, kappa, poly_gamma, &
                                    variance_tolerance, failure_count)
        Character(Len=*), Intent(In) :: label
        Integer(4), Intent(In) :: distribution
        Real(8), Intent(In) :: kappa, poly_gamma, variance_tolerance
        Integer, Intent(InOut) :: failure_count
        Integer :: i
        Real(8) :: vx, vy, vz, transformed, sum_value, sum_square
        Real(8) :: beta, vmax, transform_power, mean_value, variance_value
        Real(8) :: tangent_sum, tangent_square, tangent_mean, tangent_variance
        Real(8) :: expected_variance, largest_normal

        beta = electron_mass / (2.0D0 * KB_VD * temperature)
        sum_value = 0.0D0
        sum_square = 0.0D0
        tangent_sum = 0.0D0
        tangent_square = 0.0D0
        largest_normal = 0.0D0
        Do i = 1, sample_count
            Call SampleThermalBoundaryVelocity(distribution, electron_mass, temperature, &
                kappa, poly_gamma, vx, vy, vz)
            If (vx <= 0.0D0 .Or. .Not. ieee_is_finite(vx) .Or. &
                .Not. ieee_is_finite(vy) .Or. .Not. ieee_is_finite(vz)) Then
                failure_count = failure_count + 1
                Print*, Trim(label), ': thermal normal velocity check failed'
                Return
            End If

            Select Case (distribution)
            Case (DIST_MAXWELLIAN)
                transformed = Exp(-beta * vx * vx)
            Case (DIST_KAPPA)
                transformed = (1.0D0 + beta * vx * vx / (kappa - 1.5D0))**(1.0D0 - kappa)
            Case (DIST_POLYTROPIC)
                vmax = Sqrt(poly_gamma / ((poly_gamma - 1.0D0) * beta))
                transform_power = (poly_gamma + 1.0D0) / (2.0D0 * (poly_gamma - 1.0D0))
                transformed = Max(0.0D0, 1.0D0 - (vx / vmax)**2)**transform_power
            End Select
            sum_value = sum_value + transformed
            sum_square = sum_square + transformed * transformed
            tangent_sum = tangent_sum + vy + vz
            tangent_square = tangent_square + vy * vy + vz * vz
            largest_normal = Max(largest_normal, vx)
        End Do

        mean_value = sum_value / Dble(sample_count)
        variance_value = sum_square / Dble(sample_count) - mean_value**2
        Print '(A,2(A,F10.6))', Trim(label), ' boundary-CDF mean=', mean_value, &
              ' variance=', variance_value
        If (Abs(mean_value - 0.5D0) > 0.01D0 .Or. &
            Abs(variance_value - 1.0D0 / 12.0D0) > 0.01D0) Then
            failure_count = failure_count + 1
            Print*, Trim(label), ': thermal CDF check failed'
        End If

        expected_variance = KB_VD * temperature / electron_mass
        tangent_mean = tangent_sum / Dble(2 * sample_count)
        tangent_variance = tangent_square / Dble(2 * sample_count) - tangent_mean**2
        tangent_variance = tangent_variance / expected_variance
        Print '(A,2(A,ES12.4))', Trim(label), ' tangent mean/vth=', &
              tangent_mean / Sqrt(expected_variance), &
              ' tangent variance/expected=', tangent_variance
        Print '(A,A,ES12.4)', Trim(label), ' largest normal/vth=', &
              largest_normal / Sqrt(expected_variance)
        If (Abs(tangent_mean) / Sqrt(expected_variance) > 0.04D0 .Or. &
            Abs(tangent_variance - 1.0D0) > variance_tolerance) Then
            failure_count = failure_count + 1
            Print*, Trim(label), ': thermal tangential distribution check failed'
        End If
    End Subroutine CheckThermalBoundary

End Program VelocityDistributionSmoke
