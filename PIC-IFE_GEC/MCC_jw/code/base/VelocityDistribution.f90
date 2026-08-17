Module ModuleVelocityDistribution
    Implicit None

    ! The paper model prescribes one-dimensional component distributions.
    ! For Kappa, f(v) is proportional to
    ! [1 + v**2 / ((2*kappa-3)*kB*T/m)]**(-kappa).
    ! This is the convention consistent with the boundary formula, density
    ! closure, and <v**2> = kB*T/m used by the diagnostics.

    Integer(4), Parameter :: DIST_MAXWELLIAN = 1
    Integer(4), Parameter :: DIST_KAPPA = 2
    Integer(4), Parameter :: DIST_POLYTROPIC = 3
    Integer(4), Parameter :: BOUNDARY_THERMAL = 1
    Integer(4), Parameter :: BOUNDARY_SPECULAR = 2

    Real(8), Parameter :: PI_VD = 3.141592653589793238D0
    Real(8), Parameter :: KB_VD = 1.3807D-23

    Integer(4), Save :: ElectronDistribution = DIST_MAXWELLIAN
    Integer(4), Save :: LeftBoundaryMode = BOUNDARY_THERMAL
    Real(8), Save :: ElectronKappa = 2.0D0
    Real(8), Save :: ElectronPolytropicGamma = 2.0D0
    Logical, Save :: DistributionConfigInitialized = .False.

    Interface
        Subroutine DRandom(value)
            Real(8), Intent(Out) :: value
        End Subroutine DRandom
    End Interface

Contains

    Subroutine InitializeVelocityDistributionConfig()
        Character(Len=64) :: value
        Integer :: status, ios

        If (DistributionConfigInitialized) Return

        value = ''
        Call Get_Environment_Variable('PIC_ELECTRON_DISTRIBUTION', value, Status=status)
        If (status == 0 .And. Len_Trim(value) > 0) Then
            Call LowerCaseInPlace(value)
            Select Case (Trim(AdjustL(value)))
            Case ('maxwellian', 'maxwell', '1')
                ElectronDistribution = DIST_MAXWELLIAN
            Case ('kappa', '2')
                ElectronDistribution = DIST_KAPPA
            Case ('polytropic', 'poly', '3')
                ElectronDistribution = DIST_POLYTROPIC
            Case Default
                Print*, 'Invalid PIC_ELECTRON_DISTRIBUTION: ', Trim(value)
                Stop 2
            End Select
        End If

        value = ''
        Call Get_Environment_Variable('PIC_KAPPA', value, Status=status)
        If (status == 0 .And. Len_Trim(value) > 0) Then
            Read(value, *, Iostat=ios) ElectronKappa
            If (ios /= 0) Then
                Print*, 'Invalid PIC_KAPPA: ', Trim(value)
                Stop 2
            End If
        End If

        value = ''
        Call Get_Environment_Variable('PIC_POLYTROPIC_GAMMA', value, Status=status)
        If (status == 0 .And. Len_Trim(value) > 0) Then
            Read(value, *, Iostat=ios) ElectronPolytropicGamma
            If (ios /= 0) Then
                Print*, 'Invalid PIC_POLYTROPIC_GAMMA: ', Trim(value)
                Stop 2
            End If
        End If

        value = ''
        Call Get_Environment_Variable('PIC_LEFT_BOUNDARY', value, Status=status)
        If (status == 0 .And. Len_Trim(value) > 0) Then
            Call LowerCaseInPlace(value)
            Select Case (Trim(AdjustL(value)))
            Case ('thermal', 'thermal-reservoir', 'diffusive', '1')
                LeftBoundaryMode = BOUNDARY_THERMAL
            Case ('specular', 'energy-conserving', '2')
                LeftBoundaryMode = BOUNDARY_SPECULAR
            Case Default
                Print*, 'Invalid PIC_LEFT_BOUNDARY: ', Trim(value)
                Stop 2
            End Select
        End If

        If (ElectronKappa <= 1.5D0) Then
            Print*, 'PIC_KAPPA must be greater than 1.5'
            Stop 2
        End If
        If (ElectronPolytropicGamma <= 1.0D0 .Or. ElectronPolytropicGamma > 3.0D0) Then
            Print*, 'PIC_POLYTROPIC_GAMMA must satisfy 1 < gamma <= 3'
            Stop 2
        End If

        DistributionConfigInitialized = .True.
        Print*, 'PIC electron distribution = ', Trim(DistributionName(ElectronDistribution))
        If (ElectronDistribution == DIST_KAPPA) Print*, 'PIC kappa = ', ElectronKappa
        If (ElectronDistribution == DIST_POLYTROPIC) Then
            Print*, 'PIC polytropic gamma = ', ElectronPolytropicGamma
        End If
        Print*, 'PIC left boundary = ', Trim(BoundaryName(LeftBoundaryMode))
    End Subroutine InitializeVelocityDistributionConfig

    Subroutine SampleInitialVelocity(distribution, mass, temperature, kappa, poly_gamma, vx, vy, vz)
        Integer(4), Intent(In) :: distribution
        Real(8), Intent(In) :: mass, temperature, kappa, poly_gamma
        Real(8), Intent(Out) :: vx, vy, vz

        Call ValidateThermalInputs(mass, temperature)
        Select Case (distribution)
        Case (DIST_MAXWELLIAN)
            vx = SampleMaxwellianComponent(mass, temperature)
            vy = SampleMaxwellianComponent(mass, temperature)
            vz = SampleMaxwellianComponent(mass, temperature)
        Case (DIST_KAPPA)
            vx = SampleKappaComponent(mass, temperature, kappa)
            vy = SampleKappaComponent(mass, temperature, kappa)
            vz = SampleKappaComponent(mass, temperature, kappa)
        Case (DIST_POLYTROPIC)
            vx = SamplePolytropicComponent(mass, temperature, poly_gamma)
            vy = SamplePolytropicComponent(mass, temperature, poly_gamma)
            vz = SamplePolytropicComponent(mass, temperature, poly_gamma)
        Case Default
            Print*, 'Unknown initial velocity distribution: ', distribution
            Stop 2
        End Select
    End Subroutine SampleInitialVelocity

    Subroutine SampleThermalBoundaryVelocity(distribution, mass, temperature, kappa, poly_gamma, vx, vy, vz)
        Integer(4), Intent(In) :: distribution
        Real(8), Intent(In) :: mass, temperature, kappa, poly_gamma
        Real(8), Intent(Out) :: vx, vy, vz
        Real(8) :: beta, radial, angle, uniform_value, vmax, exponent_value

        Call ValidateThermalInputs(mass, temperature)
        beta = mass / (2.0D0 * KB_VD * temperature)

        Select Case (distribution)
        Case (DIST_MAXWELLIAN)
            uniform_value = SafeUniform()
            vx = Sqrt(-Log(uniform_value) / beta)
            radial = Sqrt(-Log(SafeUniform()) / beta)
            angle = 2.0D0 * PI_VD * SafeUniform()
            vy = radial * Cos(angle)
            vz = radial * Sin(angle)
        Case (DIST_KAPPA)
            If (kappa <= 1.5D0) Then
                Print*, 'Kappa thermal sampling requires kappa > 1.5'
                Stop 2
            End If
            uniform_value = SafeUniform()
            vx = Sqrt((kappa - 1.5D0) / beta * &
                      (uniform_value**(-1.0D0 / (kappa - 1.0D0)) - 1.0D0))
            vy = SampleKappaComponent(mass, temperature, kappa)
            vz = SampleKappaComponent(mass, temperature, kappa)
        Case (DIST_POLYTROPIC)
            If (poly_gamma <= 1.0D0 .Or. poly_gamma > 3.0D0) Then
                Print*, 'Polytropic thermal sampling requires 1 < gamma <= 3'
                Stop 2
            End If
            vmax = Sqrt(poly_gamma / ((poly_gamma - 1.0D0) * beta))
            exponent_value = 2.0D0 * (poly_gamma - 1.0D0) / (poly_gamma + 1.0D0)
            vx = vmax * Sqrt(1.0D0 - SafeUniform()**exponent_value)
            vy = SamplePolytropicComponent(mass, temperature, poly_gamma)
            vz = SamplePolytropicComponent(mass, temperature, poly_gamma)
        Case Default
            Print*, 'Unknown thermal boundary distribution: ', distribution
            Stop 2
        End Select
    End Subroutine SampleThermalBoundaryVelocity

    Function SampleMaxwellianComponent(mass, temperature) Result(value)
        Real(8), Intent(In) :: mass, temperature
        Real(8) :: value, sigma, radius, angle

        sigma = Sqrt(KB_VD * temperature / mass)
        radius = Sqrt(-2.0D0 * Log(SafeUniform()))
        angle = 2.0D0 * PI_VD * SafeUniform()
        value = sigma * radius * Cos(angle)
    End Function SampleMaxwellianComponent

    Function SampleKappaComponent(mass, temperature, kappa) Result(value)
        Real(8), Intent(In) :: mass, temperature, kappa
        Real(8) :: value, degrees_freedom, scale, chi_square

        If (kappa <= 1.5D0) Then
            Print*, 'Kappa initialization requires kappa > 1.5'
            Stop 2
        End If

        ! This Student-t construction samples the one-dimensional Kappa
        ! component exactly without truncating the suprathermal tail.
        degrees_freedom = 2.0D0 * kappa - 1.0D0
        scale = Sqrt((2.0D0 * kappa - 3.0D0) / degrees_freedom * &
                     KB_VD * temperature / mass)
        chi_square = 2.0D0 * SampleGammaShape(0.5D0 * degrees_freedom)
        value = scale * SampleStandardNormal() / Sqrt(chi_square / degrees_freedom)
    End Function SampleKappaComponent

    Function SamplePolytropicComponent(mass, temperature, poly_gamma) Result(value)
        Real(8), Intent(In) :: mass, temperature, poly_gamma
        Real(8) :: value, beta_shape, sample_left, sample_right, vmax, fraction

        If (poly_gamma <= 1.0D0 .Or. poly_gamma > 3.0D0) Then
            Print*, 'Polytropic initialization requires 1 < gamma <= 3'
            Stop 2
        End If

        beta_shape = (poly_gamma + 1.0D0) / (2.0D0 * (poly_gamma - 1.0D0))
        sample_left = SampleGammaShape(beta_shape)
        sample_right = SampleGammaShape(beta_shape)
        fraction = sample_left / (sample_left + sample_right)
        vmax = Sqrt(2.0D0 * poly_gamma / (poly_gamma - 1.0D0) * &
                     KB_VD * temperature / mass)
        value = vmax * (2.0D0 * fraction - 1.0D0)
    End Function SamplePolytropicComponent

    Function SampleStandardNormal() Result(value)
        Real(8) :: value

        value = Sqrt(-2.0D0 * Log(SafeUniform())) * &
                Cos(2.0D0 * PI_VD * SafeUniform())
    End Function SampleStandardNormal

    Function SampleGammaShape(shape) Result(value)
        Real(8), Intent(In) :: shape
        Real(8) :: value, d, c, normal_value, candidate, uniform_value

        If (shape < 1.0D0) Then
            Print*, 'Internal gamma sampler requires shape >= 1'
            Stop 2
        End If

        d = shape - 1.0D0 / 3.0D0
        c = 1.0D0 / Sqrt(9.0D0 * d)
        Do
            normal_value = SampleStandardNormal()
            candidate = 1.0D0 + c * normal_value
            If (candidate <= 0.0D0) Cycle
            candidate = candidate**3
            uniform_value = SafeUniform()
            If (uniform_value < 1.0D0 - 0.0331D0 * normal_value**4) Exit
            If (Log(uniform_value) < 0.5D0 * normal_value**2 + &
                d * (1.0D0 - candidate + Log(candidate))) Exit
        End Do
        value = d * candidate
    End Function SampleGammaShape

    Function SafeUniform() Result(value)
        Real(8) :: value

        Call DRandom(value)
        value = Max(value, Tiny(1.0D0))
        value = Min(value, 1.0D0 - Epsilon(1.0D0))
    End Function SafeUniform

    Subroutine ValidateThermalInputs(mass, temperature)
        Real(8), Intent(In) :: mass, temperature

        If (mass <= 0.0D0 .Or. temperature <= 0.0D0) Then
            Print*, 'Velocity sampling requires positive mass and temperature'
            Stop 2
        End If
    End Subroutine ValidateThermalInputs

    Function DistributionName(distribution) Result(name)
        Integer(4), Intent(In) :: distribution
        Character(Len=16) :: name

        Select Case (distribution)
        Case (DIST_MAXWELLIAN)
            name = 'maxwellian'
        Case (DIST_KAPPA)
            name = 'kappa'
        Case (DIST_POLYTROPIC)
            name = 'polytropic'
        Case Default
            name = 'unknown'
        End Select
    End Function DistributionName

    Function BoundaryName(boundary_mode) Result(name)
        Integer(4), Intent(In) :: boundary_mode
        Character(Len=16) :: name

        Select Case (boundary_mode)
        Case (BOUNDARY_THERMAL)
            name = 'thermal'
        Case (BOUNDARY_SPECULAR)
            name = 'specular'
        Case Default
            name = 'unknown'
        End Select
    End Function BoundaryName

    Subroutine LowerCaseInPlace(value)
        Character(Len=*), Intent(InOut) :: value
        Integer :: i, code

        Do i = 1, Len(value)
            code = Iachar(value(i:i))
            If (code >= Iachar('A') .And. code <= Iachar('Z')) Then
                value(i:i) = Achar(code + Iachar('a') - Iachar('A'))
            End If
        End Do
    End Subroutine LowerCaseInPlace

End Module ModuleVelocityDistribution
