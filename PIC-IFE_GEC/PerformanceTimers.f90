Module ModulePerformanceTimers
    Use Omp_Lib, Only: Omp_Get_Wtime, Omp_Get_Max_Threads
    Implicit None

    Integer, Parameter :: Perf_Move = 1
    Integer, Parameter :: Perf_Deposit = 2
    Integer, Parameter :: Perf_Field_Solve = 3
    Integer, Parameter :: Perf_Efield = 4
    Integer, Parameter :: Perf_Output = 5
    Integer, Parameter :: Perf_Count = 5

    Real(8), Save :: ComponentSeconds(Perf_Count) = 0.0D0
    Real(8), Save :: RunStart = 0.0D0
    Logical, Save :: ProfilingEnabled = .True.

Contains

    Subroutine InitializePerformanceTimers()
        Character(Len=16) :: EnvValue
        Integer :: EnvStatus

        ComponentSeconds = 0.0D0
        EnvValue = ''
        Call Get_Environment_Variable('PIC_PROFILE', EnvValue, Status=EnvStatus)
        If (EnvStatus == 0 .And. Len_Trim(EnvValue) > 0) Then
            ProfilingEnabled = Trim(EnvValue) /= '0' .And. Trim(EnvValue) /= 'false' .And. &
                               Trim(EnvValue) /= 'FALSE'
        Else
            ProfilingEnabled = .True.
        End If
        RunStart = Omp_Get_Wtime()
    End Subroutine InitializePerformanceTimers

    Real(8) Function PerformanceNow()
        PerformanceNow = Omp_Get_Wtime()
    End Function PerformanceNow

    Subroutine AddPerformanceTime(Category, StartTime)
        Integer, Intent(In) :: Category
        Real(8), Intent(In) :: StartTime

        If (.Not. ProfilingEnabled) Return
        If (Category < 1 .Or. Category > Perf_Count) Return
        ComponentSeconds(Category) = ComponentSeconds(Category) + Omp_Get_Wtime() - StartTime
    End Subroutine AddPerformanceTime

    Subroutine WritePerformanceSummary(StepCount)
        Integer, Intent(In) :: StepCount
        Integer :: UnitNumber, IoStatus, Category
        Real(8) :: TotalSeconds, AccountedSeconds, OtherSeconds, Fraction, SecondsPerStep
        Real(8) :: TotalSecondsPerStep
        Character(Len=16), Dimension(Perf_Count) :: Names

        If (.Not. ProfilingEnabled) Return

        Names = [Character(Len=16) :: 'particle_move', 'charge_deposit', 'field_solve', &
                                          'field_recovery', 'output']
        TotalSeconds = Max(0.0D0, Omp_Get_Wtime() - RunStart)
        AccountedSeconds = Sum(ComponentSeconds)
        OtherSeconds = Max(0.0D0, TotalSeconds - AccountedSeconds)

        Open(Newunit=UnitNumber, File='./OUTPUT/performance_summary.csv', Status='Replace', &
             Action='Write', Iostat=IoStatus)
        If (IoStatus /= 0) Then
            Write(*,*) 'Warning: cannot write OUTPUT/performance_summary.csv'
            Return
        End If

        Write(UnitNumber,'(A)') 'component,seconds,fraction,seconds_per_step,steps,omp_threads'
        Do Category = 1, Perf_Count
            If (TotalSeconds > 0.0D0) Then
                Fraction = ComponentSeconds(Category) / TotalSeconds
            Else
                Fraction = 0.0D0
            End If
            If (StepCount > 0) Then
                SecondsPerStep = ComponentSeconds(Category) / Real(StepCount, 8)
            Else
                SecondsPerStep = 0.0D0
            End If
            Write(UnitNumber,'(A,",",ES16.8,",",ES16.8,",",ES16.8,",",I0,",",I0)') &
                Trim(Names(Category)), ComponentSeconds(Category), Fraction, SecondsPerStep, &
                StepCount, Omp_Get_Max_Threads()
        End Do

        If (TotalSeconds > 0.0D0) Then
            Fraction = OtherSeconds / TotalSeconds
        Else
            Fraction = 0.0D0
        End If
        If (StepCount > 0) Then
            SecondsPerStep = OtherSeconds / Real(StepCount, 8)
            TotalSecondsPerStep = TotalSeconds / Real(StepCount, 8)
        Else
            SecondsPerStep = 0.0D0
            TotalSecondsPerStep = 0.0D0
        End If
        Write(UnitNumber,'(A,",",ES16.8,",",ES16.8,",",ES16.8,",",I0,",",I0)') &
            'other', OtherSeconds, Fraction, SecondsPerStep, StepCount, Omp_Get_Max_Threads()
        Write(UnitNumber,'(A,",",ES16.8,",",ES16.8,",",ES16.8,",",I0,",",I0)') &
            'total', TotalSeconds, 1.0D0, TotalSecondsPerStep, &
            StepCount, Omp_Get_Max_Threads()
        Close(UnitNumber)

        Write(*,*) 'Performance summary: OUTPUT/performance_summary.csv'
        Write(*,*) 'OpenMP max threads = ', Omp_Get_Max_Threads()
    End Subroutine WritePerformanceSummary

End Module ModulePerformanceTimers
