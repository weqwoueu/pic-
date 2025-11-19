module ParticleDistribution
   implicit none
   contains
      ! Function to compute the density function
      function density(x, n0, a, b) result(d)
         real(8), intent(in) :: x, n0, a, b
         real(8) :: d
         d = n0 * (1.0d0 - exp(-((x - a)**2 / b)))
      end function density

      ! Function to compute the cumulative distribution function (CDF)
      function cumulative_distribution(x, n0, a, b) result(cdf)
         real(8), intent(in) :: x, n0, a, b
         real(8) :: cdf
         ! Simple numerical integration for CDF
         ! Here, we'll use a trapezoidal rule for demonstration purposes
         integer :: n, i
         real(8) :: dx, sum, x0
         parameter (n = 1000)
         dx = x / real(n)
         sum = 0.5 * (density(0.0d0, n0, a, b) + density(x, n0, a, b))
         do i = 1, n - 1
            x0 = i * dx
            sum = sum + density(x0, n0, a, b)
         end do
         cdf = sum * dx
      end function cumulative_distribution
end module ParticleDistribution