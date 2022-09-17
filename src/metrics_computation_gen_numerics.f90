!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
module gen_numerics
!-------------------------------------------------------------------------------

! Module containing generic numerical subroutine, like non linear root finding,
! bisection rules and so forth

use data_types, only: double

implicit none

private
public zbren, brent

contains


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
real(kind = double) function brent(ax, bx, cx, func, tol, xmin, ierror)
!-------------------------------------------------------------------------------

  real(kind = double), intent(in) :: ax, bx, cx, tol
  real(kind = double), external :: func
  real(kind = double), intent(out) :: xmin
  integer, intent(out) :: ierror

  real(kind = double) :: a, b
  integer :: ierror2 = 0

  a = min(ax, cx)
  b = max(ax, cx)

  xmin = fmin(f, a, b, tol, ierror)
  brent = f(xmin)
  ierror = max(ierror, ierror2)

  contains

  real(kind = double) function f(x)
    real(kind = double) :: x
    integer :: ierror
    f = func(x, ierror)
    ierror2 = max(ierror2, ierror)
  end function

!-------------------------------------------------------------------------------
end function brent
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine zbren(func_value, func, x1, x2, tol, xmin, ierror)
!-------------------------------------------------------------------------------

  real(kind = double), intent(in):: func_value
  real(kind = double), external :: func
  real(kind = double), intent(in) :: x1, x2
  real(kind = double), intent(in) :: tol

  real(kind = double), intent(out) :: xmin
  integer, intent(out) :: ierror

  xmin = zeroin(f, x1, x2, tol, ierror)

  contains

  real(kind = double) function f(x)
    real(kind = double) :: x
    f = func(x) - func_value
  end function

!----------------------------------------------------------------
end subroutine zbren
!----------------------------------------------------------------
!----------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
real(kind = double) function fmin(f, lower, upper, tol, ierror)
!-------------------------------------------------------------------------------

  ! An approximation x to the point where f attains a minimum on the interval
  ! (lower,upper) is determined.
  !
  ! Input:
  !   f      function subprogram which evaluates  f(x)  for any  x in the
  !          interval  (lower,upper)
  !   lower  left endpoint of initial interval
  !   upper  right endpoint of initial interval
  !   tol    desired length of the interval of uncertainty of the final result
  !          (.ge. 0.0_double)
  !
  ! Output:
  !   fmin   abcissa approximating the point where f attains a minimum
  !   ierror error indicator
  !
  ! The method used is a combination of golden section search and successive
  ! parabolic interpolation.  Convergence is never much slower than that for a
  ! fibonacci search.  If f has a continuous second derivative which is positive
  ! at the minimum (which is not at lower or upper), then convergence is
  ! superlinear, and usually of the order of about  1.324.
  !
  ! The function f is never evaluated at two points closer together than
  ! eps*abs(fmin) + (tol/3), where eps is  approximately the square root of the
  ! relative machine precision.  If f is a unimodal function and the computed
  ! values of f are always unimodal when separated by at least
  ! eps*abs(x)+(tol/3), then fmin approximates the abcissa of the global minimum
  ! of f on the interval lower,upper with an error less than
  ! 3*eps*abs(fmin)+tol.  If f is not unimodal, then fmin may approximate a
  ! local, but perhaps non-global, minimum to the same accuracy.
  !
  ! This function subprogram is a slightly modified version of the algol 60
  ! procedure localmin given in Richard Brent, Algorithms for Minimization
  ! without Derivatives, Prentice-Hall, Inc. (1973).
  !
  ! Source: https://netlib.org/fmm/fmin.f
  !
  ! Modified:
  !
  !   Brian J Smith <brian-j-smith@uiowa.edu>
  !   2022-09-16

  real(kind = double), external :: f
  real(kind = double) :: lower, upper, tol
  integer, intent(out) :: ierror

  real(kind = double) :: a, b, c, d, e, eps, p, q, r, tol1, tol2, u, v, w, xm
  real(kind = double) :: fu, fv, fw, fx, x

  integer :: i, maxiter = 1000
  logical :: golden_section

  ierror = 1

  !
  ! is the squared inverse of the golden ratio
  !
  c = 0.5_double * (3.0_double - sqrt(5.0_double))

  !
  ! eps is approximately the square root of the relative machine precision.
  !
  eps = sqrt(epsilon(eps))

  !
  ! initialization
  !
  a = lower
  b = upper
  v = a + c * (b - a)
  w = v
  x = v
  e = 0.0_double
  fx = f(x)
  fv = fx
  fw = fx

  !
  ! main loop starts here
  !
  do i = 1, maxiter

    xm = 0.5_double * (a + b)
    tol1 = eps * abs(x) + tol / 3.0_double
    tol2 = 2.0_double * tol1

    !
    ! check stopping criterion
    !
    if (abs(x - xm) <= (tol2 - 0.5_double * (b - a))) then
      ierror = 0
      exit
    end if

    !
    ! is golden-section necessary
    !
    golden_section = abs(e) <= tol1

    if (.not. golden_section) then

      !
      ! fit parabola
      !
      r = (x - w) * (fx - fv)
      q = (x - v) * (fx - fw)
      p = (x - v) * q - (x - w) * r
      q = 2.0_double * (q - r)
      if (q > 0.0_double) p = -p
      q = abs(q)
      r = e
      e = d

      golden_section = abs(p) >= abs(0.5_double * q * r) .or. &
        p <= q * (a - x) .or. p >= q * (b - x)

    end if

    if (golden_section) then

      if (x >= xm) then
        e = a - x
      else
        e = b - x
      end if
      d = c * e

    else

      !
      ! a parabolic interpolation step
      !
      d = p / q
      u = x + d
      !
      ! f must not be evaluated too close to lower or upper
      !
      if ((u - a) < tol2 .or. (b - u) < tol2) d = sign(tol1, xm - x)

    end if

    !
    ! f must not be evaluated too close to x
    !
    if (abs(d) >= tol1) then
      u = x + d
    else
      u = x + sign(tol1, d)
    end if

    fu = f(u)

    !
    ! update  a, b, v, w, and x
    !
    if (fu <= fx) then

      if (u >= x) then
        a = x
      else
        b = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if (u < x) then
        a = u
      else
        b = u
      end if

      if (fu <= fw .or. w == x) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if (fu <= fv .or. v == x .or. v == w) then
        v = u
        fv = fu
      end if

    end if

  end do

  fmin = x

!-------------------------------------------------------------------------------
end function fmin
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
real(kind = double) function zeroin(f, lower, upper, tol, ierror)
!-------------------------------------------------------------------------------

  ! A zero of the function f(x) is computed in the interval lower,upper.
  !
  ! Input:
  !   lower  left endpoint of initial interval
  !   upper  right endpoint of initial interval
  !   f      function subprogram which evaluates f(x) for any x in the interval
  !          lower,upper
  !   tol    desired length of the interval of uncertainty of the final result
  !          (.ge.0.)
  !
  ! Output:
  !   zeroin abscissa approximating a zero of f in the interval lower,upper
  !   ierror error indicator
  !
  ! It is assumed that f(lower) and f(upper) have opposite signs.  This is
  ! checked, and an error message is printed if this is not satisfied.  Zeroin
  ! returns a zero x in the given interval lower,upper to within a tolerance
  ! 4*macheps*abs(x)+tol, where macheps is the  relative machine precision
  ! defined as the smallest representable number such that 1.+macheps .gt. 1.
  !
  ! This function subprogram is a slightly modified translation of the algol 60
  ! procedure zero given in Richard Brent, Algorithms for Minimization without
  ! Derivatives, Prentice-Hall, Inc. (1973).
  !
  ! Source: https://netlib.org/go/zeroin.f
  !
  ! Modified:
  !
  !   Brian J Smith <brian-j-smith@uiowa.edu>
  !   2022-09-16

  real(kind = double), external :: f
  real(kind = double), intent(in) :: lower, upper, tol
  integer, intent(out) :: ierror

  real(kind = double) :: a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s
  integer :: i, maxiter = 1000

  ierror = 1

  eps = epsilon(eps)

  a = lower
  b = upper
  fa = f(a)
  fb = f(b)
  fc = fb

  if (fa * fb > 0.0_double) then
    return
  end if

  do i = 1, maxiter

    if (fb * fc > 0.0_double) then
      c = a
      fc = fa
      d = b - a
      e = d
    end if

    if (abs(fc) < abs(fb)) then
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
    end if

    tol1 = 2.0_double * eps * abs(b) + 0.5_double * tol
    xm = 0.5_double * (c - b)

    if (abs(xm) <= tol1 .or. fb == 0.0_double) then
      ierror = 0
      exit
    end if

    !
    ! see if a bisection is forced
    !
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s = fb / fa
      if (a == c) then
        !
        ! linear interpolation
        !
        p = 2.0_double * xm * s
        q = 1.0_double - s
      else
        !
        ! inverse quadratic interpolation
        !
        q = fa / fc
        r = fb / fc
        p = s * (2.0_double * xm * q * (q - r) - (b - a) * (r - 1.0_double))
        q = (q - 1.0_double) * (r - 1.0_double) * (s - 1.0_double)
      end if
      if (p <= 0.0_double) then
        p = -p
      else
        q = -q
      end if
      if (2.0_double * p >= &
        min(3.0_double * xm * q - abs(tol1 * q), abs(e * q))) &
      then
        d = xm
        e = d
      else
        e = d
        d = p / q
      end if
    else
      d = xm
      e = d
    end if

    a = b
    fa = fb
    if (abs(d) > tol1) then
      b = b + d
    else
      b = b + sign(tol1, xm)
    end if
    fb = f(b)

  end do

  zeroin = b

!-------------------------------------------------------------------------------
end function zeroin
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
end module gen_numerics
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
