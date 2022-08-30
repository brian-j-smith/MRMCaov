module toms611
! Module that contains the Transactions for Mathematical Sofware procedures to
! implement the algorithm 611, known as TOMS611.
! NOTES:
!       - Some functions that did not do anything have been eliminated (timestamp for example)
!       - Very archaic syntax have been cleared, more work needs to be done
!       - All write is done to file pu, which is defined in the calling program (in element iv(21))
!       - Module created and files modified at UC in May-June 2005
!       - Note that only the procedures which are after a "public::" will be callable from other procedures
!         outside the modulus

use data_types
implicit none

private

public:: deflt ! initialization of the optimizer
public:: humsl ! implementation of Newton-Raphson, trust zone version style 611
public:: sumsl ! implementation of Quasi-Newton, trust zone version style 611


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assst( iv, liv, lv, v )
!
!*******************************************************************************
!
!! ASSST: assess a candidate step.
!
!
!   purpose
!
!        this subroutine is called by an unconstrained minimization
!     routine to assess the next candidate step.  it may recommend one
!     of several courses of action, such as accepting the step, recom-
!     puting it using the same or a new quadratic model, or halting due
!     to convergence or false convergence.  see the return code listing
!     below.
!
!  Parameters:
!
!  iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
! liv (in)  length of iv array.
!  lv (in)  length of v array.
!   v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!
!   iv values referenced
!
!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assst.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmaxs) but do not
!                       evaulate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmaxs)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precdence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 1 if v(f) has been restored and x should be
!             restored to its initial value, to 2 if x should be saved,
!             to 3 if x should be restored from the saved value, and to
!             0 otherwise.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assst sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.
!
!   v values referenced
!
! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assst returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d*step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
!             i.e., for v(nreduc) >= 0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm) > v(lmaxs), and if
!             v(preduc) <= v(sctol) * abs(v(f0)), then assst re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assst repeats this test with v(preduc) computed for
!             a step of length v(lmaxs) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             newton step.  if assst is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmaxs) for
!             use in the singular convervence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd)*step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0_double is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (in) scaled relative change in x caused by step, computed
!             (e.g.) by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
!             assst returns with iv(irc) = 8 or 9.
! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a newton step
!             (v(stppar) = 0) having v(reldx) <= v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assst returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx) <= v(xftol),
!             then assst returns with iv(irc) = 12.
!
!  notes
!
!   application and usage restrictions
!
!        this routine is called as part of the nl2sol (nonlinear
!     least-squares) package.  it may be used in any unconstrained
!     minimization solver that uses dogleg, goldfeld-quandt-trotter,
!     or levenberg-marquardt steps.
!
!   algorithm notes
!
!        see (1) for further discussion of the assessing and model
!     switching strategies.  while nl2sol considers only two models,
!     assst is designed to handle any number of models.
!
!   usage notes
!
!        on the first call of an iteration, only the i/o variables
!     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
!     v(preduc) need have been initialized.  between calls, no i/o
!     values execpt step, x, iv(model), v(f) and the stopping toler-
!     ances should be changed.
!        after a return for convergence or false convergence, one can
!     change the stopping tolerances and call assst again, in which
!     case the stopping tests will be repeated.
!
!  Reference:
!
!    J E Dennis, David Gay, and R E Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    M J D Powell,
!    A Fortran Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!   history
!
!        john dennis designed much of this routine, starting with
!     ideas in (2). roy welsch suggested the model switching strategy.
!        david gay and stephen peters cast this subroutine into a more
!     portable form (winter 1977), and david gay cast it into its
!     present form (fall 1978).
!
! NOTE: LP uchicago, august 2005. This subroutine is a complete mess, it is surprising that it works at
!       all. Some day a hero will clean it up. For now, I have only eliminated the obsolescent features
!       The whole web structure needs to be disentangled and some of the branches need to become internal
!       procedures to eliminate the webbing.
implicit none

  integer liv, lv
  integer iv(liv)
  real(kind=double) v(lv)
  logical goodx
  integer i, nfc
  real(kind=double) emax, emaxs, gts, rfac1, xmax
  real(kind=double) half, one, onep2, two
  integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0
  integer gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall
  integer nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn
  integer rdfcmx, reldx, restor, rfctol, sctol, stage, stglim
  integer stppar, switch, toobig, tuner1, tuner2, tuner3, xctol
  integer xftol, xirc

  parameter ( half=0.5_double, one=1._double, onep2=1.2_double, two=2._double)
  parameter ( irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7 )
  parameter ( radinc=8, restor=9, stage=10, stglim=11, switch=12 )
  parameter ( toobig=2, xirc=13)
  parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18 )
  parameter (f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4 )
  parameter (incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7 )
  parameter (radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32 )
  parameter (sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28 )
  parameter (xctol=33, xftol=34)
!
  nfc = iv(nfcall)
  iv(switch) = 0
  iv(restor) = 0
  rfac1 = one
  goodx = .true.
  i = iv(irc)

  ! depending on the value of i (iv(irc) ) branch out - LP
  ! ideally this step will not have the gotos anymore, but instead a number of operations followed
  ! by call to functions for operations common to the different branches
  select case(i)
  case(1)
        goto 20
  case(2)
        goto 30
  case(3:4)
        goto 10
  case(5)
        goto 40
  case(6)
        goto 280
  case(7:11)
        goto 220
  case(12)
        goto 170
  case default ! the value needs to be set, and then we return
        iv(irc) = 13
        return
  end select

!
!   initialize for new iteration
!
 10   iv(stage) = 1
  iv(radinc) = 0
  v(flstgd) = v(f0)
  if (iv(toobig) == 0) go to 110
     iv(stage) = -1
     iv(xirc) = i
     go to 60
!
!   step was recomputed with new model or smaller radius
!   first decide which
!
 20   if (iv(model) /= iv(mlstgd)) go to 30
!         old model retained, smaller radius tried
!         do not consider any more new models this iteration
     iv(stage) = iv(stglim)
     iv(radinc) = -1
     go to 110
!
!   a new model is being tried.  decide whether to keep it.
!
 30   iv(stage) = iv(stage) + 1
!
!      now we add the possibiltiy that step was recomputed with
!      the same model, perhaps because of an oversized step.
!
 40   if (iv(stage) > 0) go to 50
!
!         step was recomputed because it was too big.
!
     if (iv(toobig) /= 0) go to 60
!
!         restore iv(stage) and pick up where we left off.
!
     iv(stage) = -iv(stage)
     i = iv(xirc)

     select case(i) ! LP
     case(1)
          goto 20
     case(2)
          goto 30
     case(3:4)
          goto 110
     case(5)
          goto 70
     end select

 50   if (iv(toobig) == 0) go to 70
!
!   handle oversize step
!
  if (iv(radinc) > 0) go to 80
     iv(stage) = -iv(stage)
     iv(xirc) = iv(irc)

 60      v(radfac) = v(decfac)
     iv(radinc) = iv(radinc) - 1
     iv(irc) = 5
     iv(restor) = 1
     return

 70   if (v(f) < v(flstgd)) go to 110
!
!     the new step is a loser.  restore old model.
!
  if (iv(model) == iv(mlstgd)) go to 80
     iv(model) = iv(mlstgd)
     iv(switch) = 1
!
!      restore step, etc. only if a previous step decreased v(f).
!
 80   if (v(flstgd) >= v(f0)) go to 110
     iv(restor) = 1
     v(f) = v(flstgd)
     v(preduc) = v(plstgd)
     v(gtstep) = v(gtslst)
     if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
     v(dstnrm) = v(dstsav)
     nfc = iv(nfgcal)
     goodx = .false.
!
 110  v(fdif) = v(f0) - v(f)
  if (v(fdif) > v(tuner2) * v(preduc)) go to 140
  if(iv(radinc)>0) go to 140
!
!         no (or only a trivial) function decrease
!         so try new model or smaller radius
!
     if (v(f) < v(f0)) go to 120
          iv(mlstgd) = iv(model)
          v(flstgd) = v(f)
          v(f) = v(f0)
          iv(restor) = 1
          go to 130
 120     iv(nfgcal) = nfc
 130     iv(irc) = 1
     if (iv(stage) < iv(stglim)) go to 160
          iv(irc) = 5
          iv(radinc) = iv(radinc) - 1
          go to 160
!
!  Nontrivial function decrease achieved
!
 140  iv(nfgcal) = nfc
  rfac1 = 1.0_double
  v(dstsav) = v(dstnrm)
  if (v(fdif) > v(preduc)*v(tuner1)) go to 190
!
!   decrease was much less than predicted -- either change models
!   or accept step with decreased radius.
!
  if (iv(stage) >= iv(stglim)) go to 150
!         consider switching models
     iv(irc) = 2
     go to 160
!
!      accept step with decreased radius
!
 150  iv(irc) = 4
!
!   set v(radfac) to fletcher*s decrease factor
!
 160  iv(xirc) = iv(irc)
  emax = v(gtstep) + v(fdif)
  v(radfac) = half * rfac1

  if (emax < v(gtstep)) then
    v(radfac) = rfac1 * max (v(rdfcmn),half * v(gtstep)/emax)
  end if
!
!   do false convergence test
!
 170  if (v(reldx) <= v(xftol)) go to 180
     iv(irc) = iv(xirc)
     if (v(f) < v(f0)) go to 200
          go to 230
!
 180  iv(irc) = 12
  go to 240
!
!   handle good function decrease
!
 190  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 210
!
!      increasing radius looks worthwhile.  see if we just
!      recomputed step with a decreased radius or restored step
!      after recomputing it with a larger radius.
!
  if (iv(radinc) < 0) go to 210
  if (iv(restor) == 1) go to 210
!
!         we did not.  try a longer step unless this was a newton
!         step.
!
     v(radfac) = v(rdfcmx)
     gts = v(gtstep)
     if (v(fdif) < (half/v(radfac) - 1.0_double ) * gts) then
       v(radfac) = max (v(incfac), half*gts/(gts + v(fdif)))
     end if
     iv(irc) = 4
     if (v(stppar) .speq. 0.0_double ) go to 230
     if (v(dst0) >= 0.0_double .and. (v(dst0) < two*v(dstnrm) &
              .or. v(nreduc) < onep2*v(fdif)))  then
       go to 230
     end if
!
!  step was not a newton step.  recompute it with a larger radius.
!
          iv(irc) = 5
          iv(radinc) = iv(radinc) + 1
!
!   save values corresponding to good step
!
 200  v(flstgd) = v(f)
  iv(mlstgd) = iv(model)
  if (iv(restor) /= 1) iv(restor) = 2
  v(dstsav) = v(dstnrm)
  iv(nfgcal) = nfc
  v(plstgd) = v(preduc)
  v(gtslst) = v(gtstep)
  go to 230
!
!  Accept step with radius unchanged.
!
 210  v(radfac) = 1.0_double
  iv(irc) = 3
  go to 230
!
!   come here for a restart after convergence
!
 220  iv(irc) = iv(xirc)
  if (v(dstsav) >= 0.0_double ) go to 240
     iv(irc) = 12
     go to 240
!
!  Perform convergence tests
!
 230  iv(xirc) = iv(irc)
 240  if (iv(restor) == 1 .and. v(flstgd) < v(f0)) iv(restor) = 3
  if (abs(v(f)) < v(afctol)) iv(irc) = 10

  if (half * v(fdif) > v(preduc)) then
    return
  end if

  emax = v(rfctol) * abs(v(f0))
  emaxs = v(sctol) * abs(v(f0))
  if (v(dstnrm) > v(lmaxs) .and. v(preduc) <= emaxs) then
    iv(irc) = 11
  end if
  if (v(dst0) < 0.0_double ) go to 250
  i = 0

  if (( v(nreduc) > 0.0_double .and. v(nreduc) <= emax) .or. &
      (  (v(nreduc) .speq. 0.0_double) .and. (v(preduc) .speq. 0.0_double) )) then
    i = 2
  end if

  if (  (v(stppar) .speq. 0.0_double) .and. v(reldx) <= v(xctol) .and. goodx) then
    i = i + 1
  end if

  if (i > 0) iv(irc) = i + 6
!
!   consider recomputing step of length v(lmaxs) for singular
!   convergence test.
!
 250  if (iv(irc) > 5 .and. iv(irc) /= 12) then
     return
  end if

  if (v(dstnrm) > v(lmaxs)) go to 260
     if (v(preduc) >= emaxs) then
       return
     end if
          if (v(dst0) <= 0.0_double ) go to 270
               if (half * v(dst0) <= v(lmaxs)) then
                 return
               end if
                    go to 270
 260  if (half * v(dstnrm) <= v(lmaxs)) then
        return
      end if
  xmax = v(lmaxs) / v(dstnrm)
  if (xmax * (two - xmax) * v(preduc) >= emaxs) then
    return
  end if
 270  if (v(nreduc) < 0.0_double ) go to 290
!
!   recompute v(preduc) for use in singular convergence test
!
  v(gtslst) = v(gtstep)
  v(dstsav) = v(dstnrm)
  if (iv(irc) == 12) v(dstsav) = -v(dstsav)
  v(plstgd) = v(preduc)
  i = iv(restor)
  iv(restor) = 2
  if (i == 3) iv(restor) = 0
  iv(irc) = 6
  return
!
!  Perform singular convergence test with recomputed v(preduc)
!
 280  v(gtstep) = v(gtslst)
  v(dstnrm) = abs(v(dstsav))
  iv(irc) = iv(xirc)
  if (v(dstsav) <= 0.0_double ) iv(irc) = 12
  v(nreduc) = -v(preduc)
  v(preduc) = v(plstgd)
  iv(restor) = 3
 290  if (-v(nreduc) <= v(rfctol) * abs(v(f0))) iv(irc) = 11

end subroutine assst

subroutine dbdog(dig, lv, n, nwtstp, step, v)
!
!*******************************************************************************
!
!! DBDOG: compute double dogleg step
!
!
!  Discussion:
!
!        this subroutine computes a candidate step (for use in an uncon-
!     strained minimization code) by the double dogleg algorithm of
!     dennis and mei (ref. 1), which is a variation on powell*s dogleg
!     scheme (ref. 2, p. 95).
!
!  Parameters:
!
!    dig (input) diag(d)**-2 * g -- see algorithm notes.
!      g (input) the current gradient vector.
!     lv (input) length of v.
!      n (input) number of components in  dig, g, nwtstp,  and  step.
! nwtstp (input) negative newton step -- see algorithm notes.
!   step (output) the computed step.
!      v (i/o) values array, the following components of which are
!             used here...
! v(bias)   (input) bias for relaxed newton step, which is v(bias) of
!             the way from the full newton to the fully relaxed newton
!             step.  recommended value = 0.8 .
! v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
! v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
!             unless v(stppar) = 0 -- see algorithm notes.
! v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
! v(grdfac) (output) the coefficient of  dig  in the step returned --
!             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
! v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
!             algorithm notes.
! v(gtstep) (output) inner product between g and step.
! v(nreduc) (output) function reduction predicted for the full newton
!             step.
! v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
!             see v(grdfac) above.
! v(preduc) (output) function reduction predicted for the step returned.
! v(radius) (input) the trust region radius.  d times the step returned
!             has 2-norm v(radius) unless v(stppar) = 0.
! v(stppar) (output) code telling how step was computed... 0 means a
!             full newton step.  between 0 and 1 means v(stppar) of the
!             way from the newton to the relaxed newton step.  between
!             1 and 2 means a true double dogleg step, v(stppar) - 1 of
!             the way from the relaxed newton to the cauchy step.
!             greater than 2 means 1 / (v(stppar) - 1) times the cauchy
!             step.
!
!  notes
!
!   algorithm notes
!
!        let  g  and  h  be the current gradient and hessian approxima-
!     tion respectively and let d be the current scale vector.  this
!     routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
!     the step computed is the same one would get by replacing g and h
!     by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
!     computing step, and translating step back to the original
!     variables, i.e., premultiplying it by diag(d)**-1.
!
!  Reference:
!
!    J E Dennis, H H W Mei,
!    Two New Unconstrained Optimization Algorithms Which Use
!    Function and Gradient Values,
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    M J D Powell,
!    A Hybrid Method for Non-linear Equations,
!    in Numerical Methods for Non-linear Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
implicit none
  integer lv, n
  real(kind=double) dig(n), nwtstp(n), step(n), v(lv)
  integer i
  real(kind=double) cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm
  real(kind=double) nwtnrm, relax, rlambd, t, t1, t2
  real(kind=double) half, two
  integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep
  integer nreduc, nwtfac, preduc, radius, stppar
  parameter (half=0.5_double, two=2._double)
  parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45 )
  parameter ( gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7 )
  parameter ( radius=8, stppar=5)
!
  nwtnrm = v(dst0)
  rlambd = 1.0_double
  if (nwtnrm > 0.0_double ) rlambd = v(radius) / nwtnrm
  gnorm = v(dgnorm)
  ghinvg = two * v(nreduc)
  v(grdfac) = 0.0_double
  v(nwtfac) = 0.0_double
  if (rlambd < 1.0_double ) go to 30
!
!  The newton step is inside the trust region
!
     v(stppar) = 0.0_double
     v(dstnrm) = nwtnrm
     v(gtstep) = -ghinvg
     v(preduc) = v(nreduc)
     v(nwtfac) = -1.0_double
     do i = 1, n
            step(i) = -nwtstp(i)
     enddo
     return

 30   v(dstnrm) = v(radius)
  cfact = (gnorm / v(gthg))**2
!
!      cauchy step = -cfact * g.
!
  cnorm = gnorm * cfact
  relax = 1.0_double - v(bias) * ( 1.0_double - gnorm*cnorm/ghinvg)
  if (rlambd < relax) go to 50
!
!  Step is between relaxed newton and full newton steps
!
     v(stppar) = 1.0_double -  (rlambd - relax) / ( 1.0_double - relax)
     t = -rlambd
     v(gtstep) = t * ghinvg
     v(preduc) = rlambd * ( 1.0_double - half*rlambd) * ghinvg
     v(nwtfac) = t
     do i = 1, n
            step(i) = t * nwtstp(i)
     enddo
     return

 50   if (cnorm < v(radius)) go to 70
!
!  the cauchy step lies outside the trust region --
!  step = scaled cauchy step
!
     t = -v(radius) / gnorm
     v(grdfac) = t
     v(stppar) = 1.0_double  +  cnorm / v(radius)
     v(gtstep) = -v(radius) * gnorm
  v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
     do i = 1, n
            step(i) = t * dig(i)
     enddo
     return
!
!      compute dogleg step between cauchy and relaxed newton
!      femur = relaxed newton step minus cauchy step
!
 70   ctrnwt = cfact * relax * ghinvg / gnorm
!
!     ctrnwt = inner prod. of cauchy and relaxed newton steps,
!     scaled by gnorm**-1.
!
  t1 = ctrnwt - gnorm*cfact**2
!
!      t1 = inner prod. of femur and cauchy step, scaled by
!      gnorm**-1.
!
  t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
  t = relax * nwtnrm
  femnsq = (t/gnorm)*t - ctrnwt - t1
!
!  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
!
  t = t2 / (t1 + sqrt(t1**2 + femnsq*t2))
!
!  dogleg step  =  cauchy step  +  t * femur.
!
  t1 = (t - 1.0_double ) * cfact
  v(grdfac) = t1
  t2 = -t * relax
  v(nwtfac) = t2
  v(stppar) = two - t
  v(gtstep) = t1*gnorm**2 + t2*ghinvg
  v(preduc) = -t1*gnorm * ((t2 + 1.0_double )*gnorm) &
                  - t2 * ( 1.0_double + half*t2)*ghinvg &
                   - half * (v(gthg)*t1)**2

  do i = 1, n
    step(i) = t1*dig(i) + t2*nwtstp(i)
  end do


end subroutine dbdog
!-------------------------------------------------------


subroutine deflt(alg, iv, liv, lv, v)
!
!*******************************************************************************
!
!! DEFLT: supply default values to iv and v.
!
!   alg = 1 means regression constants.
!   alg = 2 means general unconstrained optimization constants.
!
implicit none
  integer liv, lv
!
  integer alg, iv(liv)
  real(kind=double) v(lv)
  integer miv, mv
  integer miniv(2), minv(2)
  integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits
  integer ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter
  integer nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm
  integer prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed
  integer vsave, x0prt

  parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71 )
  parameter (ierr=75, inith=25, inits=25, ipivot=76, ivneed=3 )
  parameter (lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18 )
  parameter (nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20 )
  parameter (parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57 )
  parameter (rmat=78, solprt=22, statpr=23, vneed=4, vsave=60 )
  parameter (x0prt=24)

  data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/
!
  if (alg < 1 .or. alg > 2) go to 40
  miv = miniv(alg)
  if (liv < miv) go to 20
  mv = minv(alg)
  if (lv < mv) go to 30
  call vdflt(alg, lv, v)
  iv(1) = 12
  iv(algsav) = alg
  iv(ivneed) = 0
  iv(lastiv) = miv
  iv(lastv) = mv
  iv(lmat) = mv + 1
  iv(mxfcal) = 200
  iv(mxiter) = 150
  iv(outlev) = 1
  iv(parprt) = 1
  iv(perm) = miv + 1
  iv(prunit) = 6
  iv(solprt) = 1
  iv(statpr) = 1
  iv(vneed) = 0
  iv(x0prt) = 1

  if (alg >= 2) go to 10
!
!   regression  values
!
  iv(covprt) = 3
  iv(covreq) = 1
  iv(dtype) = 1
  iv(hc) = 0
  iv(ierr) = 0
  iv(inits) = 0
  iv(ipivot) = 0
  iv(nvdflt) = 32
  iv(parsav) = 67
  iv(qrtyp) = 1
  iv(rdreq) = 3
  iv(rmat) = 0
  iv(vsave) = 58
  return
!
!   general optimization values
!
 10   iv(dtype) = 0
  iv(inith) = 1
  iv(nfcov) = 0
  iv(ngcov) = 0
  iv(nvdflt) = 25
  iv(parsav) = 47
  return

 20   iv(1) = 15
  return

 30   iv(1) = 16
  return

 40   iv(1) = 67


end subroutine deflt
!*************************************************8

real(kind=double) function dotprd(p, x, y)
!
!*******************************************************************************
!
!! DOTPRD returns the inner product of the p-vectors x and y.
!
implicit none
  integer p
!
  real(kind=double) x(p), y(p)
  integer i
  real(kind=double) t
  real(kind=double), save :: sqteta = 0.0_double
!
  dotprd = 0.0_double

  if (sqteta .speq. 0.0_double ) sqteta = rmdcon(2)
  do 20 i = 1, p
     t = max (abs(x(i)), abs(y(i)))
     if (t > 1.0_double ) go to 10
     if (t < sqteta) go to 20
     t = (x(i)/sqteta)*y(i)
     if (abs(t) < sqteta) go to 20
 10      dotprd = dotprd + x(i)*y(i)
 20   continue

end function dotprd

subroutine dupdu(d, hdiag, iv, liv, lv, n, v)
!
!*******************************************************************************
!
!! DUPDU: update scale vector d for humsl
!
!
implicit none
  integer liv, lv, n
  integer iv(liv)
  real(kind=double) d(n), hdiag(n), v(lv)
  integer dtoli, d0i, i
  real(kind=double) t, vdfac
  integer dfac, dtol, dtype, niter
!
  parameter (dfac=41, dtol=59, dtype=16, niter=31)
!
  i = iv(dtype)
  if (i == 1) go to 10
     if (iv(niter) > 0) then
       return
     end if

 10   dtoli = iv(dtol)
  d0i = dtoli + n
  vdfac = v(dfac)
  do i = 1, n
     t = max (sqrt(abs(hdiag(i))), vdfac*d(i))
     if (t < v(dtoli)) t = max (v(dtoli), v(d0i))
     d(i) = t
     dtoli = dtoli + 1
     d0i = d0i + 1
  end do

end subroutine dupdu
!************************************


subroutine gqtst(d, dig, dihdi, ka, l, p, step, v, w)
!
!*******************************************************************************
!
!! GQTST: compute goldfeld-quandt-trotter step by more-hebden technique
!   (nl2sol version 2.2), modified a la more and sorensen
!
!   purpose
!
!        given the (compactly stored) lower triangle of a scaled
!     hessian (approximation) and a nonzero scaled gradient vector,
!     this subroutine computes a goldfeld-quandt-trotter step of
!     approximate length v(radius) by the more-hebden technique.  in
!     other words, step is computed to (approximately) minimize
!     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
!     2-norm of d*step is at most (approximately) v(radius), where
!     g  is the gradient,  h  is the hessian, and  d  is a diagonal
!     scale matrix whose diagonal is stored in the parameter d.
!     (gqtst assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
!
!   parameter description
!
!     d (in)  = the scale vector, i.e. the diagonal of the scale
!              matrix  d  mentioned above under purpose.
!   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
!              step = 0  and  v(stppar) = 0  are returned.
! dihdi (in)  = lower triangle of the scaled hessian (approximation),
!              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
!              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
!    ka (i/o) = the number of hebden iterations (so far) taken to deter-
!              mine step.  ka < 0 on input means this is the first
!              attempt to determine step (for the present dig and dihdi)
!              -- ka is initialized to 0 in this case.  output with
!              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
!     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
!     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
!  step (i/o) = the step computed.
!     v (i/o) contains various constants and variables described below.
!     w (i/o) = workspace of length 4*p + 6.
!
!   entries in v
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (output) = 2-norm of d*step.
! v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
!             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
! v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
!             step returned, psi(step) will exceed its optimal value
!             by less than -v(epslon)*psi(step).  suggested value = 0.1.
! v(gtstep) (out) = inner product between g and step.
! v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
!             h only -- v(nreduc) is set to zero otherwise).
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
!             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
! v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
!             described below under algorithm notes.  if h + alpha*d**2
!             (see algorithm notes) is (nearly) singular, however,
!             then v(stppar) = -alpha.
!
!   usage notes
!
!     if it is desired to recompute step using a different value of
!     v(radius), then this routine may be restarted by calling it
!     with all parameters unchanged except v(radius).  (this explains
!     why step and w are listed as i/o).  on an initial call (one with
!     ka < 0), step and w need not be initialized and only compo-
!     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
!     v(rad0) of v must be initialized.
!
!   algorithm notes
!
!        the desired g-q-t step (ref. 2, 3, 4, 6) satisfies
!     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
!     h + alpha*d**2 is positive semidefinite.  alpha and step are
!     computed by a scheme analogous to the one described in ref. 5.
!     estimates of the smallest and largest eigenvalues of the hessian
!     are obtained from the gerschgorin circle theorem enhanced by a
!     simple form of the scaling described in ref. 7.  cases in which
!     h + alpha*d**2 is nearly (or exactly) singular are handled by
!     the technique discussed in ref. 2.  in these cases, a step of
!     (exact) length v(radius) is returned for which psi(step) exceeds
!     its optimal value by less than -v(epslon)*psi(step).  the test
!     suggested in ref. 6 for detecting the special case is performed
!     once two matrix factorizations have been done -- doing so sooner
!     seems to degrade the performance of optimization routines that
!     call this routine.
!
!  Reference:
!
!    J E Dennis, David Gay, and R E Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, pages 186-197, 1981.
!
!    S M Goldfeld, R E Quandt, H F Trotter,
!    Maximization by Quadratic Hill-climbing,
!    Econometrica,
!    Volume 34, pages 541-551, 1966.
!
!    M D Hebden,
!    An Algorithm for Minimization Using Exact Second Derivatives,
!    Report TP 515, Theoretical Physics Division,
!    AERE Harwell, Oxon., England, 1973.
!
!    J J More,
!    The Levenberg-Marquardt Algorithm, Implementation and Theory,
!    Springer Lecture Notes in Mathematics Number 630, pages 105-116,
!    edited by G A Watson,
!    Springer-Verlag, Berlin and New York, 1978.
!
!    J J More and D C Sorensen,
!    Computing a Trust Region Step,
!    Technical Report ANL-81-83,
!    Argonne National Lab, 1981.
!
!    Richard Varga,
!    Minimal Gerschgorin Sets,
!    Pacific Journal of Mathematics,
!    Volume 15, pages 719-729, 1965.
!
implicit none
  integer ka, p
  real(kind=double) d(p), dig(p), v(21), step(p)
  real(kind=double),dimension(p*(p+1)/2) :: dihdi
  real(kind=double),dimension(p*(p+1)/2) :: l
  real(kind=double),dimension(4*p+7)     :: w

  logical restrt
  integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc
  integer j, k, kalim, kamin, k1, lk0, phipin, q, q0, uk0, x
  real(kind=double) alphak, aki, akk, delta, dst, eps, gtsta, lk
  real(kind=double) oldphi
  real(kind=double) phi, phimax, phimin, psifac, rad, radsq
  real(kind=double) root, si, sk, sw, t, twopsi, t1, t2, uk, wi
  real(kind=double) big, dgxfac, epsfac, four, half, kappa, negone
  real(kind=double) one, p001, six, three, two, zero
  integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc
  integer phmnfc, phmxfc, preduc, radius, rad0

  parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4 )
  parameter ( nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8 )
  parameter ( rad0=9, stppar=5)

  parameter (epsfac=50.0_double, four=4.0_double, half=0.5_double )
  parameter ( kappa=2.0_double, negone=-1.0_double, one=1.0_double, p001=1.0e-3_double )
  parameter ( six=6.0_double, three=3.0_double, two=2.0_double, zero=0.0_double)

  save dgxfac

  data big/0._double/, dgxfac/0._double/
!
!  Store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
!
  dggdmx = p + 1
!
!  Store gerschgorin over- and underestimates of the largest
!  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
!  and w(emin) respectively.
!
  emax = dggdmx + 1
  emin = emax + 1

!
!  For use in recomputing step, the final values of lk, uk, dst,
!  and the inverse derivative of more*s phi at 0 (for pos. def.
!  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
!  respectively.
!
  lk0 = emin + 1
  phipin = lk0 + 1
  uk0 = phipin + 1
  dstsav = uk0 + 1

!
!  Store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
!
  diag0 = dstsav
  diag = diag0 + 1
!
!  Store -d*step in w(q),...,w(q0+p).
!
  q0 = diag0 + p
  q = q0 + 1

!
!  Allocate storage for scratch vector x
!
  x = q + p
  rad = v(radius)
  radsq = rad**2

!
!  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of d*step.
!
  phimax = v(phmxfc) * rad
  phimin = v(phmnfc) * rad
  psifac = two * v(epslon) / (three * (four * (v(phmnfc) + 1.0_double ) * &
                    (kappa + 1.0_double )  +  kappa  +  two) * rad**2)


!
!  OLDPHI is used to detect limits of numerical accuracy.  if
!  we recompute step and it does not change, then we accept it.
!
  oldphi = 0.0_double
  eps = v(epslon)
  irc = 0
  restrt = .false.
  kalim = ka + 50

!
!  Start or restart, depending on ka
!
  if (ka >= 0) go to 290
!
!   fresh start
!
  k = 0
  uk = negone
  ka = 0
  kalim = 50
  v(dgnorm) = v2norm(p, dig)
  v(nreduc) = 0.0_double
  v(dst0) = 0.0_double
  kamin = 3



  if (v(dgnorm) .speq. 0.0_double ) kamin = 0
!
!   store diag(dihdi) in w(diag0+1),...,w(diag0+p)
!
  j = 0
  do i = 1, p
    j = j + i
    k1 = diag0 + i
    w(k1) = dihdi(j)
  end do
!
!  determine w(dggdmx), the largest element of dihdi
!
  t1 = 0.0_double
  j = p * (p + 1) / 2
  do i = 1, j
     t = abs(dihdi(i))
     if (t1 < t) t1 = t
  end do
  w(dggdmx) = t1
!
!   try alpha = 0
!
 30   call lsqrt(1, p, l, dihdi, irc)
  if (irc == 0) go to 50
!
!         indef. h -- underestimate smallest eigenvalue, use this
!         estimate to initialize lower bound lk on alpha.
!
     j = irc*(irc+1)/2
     t = l(j)
     l(j) = 1.0_double
     do i = 1, irc
       w(i) = 0.0_double
     end do
     w(irc) = one
     call litvmu(irc, w, l, w)
     t1 = v2norm(irc, w)
     lk = -t / t1 / t1
     v(dst0) = -lk
     if (restrt) go to 210
     go to 70
!
!      positive definite h -- compute unmodified newton step.
!
 50   lk = 0.0_double
  t = lsvmin(p, l, w(q), w(q))
  if (t >= one) go to 60
     if (big <= 0.0_double ) big = rmdcon(6)
     if (v(dgnorm) >= t*t*big) go to 70
 60   call livmul(p, w(q), l, dig)
  gtsta = dotprd(p, w(q), w(q))
  v(nreduc) = half * gtsta
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  v(dst0) = dst
  phi = dst - rad
  if (phi <= phimax) go to 260
  if (restrt) go to 210
!
!  Prepare to compute gerschgorin estimates of largest (and
!  smallest) eigenvalues.
!
 70   k = 0
  do i = 1, p
     wi = 0.0_double
     im1 = i - 1
     do j = 1, im1
          k = k + 1
          t = abs(dihdi(k))
          wi = wi + t
          w(j) = w(j) + t
     end do
     w(i) = wi
     k = k + 1
  end do
!
!   (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)
!
  k = 1
  t1 = w(diag) - w(1)
  do i = 2, p
     j = diag0 + i
     t = w(j) - w(i)
     if ( t < t1 ) then
       t1 = t
       k = i
     end if
  end do

  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = k*(k-1)/2 + 1
  inc = 1
  t = 0.0_double
  do i = 1, p
     if (i >= k) inc = i
     if (i /= k) then
        aki = abs(dihdi(k1))
        si = w(i)
        j = diag0 + i
        t1 = half * (akk - w(j) + si - aki)
        t1 = t1 + sqrt(t1*t1 + sk*aki)
        if (t < t1) t = t1
     endif
     k1 = k1 + inc
  end do

  w(emin) = akk - t
  uk = v(dgnorm)/rad - w(emin)
  if (v(dgnorm) .speq. 0.0_double ) uk = uk + p001 + p001*uk
  if (uk <= 0.0_double) uk = p001
!
!   compute gerschgorin (over-)estimate of largest eigenvalue
!
  k = 1
  t1 = w(diag) + w(1)


  if (p > 1) then
      do i = 2, p
          j = diag0 + i
          t = w(j) + w(i)
          if (t <= t1) cycle
          t1 = t
          k = i
       enddo
  endif

  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = k*(k-1)/2 + 1
  inc = 1
  t = 0.0_double
  do i = 1, p
     ! Change increment if i is larger than k
     if(i >= k) inc = i
     ! update t1 and t (if appropriate)
     if(i /= k) then
           aki = abs(dihdi(k1))
           si = w(i)
           j = diag0 + i
           t1 = half * (w(j) + si - aki - akk)
           t1 = t1 + sqrt(t1*t1 + sk*aki)
           if (t < t1) t = t1
      endif
      k1 = k1 + inc
  enddo

  w(emax) = akk + t
  lk = max (lk, v(dgnorm)/rad - w(emax))

!
!      alphak = current value of alpha (see alg. notes above).  we
!      use more*s scheme for initializing it.
!
  alphak = abs(v(stppar)) * v(rad0)/rad
  if (irc /= 0) go to 210
!
!   compute l0 for positive definite h
!
  call livmul(p, w, l, w(q))
  t = v2norm(p, w)
  w(phipin) = dst / t / t
  lk = max (lk, phi*w(phipin))
!
!   safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)
!
 210  ka = ka + 1
  if (-v(dst0) >= alphak .or. alphak < lk .or. alphak >= uk) then
    alphak = uk * max (p001, sqrt(lk/uk))
  end if
  if (alphak <= 0.0_double) alphak = half * uk
  if (alphak <= 0.0_double) alphak = uk
  k = 0
  do i = 1, p
     k = k + i
     j = diag0 + i
     dihdi(k) = w(j) + alphak
  enddo
!
!   try computing cholesky decomposition
!
  call lsqrt(1, p, l, dihdi, irc)
  if (irc == 0) go to 240
!
!   (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
!   smallest eigenvalue for use in updating lk
!
  j = (irc*(irc+1))/2
  t = l(j)
  l(j) = one
  do  i = 1, irc
     w(i) = 0.0_double
  enddo
  w(irc) = one
  call litvmu(irc, w, l, w)
  t1 = v2norm(irc, w)
  lk = alphak - t/t1/t1
  v(dst0) = -lk
  go to 210
!
!   alphak makes (d**-1)*h*(d**-1) positive definite.
!   compute q = -d*step, check for convergence.
!
  240  call livmul(p, w(q), l, dig)
  gtsta = dotprd(p, w(q), w(q))
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  phi = dst - rad
  if (phi <= phimax .and. phi >= phimin) go to 270
  if (phi .speq. oldphi) go to 270
  oldphi = phi
  if (phi < 0.0_double) go to 330
!
!   unacceptable alphak -- update lk, uk, alphak
!
 250  if (ka >= kalim) go to 270
!
!      the following dmin1 is necessary because of restarts
!
  if (phi < 0.0_double) uk = min (uk, alphak)
!
!     kamin = 0 only iff the gradient vanishes
!
  if (kamin == 0) go to 210
  call livmul(p, w, l, w(q))
  t1 = v2norm(p, w)
  alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
  lk = max (lk, alphak)
  go to 210
!
!   acceptable step on first try
!
 260  alphak = 0.0_double
!
!   successful step in general.  compute step = -(d**-1)*q
!
 270  do i = 1, p
        j = q0 + i
        step(i) = -w(j)/d(i)
 enddo

  v(gtstep) = -gtsta
  v(preduc) = half * (abs(alphak)*dst*dst + gtsta)
  go to 410
!
!   restart with new radius
!
 290  if (v(dst0) <= 0.0_double .or. v(dst0) - rad > phimax) go to 310
!
!      prepare to return newton step
!
     restrt = .true.
     ka = ka + 1
     k = 0
     do i = 1, p
          k = k + i
          j = diag0 + i
          dihdi(k) = w(j)
     enddo
     uk = negone
     go to 30

 310  kamin = ka + 3
  if (v(dgnorm) .speq. 0.0_double) kamin = 0
  if (ka == 0) go to 50

  dst = w(dstsav)
  alphak = abs(v(stppar))
  phi = dst - rad
  t = v(dgnorm)/rad
  uk = t - w(emin)
  if (v(dgnorm) .speq. 0.0_double) uk = uk + p001 + p001*uk
  if (uk <= 0.0_double) uk = p001
  if (rad > v(rad0)) go to 320
!
!         smaller radius
!
     lk = 0.0_double
     if (alphak > 0.0_double) lk = w(lk0)
     lk = max (lk, t - w(emax))
     if (v(dst0) > 0.0_double) lk = max (lk, (v(dst0)-rad)*w(phipin))
     go to 250
!
!      bigger radius.
!
 320  if (alphak > 0.0_double) uk = min (uk, w(uk0))
  lk = max (zero, -v(dst0), t - w(emax))
  if (v(dst0) > 0.0_double) lk = max (lk, (v(dst0)-rad)*w(phipin))
  go to 250
!
!   decide whether to check for special case... in practice (from
!   the standpoint of the calling optimization code) it seems best
!   not to check until a few iterations have failed -- hence the
!   test on kamin below.
!
 330  delta = alphak + min (zero, v(dst0))
  twopsi = alphak*dst*dst + gtsta
  if (ka >= kamin) go to 340
!     if the test in ref. 2 is satisfied, fall through to handle
!     the special case (as soon as the more-sorensen test detects
!     it).
  if (delta >= psifac*twopsi) go to 370
!
!   check for the special case of  h + alpha*d**2  (nearly)
!   singular.  use one step of inverse power method with start
!   from lsvmin to obtain approximate eigenvector corresponding
!   to smallest eigenvalue of (d**-1)*h*(d**-1).  lsvmin returns
!   x and w with  l*w = x.
!
 340  t = lsvmin(p, l, w(x), w)
!
!      normalize w.
!
  do i = 1, p
      w(i) = t*w(i)
  enddo
!
!      complete current inv. power iter. -- replace w by (l**-t)*w.
!
  call litvmu(p, w, l, w)
  t2 = one/v2norm(p, w)
  do  i = 1, p
      w(i) = t2*w(i)
  enddo

  t = t2 * t
!
!   now w is the desired approximate (unit) eigenvector and
!   t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
!
  sw = dotprd(p, w(q), w)
  t1 = (rad + dst) * (rad - dst)
  root = sqrt(sw*sw + t1)
  if (sw < 0.0_double) root = -root
  si = t1 / (sw + root)
!
!   the actual test for the special case...
!
  if ((t2*si)**2 <= eps*(dst**2 + alphak*radsq)) go to 380
!
!   update upper bound on smallest eigenvalue (when not positive)
!   (as recommended by more and sorensen) and continue...
!
  if (v(dst0) <= 0.0_double) v(dst0) = min (v(dst0), t2**2 - alphak)
  lk = max (lk, -v(dst0))
!
!   check whether we can hope to detect the special case in
!   the available arithmetic.  accept step as it is if not.
!
!      if not yet available, obtain machine dependent value dgxfac.
 370  if (dgxfac .speq. 0.0_double) dgxfac = epsfac * rmdcon(3)

  if (delta > dgxfac*w(dggdmx)) go to 250
     go to 270
!
!   special case detected... negate alphak to indicate special case
!
 380  alphak = -alphak
  v(preduc) = half * twopsi
!
!   accept current step if adding si*w would lead to a
!   further relative reduction in psi of less than v(epslon)/3.
!
  t1 = 0.0_double
  t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
  if (t < eps*twopsi/six) go to 390
     v(preduc) = v(preduc) + t
     dst = rad
     t1 = -si
 390  do i = 1, p
      j = q0 + i
      w(j) = t1*w(i) - w(j)
      step(i) = w(j) / d(i)
   enddo

  v(gtstep) = dotprd(p, dig, w(q))
!
!   save values for use in a possible restart
!
 410  v(dstnrm) = dst
  v(stppar) = alphak
  w(lk0) = lk
  w(uk0) = uk
  v(rad0) = rad
  w(dstsav) = dst
!
!  Restore diagonal of dihdi
!
  j = 0
  do i = 1, p
     j = j + i
     k = diag0 + i
     dihdi(j) = w(k)
  end do

end subroutine gqtst
!************************************************


subroutine humit(d, fx, g, h, iv, lh, liv, lv, n, v, x)
!
!*******************************************************************************
!
!! HUMIT carries out humsl (unconstrained minimization) iterations, using
!   hessian matrix provided by the caller.
!
!  parameter usage
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! h.... lower triangle of the hessian, stored rowwise.
! iv... integer value array.
! lh... length of h = p*(p+1)/2.
! liv.. length of iv (at least 60).
! lv... length of v (at least 78 + n*(n+21)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... parameter vector.
!
!   discussion
!
!        parameters iv, n, v, and x are the same as the corresponding
!     ones to humsl (which see), except that v can be shorter (since
!     the part of v that humsl uses for storing g and h is not needed).
!     moreover, compared with humsl, iv(1) may have the two additional
!     output values 1 and 2, which are explained below, as is the use
!     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!     output value from humsl, is not referenced by humit or the
!     subroutines it calls.
!
! iv(1) = 1 means the caller should set sumitfx to f(x), the function value
!             at x, and call humit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             computed (e.g. if overflow would occur), which may happen
!             because of an oversized step.  in this case the caller
!             should set iv(toobig) = iv(2) to 1, which will cause
!             humit to ignore fx and try a smaller step.  the para-
!             meter nf that humsl passes to calcf (for possible use by
!             calcgh) is a copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient of f at
!             x, and h to the lower triangle of h(x), the hessian of f
!             at x, and call humit again, having changed none of the
!             other parameters except perhaps the scale vector d.
!                  the parameter nf that humsl passes to calcg is
!             iv(nfgcal) = iv(7).  if g(x) and h(x) cannot be evaluated,
!             then the caller may set iv(nfgcal) to 0, in which case
!             humit will return with iv(1) = 65.
!                  note -- humit overwrites h with the lower triangle
!             of  diag(d)**-1 * h(x) * diag(d)**-1.
!
implicit none
  integer lh, liv, lv, n
  integer iv(liv)
  real(kind=double):: fx
  real(kind=double), dimension(n):: d
  real(kind=double), dimension(n):: g
  real(kind=double), dimension(lh):: h
  real(kind=double), dimension(n):: x
  real(kind=double),dimension(lv):: v

  integer dg1, dummy, i, j, k, l, lstgst, nn1o2, step1
  integer temp1, w1, x01
  real(kind=double) t
  real(kind=double) one, onep2

  integer cnvcod, dg, dgnorm, dinit, dstnrm, dtinit, dtol
  integer dtype, d0init, f, f0, fdif, gtstep, incfac, irc, kagqt
  integer lmat, lmax0, lmaxs, mode, model, mxfcal, mxiter, nextv
  integer nfcall, nfgcal, ngcall, niter, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, stppar
  integer toobig, tuner4, tuner5, vneed, w, xirc, x0

  parameter (cnvcod=55, dg=37, dtol=59, dtype=16, irc=29, kagqt=33 )
  parameter ( lmat=42, mode=35, model=5, mxfcal=17, mxiter=18 )
  parameter ( nextv=47, nfcall=6, nfgcal=7, ngcall=30, niter=31 )
  parameter ( radinc=8, restor=9, step=40, stglim=11, stlstg=41 )
  parameter ( toobig=2, vneed=4, w=34, xirc=13, x0=43)
  parameter (dgnorm=1, dinit=38, dstnrm=2, dtinit=39, d0init=40 )
  parameter ( f=10, f0=13, fdif=11, gtstep=4, incfac=23, lmax0=35 )
  parameter ( lmaxs=36, preduc=7, radfac=16, radius=8, rad0=9 )
  parameter ( reldx=17, stppar=5, tuner4=29, tuner5=30)
!
  parameter (one=1._double, onep2=1.2_double)
!

  i = iv(1)
  if (i == 1) go to 30
  if (i == 2) go to 40
!
!   check validity of iv and v input values
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
   iv(vneed) = iv(vneed) + n*(n+21)/2 + 7
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2
  if (i > 12) then
    return
  end if
  nn1o2 = n * (n + 1) / 2

 if (lh >= nn1o2) then
      select case (i)
      case(1:6)
          goto 210
      case(7,9)
          goto 160
      case(8)
           goto 120
      case(10:11)
           goto 10
      case(12)
           goto 20
      end select
endif



iv(1) = 66
go to 350
!
!   storage allocation
!
 10   iv(dtol) = iv(lmat) + nn1o2
  iv(x0) = iv(dtol) + 2*n
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(dg) = iv(stlstg) + n
  iv(w) = iv(dg) + n
  iv(nextv) = iv(w) + 4*n + 7
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0_double
  v(stppar) = 0.0_double
  if (v(dinit) >= 0.0_double) call vscopy(n, d, v(dinit))
  k = iv(dtol)
  if (v(dtinit) > 0.0_double) call vscopy(n, v(k), v(dtinit))
  k = k + n
  if (v(d0init) > 0.0_double) call vscopy(n, v(k), v(d0init))
  iv(1) = 1
  return

 30   v(f) = fx
  if (iv(mode) >= 0) go to 210
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 350
!
!   make sure gradient could be computed
!
 40   if (iv(nfgcal) /= 0) go to 50
     iv(1) = 65
     go to 350
!
!   update the scale vector d
!
 50   dg1 = iv(dg)
  if (iv(dtype) <= 0) go to 70
  k = dg1
  j = 0
  do i = 1, n
     j = j + i
     v(k) = h(j)
     k = k + 1
  enddo
  call dupdu(d, v(dg1), iv, liv, lv, n, v)
!
!   compute scaled gradient and its norm
!
 70   dg1 = iv(dg)
  k = dg1
  do  i = 1, n
     v(k) = g(i) / d(i)
     k = k + 1
  enddo
  v(dgnorm) = v2norm(n, v(dg1))
!
!   compute scaled hessian
!
  k = 1
  do  i = 1, n
     t = one / d(i)
     do  j = 1, i
          h(k) = t * h(k) / d(j)
          k = k + 1
     enddo
  enddo

 !
  if (iv(cnvcod) /= 0) go to 340
  if (iv(mode) == 0) go to 300
!
!   allow first step to have scaled 2-norm at most v(lmax0)
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  main loop
!
!   print iteration summary, check iteration limit
!
 110 call itsum(d, g, iv, liv, lv, n, v, x)
 120  k = iv(niter)
  if (k < iv(mxiter)) go to 130
     iv(1) = 10
     go to 350

 130  iv(niter) = k + 1
!
!   initialize for start of next iteration
!
  dg1 = iv(dg)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1

!
!      copy x to x0
!
  call vcopy(n, v(x01), x)
!
!   update radius
!
  if (k == 0) go to 150
  step1 = iv(step)
  k = step1
  do i = 1, n
       v(k) = d(i) * v(k)
       k = k + 1
  enddo

  v(radius) = v(radfac) * v2norm(n, v(step1))
!
!   check stopx and function evaluation limit
!
 150  if (.not. stopx(dummy)) go to 170
     iv(1) = 11
     go to 180
!
!      come here when restarting after func. eval. limit or stopx.
!
 160  if (v(f) >= v(f0)) go to 170
     v(radfac) = one
     k = iv(niter)
     go to 130

 170  if (iv(nfcall) < iv(mxfcal)) go to 190

     iv(1) = 9
 180     if (v(f) >= v(f0)) go to 350
!
!         in case of stopx or function evaluation limit with
!         improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 290
!
!  Compute candidate step
!
 190  step1 = iv(step)
  dg1 = iv(dg)
  l = iv(lmat)
  w1 = iv(w)
  call gqtst(d, v(dg1), h, iv(kagqt), v(l), n, v(step1), v, v(w1))
  if (iv(irc) == 6) go to 210
!
!   check whether evaluating f(x0 + step) looks worthwhile
!
  if (v(dstnrm) <= 0.0_double) go to 210
  if (iv(irc) /= 5) go to 200
  if (v(radfac) <= one) go to 200
  if (v(preduc) <= onep2 * v(fdif)) go to 210
!
!   compute f(x0 + step)
!
 200  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 210  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 220
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))
 220  k = iv(irc)

select case(k) ! LP August 2005
case(1,5)
    goto 230
case(2:4)
    goto 260
case(6)
    goto 240
case(7:12)
    goto 250
case(13)
    goto 330
case(14)
    goto 300
end select

!      recompute step with new radius
!
 230     v(radius) = v(radfac) * v(dstnrm)
     go to 150
!
!   compute step of length v(lmaxs) for singular convergence test.
!
 240  v(radius) = v(lmaxs)
  go to 190
!
!  Convergence or false convergence
!
 250  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 340
     if (iv(xirc) == 14) go to 340
          iv(xirc) = 14
!
!  Process acceptable step.
!
 260  if (iv(irc) /= 3) go to 290
     temp1 = lstgst
!
!      prepare for gradient tests
!      set  temp1 = hessian * step + g(x0)
!                 = diag(d) * (h * step + g(x0))
!
!        use x0 vector as temporary.
!
     k = x01
     do i = 1, n
          v(k) = d(i) * v(step1)
          k = k + 1
          step1 = step1 + 1
    enddo
     call slvmul(n, v(temp1), h, v(x01))
     do i = 1, n
          v(temp1) = d(i) * v(temp1) + g(i)
          temp1 = temp1 + 1
     enddo
!
!   compute gradient and hessian
!
 290  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return

 300  iv(1) = 2
  if (iv(irc) /= 3) go to 110
!
!   set v(radfac) by gradient tests
!
  temp1 = iv(stlstg)
  step1 = iv(step)
!
!      set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
!
  k = temp1
  do i = 1, n
     v(k) = (v(k) - g(i)) / d(i)
     k = k + 1
  enddo
!
!  Do gradient tests
!
  if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) go to 320
       if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  go to 110
 320            v(radfac) = v(incfac)
            go to 110
!
!   misc. details
!
!   bad parameters to assess
!
 330  iv(1) = 64
  go to 350
!
!   print summary of final iteration and other requested items
!
 340  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 350 call itsum(d, g, iv, liv, lv, n, v, x)

end subroutine humit
!*****************************************************


subroutine humsl(n, d, x, calcf, calcgh, iv, liv, lv, v, uiparm, urparm, ufparm)
!
!*******************************************************************************
!
!! HUMSL minimizes a general unconstrained objective function using
!   (analytic) gradient and hessian provided by the caller.
!
!  discussion
!
!        this routine is like sumsl, except that the subroutine para-
!     meter calcg of sumsl (which computes the gradient of the objec-
!     tive function) is replaced by the subroutine parameter calcgh,
!     which computes both the gradient and (lower triangle of the)
!     hessian of the objective function.  the calling sequence is...
!             call calcgh(n, x, nf, g, h, uiparm, urparm, ufparm)
!     parameters n, x, nf, g, uiparm, urparm, and ufparm are the same
!     as for sumsl, while h is an array of length n*(n+1)/2 in which
!     calcgh must store the lower triangle of the hessian at x.  start-
!     ing at h(1), calcgh must store the hessian entries in the order
!     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
!        the value printed (by itsum) in the column labelled stppar
!     is the levenberg-marquardt used in computing the current step.
!     zero means a full newton step.  if the special case described in
!     ref. 1 is detected, then stppar is negated.  the value printed
!     in the column labelled npreldf is zero if the current hessian
!     is not positive definite.
!        it sometimes proves worthwhile to let d be determined from the
!     diagonal of the hessian matrix by setting iv(dtype) = 1 and
!     v(dinit) = 0.  the following iv and v components are relevant...
!
! iv(dtol)  iv(59) gives the starting subscript in v of the dtol
!             array used when d is updated.  (iv(dtol) can be
!             initialized by calling humsl with iv(1) = 13.)
! iv(dtype).... iv(16) tells how the scale vector d should be chosen.
!             iv(dtype) <= 0 means that d should not be updated, and
!             iv(dtype) >= 1 means that d should be updated as
!             described below with v(dfac).  default = 0.
! v(dfac)  v(41) and the dtol and d0 arrays (see v(dtinit) and
!             v(d0init)) are used in updating the scale vector d when
!             iv(dtype) > 0.  (d is initialized according to
!             v(dinit), described in sumsl.)  let
!                  d1(i) = max(sqrt(abs(h(i,i))), v(dfac)*d(i)),
!             where h(i,i) is the i-th diagonal element of the current
!             hessian.  if iv(dtype) = 1, then d(i) is set to d1(i)
!             unless d1(i) < dtol(i), in which case d(i) is set to
!                  max(d0(i), dtol(i)).
!             if iv(dtype) >= 2, then d is updated during the first
!             iteration as for iv(dtype) = 1 (after any initialization
!             due to v(dinit)) and is left unchanged thereafter.
!             default = 0.6.
! v(dtinit)... v(39), if positive, is the value to which all components
!             of the dtol array (see v(dfac)) are initialized.  if
!             v(dtinit) = 0, then it is assumed that the caller has
!             stored dtol in v starting at v(iv(dtol)).
!             default = 10**-6.
! v(d0init)... v(40), if positive, is the value to which all components
!             of the d0 vector (see v(dfac)) are initialized.  if
!             v(dfac) = 0, then it is assumed that the caller has
!             stored d0 in v starting at v(iv(dtol)+n).  default = 1.0.
!
!  Reference:
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, pages 186-197, 1981.
!
implicit none
  integer liv, lv, n
  integer iv(liv), uiparm(*)
  real(kind=double) d(n), x(n), v(lv), urparm(n)
!     dimension v(78 + n*(n+12))
  external ufparm
  integer g1, h1, iv1, lh, nf
  real(kind=double) f
  integer g, h, nextv, nfcall, nfgcal, toobig, vneed
!
  parameter (nextv=47, nfcall=6, nfgcal=7, g=28, h=56, toobig=2, vneed=4)
!


  lh = n * (n + 1) / 2
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+3)/2
  end if
  iv1 = iv(1)
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  h1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)
  h1 = iv(h)

 20 continue
  call humit(d, f, v(g1), v(h1), iv, lh, liv, lv, n, v, x)
  if ( iv(1)  < 2 ) then
     nf = iv(nfcall)
     call calcf(n, x, nf, f, uiparm, urparm, ufparm)
     if (nf <= 0) iv(toobig) = 1
     go to 20
  elseif (iv(1)  ==  2 ) then
     call calcgh(n, x(1:n), iv(nfgcal), v(g1:g1+n-1), v(h1:h1+n*(n+1)/2-1), uiparm, urparm, ufparm)
     go to 20
  elseif (iv(1) /= 14) then
        return
  end if

!
!   storage allocation
!
  iv(g) = iv(nextv)
  iv(h) = iv(g) + n
  iv(nextv) = iv(h) + n*(n+1)/2
  if (iv1 /= 13) go to 10

end subroutine humsl

subroutine itsum(d, g, iv, liv, lv, p, v, x)
!
!*******************************************************************************
!
!! ITSUM prints an iteration summary.
!
!
implicit none
  integer liv, lv, p
  integer iv(liv)
  real(kind=double) d(p), g(p), v(lv), x(p)
  integer alg, i, iv1, m, nf, ng, ol, pu
  character(len=4) model1(6), model2(6)
  real(kind=double) nreldf, oldf, preldf, reldf
  integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov
  integer ngcall, niter, nreduc, outlev, preduc, prntit, prunit
  integer reldx, solprt, statpr, stppar, sused, x0prt

  parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30 )
  parameter ( ngcov=53, niter=31, outlev=19, prntit=39, prunit=21 )
  parameter ( solprt=22, statpr=23, sused=64, x0prt=24)
  parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, reldx=17 )
  parameter ( stppar=5)

  data model1/'    ','    ','    ','    ','  g ','  s '/
  data model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
!


  pu = iv(prunit)
  if (pu == 0) then
    return
  end if
  iv1 = iv(1)
  if (iv1 > 62) iv1 = iv1 - 51
  ol = iv(outlev)
  alg = iv(algsav)
  if (iv1 < 2 .or. iv1 > 15) go to 370
  if (iv1 >= 12) go to 120
  if (iv1 == 2 .and. iv(niter) == 0) go to 390
  if (ol == 0) go to 120
  if (iv1 >= 10 .and. iv(prntit) == 0) go to 120
  if (iv1 > 2) go to 10
     iv(prntit) = iv(prntit) + 1
     if (iv(prntit) < iabs(ol)) then
       return
     end if
 10   nf = iv(nfcall) - iabs(iv(nfcov))
  iv(prntit) = 0
  reldf = 0.0_double
  preldf = 0.0_double
  oldf = max (abs(v(f0)), abs(v(f)))
  if (oldf <= 0.0_double) go to 20
     reldf = v(fdif) / oldf
     preldf = v(preduc) / oldf
 20   if (ol > 0) go to 60
!
!         print short summary line
!
     if (iv(needhd) == 1 .and. alg == 1) write(pu,30)
 30   format(/"   it   nf",6x,"f",7x,"reldf",3x,"preldf",3x,"reldx",2x,"model  stppar")
     if (iv(needhd) == 1 .and. alg == 2) write(pu,40)
 40   format(/"    it   nf",7x,"f",8x,"hreldf",4x,"preldf",4x,"reldx",3x,"stppar")
     iv(needhd) = 0
     if (alg == 2) go to 50
     m = iv(sused)
     write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
       model1(m), model2(m), v(stppar)
     go to 120

 50  write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx), v(stppar)
     go to 120
!
!      print long summary line
!
 60   if (iv(needhd) == 1 .and. alg == 1) write(pu,70)
 70   format(/"    it   nf",6x,"f",7x,"reldf",3x,"preldf",3x,"reldx", &
        2x,"model  stppar",2x,"d*step",2x,"npreldf")
  if (iv(needhd) == 1 .and. alg == 2) write(pu,80)
 80   format(/"    it   nf",7x,"f",8x,"reldf",4x,"preldf",4x,"reldx", &
        3x,"stppar",3x,"d*step",3x,"npreldf")
  iv(needhd) = 0
  nreldf = 0.0_double
  if (oldf > 0.0_double) nreldf = v(nreduc) / oldf
  if (alg == 2) go to 90
  m = iv(sused)
  write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
              model1(m), model2(m), v(stppar), v(dstnrm), nreldf
  go to 120

 90     write(pu,110) iv(niter), nf, v(f), reldf, preldf, &
                      v(reldx), v(stppar), v(dstnrm), nreldf
 100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
 110  format(i6,i5,d11.3,2(d10.2),3(d9.1),d10.2)

 120  if (iv(statpr) < 0) go to 430
  select case (iv1)
  case(1:2)
       goto 999
  case(3)
       goto 130
  case(4)
       goto 150
  case(5)
       goto 170
  case(6)
       goto 190
  case(7)
       goto 210
  case(8)
       goto 230
  case(9)
       goto 250
  case(10)
       goto 270
  case(11)
       goto 290
  case(12)
       goto 310
  case(13)
       goto 330
  case(14)
       goto 350
  case(15)
       goto 520
  end select

 130  write(pu,140)
 140  format(/' x-convergence')
  go to 430

 150  write(pu,160)
 160  format(/'relative function convergence')
  go to 430

 170  write(pu,180)
 180  format(/'x- and relative function convergence')
  go to 430

 190  write(pu,200)
 200  format(/'Absolute function convergence.')
  go to 430

 210  write(pu,220)
 220  format(/'Singular convergence.')
  go to 430

 230  write(pu,240)
 240  format(/'False convergence.')
  go to 430

 250  write(pu,260)
 260  format(/'Function evaluation limit.')
  go to 430

 270  write(pu,280)
 280  format(/'Iteration limit.')
  go to 430

 290  write(pu,300)
 300  format(/'STOPX')
  go to 430

 310  write(pu,320)
 320  format(/'Initial f(x) cannot be computed.')
  go to 390

 330  write(pu,340)
 340  format(/'Bad parameters to assess.')
  go to 999

 350  write(pu,360)
 360  format(/'Gradient could not be computed.')
  if (iv(niter) > 0) go to 480
  go to 390

 370  write(pu,380) iv(1)
 380  format(/'iv(1) =',i5)
  go to 999
!
!   initial call on itsum
!
 390  if (iv(x0prt) /= 0) write(pu,400) (i, x(i), d(i), i = 1, p)
 400  format(/"     i     initial x(i)",8x,"d(i)"//(1x,i5,d17.6,d14.3))
!     the following are to avoid undefined variables when the
!     function evaluation limit is 1...
!
  v(dstnrm) = 0.0_double
  v(fdif) = 0.0_double
  v(nreduc) = 0.0_double
  v(preduc) = 0.0_double
  v(reldx)  = 0.0_double
  if (iv1 >= 12) go to 999
  iv(needhd) = 0
  iv(prntit) = 0
  if (ol == 0) go to 999
  if (ol < 0 .and. alg == 1) write(pu,30)
  if (ol < 0 .and. alg == 2) write(pu,40)
  if (ol > 0 .and. alg == 1) write(pu,70)
  if (ol > 0 .and. alg == 2) write(pu,80)
  if (alg == 1) write(pu,410) v(f)
  if (alg == 2) write(pu,420) v(f)
 410  format(/"     0    1",d10.3)
!365  format(/11h     0    1,e11.3)
 420  format(/"     0    1",d11.3)
  go to 999
!
!   print various information requested on solution
!
 430  iv(needhd) = 1
  if (iv(statpr) == 0) go to 480
     oldf = max (abs(v(f0)), abs(v(f)))
     preldf = 0.0_double
     nreldf = 0.0_double
     if (oldf <= 0.0_double) go to 440
          preldf = v(preduc) / oldf
          nreldf = v(nreduc) / oldf
 440     nf = iv(nfcall) - iv(nfcov)
     ng = iv(ngcall) - iv(ngcov)
     write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
 450  format(/" function",d17.6,"   reldx",d17.3/" func. evals", &
    i8,9x,"grad. evals",i8/" preldf",d16.3,6x,"npreldf",d15.3)

     if (iv(nfcov) > 0) write(pu,460) iv(nfcov)
 460     format(/1x,i4," extra func. evals for covariance and diagnostics.")
     if (iv(ngcov) > 0) write(pu,470) iv(ngcov)
 470     format(1x,i4," extra grad. evals for covariance and diagnostics.")

 480  if (iv(solprt) == 0) go to 999
     iv(needhd) = 1
     write(pu,490)
 490  format(/"     i      final x(i)",8x,"d(i)",10x,"g(i)"/)
     do 500 i = 1, p
          write(pu,510) i, x(i), d(i), g(i)
 500          continue
 510     format(1x,i5,d16.6,2d14.3)
  go to 999

 520  write(pu,530)
 530  format(/'Inconsistent dimensions.')
 999  continue

  return
end subroutine itsum
!********************************************************************

subroutine litvmu(n, x, l, y)
!
!*******************************************************************************
!
!! LITVMU solves  (l**t)*x = y,  where  l  is an  n x n  lower triangular
!   matrix stored compactly by rows.  x and y may occupy the same
!   storage.
!
implicit none
  integer n
  real(kind=double) x(n), l(n*(n+1)/2), y(n)
  integer i, ii, ij, i0, j
  real(kind=double) xi
!
  x(1:n) = y(1:n)

  i0 = n*(n+1)/2

  do ii = 1, n
    i = n+1 - ii
    xi = x(i)/l(i0)
    x(i) = xi
    if ( i <= 1 ) then
      exit
    end if
    i0 = i0 - i
    if ( xi .spne. 0.0_double ) then
      do j = 1, i-1
        ij = i0 + j
        x(j) = x(j) - xi*l(ij)
      end do
    end if
  end do

  return
end subroutine litvmu
!*****************************************
subroutine livmul(n, x, l, y)
!
!*******************************************************************************
!
!! LIVMUL solves l*x = y, where  l  is an  n x n  lower triangular
!   matrix stored compactly by rows.  x and y may occupy the same
!   storage.
!
implicit none
  integer n
  real(kind=double) x(n), l(n*(n+1)/2), y(n)

  integer i, j, k
  real(kind=double) t
!

  do k = 1, n
    if (y(k) .spne. 0.0_double ) go to 20
    x(k) = 0.0_double
  end do
  return
 20   j = k*(k+1)/2
  x(k) = y(k) / l(j)
  if (k >= n) then
    return
  end if
  k = k + 1
  do i = k, n
     t = dotprd(i-1, l(j+1), x)
     j = j + i
     x(i) = (y(i) - t)/l(j)
  end do

  return
end subroutine livmul
!***************************************

subroutine lsqrt(n1, n, l, a, irc)
!
!*******************************************************************************
!
!! LSQRT computes rows n1 through n of the cholesky factor  l  of
!   a = l*(l**t),  where  l  and the lower triangle of  a  are both
!   stored compactly by rows (and may occupy the same storage).
!   irc = 0 means all went well.  irc = j means the leading
!   principal  j x j  submatrix of  a  is not positive definite --
!   and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
!
!   parameters
!
implicit none
  integer n1, n, irc, im1
  real(kind=double), dimension(n*(n+1)/2)   :: l
  real(kind=double), dimension(n*(n+1)/2) :: a
  integer i, ij, ik, i0, j, jk, j0, k
  real(kind=double) t, td

  i0 = n1 * (n1 - 1) / 2
  do 50 i = n1, n
     td = 0.0_double
     if (i == 1) go to 40
     j0 = 0
     im1 = i - 1
     do 30 j = 1, im1
          t = 0.0_double
          do k = 1, j-1
            ik = i0 + k
            jk = j0 + k
            t = t + l(ik)*l(jk)
          end do
          ij = i0 + j
          j0 = j0 + j
          t = (a(ij) - t) / l(j0)
          l(ij) = t
          td = td + t*t
 30           continue
 40      i0 = i0 + i
     t = a(i0) - td
     if (t <= 0.0_double) go to 60
     l(i0) = sqrt(t)
 50      continue

  irc = 0
  return

 60   l(i0) = t
  irc = i

end subroutine lsqrt
!******************************************8

real (kind=double) function lsvmin(p, l, x, y)
!
!*******************************************************************************
!
!! LSVMIN estimates smallest sing. value of packed lower triang. matrix l
!
!   parameter declarations
!
!   purpose
!
!     this function returns a good over-estimate of the smallest
!     singular value of the packed lower triangular matrix l.
!
!   parameter description
!
!  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
!  l (in)  = array holding the elements of  l  in row order, i.e.
!             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
!  x (out) if lsvmin returns a positive value, then x is a normalized
!             approximate left singular vector corresponding to the
!             smallest singular value.  this approximation may be very
!             crude.  if lsvmin returns zero, then some components of x
!             are zero and the rest retain their input values.
!  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
!             unnormalized approximate right singular vector correspond-
!             ing to the smallest singular value.  this approximation
!             may be crude.  if lsvmin returns zero, then y retains its
!             input value.  the caller may pass the same vector for x
!             and y (nonstandard fortran usage), in which case y over-
!             writes x (for nonzero lsvmin returns).
!
!   algorithm notes
!
!     the algorithm is based on (1), with the additional provision that
!     lsvmin = 0 is returned if the smallest diagonal element of l
!     (in magnitude) is not more than the unit roundoff times the
!     largest.  the algorithm uses a random number generator proposed
!     in (4), which passes the spectral test with flying colors -- see
!     (2) and (3).
!
!  Feference:
!
!    Alan Cline, Cleve Moler, G Stewart, and J H Wilkinson,
!    An estimate for the Condition Number of a Matrix,
!    Report TM-310, 1977,
!    Applied Mathematics Division,
!    Argonne National Laboratory.
!
!     (2) hoaglin, d.c. (1976), theoretical properties of congruential
!         random-number generators --  an empirical view,
!         memorandum ns-340, dept. of statistics, harvard univ.
!
!     (3) knuth, d.e. (1969), the art of computer programming, vol. 2
!         (seminumerical algorithms), addison-wesley, reading, mass.
!
!     (4) smith, c.s. (1971), multiplicative pseudo-random number
!         generators with prime modulus, j. assoc. comput. mach. 18,
!         pp. 586-593.
!
implicit none
  integer p

  real(kind=double) l(p*(p+1)/2), x(p), y(p)
!
  integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pm1
  real(kind=double) b, sminus, splus, t, xminus, xplus
  real(kind=double) half, one, r9973
!
  parameter (half=0.5_double, one=1._double, r9973=9973._double)
!
  ix = 2
  pm1 = p - 1
!
!   first check whether to return lsvmin = 0 and initialize x
!
  ii = 0
  j0 = p*pm1/2
  jj = j0 + p
  if (l(jj) .speq. 0.0_double) go to 110
  ix = mod(3432*ix, 9973)
  b = half*(one + float(ix)/r9973)
  xplus = b / l(jj)
  x(p) = xplus
  if (p <= 1) go to 60
  do i = 1, pm1
     ii = ii + i
     if (l(ii) .speq. 0.0_double) go to 110
     ji = j0 + i
     x(i) = xplus * l(ji)
  end do
!
!   solve (l**t)*x = b, where the components of b have randomly
!   chosen magnitudes in (.5,1) with signs chosen to make x large.
!
!     do j = p-1 to 1 by -1...
  do 50 jjj = 1, pm1
     j = p - jjj
!        determine x(j) in this iteration. note for i = 1,2,...,j
!        that x(i) holds the current partial sum for row i.
     ix = mod(3432*ix, 9973)
     b = half*(one + float(ix)/r9973)
     xplus = (b - x(j))
     xminus = (-b - x(j))
     splus = abs(xplus)
     sminus = abs(xminus)
     jm1 = j - 1
     j0 = j*jm1/2
     jj = j0 + j
     xplus = xplus/l(jj)
     xminus = xminus/l(jj)
     if (jm1 == 0) go to 30
     do 20 i = 1, jm1
          ji = j0 + i
          splus = splus + abs(x(i) + l(ji)*xplus)
          sminus = sminus + abs(x(i) + l(ji)*xminus)
 20           continue
 30      if (sminus > splus) xplus = xminus
     x(j) = xplus
!        update partial sums
     if (jm1 > 0) call vaxpy(jm1, x, xplus, l(j0+1), x)
 50      continue
!
!   normalize x
!
 60   t = one/v2norm(p, x)
  do i = 1, p
       x(i) = t*x(i)
  enddo
!
!   solve l*y = x and return lsvmin = 1/twonorm(y)
!
 do j = 1, p
     jm1 = j - 1
     j0 = j*jm1/2
     jj = j0 + j
     t = 0.0_double
     if (jm1 > 0) t = dotprd(jm1, l(j0+1), y)
     y(j) = (x(j) - t) / l(jj)
 enddo

  lsvmin = one/v2norm(p, y)
  return

 110  lsvmin = 0.0_double

end function lsvmin
!*****************************************************


subroutine ltvmul(n, x, l, y)
!
!*******************************************************************************
!
!! LTVMUL computes  x = (l**t)*y, where  l  is an  n x n  lower
!   triangular matrix stored compactly by rows.  x and y may
!   occupy the same storage.
!
implicit none
  integer n
  real(kind=double) x(n), l(n*(n+1)/2), y(n)

  integer i, ij, i0, j
  real(kind=double) yi
!
  i0 = 0
  do 20 i = 1, n
     yi = y(i)
     x(i) = 0.0_double
     do 10 j = 1, i
          ij = i0 + j
          x(j) = x(j) + yi*l(ij)
 10           continue
     i0 = i0 + i
 20      continue

  return
end subroutine ltvmul
!**************************************************

subroutine lupdat(beta, gamma, l, lambda, lplus, n, w, z)
!
!*******************************************************************************
!
!! LUPDAT computes lplus = secant update of l
!
!
!  parameter usage
!
!   beta = scratch vector.
!  gamma = scratch vector.
!      l (input) lower triangular matrix, stored rowwise.
! lambda = scratch vector.
!  lplus (output) lower triangular matrix, stored rowwise, which may
!             occupy the same storage as  l.
!      n (input) length of vector parameters and order of matrices.
!      w (input, destroyed on output) right singular vector of rank 1
!             correction to  l.
!      z (input, destroyed on output) left singular vector of rank 1
!             correction to  l.
!
!  notes
!
!   application and usage restrictions
!
!        this routine updates the cholesky factor  l  of a symmetric
!     positive definite matrix to which a secant update is being
!     applied -- it computes a cholesky factor  lplus  of
!     l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
!     and  z  have been chosen so that the updated matrix is strictly
!     positive definite.
!
!   algorithm notes
!
!        this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
!     to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
!     is an orthogonal matrix that makes the result lower triangular.
!        lplus may have some negative diagonal elements.
!
!  Reference:
!
!    D Goldfarb,
!    Factorized Variable Metric Methods for Unconstrained Optimization,
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
implicit none
  integer n
  real(kind=double) beta(n), gamma(n), l(n*(n+1)/2), lambda(n), lplus(n*(n+1)/2), w(n), z(n)

  integer i, ij, j, jj, jp1, k, nm1, np1
  real(kind=double) a, b, bj, eta, gj, lj, lij, ljj, nu, s, theta, wj, zj
  real(kind=double), parameter :: one = 1._double

  nu = one
  eta = 0.0_double
  if (n <= 1) go to 30
  nm1 = n - 1
!
!   temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
!   lambda(j).
!
  s = 0.0_double
  do 10 i = 1, nm1
     j = n - i
     s = s + w(j+1)**2
     lambda(j) = s
 10      continue
!
!   compute lambda, gamma, and beta by goldfarb*s recurrence 3.
!
  do 20 j = 1, nm1
     wj = w(j)
     a = nu*z(j) - eta*wj
     theta = one + a*wj
     s = a*lambda(j)
     lj = sqrt(theta**2 + a*s)
     if (theta > 0.0_double) lj = -lj
     lambda(j) = lj
     b = theta*wj + s
     gamma(j) = b * nu / lj
     beta(j) = (a - b*eta) / lj
     nu = -nu / lj
     eta = -(eta + (a**2)/(theta - lj)) / lj
 20      continue
 30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
!
!   update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
!
  np1 = n + 1
  jj = n * (n + 1) / 2
  do 60 k = 1, n
     j = np1 - k
     lj = lambda(j)
     ljj = l(jj)
     lplus(jj) = lj * ljj
     wj = w(j)
     w(j) = ljj * wj
     zj = z(j)
     z(j) = ljj * zj
     if (k == 1) go to 50
     bj = beta(j)
     gj = gamma(j)
     ij = jj + j
     jp1 = j + 1
     do 40 i = jp1, n
          lij = l(ij)
          lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
          w(i) = w(i) + lij*wj
          z(i) = z(i) + lij*zj
          ij = ij + i
 40           continue
 50      jj = jj - j
 60      continue

end subroutine lupdat
!**********************************************


subroutine lvmul(n, x, l, y)
!
!*******************************************************************************
!
!! LVMUL computes x = l*y, where  l  is an  n x n  lower triangular
!   matrix stored compactly by rows.  x and y may occupy the same
!   storage.
!
implicit none
  integer n
  real(kind=double) x(n), l(n*(n+1)/2), y(n)

  integer i, ii, ij, i0, j, np1
  real(kind=double) t

  np1 = n + 1
  i0 = n*(n+1)/2
  do 20 ii = 1, n
     i = np1 - ii
     i0 = i0 - i
     t = 0.0_double
     do 10 j = 1, i
          ij = i0 + j
          t = t + l(ij)*y(j)
 10           continue
     x(i) = t
 20      continue

end subroutine lvmul
!**************************************************************

subroutine parck(alg, d, iv, liv, lv, n, v)
!
!*******************************************************************************
!
!! PARCK checks parameters, prints changed values.
!
!   alg = 1 for regression, alg = 2 for general unconstrained opt.
!
implicit none
  integer alg, liv, lv, n
  integer iv(liv)
  real(kind=double) d(n), v(lv)
  integer max0
  integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
  integer ijmp, jlim(2), miniv(2), ndflt(2)
  character(len=1) varnm(2), sh(2)
  character(len=4) cngd(3), dflt(3), vn(2,34), which(3)
  real(kind=double) big, machep, tiny, vk, vm(34), vx(34)
  integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed
  integer lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn
  integer parprt, parsav, perm, prunit, vneed
!
  parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19 )
  parameter ( inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42 )
  parameter ( nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20 )
  parameter ( parsav=49, perm=58, prunit=21, vneed=4)
  save big, machep, tiny

  data big/0._double/, machep/-1._double/, tiny/1._double/

     data vn(1,1),vn(2,1)/'epsl','on..'/
     data vn(1,2),vn(2,2)/'phmn','fc..'/
     data vn(1,3),vn(2,3)/'phmx','fc..'/
     data vn(1,4),vn(2,4)/'decf','ac..'/
     data vn(1,5),vn(2,5)/'incf','ac..'/
     data vn(1,6),vn(2,6)/'rdfc','mn..'/
     data vn(1,7),vn(2,7)/'rdfc','mx..'/
     data vn(1,8),vn(2,8)/'tune','r1..'/
     data vn(1,9),vn(2,9)/'tune','r2..'/
     data vn(1,10),vn(2,10)/'tune','r3..'/
     data vn(1,11),vn(2,11)/'tune','r4..'/
     data vn(1,12),vn(2,12)/'tune','r5..'/
     data vn(1,13),vn(2,13)/'afct','ol..'/
     data vn(1,14),vn(2,14)/'rfct','ol..'/
     data vn(1,15),vn(2,15)/'xcto','l...'/
     data vn(1,16),vn(2,16)/'xfto','l...'/
     data vn(1,17),vn(2,17)/'lmax','0...'/
     data vn(1,18),vn(2,18)/'lmax','s...'/
     data vn(1,19),vn(2,19)/'scto','l...'/
     data vn(1,20),vn(2,20)/'dini','t...'/
     data vn(1,21),vn(2,21)/'dtin','it..'/
     data vn(1,22),vn(2,22)/'d0in','it..'/
     data vn(1,23),vn(2,23)/'dfac','....'/
     data vn(1,24),vn(2,24)/'dltf','dc..'/
     data vn(1,25),vn(2,25)/'dltf','dj..'/
     data vn(1,26),vn(2,26)/'delt','a0..'/
     data vn(1,27),vn(2,27)/'fuzz','....'/
     data vn(1,28),vn(2,28)/'rlim','it..'/
     data vn(1,29),vn(2,29)/'cosm','in..'/
     data vn(1,30),vn(2,30)/'hube','rc..'/
     data vn(1,31),vn(2,31)/'rspt','ol..'/
     data vn(1,32),vn(2,32)/'sigm','in..'/
     data vn(1,33),vn(2,33)/'eta0','....'/
     data vn(1,34),vn(2,34)/'bias','....'/

  data vm(1)/1.0e-3_double/, vm(2)/-0.99_double/, vm(3)/1.0e-3_double/, vm(4)/1.0e-2_double/
  data vm(5)/1.2_double/, vm(6)/1.e-2_double/, vm(7)/1.2_double/, vm(8)/0._double/
  data vm(9)/0._double/, vm(10)/1.e-3_double/, vm(11)/-1._double/, vm(13)/0._double/
  data vm(15)/0._double/, vm(16)/0._double/, vm(19)/0._double/, vm(20)/-10._double/
  data vm(21)/0._double/, vm(22)/0._double/, vm(23)/0._double/, vm(27)/1.01_double/
  data vm(28)/1.e+10_double/, vm(30)/0._double/, vm(31)/0._double/, vm(32)/0._double/
  data vm(34)/0._double/

  data vx(1)/0.9_double/, vx(2)/-1.e-3_double/, vx(3)/1.e+1_double/, vx(4)/0.8_double/
  data vx(5)/1.e+2_double/, vx(6)/0.8_double/, vx(7)/1.e+2_double/, vx(8)/0.5_double/
  data vx(9)/0.5_double/, vx(10)/1._double/, vx(11)/1._double/, vx(14)/0.1_double/
  data vx(15)/1._double/, vx(16)/1._double/, vx(19)/1._double/, vx(23)/1._double/
  data vx(24)/1._double/, vx(25)/1._double/, vx(26)/1._double/, vx(27)/1.e+10_double/
  data vx(29)/1._double/, vx(31)/1._double/, vx(32)/1._double/, vx(33)/1._double/
  data vx(34)/1._double/

  data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
  data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
  data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
  data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
  data miniv(1)/80/, miniv(2)/59/
!
  pu = 0
  if (prunit <= liv) pu = iv(prunit)
  if (alg < 1 .or. alg > 2) go to 340
  if (iv(1) == 0) call deflt(alg, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 /= 13 .and. iv1 /= 12) go to 10
  miv1 = miniv(alg)
  if (perm <= liv) miv1 = max0(miv1, iv(perm) - 1)
  if (ivneed <= liv) miv2 = miv1 + max0(iv(ivneed), 0)
  if (lastiv <= liv) iv(lastiv) = miv2
  if (liv < miv1) go to 300
  iv(ivneed) = 0
  iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
  iv(vneed) = 0
  if (liv < miv2) go to 300
  if (lv < iv(lastv)) go to 320
 10   if (alg == iv(algsav)) go to 30
  if (pu /= 0) write(pu,20) alg, iv(algsav)
 20 format(/" the first parameter to deflt should be",i3," rather than",i3)
     iv(1) = 82
     return
 30   if (iv1 < 12 .or. iv1 > 14) go to 60
     if (n >= 1) go to 50
          iv(1) = 81
          if (pu == 0) then
            return
          end if
          write(pu,40) varnm(alg), n
 40           format(/" /// bad",a1," =",i5)
          return
 50      if (iv1 /= 14) iv(nextiv) = iv(perm)
     if (iv1 /= 14) iv(nextv) = iv(lmat)
     if (iv1 == 13) then
       return
     end if
     k = iv(parsav) - epslon
     call vdflt(alg, lv-k, v(k+1))
     iv(dtype0) = 2 - alg
     iv(oldn) = n
     which(1) = dflt(1)
     which(2) = dflt(2)
     which(3) = dflt(3)
     go to 110
 60   if (n == iv(oldn)) go to 80
     iv(1) = 17
     if (pu == 0) then
       return
     end if
     write(pu,70) varnm(alg), iv(oldn), n
 70      format(/" /// ",1a1," changed from",i5," to" ,i5)
     return

 80   if (iv1 <= 11 .and. iv1 >= 1) go to 100
     iv(1) = 80
     if (pu /= 0) write(pu,90) iv1
 90      format(/" ///  iv(1) =",i5," should be between 0 and 14.")
     return

 100  which(1) = cngd(1)
  which(2) = cngd(2)
  which(3) = cngd(3)

 110  if (iv1 == 14) iv1 = 12
  if (big > tiny) go to 120
     tiny = rmdcon(1)
     machep = rmdcon(3)
     big = rmdcon(6)
     vm(12) = machep
     vx(12) = big
     vx(13) = big
     vm(14) = machep
     vm(17) = tiny
     vx(17) = big
     vm(18) = tiny
     vx(18) = big
     vx(20) = big
     vx(21) = big
     vx(22) = big
     vm(24) = machep
     vm(25) = machep
     vm(26) = machep
     vx(28) = rmdcon(5)
     vm(29) = machep
     vx(30) = big
     vm(33) = machep
 120  m = 0
  i = 1
  j = jlim(alg)
  k = epslon
  ndfalt = ndflt(alg)
  do 150 l = 1, ndfalt
    vk = v(k)
    if (vk >= vm(i) .and. vk <= vx(i)) go to 140
      m = k
      if (pu /= 0) write(pu,130) vn(1,i), vn(2,i), k, vk,vm(i), vx(i)
 130  format(/" ///  ",2a4,".. v(",i2,") =",d11.3," should", &
      " be between",d11.3," and",d11.3)
 140  k = k + 1
     i = i + 1
     if (i == j) i = ijmp
 150     continue

  if (iv(nvdflt) == ndfalt) go to 170
     iv(1) = 51
     if (pu == 0) then
       return
     end if
     write(pu,160) iv(nvdflt), ndfalt
 160     format(" iv(nvdflt) =",i5," rather than ",i5)
     return
 170  if ((iv(dtype) > 0 .or. v(dinit) > 0.0_double) .and. iv1 == 12) then
             go to 200
  end if
  do 190 i = 1, n
     if (d(i) > 0.0_double) go to 190
          m = 18
          if (pu /= 0) write(pu,180) i, d(i)
 180     format(/" ///  d(",i3,") =",d11.3," should be positive")
 190     continue
 200  if (m == 0) go to 210
     iv(1) = m
     return

 210  if (pu == 0 .or. iv(parprt) == 0) then
        return
      end if
  if (iv1 /= 12 .or. iv(inits) == alg-1) go to 230
     m = 1
     write(pu,220) sh(alg), iv(inits)
 220 format(/" nondefault values...."/" init",a1,"      iv(25) =",i3)
 230  if (iv(dtype) == iv(dtype0)) go to 250
     if (m == 0) write(pu,260) which
     m = 1
     write(pu,240) iv(dtype)
 240     format(" dtype      iv(16) =",i3)
 250  i = 1
  j = jlim(alg)
  k = epslon
  l = iv(parsav)
  ndfalt = ndflt(alg)
  do 290 ii = 1, ndfalt
     if (v(k) .speq. v(l)) go to 280
          if (m == 0) write(pu,260) which
 260          format(/1x,3a4,"alues...."/)
          m = 1
          write(pu,270) vn(1,i), vn(2,i), k, v(k)
 270          format(1x,2a4,".. v(",i2,") =",d15.7)
 280     k = k + 1
     l = l + 1
     i = i + 1
     if (i == j) i = ijmp
 290     continue

  iv(dtype0) = iv(dtype)
  parsv1 = iv(parsav)
  call vcopy(iv(nvdflt), v(parsv1), v(epslon))
  return

 300  iv(1) = 15
  if (pu == 0) then
    return
  end if
  write(pu,310) liv, miv2
 310  format(/" /// liv =",i5," must be at least",i5)
  if (liv < miv1) then
    return
  end if
  if (lv < iv(lastv)) go to 320
  return

 320  iv(1) = 16
  if (pu == 0) then
    return
  end if
  write(pu,330) lv, iv(lastv)
 330  format(/" /// lv =",i5," must be at least",i5)
  return

 340  iv(1) = 67
  if (pu == 0) then
    return
  end if
  write(pu,350) alg
 350  format(/" /// alg =",i5," must be 1 or 2")

end subroutine parck
!*******************************************************************************

real(kind=double) function reldst (p, d, x, x0)
!
!*******************************************************************************
!
!! RELDST compute and return relative difference between x and x0
!   nl2sol version 2.2
!
implicit none
  integer p
  real(kind=double) d(p), x(p), x0(p)
  integer i
  real(kind=double) emax, t, xmax

  emax = 0.0_double
  xmax = 0.0_double
  do i = 1, p
     t = abs(d(i) * (x(i) - x0(i)))
     if (emax < t) emax = t
     t = d(i) * (abs(x(i)) + abs(x0(i)))
     if (xmax < t) xmax = t
  end do
  reldst = 0.0_double
  if (xmax > 0.0_double) reldst = emax / xmax

end function reldst
!********************************************************************

real(kind=double) function rmdcon(k)
!
!*******************************************************************************
!
!! RMDCON return machine dependent constants used by nl2sol
!
! comments below contain data statements for various machines.
! to convert to another machine, place a c in column 1 of the
! data statement line(s) that correspond to the current machine
! and remove the c from column 1 of the data statement line(s)
! that correspond to the new machine.
!
!   the constant returned depends on k...
!
!         k = 1... smallest pos. eta such that -eta exists.
!         k = 2... square root of eta.
!         k = 3... unit roundoff = smallest pos. no. machep such
!                  that 1 + machep > 1 .and. 1 - machep < 1.
!         k = 4... square root of machep.
!         k = 5... square root of big (see k = 6).
!         k = 6... largest machine no. big such that -big exists.
!
  implicit none
  integer k
  real(kind=double) big, eta, machep
  integer bigi(4), etai(4), machei(4)
  equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
!
!
!
  data big/1.7e+38_double/, eta/2.939e-39_double/, machep/1.3877788e-17_double/
!
!  body -
!
select case (k)
 case(1)
     rmdcon = eta
 case(2)
     rmdcon = sqrt(256._double*eta)/16._double
 case(3)
     rmdcon = machep
 case(4)
     rmdcon = sqrt(machep)
 case(5)
     rmdcon = sqrt(big/256._double)*16._double
 case(6:)
     rmdcon = big
 case default ! I doubt that this was intentional,but this is the way the old code worked LP 2005
     rmdcon = eta
 end select


end function rmdcon
!*********************************************************

subroutine sgrad2 (alpha, d, eta0, fx, g, irc, n, w, x)
!
!*******************************************************************************
!
!! SGRAD2 computes finite difference gradient by stewart's scheme
!
!      parameters
!
!
!      purpose
!
!        this subroutine uses an embellished form of the finite-differ-
!     ence scheme proposed by stewart (ref. 1) to approximate the
!     gradient of the function f(x), whose values are supplied by
!     reverse communication.
!
!      parameter description
!
!  alpha in  (approximate) diagonal elements of the hessian of f(x).
!      d in  scale vector such that d(i)*x(i), i = 1,...,n, are in
!             comparable units.
!   eta0 in  estimated bound on relative error in the function value...
!             (true value) = (computed value)*(1+e),   where
!             abs(e) <= eta0.
!     fx i/o on input,  fx  must be the computed value of f(x).  on
!             output with irc = 0, fx has been restored to its original
!             value, the one it had when sgrad2 was last called with
!             irc = 0.
!      g i/o on input with irc = 0, g should contain an approximation
!             to the gradient of f near x, e.g., the gradient at the
!             previous iterate.  when sgrad2 returns with irc = 0, g is
!             the desired finite-difference approximation to the
!             gradient at x.
!    irc i/o input/return code... before the very first call on sgrad2,
!             the caller must set irc to 0.  whenever sgrad2 returns a
!             nonzero value for irc, it has perturbed some component of
!             x... the caller should evaluate f(x) and call sgrad2
!             again with fx = f(x).
!      n in  the number of variables (components of x) on which f
!             depends.
!      x i/o on input with irc = 0, x is the point at which the
!             gradient of f is desired.  on output with irc nonzero, x
!             is the point at which f should be evaluated.  on output
!             with irc = 0, x has been restored to its original value
!             (the one it had when sgrad2 was last called with irc = 0)
!             and g contains the desired gradient approximation.
!      w i/o work vector of length 6 in which sgrad2 saves certain
!             quantities while the caller is evaluating f(x) at a
!             perturbed x.
!
!      application and usage restrictions
!
!        this routine is intended for use with quasi-newton routines
!     for unconstrained minimization (in which case  alpha  comes from
!     the diagonal of the quasi-newton hessian approximation).
!
!      algorithm notes
!
!        this code departs from the scheme proposed by stewart (ref. 1)
!     in its guarding against overly large or small step sizes and its
!     handling of special cases (such as zero components of alpha or g).
!
!  Reference:
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference
!    Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages. 72-83, 1967.
!
implicit none
  integer irc, n
  real(kind=double) alpha(n), d(n), eta0, fx, g(n), w(6), x(n)

  integer fh, fx0, hsave, i, xisave
  real(kind=double) aai, afx, afxeta, agi, alphai, axi, axibar
  real(kind=double) discon, eta, gi, h, hmin
  real(kind=double) c2000, four, hmax0, hmin0, h0, machep, one, p002
  real(kind=double) three, two

  parameter (c2000=2.0e+3_double, four=4.0_double, hmax0=0.02_double, hmin0=5.0e+1_double )
  parameter ( one=1.0_double, p002=0.002_double, three=3.0_double )
  parameter ( two=2.0_double)

  parameter (fh=3, fx0=4, hsave=5, xisave=6)
!
 if(irc == 0) then
       goto 100
 elseif(irc < 0) then
       goto 140
 else
       goto 210
 endif

!
!      fresh start -- get machine-dependent constants
!
!     store machep in w(1) and h0 in w(2), where machep is the unit
!     roundoff (the smallest positive number such that
!     1 + machep > 1  and  1 - machep < 1),  and  h0 is the
!     square-root of machep.
!
 100  w(1) = rmdcon(3)
  w(2) = sqrt(w(1))
!
  w(fx0) = fx
!
!      increment  i  and start computing  g(i)
!
 110  i = iabs(irc) + 1
  if (i > n) go to 300
     irc = i
     afx = abs(w(fx0))
     machep = w(1)
     h0 = w(2)
     hmin = hmin0 * machep
     w(xisave) = x(i)
     axi = abs(x(i))
     axibar = max (axi, one/d(i))
     gi = g(i)
     agi = abs(gi)
     eta = abs(eta0)
     if (afx > 0.0_double) eta = max (eta, agi*axi*machep/afx)
     alphai = alpha(i)
     if (alphai .speq. 0.0_double) go to 170
     if (  (gi .speq. 0.0_double) .or. (fx .speq. 0.0_double)  ) go to 180
     afxeta = afx*eta
     aai = abs(alphai)
!
!  compute h = stewart*s forward-difference step size.
!
     if (gi**2 <= afxeta*aai) go to 120
          h = two*sqrt(afxeta/aai)
          h = h*(one - aai*h/(three*aai*h + four*agi))
          go to 130
 120     h = two*(afxeta*agi/(aai**2))**(one/three)
     h = h*(one - two*agi/(three*aai*h + four*agi))
!
!         ensure that  h  is not insignificantly small
!
 130     h = max (h, hmin*axibar)
!
!        use forward difference if bound on truncation error is at
!        most 10**-3.
!
     if (aai*h <= p002*agi) go to 160
!
!        compute h = stewart*s step for central difference.
!
     discon = c2000*afxeta
     h = discon/(agi + sqrt(gi**2 + aai*discon))
!
!         ensure that  h  is neither too small nor too big
!
     h = max (h, hmin*axibar)
     if (h >= hmax0*axibar) h = axibar * h0**(two/three)
!
!  compute central difference
!
     irc = -i
     go to 200

 140     h = -w(hsave)
     i = iabs(irc)
     if (h > 0.0_double) go to 150
     w(fh) = fx
     go to 200

 150     g(i) = (w(fh) - fx) / (two * h)
     x(i) = w(xisave)
     go to 110
!
!  Compute forward differences in various cases
!
 160     if (h >= hmax0*axibar) h = h0 * axibar
     if (alphai*gi < 0.0_double) h = -h
     go to 200
 170     h = axibar
     go to 200
 180     h = h0 * axibar

 200     x(i) = w(xisave) + h
     w(hsave) = h
     return
!
!      compute actual forward difference
!
 210     g(irc) = (fx - w(fx0)) / w(hsave)
     x(irc) = w(xisave)
     go to 110
!
!  Restore fx and indicate that g has been computed
!
 300  fx = w(fx0)
  irc = 0

end subroutine sgrad2
!*****************************************************************************

subroutine slvmul(p, y, s, x)
!
!*******************************************************************************
!
!! SLVMUL sets y = s * x,  s = p x p symmetric matrix.
!   lower triangle of  s  stored rowwise.
!
implicit none
  integer p
  real(kind=double) s(p*(p+1)/2), x(p), y(p)
  integer i, im1, j, k
  real(kind=double) xi
!
  j = 1
  do 10 i = 1, p
     y(i) = dotprd(i, s(j), x)
     j = j + i
 10      continue

  if (p <= 1) then
    return
  end if
  j = 1
  do 40 i = 2, p
     xi = x(i)
     im1 = i - 1
     j = j + 1
     do k = 1, im1
          y(k) = y(k) + s(j)*xi
          j = j + 1
     end do
 40      continue

end subroutine slvmul
!**************************************************************

subroutine smsno ( n, d, x, calcf, iv, liv, lv, v, uiparm, urparm, ufparm )
!
!*******************************************************************************
!
!! SMSNO minimizes a general unconstrained objective function using
!   finite-difference gradients and secant hessian approximations.
!
!
!   purpose
!
!        this routine interacts with SNOIT in an attempt
!     to find an n-vector  x*  that minimizes the (unconstrained)
!     objective function computed by  calcf.  (often the  x*  found is
!     a local minimizer rather than a global one.)
!
!  Parameters:
!
!        the parameters for smsno are the same as those for sumsl
!     (which see), except that calcg is omitted.  instead of calling
!     calcg to obtain the gradient of the objective function at x,
!     smsno calls sgrad2, which computes an approximation to the
!     gradient by finite (forward and central) differences using the
!     method of ref. 1.  the following input component is of interest
!     in this regard (and is not described in sumsl).
!
! v(eta0)  v(42) is an estimated bound on the relative error in the
!             objective function value computed by calcf...
!                  (true value) = (computed value) * (1 + e),
!             where abs(e) <= v(eta0).  default = machep * 10**3,
!             where machep is the unit roundoff.
!
!        the output values iv(nfcall) and iv(ngcall) have different
!     meanings for smsno than for sumsl...
!
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations) excluding those made only for
!             computing gradients.  the input value iv(mxfcal) is a
!             limit on iv(nfcall).
! iv(ngcall)... iv(30) is the number of function evaluations made only
!             for computing gradients.  the total number of function
!             evaluations is thus  iv(nfcall) + iv(ngcall).
!
!  Reference:
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference
!    Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages. 72-83, 1967.
!
implicit none
  integer n, liv, lv
  integer iv(liv), uiparm(*)
  real(kind=double) d(n), x(n), v(lv), urparm(n)
!     dimension v(77 + n*(n+17)/2), uiparm(*), urparm(*)
  external calcf, ufparm
  integer nf
  real(kind=double) fx
  integer nfcall, toobig
  parameter (nfcall=6, toobig=2)
!
  do

    call snoit ( d, fx, iv, liv, lv, n, v, x )

    if ( iv(1) > 2 ) then
      exit
    end if

    nf = iv(nfcall)

    call calcf ( n, x, nf, fx, uiparm, urparm, ufparm )

    if ( nf <= 0 ) then
      iv(toobig) = 1
    end if

  end do

end subroutine smsno
!******************************************************************

subroutine snoit ( d, fx, iv, liv, lv, n, v, x )
!
!*******************************************************************************
!
!! SNOIT iteration driver for smsno...
!   minimize general unconstrained objective function using
!   finite-difference gradients and secant hessian approximations.
!
!   purpose
!
!        this routine interacts with subroutine  sumit  in an attempt
!     to find an n-vector  x*  that minimizes the (unconstrained)
!     objective function  fx = f(x)  computed by the caller.  (often
!     the  x*  found is a local minimizer rather than a global one.)
!
!   parameters
!
!        the parameters for snoit are the same as those for sumsl
!     (which see), except that calcf, calcg, uiparm, urparm, and ufparm
!     are omitted, and a parameter  fx  for the objective function
!     value at x is added.  instead of calling calcg to obtain the
!     gradient of the objective function at x, snoit calls sgrad2,
!     which computes an approximation to the gradient by finite
!     (forward and central) differences using the method of ref. 1.
!     the following input component is of interest in this regard
!     (and is not described in sumsl).
!
! v(eta0)  v(42) is an estimated bound on the relative error in the
!             objective function value computed by calcf...
!                  (true value) = (computed value) * (1 + e),
!             where abs(e) <= v(eta0).  default = machep * 10**3,
!             where machep is the unit roundoff.
!
!        the output values iv(nfcall) and iv(ngcall) have different
!     meanings for smsno than for sumsl...
!
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations) excluding those made only for
!             computing gradients.  the input value iv(mxfcal) is a
!             limit on iv(nfcall).
! iv(ngcall)... iv(30) is the number of function evaluations made only
!             for computing gradients.  the total number of function
!             evaluations is thus  iv(nfcall) + iv(ngcall).
!
!  Reference:
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference
!    Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages. 72-83, 1967.
!
implicit none
  integer liv, lv, n
  integer iv(liv)
  real(kind=double) d(n), fx, x(n), v(lv)
!     dimension v(77 + n*(n+17)/2)
!

  integer alpha, g1, i, iv1, j, k, w
  real(kind=double)  zero

  integer eta0, f, g, lmat, nextv, nfgcal, ngcall
  integer niter, sgirc, toobig, vneed

  parameter ( eta0=42, f=10, g=28, lmat=42, nextv=47 )
  parameter ( nfgcal=7, ngcall=30, niter=31, sgirc=57 )
  parameter ( toobig=2, vneed=4)
!
  parameter ( zero=0._double)
!
  iv1 = iv(1)
  if (iv1 == 1) go to 10
  if (iv1 == 2) go to 50
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + 2*n + 6
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)

 20   call sumit(d, fx, v(g1), iv, liv, lv, n, v, x)
  ! Note that this function should be rewritten without gotos - LP august 2005
  if (iv(1) == 2) then
    goto 30
  elseif(iv(1) < 2) then
    goto 999
  else
    goto 70
  endif
!
!   compute gradient
!
 30   if (iv(niter) == 0) call vscopy(n, v(g1), zero)
  j = iv(lmat)
  k = g1 - n
  do i = 1, n
     v(k) = dotprd(i, v(j), v(j))
     k = k + 1
     j = j + i
  end do
!
!      undo increment of iv(ngcall) done by sumit
!
  iv(ngcall) = iv(ngcall) - 1
!
!      store return code from sgrad2 in iv(sgirc)
!
  iv(sgirc) = 0
!
!      x may have been restored, so copy back fx...
!
  fx = v(f)
  go to 60
!
!      gradient loop
!
 50   if (iv(toobig) == 0) go to 60
  iv(nfgcal) = 0
  go to 10

 60   g1 = iv(g)
  alpha = g1 - n
  w = alpha - 6
  call sgrad2(v(alpha), d, v(eta0), fx, v(g1), iv(sgirc), n, v(w),x)
  if (iv(sgirc) == 0) go to 10
     iv(ngcall) = iv(ngcall) + 1
     return

 70   if (iv(1) /= 14) then
        return
      end if
!
!   storage allocation
!
  iv(g) = iv(nextv) + n + 6
  iv(nextv) = iv(g) + n
  if (iv1 /= 13) go to 10

 999  continue

  return
end subroutine snoit
!******************************************************

logical function stopx(idummy)
!
!*******************************************************************************
!
!! STOPX ???
!
!    **parameters...
implicit none
  integer, optional:: idummy ! defined as optional in june 2005 - LP, Uchicago
!
!
!    **purpose...
!     this function may serve as the stopx (asynchronous interruption)
!     function for the nl2sol (nonlinear least-squares) package at
!     those installations which do not wish to implement a
!     dynamic stopx.
!
!    **algorithm notes...
!     at installations where the nl2sol system is used
!     interactively, this dummy stopx should be replaced by a
!     function that returns .true. if and only if the interrupt
!     (break) key has been pressed since the last call on stopx.
!
!
  stopx = .false.

end function stopx

subroutine sumit(d, fx, g, iv, liv, lv, n, v, x)
!
!*******************************************************************************
!
!! SUMIT carry out  (unconstrained minimization) iterations, using
!   double-dogleg/bfgs steps.
!
!  parameter usage
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! iv... integer value array.
! liv.. length of iv (at least 60).
! lv... length of v (at least 71 + n*(n+13)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... vector of parameters to be optimized.
!
!   discussion
!
!        parameters iv, n, v, and x are the same as the corresponding
!     ones to sumsl (which see), except that v can be shorter (since
!     the part of v that sumsl uses for storing g is not needed).
!     moreover, compared with sumsl, iv(1) may have the two additional
!     output values 1 and 2, which are explained below, as is the use
!     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!     output value from sumsl (and smsno), is not referenced by
!     sumit or the subroutines it calls.
!        fx and g need not have been initialized when sumit is called
!     with iv(1) = 12, 13, or 14.
!
! iv(1) = 1 means the caller should set fx to f(x), the function value
!             at x, and call sumit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             (e.g. if overflow would occur), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
!             nore fx and try a smaller step.  the parameter nf that
!             sumsl passes to calcf (for possible use by calcg) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient vector
!             of f at x, and call sumit again, having changed none of
!             the other parameters except possibly the scale vector d
!             when iv(dtype) = 0.  the parameter nf that sumsl passes
!             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
!             evaluated, then the caller may set iv(nfgcal) to 0, in
!             which case sumit will return with iv(1) = 65.
!
implicit none
  integer liv, lv, n
!
  integer iv(liv)
  real(kind=double) d(n), fx, g(n), v(lv), x(n)
  integer dg1, dummy, g01, i, k, l, lstgst, nwtst1, step1
  integer        temp1, w, x01, z
  real(kind=double) t
  real(kind=double) half, negone, one, onep2, zero

  integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif
  integer gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0
  integer lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, nfgcal
  integer ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, toobig
  integer tuner4, tuner5, vneed, xirc, x0

  parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33 )
  parameter ( mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6 )
  parameter ( nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8 )
  parameter ( restor=9, step=40, stglim=11, stlstg=41, toobig=2 )
  parameter ( vneed=4, xirc=13, x0=43)

  parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13 )
  parameter ( fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42 )
  parameter ( lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7 )
  parameter ( radfac=16, radius=8, rad0=9, reldx=17, tuner4=29 )
  parameter ( tuner5=30)

  parameter (half=0.5_double, negone=-1._double, one=1._double, onep2=1.2_double, zero=0._double)
!
  i = iv(1)
  if (i == 1) go to 50
  if (i == 2) go to 60
!
!   check validity of iv and v input values
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+13)/2
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2

  select case (i)
  case(1:6)
    goto 180
  case(7,9)
    goto 120
  case(8)
    goto 90
  case(10:11)
    goto 10
  case(12)
    goto 20
  case default
    return
  end select

!   storage allocation
!
10    l = iv(lmat)
  iv(x0) = l + n*(n+1)/2
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(g0) = iv(stlstg) + n
  iv(nwtstp) = iv(g0) + n
  iv(dg) = iv(nwtstp) + n
  iv(nextv) = iv(dg) + n
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0_double
  if (v(dinit) >= 0.0_double) call vscopy(n, d, v(dinit))
  if (iv(inith) /= 1) go to 40
!
!      set the initial hessian approximation to diag(d)**-2
!
     l = iv(lmat)
     call vscopy(n*(n+1)/2, v(l), zero)
     k = l - 1
     do 30 i = 1, n
          k = k + i
          t = d(i)
          if (t <= 0.0_double) t = one
          v(k) = t
 30           continue
!
!   compute initial function value
!
 40   iv(1) = 1
  return
!
 50   v(f) = fx
  if (iv(mode) >= 0) go to 180
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 300
!
!   make sure gradient could be computed
!
 60   if (iv(nfgcal) /= 0) go to 70
     iv(1) = 65
     go to 300

 70   dg1 = iv(dg)
  call vvmulp(n, v(dg1), g, d, -1)
  v(dgnorm) = v2norm(n, v(dg1))

  if (iv(cnvcod) /= 0) go to 290
  if (iv(mode) == 0) go to 250
!
!   allow first step to have scaled 2-norm at most v(lmax0)
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  main loop
!
!   print iteration summary, check iteration limit
!
 80    call itsum(d, g, iv, liv, lv, n, v, x)
 90   k = iv(niter)
  if (k < iv(mxiter)) go to 100
     iv(1) = 10
     go to 300
!
!   update radius
!
 100  iv(niter) = k + 1
  if(k>0)v(radius) = v(radfac) * v(dstnrm)
!
!   initialize for start of next iteration
!
  g01 = iv(g0)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1
!
!      copy x to x0, g to g0
!
  call vcopy(n, v(x01), x)
  call vcopy(n, v(g01), g)
!
!   check stopx and function evaluation limit
!
 110  if (.not. stopx(dummy)) go to 130
     iv(1) = 11
     go to 140
!
!      come here when restarting after func. eval. limit or stopx.
!
 120  if (v(f) >= v(f0)) go to 130
     v(radfac) = one
     k = iv(niter)
     go to 100

 130  if (iv(nfcall) < iv(mxfcal)) go to 150
     iv(1) = 9
 140     if (v(f) >= v(f0)) go to 300
!
!         in case of stopx or function evaluation limit with
!         improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 240
!
!  Compute candidate step
!
 150  step1 = iv(step)
  dg1 = iv(dg)
  nwtst1 = iv(nwtstp)
  if (iv(kagqt) >= 0) go to 160
     l = iv(lmat)
     call livmul(n, v(nwtst1), v(l), g)
     v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
     call litvmu(n, v(nwtst1), v(l), v(nwtst1))
     call vvmulp(n, v(step1), v(nwtst1), d, 1)
     v(dst0) = v2norm(n, v(step1))
     call vvmulp(n, v(dg1), v(dg1), d, -1)
     call ltvmul(n, v(step1), v(l), v(dg1))
     v(gthg) = v2norm(n, v(step1))
     iv(kagqt) = 0
 160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
  if (iv(irc) == 6) go to 180
!
!   check whether evaluating f(x0 + step) looks worthwhile
!
  if (v(dstnrm) <= 0.0_double) go to 180
  if (iv(irc) /= 5) go to 170
  if (v(radfac) <= one) go to 170
  if (v(preduc) <= onep2 * v(fdif)) go to 180
!
!  Compute f(x0 + step)
!
 170  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 180  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 190
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))

 190  k = iv(irc)

  select case (k)
  case(1,5)
     goto 200
  case(2:4)
     goto 230
  case(6)
     goto 210
  case(7:12)
      goto 220
  case(13)
       goto 280
  case(14)
       goto 250
  end select


!      recompute step with changed radius
!
 200     v(radius) = v(radfac) * v(dstnrm)
     go to 110
!
!   compute step of length v(lmaxs) for singular convergence test.
!
 210  v(radius) = v(lmaxs)
  go to 150
!
!   convergence or false convergence
!
 220  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 290
     if (iv(xirc) == 14) go to 290
          iv(xirc) = 14
!
!  Process acceptable step.
!
 230  if (iv(irc) /= 3) go to 240
     step1 = iv(step)
     temp1 = iv(stlstg)
!
!      set  temp1 = hessian * step  for use in gradient tests
!
     l = iv(lmat)
     call ltvmul(n, v(temp1), v(l), v(step1))
     call lvmul(n, v(temp1), v(l), v(temp1))
!
!   compute gradient
!
 240  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return
!
!   initializations -- g0 = g - g0, etc.
!
 250  g01 = iv(g0)
  call vaxpy(n, v(g01), negone, v(g01), g)
  step1 = iv(step)
  temp1 = iv(stlstg)
  if (iv(irc) /= 3) go to 270
!
!   set v(radfac) by gradient tests
!
!  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
!
     call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
     call vvmulp(n, v(temp1), v(temp1), d, -1)
!
!  Do gradient tests
!
     if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) then
       go to 260
     end if

     if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  then
       go to 270
     end if

 260               v(radfac) = v(incfac)
!
!   update h, loop
!
 270  w = iv(nwtstp)
  z = iv(x0)
  l = iv(lmat)
  call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
!
!  Use the n-vectors starting at v(step1) and v(g01) for scratch.
!
  call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), v(z))
  iv(1) = 2
  go to 80
!
!   misc. details
!
!   bad parameters to assess
!
 280  iv(1) = 64
  go to 300
!
!  Print summary of final iteration and other requested items
!
 290  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 300    call itsum(d, g, iv, liv, lv, n, v, x)

end subroutine sumit
!*********************************************************8

subroutine sumsl(n, d, x, calcf, calcg, iv, liv, lv, v, uiparm, urparm, ufparm)
!
!*******************************************************************************
!
!! SUMSL minimizes a general unconstrained objective function using
!   analytic gradient and hessian approx. from secant update
!
!   purpose
!
!        this routine interacts with subroutine  sumit  in an attempt
!     to find an n-vector  x*  that minimizes the (unconstrained)
!     objective function computed by  calcf.  (often the  x*  found is
!     a local minimizer rather than a global one.)
!
!  parameter usage
!
! n  (input) the number of variables on which  f  depends, i.e.,
!                  the number of components in  x.
! d  (input/output) a scale vector such that  d(i)*x(i),
!                  i = 1,2,...,n,  are all in comparable units.
!                  d can strongly affect the behavior of sumsl.
!                  finding the best choice of d is generally a trial-
!                  and-error process.  choosing d so that d(i)*x(i)
!                  has about the same value for all i often works well.
!                  the defaults provided by subroutine deflt (see iv
!                  below) require the caller to supply d.
! x........ (input/output) before (initially) calling sumsl, the call-
!                  er should set  x  to an initial guess at  x*.  when
!                  sumsl returns,  x  contains the best point so far
!                  found, i.e., the one that gives the least value so
!                  far seen for  f(x).
! calcf.... (input) a subroutine that, given x, computes f(x).  calcf
!                  must be declared external in the calling program.
!                  it is invoked by
!                       call calcf(n, x, nf, f, uiparm, urparm, ufparm)
!                  when calcf is called, nf is the invocation
!                  count for calcf.  nf is included for possible use
!                  with calcg.  if x is out of bounds (e.g., if it
!                  would cause overflow in computing f(x)), then calcf
!                  should set nf to 0.  this will cause a shorter step
!                  to be attempted.  (if x is in bounds, then calcf
!                  should not change nf.)  the other parameters are as
!                  described above and below.  calcf should not change
!                  n, p, or x.
! calcg.... (input) a subroutine that, given x, computes g(x), the gra-
!                  dient of f at x.  calcg must be declared external in
!                  the calling program.  it is invoked by
!                       call calcg(n, x, nf, g, uiparm, urparm, ufaprm)
!                  when calcg is called, nf is the invocation
!                  count for calcf at the time f(x) was evaluated.  the
!                  x passed to calcg is usually the one passed to calcf
!                  on either its most recent invocation or the one
!                  prior to it.  if calcf saves intermediate results
!                  for use by calcg, then it is possible to tell from
!                  nf whether they are valid for the current x (or
!                  which copy is valid if two copies are kept).  if g
!                  cannot be computed at x, then calcg should set nf to
!                  0.  in this case, sumsl will return with iv(1) = 65.
!                  (if g can be computed at x, then calcg should not
!                  changed nf.)  the other parameters to calcg are as
!                  described above and below.  calcg should not change
!                  n or x.
! iv....... (input/output) an integer value array of length liv (see
!                  below) that helps control the sumsl algorithm and
!                  that is used to store various intermediate quanti-
!                  ties.  of particular interest are the initialization/
!                  return code iv(1) and the entries in iv that control
!                  printing and limit the number of iterations and func-
!                  tion evaluations.  see the section on iv input
!                  values below.
! liv...... (input) length of iv array.  must be at least 60.  if liv
!                  is too small, then sumsl returns with iv(1) = 15.
!                  when sumsl returns, the smallest allowed value of
!                  liv is stored in iv(lastiv) -- see the section on
!                  iv output values below.  (this is intended for use
!                  with extensions of sumsl that handle constraints.)
! lv....... (input) length of v array.  must be at least 71+n*(n+15)/2.
!                  (at least 77+n*(n+17)/2 for smsno, at least
!                  78+n*(n+12) for humsl).  if lv is too small, then
!                  sumsl returns with iv(1) = 16.  when sumsl returns,
!                  the smallest allowed value of lv is stored in
!                  iv(lastv) -- see the section on iv output values
!                  below.
! v........ (input/output) a floating-point value array of length lv
!                  (see below) that helps control the sumsl algorithm
!                  and that is used to store various intermediate
!                  quantities.  of particular interest are the entries
!                  in v that limit the length of the first step
!                  attempted (lmax0) and specify convergence tolerances
!                  (afctol, lmaxs, rfctol, sctol, xctol, xftol).
! uiparm... (input) user integer parameter array passed without change
!                  to calcf and calcg.
! urparm... (input) user floating-point parameter array passed without
!                  change to calcf and calcg.
! ufparm... (input) user external subroutine or function passed without
!                  change to calcf and calcg.
!
!   iv input values (from subroutine deflt)
!
! iv(1)...  on input, iv(1) should have a value between 0 and 14......
!             0 and 12 mean this is a fresh start.  0 means that
!                  deflt(2, iv, liv, lv, v)
!             is to be called to provide all default values to iv and
!             v.  12 (the value that deflt assigns to iv(1)) means the
!             caller has already called deflt and has possibly changed
!             some iv and/or v entries to non-default values.
!             13 means deflt has been called and that sumsl (and
!             sumit) should only do their storage allocation.  that is,
!             they should set the output components of iv that tell
!             where various subarrays arrays of v begin, such as iv(g)
!             (and, for humsl and humit only, iv(dtol)), and return.
!             14 means that a storage has been allocated (by a call
!             with iv(1) = 13) and that the algorithm should be
!             started.  when called with iv(1) = 13, sumsl returns
!             iv(1) = 14 unless liv or lv is too small (or n is not
!             positive).  default = 12.
! iv(inith).... iv(25) tells whether the hessian approximation h should
!             be initialized.  1 (the default) means sumit should
!             initialize h to the diagonal matrix whose i-th diagonal
!             element is d(i)**2.  0 means the caller has supplied a
!             cholesky factor  l  of the initial hessian approximation
!             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42))
!             (and stored compactly by rows).  note that iv(lmat) may
!             be initialized by calling sumsl with iv(1) = 13 (see
!             the iv(1) discussion above).  default = 1.
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcf) allowed.  if this number does not suf-
!             fice, then sumsl returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then sumsl returns with
!             iv(1) = 10.  default = 150.
! iv(outlev)... iv(19) controls the number and length of iteration sum-
!             mary lines printed (by itsum).  iv(outlev) = 0 means do
!             not print any summary lines.  otherwise, print a summary
!             line after each abs(iv(outlev)) iterations.  if iv(outlev)
!             is positive, then summary lines of length 78 (plus carri-
!             age control) are printed, including the following...  the
!             iteration and function evaluation counts, f = the current
!             function value, relative difference in function values
!             achieved by the latest step (i.e., reldf = (f0-v(f))/f01,
!             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and
!             v(f0) is the function value from the previous itera-
!             tion), the relative function reduction predicted for the
!             step just taken (i.e., preldf = v(preduc) / f01, where
!             v(preduc) is described below), the scaled relative change
!             in x (see v(reldx) below), the step parameter for the
!             step just taken (stppar = 0 means a full newton step,
!             between 0 and 1 means a relaxed newton step, between 1
!             and 2 means a double dogleg step, greater than 2 means
!             a scaled down cauchy step -- see subroutine dbldog), the
!             2-norm of the scale vector d times the step just taken
!             (see v(dstnrm) below), and npreldf, i.e.,
!             v(nreduc)/f01, where v(nreduc) is described below -- if
!             npreldf is positive, then it is the relative function
!             reduction predicted for a newton step (one with
!             stppar = 0).  if npreldf is negative, then it is the
!             negative of the relative function reduction predicted
!             for a step computed with step bound v(lmaxs) for use in
!             testing for singular convergence.
!                  if iv(outlev) is negative, then lines of length 50
!             are printed, including only the first 6 items listed
!             above (through reldx).
!             default = 1.
! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
!             fresh start or any changed v values on a restart.
!             iv(parprt) = 0 means skip this printing.  default = 1.
! iv(prunit)... iv(21) is the output unit number on which all printing
!             is done.  iv(prunit) = 0 means suppress all printing.
!             default = standard output unit (unit 6 on most systems).
! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
!             well as the gradient and the scale vector d).
!             iv(solprt) = 0 means skip this printing.  default = 1.
! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
!             ing.  these consist of the function value, the scaled
!             relative change in x caused by the most recent step (see
!             v(reldx) below), the number of function and gradient
!             evaluations (calls on calcf and calcg), and the relative
!             function reductions predicted for the last step taken and
!             for a newton step (or perhaps a step bounded by v(lmaxs)
!             -- see the descriptions of preldf and npreldf under
!             iv(outlev) above).
!             iv(statpr) = 0 means skip this printing.
!             iv(statpr) = -1 means skip this printing as well as that
!             of the one-line termination reason message.  default = 1.
! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
!             (on a fresh start only).  iv(x0prt) = 0 means skip this
!             printing.  default = 1.
!
!   (selected) iv output values --
!
! iv(1)........ on output, iv(1) is a return code....
!             3 = x-convergence.  the scaled relative difference (see
!                  v(reldx)) between the current parameter vector x and
!                  a locally optimal parameter vector is very likely at
!                  most v(xctol).
!             4 = relative function convergence.  the relative differ-
!                  ence between the current function value and its lo-
!                  cally optimal value is very likely at most v(rfctol).
!             5 = both x- and relative function convergence (i.e., the
!                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
!             6 = absolute function convergence.  the current function
!                  value is at most v(afctol) in absolute value.
!             7 = singular convergence.  the hessian near the current
!                  iterate appears to be singular or nearly so, and a
!                  step of length at most v(lmaxs) is unlikely to yield
!                  a relative function decrease of more than v(sctol).
!             8 = false convergence.  the iterates appear to be converg-
!                  ing to a noncritical point.  this may mean that the
!                  convergence tolerances (v(afctol), v(rfctol),
!                  v(xctol)) are too small for the accuracy to which
!                  the function and gradient are being computed, that
!                  there is an error in computing the gradient, or that
!                  the function or gradient is discontinuous near x.
!             9 = function evaluation limit reached without other con-
!                  vergence (see iv(mxfcal)).
!            10 = iteration limit reached without other convergence
!                  (see iv(mxiter)).
!            11 = stopx returned .true. (external interrupt).  see the
!                  usage notes below.
!            14 = storage has been allocated (after a call with
!                  iv(1) = 13).
!            17 = restart attempted with n changed
!            18 = d has a negative component and iv(dtype) <= 0.
!            19...43 = v(iv(1)) is out of range.
!            63 = f(x) cannot be computed at the initial x.
!            64 = bad parameters passed to assess (which should not
!                  occur).
!            65 = the gradient could not be computed at x (see calcg
!                  above).
!            67 = bad first parameter to deflt.
!            80 = iv(1) was out of range.
!            81 = n is not positive.
! iv(g)........ iv(28) is the starting subscript in v of the current
!             gradient vector (the one corresponding to x).
! iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is
!             only set if liv is at least 44.)
! iv(lastv).... iv(45) is the least acceptable value of lv.  (it is
!             only set if liv is large enough, at least iv(lastiv).)
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations).
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcg).
! iv(niter).... iv(31) is the number of iterations performed.
!
!   (selected) v input values (from subroutine deflt)
!
! v(bias)..... v(43) is the bias parameter used in subroutine dbldog --
!             see that subroutine for details.  default = 0.8.
! v(afctol)... v(31) is the absolute function convergence tolerance.
!             if sumsl finds a point where the function value is less
!             than v(afctol) in absolute value, and if sumsl does not
!             return with iv(1) = 3, 4, or 5, then it returns with
!             iv(1) = 6.  this test can be turned off by setting
!             v(afctol) to zero.  default = max(10**-20, machep**2),
!             where machep is the unit roundoff.
! v(dinit).... v(38), if nonnegative, is the value to which the scale
!             vector d is initialized.  default = -1.
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that sumsl attempts.  this parameter can
!             markedly affect the performance of sumsl.
! v(lmaxs).... v(36) is used in testing for singular convergence -- if
!             the function reduction predicted for a step of length
!             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where
!             f0  is the function value at the start of the current
!             iteration, and if sumsl does not return with iv(1) = 3,
!             4, 5, or 6, then it returns with iv(1) = 7.  default = 1.
! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0)
!             at the start of the current iteration, where  f0  is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then sumsl returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
! v(sctol).... v(37) is the singular convergence tolerance -- see the
!             description of v(lmaxs) above.
! v(tuner1)... v(26) helps decide when to check for false convergence.
!             this is done if the actual function decrease from the
!             current step is no more than v(tuner1) times its predict-
!             ed value.  default = 0.1.
! v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
!             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then sumsl returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.)
!             default = machep**0.5, where machep is the unit roundoff.
! v(xftol).... v(34) is the false convergence tolerance.  if a step is
!             tried that gives no more than v(tuner1) times the predict-
!             ed function decrease and that has v(reldx) <= v(xftol),
!             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or
!             7, then it returns with iv(1) = 8.  (see the description
!             of v(reldx) below.)  default = 100*machep, where
!             machep is the unit roundoff.
! v(*)........ deflt supplies to v a number of tuning constants, with
!             which it should ordinarily be unnecessary to tinker.  see
!             section 17 of version 2.2 of the nl2sol usage summary
!             (i.e., the appendix to ref. 1) for details on v(i),
!             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx,
!             tuner2, tuner3, tuner4, tuner5.
!
!   (selected) v output values
!
! v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the
!             most recently computed gradient.
! v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the
!             current step.
! v(f)........ v(10) is the current function value.
! v(f0)....... v(13) is the function value at the start of the current
!             iteration.
! v(nreduc)... v(6), if positive, is the maximum function reduction
!             possible according to the current model, i.e., the func-
!             tion reduction predicted for a newton step (i.e.,
!             step = -h**-1 * g,  where  g  is the current gradient and
!             h is the current hessian approximation).
!                  if v(nreduc) is negative, then it is the negative of
!             the function reduction predicted for a step computed with
!             a step bound of v(lmaxs) for use in testing for singular
!             convergence.
! v(preduc)... v(7) is the function reduction predicted (by the current
!             quadratic model) for the current step.  this (divided by
!             v(f0)) is used in testing for relative function
!             convergence.
! v(reldx).... v(17) is the scaled relative change in x caused by the
!             current step, computed as
!                  max(abs(d(i)*(x(i)-x0(i)), 1 <= i <= p) /
!                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p),
!             where x = x0 + step.
!
!  notes
!
!   algorithm notes
!
!        this routine uses a hessian approximation computed from the
!     bfgs update (see ref 3).  only a cholesky factor of the hessian
!     approximation is stored, and this is updated using ideas from
!     ref. 4.  steps are computed by the double dogleg scheme described
!     in ref. 2.  the steps are assessed as in ref. 1.
!
!   usage notes
!
!        after a return with iv(1) <= 11, it is possible to restart,
!     i.e., to change some of the iv and v input values described above
!     and continue the algorithm from the point where it was interrupt-
!     ed.  iv(1) should not be changed, nor should any entries of iv
!     and v other than the input values (those supplied by deflt).
!        those who do not wish to write a calcg which computes the
!     gradient analytically should call smsno rather than sumsl.
!     smsno uses finite differences to compute an approximate gradient.
!        those who would prefer to provide f and g (the function and
!     gradient) by reverse communication rather than by writing subrou-
!     tines calcf and calcg may call on sumit directly.  see the com-
!     ments at the beginning of sumit.
!        those who use sumsl interactively may wish to supply their
!     own stopx function, which should return .true. if the break key
!     has been pressed since stopx was last invoked.  this makes it
!     possible to externally interrupt sumsl (which will return with
!     iv(1) = 11 if stopx returns .true.).
!        storage for g is allocated at the end of v.  thus the caller
!     may make v longer than specified above and may allow calcg to use
!     elements of g beyond the first n as scratch storage.
!
!   portability notes
!
!        the sumsl distribution tape contains both single- and double-
!     precision versions of the sumsl source code, so it should be un-
!     necessary to change precisions.
!        only the function rmdcon contains machine-dependent
!     constants.  to change from one machine to another, it should
!     suffice to change the (few) relevant lines in these functions.
!        intrinsic functions are explicitly declared.  on certain com-
!     puters (e.g. univac), it may be necessary to comment out these
!     declarations.  so that this may be done automatically by a simple
!     program, such declarations are preceded by a comment having c/+
!     in columns 1-3 and blanks in columns 4-72 and are followed by
!     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
!        the sumsl source code is expressed in 1966 ansi standard
!     fortran.  it may be converted to fortran 77 by commenting out all
!     lines that fall between a line having c/6 in columns 1-3 and a
!     line having c/7 in columns 1-3 and by removing (i.e., replacing
!     by a blank) the c in column 1 of the lines that follow the c/7
!     line and precede a line having c/ in columns 1-2 and blanks in
!     columns 3-72.  these changes convert some data statements into
!     parameter statements, convert some variables from real to
!     character*4, and make the data statements that initialize these
!     variables use character strings delimited by primes instead
!     of hollerith constants.  (such variables and data statements
!     appear only in modules itsum and parck.  parameter statements
!     appear nearly everywhere.)  these changes also add save state-
!     ments for variables given machine-dependent constants by rmdcon.
!
!   references
!
!    J E Dennis, David Gay, and R E Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    J E Dennis, H H W Mei,
!    Two New Unconstrained Optimization Algorithms Which Use
!    Function and Gradient Values,
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    J E Dennis, J J More,
!    Quasi-Newton Methods, Motivation and Theory,
!    SIAM Review,
!    Volume 19, pages 46-89, 1977.
!
!    D Goldfarb,
!    Factorized Variable Metric Methods for Unconstrained Optimization,
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
implicit none
  integer n, liv, lv
  integer iv(liv), uiparm(1)
  real(kind=double) d(n), x(n), urparm(n),v
     dimension v(71 + n*(n+15)/2) !, uiparm(*), urparm(*)

  integer g1, iv1, nf
  real(kind=double) f,ufparm
  integer nextv, nfcall, nfgcal, g, toobig, vneed

  parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, vneed=4)
!
  external ufparm

!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + n
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)

 20   call sumit(d, f, v(g1), iv, liv, lv, n, v, x)

 if( iv(1) == 2) then
       goto 40
 elseif(iv(1) < 2) then
      goto 30
 else
      goto 50
 endif

 30   nf = iv(nfcall)
  call calcf(n, x, nf, f, uiparm, urparm, ufparm)
  if (nf <= 0) iv(toobig) = 1
  go to 20
 40   call calcg(n, x, iv(nfgcal), v(g1), uiparm, urparm, ufparm)
  go to 20

 50   if (iv(1) /= 14) then
        return
      end if
!
!   storage allocation
!
  iv(g) = iv(nextv)
  iv(nextv) = iv(g) + n
  if (iv1 /= 13) go to 10

!******************************************************************************
end subroutine sumsl
!******************************************************************************


real(kind=double) function v2norm(p, x)
!
!*******************************************************************************
!
!! V2NORM returns the 2-norm of the p-vector x, taking
!   care to avoid the most likely underflows.
!
implicit none
  integer p
!
  real(kind=double) x(p)
  integer i, j
  real(kind=double)  r, scale
  real(kind=double), save :: sqteta = 0.0_double
  real(kind=double) t, xi

!
  v2norm = 0.0_double

  if (p <= 0 ) then
    return
  end if

  if ( all ( x(1:p) .speq. 0.0_double ) ) then
    return
  end if

  scale = 0.0_double
  do i = 1, p
    if ( x(i) .spne. 0.0_double ) then
      scale = abs(x(i))
      exit
    end if
  end do

  if ( scale .speq. 0.0_double ) then
    return
  end if

  if (i < p) go to 40
     v2norm = scale
     return
 40   t = 1.0_double
  if (sqteta .speq. 0.0_double) sqteta = rmdcon(2)
!
!      sqteta is (slightly larger than) the square root of the
!      smallest positive floating point number on the machine.
!      the tests involving sqteta are done to prevent underflows.
!
  j = i + 1
  do i = j, p
    xi = abs(x(i))
    if (xi <= scale) then
      r = xi / scale
      if (r > sqteta) t = t + r*r
    else
      r = scale / xi
      if (r <= sqteta) r = 0.0_double
      t = 1.0_double  +  t * r*r
      scale = xi
    end if
  end do

  v2norm = scale * sqrt(t)

end function v2norm
!********************************************************************************


subroutine vaxpy(p, w, a, x, y)
!
!*******************************************************************************
!
!! VAXPY sets w = a*x + y  --  w, x, y = p-vectors, a = scalar
!
  implicit none
!
  integer p
  real(kind=double) a, w(p), x(p), y(p)
  integer i
!
  do i = 1, p
    w(i) = a*x(i) + y(i)
  end do

end subroutine vaxpy
!*******************************************************************************

subroutine vcopy(p, y, x)
!
!*******************************************************************************
!
!! VCOPY sets y = x, where x and y are p-vectors
!
  implicit none
!
  integer p
  real(kind=double) x(p), y(p)

  y(1:p) = x(1:p)

  return
end subroutine vcopy
!*******************************************************************************


subroutine vdflt(alg, lv, v)
!
!*******************************************************************************
!
!! VDFLT supplies default values to V.
!
!   alg = 1 means regression constants.
!   alg = 2 means general unconstrained optimization constants.
!
  implicit none
!
  integer alg, lv
  real(kind=double) v(lv)
  real(kind=double) machep, mepcrt, one, sqteps, three
  integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc
  integer dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc
  integer incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx
  integer rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2
  integer tuner3, tuner4, tuner5, xctol, xftol

  parameter (one=1._double, three=3._double)

  parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44 )
  parameter ( dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39 )
  parameter ( d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48 )
  parameter ( incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21 )
  parameter ( rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49 )
  parameter ( sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28 )
  parameter ( tuner4=29, tuner5=30, xctol=33, xftol=34)
!
  machep = rmdcon(3)
  v(afctol) = 1.e-20_double
  if (machep > 1.e-10_double) v(afctol) = machep**2
  v(decfac) = 0.5_double
  sqteps = rmdcon(4)
  v(dfac) = 0.6_double
  v(delta0) = sqteps
  v(dtinit) = 1.e-6_double
  mepcrt = machep ** (one/three)
  v(d0init) = 1._double
  v(epslon) = 0.1_double
  v(incfac) = 2._double
  v(lmax0) = 1._double
  v(lmaxs) = 1._double
  v(phmnfc) = -0.1_double
  v(phmxfc) = 0.1_double
  v(rdfcmn) = 0.1_double
  v(rdfcmx) = 4._double
  v(rfctol) = max (1.e-12_double, mepcrt**2) ! reduced in size from 10-10, LP UC Dec. 2004
  v(sctol) = v(rfctol)
  v(tuner1) = 0.1_double
  v(tuner2) = 1.e-4_double
  v(tuner3) = 0.75_double
  v(tuner4) = 0.5_double
  v(tuner5) = 0.75_double
  v(xctol) = sqteps
  v(xftol) = 1.e+2_double * machep

  if ( alg < 2 ) then
    v(cosmin) = max (1.e-6_double, 1.e+2_double * machep)
    v(dinit) = 0._double
    v(dltfdc) = mepcrt
    v(dltfdj) = sqteps
    v(fuzz) = 1.5_double
    v(huberc) = 0.7_double
    v(rlimit) = rmdcon(5)
    v(rsptol) = 1.e-3_double
    v(sigmin) = 1.e-4_double
  else
    v(bias) = 0.8_double
    v(dinit) = -1.0_double
    v(eta0) = 1.0e+3_double * machep
  end if

end subroutine vdflt
!*******************************************************************************

subroutine vscopy(p, y, s)
!
!*******************************************************************************
!
!! VSCOPY sets p-vector y to scalar s
!
  implicit none
!
  integer p
  real(kind=double) s, y(p)
!
  y(1:p) = s

  return
end subroutine vscopy
!*******************************************************************************

subroutine vvmulp(n, x, y, z, k)
!
!*******************************************************************************
!
!! VVMULP sets x(i) = y(i) * z(i)**k, 1 <= i <= n (for k = 1 or -1)
!
  implicit none
!
  integer n, k
  real(kind=double) x(n), y(n), z(n)
  integer i
!
  if ( k < 0 ) then
    do i = 1, n
      x(i) = y(i) / z(i)
    end do
  else
    do i = 1, n
      x(i) = y(i) * z(i)
    end do
  end if

  return
end subroutine vvmulp

subroutine wzbfgs (l, n, s, w, y, z)
!
!*******************************************************************************
!
!! WZBFGS compute  y  and  z  for  lupdat  corresponding to bfgs update.
!
!
!  parameter usage
!
! l (i/o) cholesky factor of hessian, a lower triang. matrix stored
!             compactly by rows.
! n (input) order of  l  and length of  s,  w,  y,  z.
! s (input) the step just taken.
! w (output) right singular vector of rank 1 correction to l.
! y (input) change in gradients corresponding to s.
! z (output) left singular vector of rank 1 correction to l.
!
!  notes
!
!   algorithm notes
!
!        when  s  is computed in certain ways, e.g. by  gqtstp  or
!     dbldog,  it is possible to save n**2/2 operations since  (l**t)*s
!     or  l*(l**t)*s is then known.
!        if the bfgs update to l*(l**t) would reduce its determinant to
!     less than eps times its old value, then this routine in effect
!     replaces  y  by  theta*y + (1 - theta)*l*(l**t)*s,  where  theta
!     (between 0 and 1) is chosen to make the reduction factor = eps.
!
  implicit none
!
  integer n
!
  real(kind=double) l(n*(n+1)/2), s(n), w(n), y(n), z(n)
!

  integer i
  real(kind=double) cs, cy, eps, epsrt, shs, ys, theta
  parameter (eps=0.1_double)
!
  call ltvmul(n, w, l, s)
  shs = dotprd(n, w, w)
  ys = dotprd(n, y, s)

  if ( ys < eps*shs) then
     theta = ( 1.0_double - eps) * shs / (shs - ys)
     epsrt = sqrt(eps)
     cy = theta / (shs * epsrt)
     cs = ( 1.0_double + (theta-1.0_double)/epsrt) / shs
  else
    cy = 1.0_double / (sqrt(ys) * sqrt(shs))
    cs = 1.0_double / shs
  end if

  call livmul(n, z, l, y)

  do i = 1, n
    z(i) = cy * z(i) - cs * w(i)
  end do

end subroutine wzbfgs

end module toms611
