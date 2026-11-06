test_that("exact normalisers agree with independent radial and angular quadrature", {
  b <- rbind(step=c(-3,-1.2),log_step=c(-.6,1.1),cos_angle=c(2,-.4),
    cos_diff_theta_to_target=c(2,-1.3))
  psi <- c(0,pi/3,pi,-pi/2)
  active <- c(TRUE,TRUE,TRUE,FALSE)
  z <- ssf_log_normalizer(b,"area",psi,active)
  numerical <- sapply(1:2,function(k) vapply(seq_along(psi),function(i) {
    radial <- integrate(function(r) exp(b[1,k]*r+(b[2,k]+1)*log(r)),0,Inf,
      rel.tol=1e-10)$value
    angular <- integrate(function(phi) exp(b[3,k]*cos(phi)+
      active[i]*b[4,k]*cos(phi-psi[i])),-pi,pi,rel.tol=1e-10)$value
    log(radial*angular)
  },numeric(1)))
  expect_equal(unname(z),numerical,tolerance=1e-9)
  high <- b[,1,drop=FALSE]; high[3,1] <- 1000
  expect_true(all(is.finite(ssf_log_normalizer(high,"area",psi))))
  bad <- b; bad[1,1] <- 0
  expect_error(ssf_log_normalizer(bad,"area"),"positive gamma")
})

test_that("exact emissions match gamma and conditional circular laws in both measures", {
  d <- prepare_ssf_steps(data.frame(ID=c("a","a","b","b"),step=c(.1,.3,.8,.2),
    angle=c(.2,-1,2,.7),relative_target_bearing=c(.1,1,-2,0),
    target_active=c(TRUE,TRUE,TRUE,FALSE)),"area")
  b <- rbind(step=c(-3,-5),log_step=c(.2,-.8),cos_angle=c(1.5,-.4),
    cos_diff_theta_to_target=c(.8,1.1))
  f <- ~ step+log_step+cos_angle+cos_diff_theta_to_target
  exact <- ssf_log_emissions(d,f,b,normalization="exact")
  independent <- sapply(1:2,function(k) {
    a <- b[3,k]; target <- b[4,k]*d$target_active
    vx <- a+target*cos(d$relative_target_bearing)
    vy <- target*sin(d$relative_target_bearing)
    concentration <- sqrt(vx^2+vy^2)
    direction <- atan2(vy,vx)
    dgamma(d$step,shape=b[2,k]+2,rate=-b[1,k],log=TRUE)-log(d$step)+
      concentration*(cos(d$angle-direction)-1)-log(2*pi)-
      log(besselI(concentration,0,expon.scaled=TRUE))
  })
  expect_equal(unname(exact),independent,tolerance=1e-12)
  radial <- convert_density_measure(d,"step_angle")
  br <- convert_movement_coefficients(b,"area","step_angle")
  er <- ssf_log_emissions(radial,f,br,normalization="exact")
  expect_equal(er,exact+log(d$step),tolerance=1e-12)
  P <- matrix(c(.8,.2,.3,.7),2,byrow=TRUE)
  a <- decode_states(exact,P,d$ID); r <- decode_states(er,P,d$ID)
  expect_equal(r$posterior,a$posterior,tolerance=1e-12)
  expect_identical(r$state,a$state)
  expect_equal(r$log_likelihood-a$log_likelihood,sum(log(d$step)),tolerance=1e-12)
})

test_that("exact compiled score agrees with independent finite differences", {
  set.seed(1704)
  d <- prepare_ssf_steps(data.frame(ID=rep(c("a","b","c"),c(3,1,5)),
    step=rgamma(9,2,4),angle=runif(9,-pi,pi),
    relative_target_bearing=c(0,pi,pi/2,0,runif(5,-pi,pi)),
    target_active=c(TRUE,TRUE,TRUE,FALSE,rep(TRUE,5))),"area")
  for (target in c(FALSE,TRUE)) for (k in c(1L,2L,3L,6L)) {
    f <- if (target) ~ cos_angle+step+cos_diff_theta_to_target+log_step else
      ~ cos_angle+step+log_step
    design <- behavSSF:::.design_exact(d,f); p <- ncol(design$X)
    theta <- rnorm(p*k+k*(k-1),sd=.4)
    objective <- behavSSF:::.native_exact_objective(design,k)
    reference <- function(v) {
      natural <- behavSSF:::.exact_transform(v,design,k)
      z <- behavSSF:::.unpack(natural,p,k,colnames(design$X))
      -decode_states(behavSSF:::.exact_emissions(design,z$coefficients),
        z$transition,design$segment)$log_likelihood
    }
    finite_difference <- vapply(seq_along(theta),function(i) {
      plus <- minus <- theta
      plus[i] <- plus[i]+1e-6;minus[i] <- minus[i]-1e-6
      (reference(plus)-reference(minus))/2e-6
    },numeric(1))
    expect_equal(objective$fn(theta),reference(theta),tolerance=1e-11)
    expect_equal(objective$gr(theta),finite_difference,tolerance=2e-6)
    # Invalid transformed gamma parameters must reject the trial, then recover.
    bad <- theta; bad[match("step",colnames(design$X))] <- 1000
    expect_equal(objective$fn(bad),1e100)
    expect_equal(objective$fn(theta),reference(theta),tolerance=1e-11)
  }
  # Cancellation gives a zero resultant: its score must remain finite.
  design <- behavSSF:::.design_exact(d,~ step+log_step+cos_angle+cos_diff_theta_to_target)
  b <- c(-4,0,2,-2)
  v <- behavSSF:::.exact_transform(b,design,1L,inverse=TRUE)
  expect_true(all(is.finite(behavSSF:::.native_exact_objective(design,1L)$gr(v))))
  expect_equal(behavSSF:::.native_exact_objective(design,1L)$fn(numeric()),1e100)
})

test_that("one-state exact fit and natural covariance agree with independent MLE equations", {
  set.seed(42)
  d <- prepare_ssf_steps(data.frame(ID=1,step=rgamma(200,shape=1.4,rate=5),
    angle=runif(200,-1.5,1.5)),"area")
  shape <- uniroot(function(a) log(a)-digamma(a)-
    (log(mean(d$step))-mean(log(d$step))),c(.01,100),tol=1e-11)$root
  angular_moment <- mean(cos(d$angle))
  concentration <- uniroot(function(k) besselI(k,1,expon.scaled=TRUE)/
    besselI(k,0,expon.scaled=TRUE)-angular_moment,c(0,100),tol=1e-11)$root
  expected <- c(-shape/mean(d$step),shape-2,concentration)
  fit <- fit_ssf(d,n_states=1,normalization="exact",long_iter=500,seed=19)
  expect_equal(unname(c(fit$coefficients)),expected,tolerance=2e-5)
  expect_true(fit$diagnostics$converged)
  expect_lt(fit$diagnostics$max_abs_gradient,.01)
  expect_equal(fit$hessian_parameterization,"log_rate_log_shape")
  independent <- function(v) {
    a <- v[2]+2;rate <- -v[1]
    -sum(dgamma(d$step,a,rate,log=TRUE)-log(d$step)+v[3]*cos(d$angle)-
      log(2*pi)-log(besselI(abs(v[3]),0,expon.scaled=TRUE))-abs(v[3]))
  }
  h <- optimHess(c(fit$coefficients),independent)
  expect_equal(fit$covariance,solve(h),tolerance=2e-4)
  r <- fit_ssf(d,n_states=1,normalization="exact",engine="R",long_iter=500,seed=19)
  expect_equal(c(r$coefficients),c(fit$coefficients),tolerance=2e-4)
  expect_equal(r$log_likelihood,fit$log_likelihood,tolerance=1e-8)
})

test_that("an exact fit on a long series passes the unscaled score check", {
  set.seed(4)
  n <- 16434L
  step <- rgamma(n, shape = .72, rate = 1.36)
  angle <- rnorm(n, sd = 1.32)
  steps <- prepare_ssf_steps(data.frame(ID = 1L, step = step,
    angle = atan2(sin(angle), cos(angle))), "area")
  fit <- fit_ssf(steps, n_states = 1L, normalization = "exact", seed = 19)
  shape <- uniroot(function(s) log(s) - digamma(s) -
    (log(mean(step)) - mean(log(step))), c(.01, 100), tol = 1e-12)$root
  kappa <- uniroot(function(k) besselI(k, 1, expon.scaled = TRUE) /
    besselI(k, 0, expon.scaled = TRUE) - mean(cos(angle)), c(0, 100), tol = 1e-12)$root
  expected <- c(-shape / mean(step), shape - 2, kappa)
  expect_true(assess_taxonomy(fit, "Other", "area")$accepted)
  expect_lt(fit$diagnostics$max_abs_gradient, 1e-3)
  expect_equal(c(coef(fit)), expected, tolerance = 1e-6)
  expect_lt(fit$diagnostics$objective_agreement, 1e-8)
  expect_lte(fit$diagnostics$max_abs_gradient, fit$diagnostics$score_refinement$before)
})

test_that("exact input contracts reject controls, missing target context and unsupported terms", {
  d <- data.frame(ID=c("a","a","b"),step=c(.1,.4,.2),angle=c(.2,1,-.4),
    previous_bearing=c(0,.2,1),target_bearing=c(1,2,-1))
  a <- prepare_ssf_steps(d,"area")
  expect_equal(a$cos_diff_theta_to_target,cos(d$previous_bearing+d$angle-d$target_bearing))
  expect_error(prepare_ssf_steps(transform(d,obs=c(1,0,1)),"area"),"observed rows only")
  expect_error(prepare_ssf_steps(transform(d,cos_angle=0),"area"),"disagrees")
  expect_error(prepare_ssf_steps(d,"area",expected_interval=1800),"time column")
  timed <- transform(d,time=c(0,3600,0))
  expect_error(prepare_ssf_steps(timed,"area",time="time",expected_interval=1800),"time gap")
  base <- prepare_ssf_steps(d[c("ID","step","angle")],"area")
  expect_error(fit_ssf(base,~step+log_step+cos_angle+cos_diff_theta_to_target,
    normalization="exact"),"target bearing")
  expect_error(fit_ssf(a,~step+log_step+cos_angle+I(step^2),normalization="exact"),"supports")
  expect_error(fit_ssf(a,~step+log_step+cos_angle+offset(step),normalization="exact"),"offsets")
  expect_error(fit_ssf(d,normalization="exact"),"prepare_ssf_steps")
})
