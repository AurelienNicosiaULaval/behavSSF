test_that("score and curvature screens are invariant to diagonal parameter units", {
  h <- matrix(c(4, 1, 1, 9), 2)
  gradient <- c(.02, -.03)
  scale <- c(1000, .001)
  changed <- h * outer(scale, scale)
  expect_equal(.scaled_ssf_score(gradient, h),
    .scaled_ssf_score(gradient * scale, changed), tolerance = 1e-14)
  expect_equal(.ssf_curvature(h), .ssf_curvature(changed), tolerance = 1e-13)
  expect_true(.positive_ssf_hessian(changed))
  expect_false(.positive_ssf_hessian(matrix(c(1,2,2,1),2)))
})

test_that("finite radial support validates positive linear length coefficients", {
  b <- rbind(step = c(.6, -2), log_step = c(0, -.5), cos_angle = c(1,.2))
  full <- movement_metrics(b, "area")
  disk <- movement_metrics(b, "area", step_max = 3.8)
  expect_false(full$valid[1])
  expect_true(all(disk$valid))
  for (i in 1:2) {
    mass <- function(r) r^(b["log_step",i]+1)*exp(b["step",i]*r)
    z <- integrate(mass,0,3.8,rel.tol=1e-11)$value
    mu <- integrate(function(r) r*mass(r),0,3.8,rel.tol=1e-11)$value/z
    second <- integrate(function(r) r^2*mass(r),0,3.8,rel.tol=1e-11)$value/z
    expect_equal(disk$mean_step[i],mu,tolerance=1e-8)
    expect_equal(disk$sd_step[i],sqrt(second-mu^2),tolerance=1e-8)
  }
  metres <- b; metres["step",] <- b["step",]/1000
  m <- movement_metrics(metres,"area",step_max=3800)
  expect_equal(m$mean_step/1000,disk$mean_step,tolerance=1e-10)
  expect_equal(m$sd_step/1000,disk$sd_step,tolerance=1e-10)
  data <- data.frame(ID=1,stratum=1:2,obs=1,w=1,step=c(.4,4))
  expect_error(prepare_ssf_data(rbind(data,data),step_max=3.8),"stratum")
  data <- data.frame(ID=1,stratum=rep(1:2,each=2),obs=rep(1:0,2),w=1,
    step=c(.4,.8,4,.6))
  expect_error(prepare_ssf_data(data,step_max=3.8),"step_max")
})

test_that("bounded taxonomy ordering runs and covariance contrasts use named states", {
  b <- rbind(step=c(.6,-2),log_step=c(0,-.5),cos_angle=c(1,.2))
  covariance <- diag(.04,6);covariance[3,6]<-covariance[6,3]<-.01
  fit <- structure(list(coefficients=b,n_states=2,step_max=3.8,
    density_measure="area",state=rep(1:2,each=20),transition=matrix(c(.8,.2,.2,.8),2),
    covariance=covariance,diagnostics=list(converged=TRUE,hessian_positive=TRUE,
      max_abs_gradient=0)),class="ssf_fit")
  taxonomy <- ssf_taxonomy("bounded",c("Encamped","Travelling"),
    rules=function(s) c(persistent=ssf_contrast(s,"cos_angle","Encamped","Travelling")$lower>0))
  a <- assess_taxonomy(fit,taxonomy,"area")
  expect_true(a$accepted)
  contrast <- ssf_contrast(a$context,"cos_angle","Encamped","Travelling")
  expect_equal(contrast$estimate,.8)
  expect_equal(contrast$std_error,sqrt(.04+.04-2*.01))
  expect_true(!is.null(rownames(a$context$covariance)))
  permutation <- c(4:6,1:3)
  fit$coefficients<-b[,2:1];fit$state<-3-fit$state
  fit$covariance<-covariance[permutation,permutation]
  other<-assess_taxonomy(fit,taxonomy,"area")
  expect_equal(ssf_contrast(other$context,"cos_angle","Encamped","Travelling"),contrast)
  fit$covariance<-NULL
  missing<-assess_taxonomy(fit,taxonomy,"area")
  expect_false(missing$accepted)
  expect_true(is.na(ssf_contrast(missing$context,"cos_angle","Encamped","Travelling")$lower))
})

test_that("common proposal constants leave conditional choices unchanged", {
  d<-prepare_ssf_data(data.frame(ID=1,stratum=rep(1:4,each=3),obs=rep(c(1,0,0),4),
    step=seq(.1,1.2,length.out=12),cos_angle=cos(1:12),w=1),density_measure="area")
  b<-rbind(step=c(-2,-1),log_step=c(0,.2),cos_angle=c(.2,1))
  f<-~step+log_step+cos_angle
  a<-ssf_log_emissions(d,f,b)
  probability<-rep(c(.8,.9,.95,.99),each=3)
  d$w<-d$w/probability
  c<-ssf_log_emissions(d,f,b)
  expect_equal(c,a-log(c(.8,.9,.95,.99)),tolerance=1e-12)
  expect_equal(c-log(2*d$w[d$obs==1]),a-log(2),tolerance=1e-12)
})

test_that("finite radial support reaches the fitted object and its summary", {
  d <- prepare_ssf_data(data.frame(ID=1,stratum=rep(1:10,each=2),
    obs=rep(1:0,10),step=seq(.1,2,length.out=20),cos_angle=cos(1:20),w=1),
    density_measure="area",step_max=3.8)
  f <- fit_ssf(d,n_states=1L,n_starts=1L,seed=1L)
  expect_equal(f$step_max,3.8)
  expect_equal(summary(f)$states$mean_step,
    movement_metrics(coef(f),"area",step_max=3.8)$mean_step)
})
