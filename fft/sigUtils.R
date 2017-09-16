prep.signals <- function(data.sets, threshold=0.01, dir='./', verbose=0) {
  # massage original data and QA
  z <- "time..sec.,X.vibration..m.s.2.,Y.vibration..m.s.2.,Z.vibration..m.s.2.,Start.Time,Start.Date,Duration..sec.,Data.Rate..Hz."
  expected <- strsplit(z, ',')[[1]]
  lst <- list()
  err <- list()
  for (dn in data.sets) {
    fn <- paste(dir, sprintf('%s.csv', dn), sep='/')
    if (verbose > 0) {
       cat('INFO: prep.signals: dn=', dn, 'fn=', fn, '\n')
    }
    x <- read.csv(fn)
    if (!all(names(x) == expected)) {
      if (verbose > 0) {
        cat('ERROR: expected=', expected, '\n')
        cat('ERROR: dn, names(x)=', dn, names(x), '\n')
      }
      err[length(err) + 1] <- paste('ERROR', dn, 'unexpected column names', sep=': ')
    }
    x <- x[,1:4]
    names(x) <- c('t','x','y','z')
    dt <- diff(x$t)
    v <- sd(dt)/mean(dt)
    v <- max(abs(dt-mean(dt)))/mean(dt)
    if (v > threshold) {
      msg <- sprintf('ERROR:%s: delta t variance=%g above the threshold=%g', dn, v, threshold)
      err[length(err) + 1] <- msg
    }
    x$norm <- sqrt(apply(x[,2:4]^2, 1, sum))
    lst[[dn]] <- x
  }
  return (list(err=err, data=lst))
}

plt.signals <- function(signals) {
  for (name in names(signals)) {
    t <- signals[[name]]$t
    v <- signals[[name]]$norm
    plot(t, v, type='l', main=name)
  }
}

plt.cdfs <- function(signals, xlab='amplitude', value='norm', titl='signals') {
  cdfs <- list()
  mx <- -Inf
  for  (name in names(signals)) mx <- max(mx, signals[[name]][[value]])
  plot(1, 1, xlim=c(0, mx), ylim=c(0,1), type='n', xlab=xlab, ylab='CDF', main=titl)
  cols <- rainbow(length(signals))
  names(cols) <- names(signals)
  lg <- NULL
  for (name in names(signals)) {
    x <- sort(signals[[name]][[value]])
    lg <- c(lg, sprintf('%s: median=%s', name, signif(median(x),2)))
    y <- (1:length(x))/length(x)
    lines(x, y, type='l', col=cols[name])
  }
  legend('bottomright', lg, col=cols, lwd=1)
}

get.fft <- function(y,             # signal
                    t,             # sampling times
                    tag=NULL,
                    return.half=TRUE,
                    mean.rm=TRUE) # keep/remove mean of the signal
{
  dt <- t[2] - t[1]
  df <- 1/(max(t) - min(t))
  fnyq <- 1/(2*dt)           # Nyquist frequency
  f <- seq(0, 2*fnyq, by=df) # frequencies for FFT output (in Hz)
  if (length(f) < length(t)) {
    cat(sprintf('WARNING: length(f)[=%s] < length(t)[=%s]\n', length(f), length(t)))
    # if (abs(length(f) - length(t)))
    mx <- max(f)
    f[length(f) + 1] <- mx + df
  }
  cat(sprintf('get.fft: %s: dt=%s, fnyq=%s, max(f)=%s\n', tag,
      signif(dt,2), signif(fnyq,3), signif(max(f))))
  if (mean.rm) y <- y - mean(y)
  ft <- 2 * fft(y)/length(y) # normalized fourier transform
  pw <- abs(ft) # power spectrum
  if (return.half) {
     i <- 1:as.integer(length(f)/2)
     f <- f[i]
     pw <- pw[i]
  }
  return (list(f=f, ft=ft, pw=pw))
}

plt.power <- function(y,        # signal
                      t,        # sampling times
                      pw,       # power spectrum
                      f,        # frequencies
                      tag=NULL,
                      show='',  # what to display
                      eps=1e-6, # small number to ignore floating-point errors
                      ext=1.1)  # how much to extend the plot beyond the last maximum
  {
    nt <- unique(c(length(y), length(t)))
    nf <- unique(c(length(f), length(pw)))
    if (length(nt) != 1 | length(nf) != 1) {
       cat('plt.power: length(y), length(t)=', length(y), length(t), '\n')
       cat('plt.power: length(pw), length(f)=', length(pw), length(f), '\n')
       msg <- sprintf('%s: y,t and pw,f lengths should be respectively equal', tag)
       stop(msg)
    }
    if (show == 'all') {
      mx <- nf
    } else {
      z <- diff(pw, differences=2)
      w <- which(z < -eps)
      mx <- min(max(w) * ext, nf) # index for the highest frequency maximum
    }
    idx <- 1:mx
    par(mfrow=c(1,2))
    plot(t, y, type='l', xlab='time [s]', ylab='signal value',
         main=sprintf('%s: original signal', tag))
    plot(f[idx], pw[idx], type='l', xlab='freq [Hz]', ylab='power',
         main=sprintf('%s: FFT output', tag))
    par(mfrow=c(1,1))
}
