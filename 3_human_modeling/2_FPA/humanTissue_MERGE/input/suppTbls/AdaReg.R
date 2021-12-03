#' @title AdaReg
#' @description To robustly estimate linear regression coefficients and detect outliers.
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param X The design matrix in linear regression.
#' @param y The response variable in linear regression.
#' @param gam.seq A sequence of gamma's (default: seq(0, 3, by=.1)).
#' @param var.gp.id A vector of groups indicating the samples in the same variance group (default: all the samples have the same variance).
#' @param tol Tolerance for interations (default: 10^(-4)).
#' @param step Step limit (default: 50).
#' @return Estimated coef and variance, and fitting info for the residuals.
#' @export

AdaReg = function(X, y, gam.seq = seq(0, 3, by=.1), var.gp.id=NULL, mu.fix.value=NULL, var.fix.value=NULL, tol=10^(-4), step=50) {
  
  # initialization
  fit = lm(y ~ X - 1)
  beta.0 = coef(fit)
  res.0 = y - X %*% beta.0  
  if(is.null(var.gp.id)) {
    var.gp.id=rep(1, length(y))
  }
  var.gp.nm = sort(unique(var.gp.id))
  var.gp.len = table(var.gp.id)[var.gp.nm]
  var.sig.gp.0 = numeric(length(var.gp.nm))
  if (length(var.gp.nm) == 1) {
    var.sig.gp.0 = var(res.0)
  } else {
    for (k in 1:length(var.gp.nm)) {
      var.sig.gp.0[k] = var(res.0[var.gp.id == var.gp.nm[k]])
    }	
  }
  #var.sig.0 = var(res.0)
  
  # iterations
  int = 1
  flag = FALSE
  mu.res.int = c()
  sd.res.int = c()
  beta.int = beta.0
  var.sig.gp.int = var.sig.gp.0
  diff.par.int = c()
  while (flag == FALSE) {
    res.rob.fit.0 = adapt.gam.rob.fit.fn(res.0, gam.seq, step=step, mu.fix=mu.fix.value, var.fix=var.fix.value, bin.num=NULL)
    if (is.null(res.rob.fit.0) | sum(is.na(res.rob.fit.0$est.hat)) > 0) {
      flag = TRUE
      break
    } else {
      w = res.rob.fit.0$x.w
      W.0 = diag(c(w/rep(var.sig.gp.0, var.gp.len)))
      A = t(X) %*% W.0 %*% X
      if (rcond(A) < 10^(-4)) {
        res.rob.fit.0 = NULL
        diff.par = NA
        diff.par.int = c( diff.par.int, diff.par)
        break;
      } else {
        beta.1 = solve(t(X) %*% W.0 %*% X) %*% (t(X) %*% W.0 %*% y)
        
        var.sig.gp.1 = numeric(length(var.gp.nm))
        for (k in 1:length(var.gp.nm)) {
          var.sig.gp.1[k] = (1+res.rob.fit.0$est.hat["gam.sel"])*(sum(w[var.gp.id == var.gp.nm[k]] * (y[var.gp.id == var.gp.nm[k]] - X[var.gp.id == var.gp.nm[k],] %*% beta.1)^2))/sum(w[var.gp.id == var.gp.nm[k]])
        }
        res.1 = y - X %*% beta.1
        
        if (sum(is.na(var.sig.gp.1)) > 0) {
          res.rob.fit.0 = NULL
          diff.par = NA
          diff.par.int = c( diff.par.int, diff.par)
          break;
        } else {
          if (min(var.sig.gp.1) < 10^(-4) ) {
            res.rob.fit.0 = NULL
            diff.par = NA
            diff.par.int = c( diff.par.int, diff.par)
            break;
          }        	
        }
        
        
        
        diff.par = sum( sum(abs(beta.0 - beta.1)) + abs(abs(var.sig.gp.0- var.sig.gp.1)))
        diff.par.int = c(diff.par.int, diff.par)
        
        if ( diff.par < tol | int > step)  {  
          flag = TRUE
          break
        }else {
          res.0 = res.1
          beta.0 = beta.1
          var.sig.gp.0 = var.sig.gp.1
          mu.res.int = c(mu.res.int, res.rob.fit.0$est.hat["mu0.hat"])
          sd.res.int = c(sd.res.int, res.rob.fit.0$est.hat["sd0.hat"])
          beta.int = cbind(beta.int, beta.0)
          var.sig.gp.int = cbind(var.sig.gp.int, var.sig.gp.0)
          int = int + 1
        }
      }
    }
    
    
  }
  
  # diagnosis
  if(!is.null(res.rob.fit.0) & sum(is.na(res.rob.fit.0$est.hat)) ==	 0) {
    if (abs(res.rob.fit.0$est.hat["mu0.hat"]) >= 1) {
      print("need further check on the fitting")
    }    	
  }
  
  
  return(list(beta.rob.fit = c(beta.0), var.sig.gp.fit=c(var.sig.gp.0), x.res=c(res.0), res.info = res.rob.fit.0$est.hat, x.w=  c(res.rob.fit.0$x.w), mu.res.int = mu.res.int, sd.res.int=sd.res.int, beta.int=beta.int, var.sig.gp.int=var.sig.gp.int, diff.par.int=diff.par.int ))
  
}


# ------------------------------------------------------------------------------
# data-adaptive selection procedure
adapt.gam.rob.fit.fn = function (x.00, gam.seq, step=50, mu.fix=NULL, var.fix=NULL, bin.num=NULL) {
  
  x.0 = x.00[!is.na(x.00)]
  nm = c("mu0.hat", "sd0.hat", "pi0.hat", 'efdr0.hat')
  par.hat = matrix(NA, length(gam.seq), length(nm))
  rownames(par.hat) = gam.seq
  colnames(par.hat) = nm
  
  for (i in 1:length(gam.seq)) {
    gam = gam.seq[i]
    x.mu = ifelse(is.null(mu.fix), mean(x.0), mu.fix)
    x.var = ifelse(is.null(var.fix), var(x.0), var.fix)
    result = est.fn(x.0, x.mu, x.var, gam, fix.mu=!is.null(mu.fix), fix.var=!is.null(var.fix), step=step)		
    mu.hat = result$mu.est
    var.hat = result$var.est
    if (!is.na(var.hat)) {
      est.result = efdr.0.fn(x.0, mu.hat, var.hat, gam, bin.num)
      par.hat[i, ] = est.result[colnames(par.hat)]
    }
    
  }
  crt.hat.0 = abs(pmin(par.hat[,'efdr0.hat'], 10) - 1)
  ind.comp = !is.na(crt.hat.0)
  gam.comp = gam.seq[ind.comp]
  if (length(gam.comp)  == 0 ) {
    est.hat=NA
    x.n=NA
    w=NA
    gam.comp=NA
    crt.hat.0=NA
    para.hat.mx=NA
  } else {
    crt.hat.0 = crt.hat.0[ind.comp]
    par.hat = matrix(par.hat[ind.comp, ], ncol=length(nm))
    colnames(par.hat) = nm
    rownames(par.hat) = gam.comp
    
    crt.est.0 = pmin(par.hat[,'efdr0.hat'], 10)
    gam.sel = gam.comp[which.min(crt.hat.0)]
    
    gam.sel.char = as.character(gam.sel)
    est.hat = c(length(x.0), gam.sel, par.hat[gam.sel.char, c("mu0.hat", "sd0.hat", "pi0.hat", "efdr0.hat")])
    names(est.hat)[1:2] = c("n.observed", "gam.sel")
    
    w.nu = dnorm(x.00, est.hat["mu0.hat"], est.hat["sd0.hat"])^est.hat["gam.sel"]
    w.nu[is.na(x.00)] = NA
    w = w.nu/sum(w.nu, na.rm=TRUE)
  }
  
  return(list(est.hat=est.hat, x.w=w, gam.comp=gam.comp, crt.hat.0=crt.hat.0,  para.hat.mx=par.hat )) 
  
}

# ------------------------------------------------------------------------------
# expected of fdr criterion
efdr.0.fn = function (x, mu.hat, var.hat, gam, bin.num=NULL) {
  x = x[!is.na(x)]
  den.fit = dnorm(x, mu.hat, sqrt(var.hat))
  frac.hat= mean(den.fit^gam)*sqrt(2*pi*var.hat)^gam * sqrt(1 + gam)
  
  my.hist = bk.cnt.fn(x, bin.num)
  bin.bk = my.hist[[1]]
  cnt.bk = my.hist[[2]]
  
  p0.hat.bin  = numeric(length(bin.bk)-1)
  p0.hat.bin[1] = pnorm(bin.bk[2], mu.hat, sqrt(var.hat))
  p0.hat.bin[length(bin.bk)-1] = 1 - pnorm(bin.bk[length(bin.bk)-1], mu.hat, sqrt(var.hat))
  for ( j in 3:(length(bin.bk)-1)) {
    p0.hat.bin[j-1] = pnorm(bin.bk[j], mu.hat, sqrt(var.hat)) - pnorm(bin.bk[j-1], mu.hat, sqrt(var.hat))
  }
  p.hat.bin =  cnt.bk/sum(cnt.bk)
  null.efdr.hat = min(1,frac.hat) * sum(p0.hat.bin^2/p.hat.bin)	
  est.sum = c(gam, mu.hat, sqrt(var.hat),  frac.hat, null.efdr.hat)
  names(est.sum) = c("gamma", "mu0.hat", "sd0.hat", "pi0.hat",  'efdr0.hat')
  return(est.sum)
}



# ------------------------------------------------------------------------------
# estimation under a fixed gamma
est.fn = function(x, mu.0, var.0, gam,  tol=10^(-4), step=step, fix.mu=FALSE, fix.var=FALSE) {
  x = x[!is.na(x)]
  n = length(x)
  dum = dnorm(x, mu.0, sqrt(var.0))^gam
  w.0 = dum/sum(dum)
  
  int = 1
  flag = FALSE
  diff.par.int = c()
  mu.int = mu.0
  var.int = var.0
  while ( flag == FALSE) {
    mu.1 = ifelse(fix.mu, mu.0, sum(w.0 * x))
    var.1 =  ifelse(fix.var, var.0, (1+gam) * sum(w.0 * (x-mu.1)^2)) 
    if (var.1 < 10^(-4)) {
      mu.0 = NA
      var.0 = NA
      diff.par = NA
      diff.par.int = c( diff.par.int, diff.par)
      break;
    }
    
    diff.par = abs(mu.1 - mu.0) + abs(sqrt(var.1) - sqrt(var.0))
    diff.par.int = c( diff.par.int, diff.par)
    if (diff.par < tol | int > step ){  
      flag = TRUE
      break;
    } else {
      mu.0 = mu.1
      var.0 = var.1
      dum = dnorm(x, mu.0, sqrt(var.0))^gam
      w.0 = dum/sum(dum)
      mu.int = c(mu.int, mu.0)
      var.int = c(var.int, var.0)
      int = int + 1
    }
  }
  return(list(mu.est=mu.0, var.est=var.0, w=w.0, diff.par.est=diff.par.int))
}


# ------------------------------------------------------------------------------
# to merge intervals s.t. each bin has positive number of data points
bk.cnt.fn = function (x, bin.num=NULL) {
  if (is.null(bin.num)) {
    if (length(x) > 1000) {
      bk.num = 20
    } 
    if (length(x) <= 1000 & length(x) > 500) {
      bk.num = 10
    } 
    if (length(x) <= 500){
      bk.num= 5
    }
  } else {
    bk.num = bin.num
  }
  h = hist(x, breaks=bk.num, plot=FALSE)
  # h$counts
  # h$breaks
  ind.zero = (1:length(h$counts))[h$counts == 0]
  if (length(ind.zero) != 0) {
    bk.start = h$breaks[ind.zero]
    bk.end = h$breaks[ind.zero + 2]
    bk.update = c(bk.start[1])
    if (length(bk.start) > 1) {
      for (i in 1: (length(bk.start)-1) ) {
        if (bk.end[i] < bk.start[i+1]) bk.update = c(bk.update, bk.end[i])					      
      }
    }
    bk.update = c(bk.update,  bk.end[length(bk.end)] )
    bk.pts = sort(unique(c(h$breaks[-c(ind.zero, ind.zero+1)], bk.update)))
    cnt = unname(table(cut(x, breaks=bk.pts)))
  } else {
    bk.pts = h$breaks
    cnt = h$counts
  }
  return( list( bk.pts, cnt))
}
