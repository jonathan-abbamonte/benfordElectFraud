Benford_calc<-function(theta, x, alpha){

  #num=prod(theta^x)*prod(gamma(alpha)*gamma(sum(alpha+x)))
  #print(prod(theta^x))
  #print(prod(gamma(alpha))*gamma(sum(alpha+x)))
  #denom=gamma(sum(alpha))*prod(gamma(alpha+x))

  prenum_fac1 = sum(x*log(theta))
  prenum_fac2 = sum(lgamma(alpha))
  prenum_fac3 = lgamma(sum(alpha+x))
  predenom = lgamma(sum(alpha))+sum(lgamma((alpha+x)))

  # print(prenum_fac1)
  # print(prenum_fac2)
  # print(prenum_fac3)
  # print(predenom)
  BF = (prenum_fac1+prenum_fac2+prenum_fac3-predenom)
  #print(num)
  #print(denom)
  #BF=num/denom
  return(BF=BF)
}
