function fnk, n, k
  ; Computes function Fn(K), undulator brightness of odd harmonic N for
  ; deflection paramter K
  term1 = k^2*n^2/(1.+K^2/2.)^2
  term2 = n*k^2/(4.*(1+k^2/2.))
  return, term1*(beselj(term2, (n-1)/2.) - beselj(term2, (n+1)/2.))^2
end
