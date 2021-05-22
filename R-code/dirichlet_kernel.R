### Dirichlet's kernel as per Equation (17c) of SAUTS

dirichlet_kernel <- function(freq,N) if(isTRUE(all.equal(freq,0))) (-1)^((N-1)*freq) else sin(N*pi*freq)/(N*sin(pi*freq))

