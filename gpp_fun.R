gpp_fun = function(beta0,alpha,gamma_1,par_0,par_z,temp_0,temp_z,temp_ref){
  beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_0 + par_z))
}

gpp_grad_fun <- as.formula(paste0("~", capture.output(gpp_fun)[2])) %>%
  deriv(c("beta0","alpha","gamma_1","par_0","par_z","temp_0","temp_z","temp_ref"), 
  function(beta0,alpha,gamma_1,par_0,par_z,temp_0,temp_z,temp_ref){})

test <- as.formula(paste0("~", capture.output(gpp_fun)[2])) %>%
  deriv(c("beta0","alpha","gamma_1","par_0","par_z","temp_0","temp_z","temp_ref"),
        function.arg = T)

gpp_grad_fun(beta0=1,alpha=1,gamma_1=1,par_0=1,par_z=1,temp_0=1,temp_z=1,temp_ref=1)
gpp_grad_fun(1,1,1,1,1,1,1,1)



gpp_fun = function(beta0,alpha,gamma_1,par,temp,temp_ref){
  beta0*gamma_1^(temp-temp_ref)*tanh((alpha/(beta0*gamma_1^(temp-temp_ref)))*par)
}

t = parse(text = 
"f = function(c,d){
return(c+d)
}")

f2 = function(b,c){eval(t)}

pars = expression("c, d")
f = function(pars){
  c+d
  }
f(eval(2,3))

f2 = function(c,d){eval(parse(text=capture.output(f)[2]))}
f2(2,3)
