*gmm code for gmm_fem_bach.dta (returning for Bachelor's Degree) (first dataset - ch1par_female.dta)
/*keep tage sqage marry efnp faminc metro rfownkid firmsize rfoklt18 indcode
*gen constant=1
rename rfoklt18 child18
tab indcode, gen(ind)
bigfirm=firmsize==3
log using gmm_fem_bach, replace
*/
tomata
mata:
		void i_mom_bd_male(todo,init1,crit,g,H)
		{
		external X,Z,y1,y2
		xc=cols(X)
		zc=cols(Z)
		suc=xc+zc
		sxz=suc+zc
		b=init1[1,1..xc]
		c=init1[1,(xc+1)..(suc+1)]
		c1=init1[1,(suc+2)..(sxz+2)]
		mu1=X*b'
		phi=normalden(mu1) 
		Phi=normal(mu1) 
		mils=phi:/Phi 
		imils=phi:/(1:-Phi) 
		m3=(y1:*mils)-((1:-y1):*imils)   
		m4=(X'*m3):/rows(X)
		q=Z,mils                         
		m5=(q'*(y1:*(y2:-q*c'))):/rows(Z)
		q1=Z,imils
		m7=(q1'*((1:-y1):*(y2:-q1*c1'))):/rows(Z)
		m6=m4\m5\m7                         
		crit=(m6'*m6)                     
		}
		

y1=trbd
y2=lwr
X=tage, marry, efnp, child18, faminc, constant
Z=tage, sqage, marry, metro, bigfirm, constant
xc=cols(X)
zc=cols(Z)
suc=xc+zc
sxz=suc+zc
b=J(1,xc,0)
c=J(1,(zc+1),0)
c1=J(1,(zc+1),0)
init1=b,c,c1
S=optimize_init()
optimize_init_evaluator(S, &i_mom_bd_male())
optimize_init_which(S,"min")
optimize_init_evaluatortype(S,"d0")
optimize_init_conv_maxiter(S, 90000)
optimize_init_params(S,init1)
p=optimize(S)
p

*step 2. Getting the W
init1hat=p
bhat=p[1,1..xc]
chat=p[1,(xc+1)..(suc+1)]
chat1=p[1,(suc+2)..(sxz+2)]
mu1hat=X*bhat'
phihat=normalden(mu1hat) 
Phihat=normal(mu1hat) 
milshat=phihat:/Phihat
imilshat=phihat:/(1:-Phihat)
m3hat=(y1:*milshat)-((1:-y1):*imilshat)
m3star=(X:*m3hat)
q12=Z,milshat
qstar=(q12:*(y1:*(y2:-q12*chat')))
q13=Z,imilshat
qthree=(q13:*((1:-y1):*(y2:-q13*chat1')))
bigm=m3star,qstar,qthree
bigw=bigm'bigm
invbigw=invsym(bigw):/rows(X)

*step3: Weighted GMM
void i_momstep5(todo,init1,crit1,g,H)
		{
		external X,Z,y1,y2,invbigw
		xc=cols(X)
		zc=cols(Z)
		suc=xc+zc
		sxz=suc+zc
		b=init1[1,1..xc]
		c=init1[1,(xc+1)..(suc+1)]
		c1=init1[1,(suc+2)..(sxz+2)]
		mu1=X*b'
		phi=normalden(mu1) 
		Phi=normal(mu1) 
		mils=phi:/Phi 
		imils=phi:/(1:-Phi) 
		m3=(y1:*mils)-((1:-y1):*imils)   
		m4=(X'*m3):/rows(X)
		q=Z,mils                         
		m5=(q'*(y1:*(y2:-q*c'))):/rows(Z)
		q1=Z,imils
		m7=(q1'*((1:-y1):*(y2:-q1*c1'))):/rows(Z)
		m6=m4\m5\m7                         
		crit1=(m6'*invbigw*m6)                     
		}   
		
beta=J(1,xc,0)
ceta=J(1,(zc+1),0)
zeta=J(1,(zc+1),0)
init1=beta,ceta,zeta
fv=optimize_init()
optimize_init_evaluator(fv,&i_momstep5())
optimize_init_which(fv,"min")
optimize_init_evaluatortype(fv,"d0")
optimize_init_conv_maxiter(fv, 90000)
optimize_init_params(fv,init1)
p3step=optimize(fv)
p3step


* gradient
*variance = 1/n*inverse((T'*Psi*T))
  *First calculating Psi in the Final Var-covar matrix
init1hat2=p3step
bhat1=p3step[1,1..xc]
chat12=p3step[1,(xc+1)..(suc+1)]
chat13=p3step[1,(suc+2)..(sxz+2)]
mu1hat1=X*bhat1'
phihat1=normalden(mu1hat1) 
Phihat1=normal(mu1hat1) 
milshat1=phihat1:/Phihat1
imilshat1=phihat1:/(1:-Phihat1)
m3hat1=(y1:*milshat1)-((1:-y1):*imilshat1)
m3star1=(X:*m3hat1)
q3=Z,milshat1
qstar1=(q3:*(y1:*(y2:-q3*chat12')))
q32=Z,imilshat1
qstar2=(q32:*((1:-y1):*(y2:-q32*chat13')))
bigm1=m3star1,qstar1,qstar2
bigw1=(bigm1'bigm1):/rows(X)
invbigw1=invsym(bigw1)


*Calculating change in moments to change in paramter value (T above)
*First val

/*m4fir=(X:*m3hat1)
qfir=Z,milshat1                         
m5fir=(qfir:*(y1:*(y2:-qfir*chat1')))
qfir2=Z,imilshat1
m5fir2=(qfir2:*((1:-y1):*(y2:-qfir2*chat12')))
m6fir=m4fir,m5fir,m5fir2
*/
* Approach Three
diff=0.000001
mom=J(cols(p),cols(p),0)
ones=J(1,rows(X),1)
for(i=1;i<=cols(p3step);i++)
{
init13=p3step
init13[i]=p3step[i]*(1+diff)
est1=init13[i]-p3step[i]
b1=init13[1,1..xc]
c1=init13[1,(xc+1)..(suc+1)]
c2=init13[1,(suc+2)..(sxz+2)]
mu11d=X*b1'
phi1d=normalden(mu11d) 
Phi1d=normal(mu11d) 
mils1d=phi1d:/Phi1d
imils1=phi1d:/(1:-Phi1d)
momd1=(y1:*mils1d)-((1:-y1):*imils1)   
momd2=(X:*momd1)
qd1=Z,mils1d                         
momd3=qd1:*(y1:*(y2:-qd1*c1'))
qd2=Z,imils1
mumd3=qd2:*((1:-y1):*(y2:-qd2*c2'))
momd4=momd2,momd3,mumd3
momd5=(((ones*momd4):/rows(X))-((ones*bigm1):/rows(X))):/est1
(ones*momd4):/rows(X)-(ones*bigm1):/rows(X)
mom[.,i]=momd5'
}

cormat=invsym(mom'*invbigw1*mom):/rows(X)
cormat
cormat2=invsym(mom*invbigw1*mom'):/rows(X)
sebdfemale=diagonal(cormat)
sebdfemale2=diagonal(cormat2)
*** End of Code
lwret1=Z[.,1..cols(Z)]*p[1,(xc+1)..(suc)]'
lwnret1=Z[.,1..cols(Z)]*p[1,(suc+2)..(sxz+1)]'
lwdiff1=lwret1-lwnret1

lwret=Z[.,1..cols(Z)]*p3step[1,(xc+1)..(suc)]'
lwnret=Z[.,1..cols(Z)]*p3step[1,(suc+2)..(sxz+1)]'
lwdiff=lwret-lwnret
	
** store results
st_matrix("p_st",p)
st_matrix("p2_st",p3step)
st_matrix("cov_st",cormat2)
st_matrix("se_st",sebdfemale2)
st_matrix("cov_st_second",cormat)
st_matrix("se_st_second",sebdfemale)
st_addvar("double","lwret")
st_store(.,"lwret",lwret)
st_addvar("double","lwnret")
st_store(.,"lwnret",lwnret)
st_addvar("double","lwdiffp3")
st_store(.,"lwdiffp3",lwdiff)
st_addvar("double","lwdiffp")
st_store(.,"lwdiffp",lwdiff1)


mata matsave gmm_bd_female *  //saves all mata workspace as gmm_ad_male.mmat
mata matuse gmm_bd_female // retrieves for future work
end
save gmm_male_bach_1, replace
