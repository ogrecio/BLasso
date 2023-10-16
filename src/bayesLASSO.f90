module variables
 implicit none

! Número de ANIMALES, DATOS Y EFECTOS
   integer (kind=4):: n_animal,n_datos,nf1,maxrec,maxrec_c,n_efectos,n_rand,n_cov,ntrait,n_cat,n_snp
   character (len=14), allocatable::id(:)
   character (len=14)::phen_id

!********************************
!**********PARAMETROS ***********
!********************************
     PARAMETER (ntrait=1)        !Numero de caracteres
 !   PARAMETER (maxrec=1370000) !Número de elementos NO nulos
!********************************

 !  common /genea/ padre,madre,cod
 !  common /gen1/ d,y,valor
   integer, allocatable:: padre (:), madre (:),codcontrol(:),ani(:)
   integer, allocatable:: cen(:) !maxan x número de efectos
   real*8, allocatable:: valor(:,:),y(:),y1(:),y2(:)
real*8, allocatable:: t(:),ln_t(:),t_c(:),snp(:),snp2(:)
   real*8, allocatable:: d(:) !diagonal de A
   integer, allocatable:: nf(:),nf_rand(:),vdf(:) !número de niveles de ef.fijos y ef.aleatorios

! ECUACIONES
!   common /matriz / irow,icol,xval
   integer:: ii,neq !número de ecuaciones del MME
!    parameter (neq=nf1+nf2+nf3+n_reb+n_ani_dat+n_animal) !Número máximo de ecuaciones
!   common /equs/ sol1,sol2,rhs
   real*8, allocatable:: error (:),solf(:),vep1(:),vep2(:),rhs (:),sol(:),sol_efectos(:)
   real*8, allocatable:: gebv(:),xpx(:)
   real*8:: mean,muf


!GIBBS SAMPLING
   integer :: umbral,sv,niter,iter,nciclos,n_ciclos,burnin,lag,lag_t,nm
    !Si es un caracter umbral es 1. Diferente de 1 para caracteres continuos.
	!sv(=1para metropolis hasting;=2para wishart invertida)
   real*8::alfa,h,r,y_mean,denomh,alfa_old,llhood,lambda,scale
   real*8::mu,ve,ve_c,vg,vg_c,sd_st,apriorig(1,2),apriorie(1,2)
   real*8,allocatable::va(:),va_c(:),inivar(:),apriorir(:,:),vs(:),aprioris(:,:),inv_tau2(:),tau_F(:)

!OTROS (contadores,semillas, formatos, etc)
	integer :: io,contador,idum,horas,minutos,i,j,posj,k,q,df,n_df
	character (len=8) ::fecha
	character (len=6) ::name
	character (len=60) ::fdatos,ftest,fpedigri,fpar
	real*8 :: tiempo1,tiempo2,tiempo
	real*8::den1,den0,ratio,sc_S,sc_ep,sc_p,sc_e,sc_r,sc_t0,sc_t1,temp
	real*8::VPS,VGS,VEpiS,VES,VRS,VPf,VGf,VEf,VRf,hf,hs,rf,rs
    real*8::x2,x1,yy,cat_df
    real*8::te_ve,mcerror_ve,mu_phen,var_phen
    real*8::kk1,kk2,kk3,kk4,kk5,kk6
    integer:: tot,acc,cont


 TYPE typephen
    character (len=1):: tipoprueba
    real*8:: fiab
    real*8::dyd
    real*8::de
 END TYPE typephen

type (typephen):: phen(26)

end module variables

MODULE inicial
!*****************************************************************************
!*En este archivo se marcan los valores iniciales, la longitud de la cadena, *
!*el periodo de quemado, el 'lag', el algoritmo para muestrear las varianzas,*
!* y se determina que caracteres son umbrales                                *
!*****************************************************************************
USE variables
CONTAINS

SUBROUTINE inic
integer :: pf

ii=1          !Número de cadenas que se lanzan
niter=1       !Numero de iteraciones por cada vez que se resuelve por Gauss-Seidel
pf=1 !read *,pf
  print *,'Enter name of PARAMETER file'
fpar='param.bl' ! read *,fpar
 open (22,file=fpar, form='formatted', status='old')
!***INTRODUZCA EL NOMBRE DE LOS FICHEROS**********
!Data file:
read (22,*)   !Saltar linea
read (22,*) fdatos
print *,'DATA FILE=', fdatos
 open (10,file=fdatos, form='formatted', status='old')
!testing file:
read (22,*) ftest
print *,'Testing FILE=', ftest
 open (11,file=ftest, form='formatted', status='old')
!****DETERMINACION DEL NUMERO DE CARACTERES, EFECTOS Y NIVELES******
read (22,*)   !Saltar linea
read (22,*) n_efectos
print *,'# EFFECTS=', n_efectos
!read (22,*) n_cov
n_cov=n_efectos
print *,'# COVARIATES=', n_cov
n_rand=0 !This version does not allowed effects with a covariance structure
ALLOCATE (nf(n_efectos),nf_rand(n_rand-1),va(n_rand-1),va_c(n_rand-1),inivar(n_rand),&
          apriorir(n_rand-1,2),aprioris(n_cov,2),inv_tau2(n_cov),tau_F(n_cov))
tau=0.01d0
tau_F=0.d0
read (22,*)   !Saltar linea
nf(1)=0
DO k=1,n_efectos-1
  IF (k.le.n_cov) THEN
   nf1=1
  ELSE
    read (22,*) nf1
	print *,'levels for effect ',k,'=', nf1
	 IF (k.ge.n_efectos-n_rand+1) THEN
	   q=q+1
	   nf_rand(q)=nf1
	 END IF
  END IF
  nf(k+1)=nf(k)+nf1
END DO

!****INFORMACION A PRIORI********
read (22,*)   !Saltar linea
do i=1,n_cov
  aprioris(i,1)=0.000001d0
  aprioris(i,2)=1.d0
end do

   read (22,*) apriorie(1,1)
   if (apriorie(1,1).eq.0) then
      apriorie(1,2)=0
   else
      read (22,*) apriorie(1,2)
   end if
   print *,'residual apriori variance=',apriorie(1,1)
   print *,'df=',apriorie(1,2)

!****DETERMINACIÓN Y NOMBRAMIENTO DE LOS CARACTERES********
read (22,*)   !Saltar linea
!name
   read (22,*) name
!umbral/continuo
   read (22,*) umbral
   if (umbral.eq.1) then
      read (22,*) n_cat !#categories
	  print *,name,' is a threshold trait with ',n_cat,' categories'
   else
	  print *,name,' is a linear trait'
   end if
!missing values are not allowed in this version
   !read (22,*) df
   df=0
   if (df.eq.1) then
      read (22,*) cat_df !identify value
	  print *,cat_df,' will be consider as a missing value in data file'
   end if


!***PARAMETROS PARA IMPLEMENTAR GIBBS**********
read (22,*)   !Saltar linea
 read (22,*) n_ciclos
 print *, 'Lenght of chain = ',n_ciclos
 read (22,*) burnin
 print *, 'Burn-in = ',burnin
read (22,*) lag
 print *,'Lag = ',lag
 sv=2 !This version is optimized for sampling variances with an inverted Wishart prior
 if (sv.eq.2) then
   print *,'Variances will be sampled with an inversed-Wishart'
 else if (sv.eq.1) then
   print *,'Variances will be sampled with Metropolis-Hasting Algorithm'
 else
   print *,'error in parameter file at sampling method'
   stop
 end if
!VARIANZAS INICIALES
read (22,*) ve  !residual variance
read (22,*) lambda  !lambda initial value
 vg=1
print *,'VARIANCES USED IN FIRST ITERATION'
print *,'residual=',ve

goto 28

28 continue
!PREPARACION DE SALIDA (no modificar)
 308 format (a6)
  write (34,308) name
  write (31,308) name
  write (34,*) '  LAMBDA      MU        SCALE'
  write (31,*) ' VE       '

close (22)

END SUBROUTINE inic
END MODULE inicial

module wish
   use variables
contains

!algoritmo WISHART INVERTIDA
subroutine wishart_inv
integer:: xx(n_efectos) !nf=numero de efectos fijos + aleatorios + covariables
real*8::rate,temp,inv_gauss


!Sample lambda parameter from a gamma distribution
rate=0.d0
do j=1,n_cov
    rate=rate+0.5d0/inv_tau2(j)
enddo
call gamma2(n_cov+1.d0,rate+1.7d0,x2,lambda)
lambda=sqrt(lambda)

write (34,*) lambda,mu,scale

!   call scS
!   call wishart(ntrait,sc_S+170.d0*5,(n_datos)+5.d0,vg)

!Sample Residual variance
sc_S=0.d0
do i=1,n_datos
    sc_S=sc_S+error(i)*error(i)
enddo
call wishart(ntrait,sc_S+1.d0,(n_datos)+3.d0,ve)

sc_S=0.d0
    do i=1,n_datos
       sc_S=sc_S+0.5*(error(i)*error(i) +sol(j)*sol(j)*(inv_tau2(j)))
    enddo
call gamma(0.5d0*(n_datos-1)+0.5d0*n_cov,sc_S,x1,ve)
ve = 1/ve


end subroutine wishart_inv
!____NON PARAMETRIC SUM OF SQUARES for additive genetic______

subroutine scS
  integer::i,j
  real*8::tt(n_datos)

  sc_S=0.d0

do i=1,n_datos
  tt(i)=dot_product(sol(:),valor(:,i))
enddo

sc_S=dot_product(tt(:),sol(:))

return
end subroutine scS
! ********************************************************************************
  subroutine wishart(ntrait,sc,ne,V)
    character (len=6) :: form*6
    integer :: i, j,ntrait
    real*8 :: ne,alfa, beta,v,sc,u
    real*8 :: se(ntrait,ntrait),t(ntrait,ntrait), a (ntrait,ntrait)
    real*8 :: var(ntrait,ntrait),l(ntrait,ntrait)
   se(1,1)=1/sc
    t = 0.
    beta=0.5
    do i=1,ntrait
       alfa=(ne-dble(i)+1.)/2.
       if (alfa .le. 0) then
          print *,'nE:', ne, ntrait
       end if
       call gamma(alfa, beta, x1, u) ! n-2,0.5,x1,u
!      v=sc/u; goto 2
       t(i,i) = sqrt(u)
       do  j = 1, i-1  ! lower
          u =  gasdev(x1)
          t(j,i) = u
       end do
    end do

    !     print *,t(1,1),t(1,2),t(2,1),t(2,2)
    a = 0.d0; l =  0.d0
    a = matmul(t, transpose(t))
    form = 'lower'

    call cholesky(se, l(1:ntrait,1:ntrait), ntrait, 'upper')
    !     print *,se(1,1),se(1,2),se(2,1),se(2,2)
    if(ntrait .ne. size(se, 1)) then
       print *, 'Check rank in wish: ', ntrait ,  size(se, 1)
    end if

    var = matmul(matmul(transpose(l),a), l)
    v=1/var(1,1)
    return
  end subroutine wishart

  subroutine gamma(alfa,beta,x1,u)
    real*8, intent(in)  :: alfa ,beta
    real*8, intent(out) :: u
    real*8 :: xmed, ref, u1, chiot, chitt, o, t, x1, unif
    ! ................................................................................
    !  print *,'alfa,beta',alfa,beta
    o=(alfa-1)/beta
    chiot=(alfa-1.)*log(o)-beta*o
1   u1= unif(x1)
    xmed=(alfa/beta)-(4*sqrt(alfa/(beta**2.)))
    t=u1*(8*sqrt(alfa/(beta**2.)))+xmed
    if (t.lt.0.) then
       !print *,'GAMMA t',u1,x,t
       goto 1
    end if
    chitt = (alfa-1)*log(t)-beta*t
    ref   = exp(chitt-chiot)
    u1= unif(x1)
    if (u1.lt.ref) then
       u = t
    else
       !print *,'GAMMA ',u1,x,ref
       goto 1
    endif
    return
  end subroutine gamma
end module wish


module solve
  use variables
contains
!_______SUB SEIDEL________________________________________________
!Resuelve las ecuaciones por Gauss-Seidel
subroutine seidel (niter,icadena)
   integer:: i,j, icadena,idum !,neq
   real*8:: sum,rhs,lhs,temp2j,ruido,temp2,var_beta,nu

mean=0.d0
do i=1,n_datos
   error(i)=error(i)+mu
   mean=mean+error(i)
enddo
mu=mean/float(n_datos)+xnormal(x1)*sqrt(ve/float(n_datos))
do i=1,n_datos
   error(i)=error(i)-mu
enddo
if (nciclos.eq.1) then
!setup xpx=diag(X«X) AL
    xpx=0
    do i=1,n_cov
        xpx(i)=dot_product(valor(:,i),valor(:,i))
    enddo
 !   do i=n_cov+1,n_cov+n_animal
 !       xpx(i)=valor(i,i)*valor(i,i)
 !   enddo
end if

scale=lambda*lambda
do j=1,n_cov
    !if (abs(sol(j)).lt.0.0000001) sol(j)=0.0001d0 !print *,'beta ',k,' lower than 1.E-08'
    nu=sqrt(ve*scale / (sol(j)*sol(j)) )
    inv_tau2(j)=inv_gauss (nu,scale,x1) 
    tau_F(j)=tau_F(j)+inv_tau2(j)

    temp=0.d0
    do i=1,nlines
        error(i)=error(i)+sol(j)*valor(i,j)
        temp=temp+error(i)*valor(i,j)
    enddo

    temp=temp/(xpx(j)+(1.d0/inv_tau2(j)))
    var_beta=ve/(xpx(j)+(1.d0/inv_tau2(j)))

    sol(j)=xnormal(x1)*sqrt(var_beta)+temp 
    do i=1,nlines
        error(i)=error(i)-sol(j)*valor(i,j)
    enddo
enddo

!!!!!*****tau(j)=ve/vs(j)   !tau
return
end subroutine seidel

end module solve


module threshold_model
   use variables
contains

!algoritmh to sample thresholds
subroutine thresholds

t(1)=0.d0
!Samplear umbrales en caso de más de 2 categorías
if (n_cat.gt.2) then
  call sct(t,sc_t0)
  do k=1,10
     ln_t(1)=0.d0
     t_c(1)=0.d0
     do q=2,n_cat-1
        ln_t(q)=dlog(t(q)-t(q-1))
        ln_t(q)=ln_t(q)+xnormal(x1)*sd_st !0.012d0    !esta desv.st. se puede aumentar o disminuir 0.001, cada 100 iteraciones,
                                             ! si el reject rate, es muy baja o muy alta, respectvmnt.(normal 60-80%)
        t_c(q)=exp(ln_t(q))+t_c(q-1)
     end do

     call sct(t_c,sc_t1)
     if ((sc_t1-sc_t0).gt.2.0d0) then
        ratio=2.0d0
     else if ((sc_t1-sc_t0).lt.-50.0d0) then
        ratio=0.0d0
     else
        ratio= exp(sc_t1-sc_t0)
     end if
     yy=unif(x1)

     if (ratio.gt.yy) then
        sc_t0=sc_t1
        t=t_c
        acc=acc+1
     end if
     tot=tot+1
     cont=cont+1
   end do
   write (*,'(a13,i10,a3,f4.2)') 'reject rate (',tot,')= ',1-acc/float(tot)
if (cont.eq.100) then
   if (1-acc/float(tot).le.0.60) then
       sd_st=sd_st+0.01d0
   else if (1-acc/float(tot).gt.0.75) then
       sd_st=sd_st-0.01d0
   end if
   tot=0
   acc=0
   cont=0
end if
   print *,t(1:n_cat-1)
   write (50,'(15f12.9)') t(1:n_cat-1)
end if
end subroutine thresholds

subroutine liabilities


!Muestrear liabilities y datos censurados
  DO i=1,n_datos
     y_mean=y(i)-error(i)

     if (n_cat.eq.2) then
         if (y1(i).eq.2) then
            y(i)=trun(t(1),y_mean+4*(sqrt(ve)),y_mean,ve,x1)  !muestrea a la derecha del umbral
         else if (y1(i).eq.1) then
           y(i)=trun(y_mean-4*(sqrt(ve)),t(1),y_mean,ve,x1)
         end if

     else if (n_cat.gt.2) then !OJO!! falta programar umbrales de +d2 categorias sin datos cens (se puede poner la ultima categoria como la penultima mas censurada)
             if (y1(i).eq.1) then
                y(i)=trun(y_mean-4*sqrt(VE),t(1),y_mean,ve,x1)
             else if(y1(i).eq.n_cat) then
                y(i)=trun(t(y1(i)),y_mean+4*dsqrt(ve),y_mean,ve,x1)
             else
                y(i)=trun(t(y1(i)-1),t(y1(i)),y_mean,ve,x1)
             end if
      end if
      error(i)=y(i)-y_mean
  END DO

return
end subroutine liabilities


subroutine sct(th,sc_t)
real*8::th(n_cat-1),sc_t,pdf1,pdf2,p,q
  sc_t=0.d0
  DO i=1,n_datos
   temp=y(i)-error(i)
        if (y1(i).eq.1) then
          call normp((th(1)-temp)/dsqrt(ve),p,q,pdf1)
          sc_t=sc_t+log(pdf1)
        else if (y1(i).eq.n_cat) then
         call normp((th(y1(i)-1)-temp)/dsqrt(ve),p,q,pdf1)
          sc_t=sc_t+log(1-pdf1)
        else
          !call normp((th(y1(i))-temp)/dsqrt(ve),p,q,pdf1)
          !call normp((th(y1(i)-1)-temp)/dsqrt(ve),p,q,pdf2)

          !sc_t=sc_t+log(pdf1-pdf2)
          !write(*,*) pdf1,pdf2
        end if
  END DO
write (*,*) 'sct',sc_t
return
end subroutine sct



end module threshold_model


module salida
use variables
contains
subroutine predictions
        integer:: i,j,n_yng
        character(len=15)::id_yng
	real*8, allocatable :: yng_snp(:)
	real*8:: y_est,y_real

 open (79,file='testing.pred.txt', form='formatted')
 n_yng=0
 DO  !Lee el número de lineas en el archivo de genealogia
  read (11,*,iostat=io)
  IF(io.ne.0) EXIT
  n_yng=n_yng+1
 END DO
 rewind (11)
 allocate (yng_snp(n_cov))
 write (79,*) 'SampleID  ACTUAL_PHENTP GEBV_GRS'

 do n=1,n_yng
    read (11,*) y_real,id_yng,yng_snp(1:n_cov)
    y_est=0.d0
    do i=1,n_cov
       y_est=y_est+yng_snp(i)*sol_efectos(i)
    enddo
    write (79,'(a15,2f15.7)') id_yng,y_real,y_est
 enddo

end subroutine predictions
end module salida

program main_gibbs
!Bayesian Analysis applied to Animal Models
!PROGRAMA PRINCIPAL

  !Modulos usados:
  use variables
  use inicial
  use wish
  use solve
  use salida
  use threshold_model

real*8:: gasdev,unif
character(50):: name1

  !Semillas
!  call random_seed(x1)
!  call random_seed(x2)
   x1 = 0.7283928517d+10
   x2 = 0.7283928517d+10
  idum=567*7345




  !Ficheros
  open (29,file='SOL_SNPs.trait', form='formatted')
  open (31,file='VAR_DISTRIB.trait', form='formatted')
  open (33,file='GIBBS_STAT.trait', form='formatted')
  open (34,file='CONVERGENCE.trait', form='formatted')
  open (37,file='Log_LIKELIHOOD.trait', form='formatted')
  open (50,file='thresholds.trait', form='formatted')

!READ PARAMETER FILE
  CALL inic

!************************DATA FILES************************************
! 100 FORMAT (f12.4,f12.0,f8.4,5i8)
! 220 format (f8.0,f12.5,5f3.0,6i8)
! 210 format(10x,2i10,f10.5)
!DATOS
DO  !number of rows in data file
 read (10,*,iostat=io)
 IF(io.ne.0) EXIT
 n_datos=n_datos+1
 n_animal=n_animal+1
END DO
PRINT '(a30,a13,a1,i6)','Number of data in file ',fdatos,'=',n_datos
io=0
PRINT '(a35,i6)','Number of animals =',n_animal
rewind (10)

neq=(nf(n_efectos)+1)

ALLOCATE (id(n_datos),padre(0:n_animal), madre(0:n_animal),codcontrol(n_animal),&
   y1(n_datos),y2(n_datos),y(n_datos),cen(n_datos),d(0:n_animal),sol(neq),error(n_datos),vs(n_cov),&
   solf(n_animal),gebv(n_animal),vep1(n_animal),vep2(n_animal),sol_efectos(neq),rhs(neq),&
   valor(n_datos,n_efectos),t(n_cat-1),ln_t(n_cat-1),t_c(n_cat-1),ani(n_datos),xpx(n_cov))

write(*,*) 'Reading data'

DO i=1,n_datos
        read (10,*,iostat=io) y(i),id(i),valor(i,1:n_cov)
        if (mod(i,500).eq.0) write (*,*) i
        y1(i)=y(i)
        y2(i)=y(i)
        ani(i)=i
        cen(i)=0
        IF (io.ne.0) EXIT
END DO
io=0


  !call ginv1(valor,n_datos,n_datos,0.000001d0,irango)


  print *,'****READING DATA FILE****'
  DO i=1,3
    print *,''
	print '(a7,i2,a13,a13)','Linea: ',i,' del fichero ',fdatos
    print '(a6,a1,f8.2,1x,a11,i2)',name,'=',y(i), 'censurado =',cen(i)
	DO k=1,10
	   IF (k.le.n_cov) THEN
	       print '(a7,i5,a1,f8.2)','Efecto ',k,'=',valor(i,k)
	   END IF
	END DO
  END DO
!PAUSE
CLOSE (10)
! 'DATA AUGMENTATION' IN CASE OF MISSING DATA
     DO i=1,n_datos
     END DO

   k=0;n_df=0
   IF (df.eq.1) THEN
     DO i=1,n_datos
	IF (y2(i).eq.cat_df) then
	   n_df=n_df+1
	END IF
     END DO
     ALLOCATE (vdf(n_df))
     vdf=0;k=0
     DO i=1,n_datos
	IF (y2(i).eq.cat_df) THEN
	    k=k+1
	    vdf(k)=i
	END IF
     END DO
   END IF

call cpu_time(tiempo1)
print *,'    START GIBBS'
gebv=0.d0;sol=0.01d0;temp=0.d0;sol_efectos=0.d0;tau=0.d0;nm=0
lambda=0;mu=0.d0;error=y;mean=0.d0;vep1=0.d0;vep2=0.d0;muf=0.d0
sd_st=0.0001d0;tot=0;acc=0;cont=0
!********************COMIENZA GIBBS***********************
!if (umbral.eq.1) error=error*100
do i=1,n_cov
  sol(i)=unif(x1)*0.1
end do
IF (umbral.eq.1) THEN
t(1)=0
do i=2,n_cat-1
   t(i)=t(i-1)+unif(x1)
end do
ENDIF
do nciclos=1,n_ciclos
 ii=1
if (mod(nciclos,10).eq.0) print '(a16,i9,a4,f16.4,a8,f16.4))','Gibbs iteration ',nciclos,';ve=',ve,';lambda=',lambda



 !RESUELVE POR GAUSS-SEIDEL
   call seidel (niter,ii)
   call wishart_inv
   if (umbral.eq.1) then
      ve=1.d0
      call thresholds
      call liabilities
   end if

!_____SE ACUMULAN LOS VALORES ADITIVOS Y DE LAS VARIANZAS______

   write (31,340) ve

   if (nciclos.gt.burnin) then
       nm=nm+1
       muf=muf+mu
       do i=1,neq
	     sol_efectos(i)=sol_efectos(i)+sol(i)
       end do
       if (mod(nciclos,500).eq.0) then
          write(*,*) '  RUNNING MEAN= ',muf/float(nm)
          open (32,file='GEBV_GRS.txt', form='formatted')
          write (32,*) 'SampleID ObservedPhenotype GEBV'
		gebv=0.d0
		do i=1,n_datos
			do j=1,n_cov
				gebv(i)=gebv(i)+(sol_efectos(j)/float(nm))*valor(i,j)
			enddo
			write (32,'(a15,3f16.8)') id(i),y2(i),gebv(i)+muf/float(nm)
		end do
		close(32)
	  endif
    end if

 340 format (f16.5,5x,f16.4,5x,f16.5)
end do

!______Final de Gibbs-Sampling. Creación de ficheros de salida____________________
print *,'writting solutions and vep'
 VEf = VES/float(nm)
 hf=hs/float(nm)
 muf=muf/float(nm)
 do i=1,neq
  sol_efectos(i)=sol_efectos(i)/float(nm)
  tau_F(i)=tau_F(i)/float(nm)
 end do

 do i=1,n_efectos
  write (29,'(a6,3x,i7,3x,f15.8,3x,f15.8)') 'Alpha',i,sol_efectos(i),tau_F(i)
 end do


 open (32,file='GEBV_GRS.txt', form='formatted')
 write (32,*) 'SampleID ObservedPhenotype GEBV'

 gebv=0.d0
 do i=1,n_datos
  do j=1,n_cov
     gebv(i)=gebv(i)+sol_efectos(j)*valor(i,j)
  enddo
  write (32,'(a15,2f16.8,2f24.8)') id(i),y2(i),gebv(i)
end do

 290 format (i6,2x,f16.8)
 320 format (i6,1x,a14,2(2x,f16.8))

call predictions

print *,'FINISH'
print *,muf
print *,'Statistics in file GIBBS_STAT'
print *,'Genetic Parameters Samples in file VARIANZAS_DISTRIB'
print *,'Estimated breeding values and VEP in file VALORES_GENETICOS'
print *,'Average values for systematic effects in file SOLUCION_EFECTOS'
print *,'Samples for convergence test in file CONVERGENCE'
stop




77 continue
write (*,*) 'ERROR, animal ',i,' ',id(i),' not found in the phenotype file'
write (*,*) 'All genotyped animals in the training set must have phenotypic record'
write (*,*) 'and must appear in the same order'
stop


end program main_gibbs

      subroutine cholesky(c, chol, ncar,form)
  !calcula el factor de cholesky (chol) de la matriz C
      implicit none
      character(len=5):: form !'upper' para q saque el triangulo superior
      integer:: i,j,k
      integer:: ncar
      real*8, intent (out):: chol(ncar,ncar)
      real*8::c(ncar,ncar)

      do i=1,ncar
        do j=1,ncar
             if (j.gt.i) then
                  chol(i,j)=0
                  cycle
             end if
             if (i.eq.j) then
                  chol(i,j)=c(i,j)
                  if (i.gt.1) then
                        do k=1,i-1
                             chol(i,j)=chol(i,j)-chol(i,k)*chol(i,k)
                        end do
                  end if
                  chol(i,j)=sqrt(chol(i,j))
              else
                  chol(i,j)=c(i,j)
                  if (j.gt.1) then
                        do k=1,j-1
                              chol(i,j)=chol(i,j)-chol(j,k)*chol(i,k)
                        end do
                  end if
                  chol(i,j)=chol(i,j)/chol(j,j)
              end if
        end do
      end do
    !Para que saque el triangulo superior
      if (form.eq.'upper') then
       chol=transpose(chol)
      end if
      return
      end subroutine cholesky


      function alnorm ( x, upper )
!*********************************************
!* ALNORM computes the cumulative density of *
!*    the standard normal distribution.      *
!*********************************************
!
!  Reference:
!
!    I Hill,
!    The Normal Integral,
!    Algorithm AS 66,
!    Applied Statistics,
!    Volume 22, Number 3, pages 424-427, 1973.
!
!  Modified:
!
!    28 March 1999
!
!  Parameters:
!
!    Input, real X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
       real*8 alnorm
       logical up
       logical upper
       real*8 x
       real*8 y
       real*8 z
       real*8 a1, a2, a3, b1, b2, c1, c2, c3, c4, c5, c6
       real*8 con, d1, d2, d3, d4, d5, ltone, p, q, r, utzero
       a1 = 5.75885480458
       a2 = 2.62433121679
       a3 = 5.92885724438
       b1 = -29.8213557807
       b2 = 48.6959930692
       c1 = -0.000000038052
       c2 = 0.000398064794
       c3 = -0.151679116635
       c4 = 4.8385912808
       c5 = 0.742380924027
       c6 = 3.99019417011
       con = 1.28
       d1 = 1.00000615302
       d2 = 1.98615381364
       d3 = 5.29330324926
       d4 = -15.1508972451
       d5 = 30.789933034
       ltone = 7.0
       p = 0.398942280444
       q = 0.39990348504
       r = 0.398942280385
       utzero = 18.66
       up = upper
       z = x

       if ( z .lt. 0.0 ) then
         up = .not. up
         z = - z
       end if
!
!  TAKE ANOTHER LOOK AT THIS SET OF CONDITIONS.
!
       if(z.gt.ltone.and.((.not.up).or.z.gt.utzero))then
         if ( up ) then
           alnorm = 0.0
         else
           alnorm = 1.0
         end if
         return
       end if
       y = 0.5 * z**2
       if ( z .le. con ) then
         alnorm=0.5-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
       else
         alnorm=r*exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4&
               /(z+c5+d5/(z+c6))))))
       end if
       if ( .not. up ) then
         alnorm = 1.0 - alnorm
       end if
       return
       end

       function ppnd ( p, ifault )
!********************************************************
!* PPND produces the normal deviate value corresponding *
!*          to lower tail area = P.                     *
!********************************************************
!  Reference:
!    J Beasley and S Springer,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 111,
!    Applied Statistics,
!    Volume 26, Number ?, pages 118-121, 1977.
!
!  Modified:
!    28 March 1999
!
!  Parameters:
!    Input:real P, value of cumulative probability densitity function.
!    0 < P < 1.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, P <= 0 or P >= 1.  PPND is returned as 0.
!
!    Output, real PPND, the normal deviate value with the property that
!    the probability of a standard normal deviate being less than or
!    equal to PPND is P.
!
       integer ifault
       real*8 p
       real*8 ppnd
       real*8 r
       real*8 a0, a1, a2, a3, b1, b2, b3, b4
       real*8 c0, c1, c2, c3, d1, d2, split
       a0 = 2.50662823884
       a1 = -18.61500062529
       a2 = 41.39119773534
       a3 = -25.44106049637
       b1 = -8.47351093090
       b2 = 23.08336743743
       b3 = -21.06224101826
       b4 = 3.13082909833
       c0 = -2.78718931138
       c1 = -2.29796479134
       c2 = 4.85014127135
       c3 = 2.32121276858
       d1 = 3.54388924762
       d2 = 1.63706781897
       split = 0.42

       ifault = 0
       if (abs(p-0.5) .le. split) then
         r = (p-0.5)**2
         ppnd = (p-0.5)*(((a3*r+a2)*r+a1)*r+a0)/((((&
              b4*r+b3)*r+b2)*r+b1)*r+1.0)
       else if ( p .gt. 0.0 .and. p .lt. 1.0 ) then
         if ( p .gt. 0.5 ) then
           r = sqrt(-log(1.0-p))
         else
           r = sqrt(-log(p))
         end if
         ppnd=(((c3*r+c2)*r+c1)*r+c0)/((d2*r&
             +d1)*r+1.0)
         if (p .lt. 0.5) then
           ppnd = -ppnd
         end if
       else

         ifault = 1
         ppnd = 0.0
       !  write ( *, * ) ' '
       !  write ( *, * ) 'PPND - Warning!'
       !  write ( *, * ) '  P <= 0 or P >=1.'
       !  write ( *, * ) '  PPND value would be infinite.'
       end if
       return
       end



subroutine normp ( z, p, q, pdf )

!*****************************************************************************80
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!  Discussion:
!
!    This is algorithm 5666 from Hart, et al.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, divides the real line into two
!    semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, real ( kind = 8 ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real ( kind = 8 ) PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  real*8:: cutoff = 7.071D0
  real*8:: expntl
  real*8:: p
  real*8:: p0 = 220.2068679123761D0
  real*8:: p1 = 221.2135961699311D0
  real*8:: p2 = 112.0792914978709D0
  real*8:: p3 = 33.91286607838300D0
  real*8:: p4 = 6.373962203531650D0
  real*8:: p5 = 0.7003830644436881D0
  real*8:: p6 = 0.03526249659989109D0
  real*8:: pdf
  real*8:: q
  real*8:: q0 = 440.4137358247522D0
  real*8:: q1 = 793.8265125199484D0
  real*8:: q2 = 637.3336333788311D0
  real*8:: q3 = 296.5642487796737D0
  real*8:: q4 = 86.78073220294608D0
  real*8:: q5 = 16.06417757920695D0
  real*8:: q6 = 1.755667163182642D0
  real*8:: q7 = 0.08838834764831844D0
  real*8:: root2pi = 2.506628274631001D0
  real*8:: z
  real*8:: zabs

  zabs = abs ( z )
!
!  37 < |Z|.
!
  if ( 37.0D0 < zabs ) then

    pdf = 0.0D0
    p = 0.0D0
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5D0 * zabs * zabs )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
          p6   * zabs &
        + p5 ) * zabs &
        + p4 ) * zabs &
        + p3 ) * zabs &
        + p2 ) * zabs &
        + p1 ) * zabs &
        + p0 ) / ((((((( &
          q7   * zabs &
        + q6 ) * zabs &
        + q5 ) * zabs &
        + q4 ) * zabs &
        + q3 ) * zabs &
        + q2 ) * zabs &
        + q1 ) * zabs &
      + q0 )
!
!  CUTOFF <= |Z|.
!
    else

      p = pdf / ( &
        zabs + 1.0D0 / ( &
        zabs + 2.0D0 / ( &
        zabs + 3.0D0 / ( &
        zabs + 4.0D0 / ( &
        zabs + 0.65D0 )))))

    end if

  end if

  if ( z < 0.0D0 ) then
    q = 1.0D0 - p
  else
    q = p
    p = 1.0D0 - q
  end if

  return
end



!___________SUB LUNIF___________________________
!Generacion de un número uniforme,s1 la semilla !
      function lunif (s1,ll)                    !
      implicit doubleprecision (a-h,o-z)        !
      doubleprecision s1,unif                   !
      lunif=int(unif(s1)*ll)+1                  !
      return                                    !
      end                                       !
!_______________________________________________!

!______________SUB UNIF_________________________
!Generacion de un número uniforme,s1 la semilla !
      function unif (s1)                        !
      implicit doubleprecision (a-h,o-z)        !
      real*8 s1,unif                   !
      s1 = mod (s1*16807.0d0,21477483647.0d0)   !
      unif = s1 / 21477483647.0d0               !
   !call random_number(unif)
   !lsol=(int((s1*16807.0d0)/2147483647.0d0))   !
   !s1=s1*16807.0d0-lsol*21477483647.0d0        !
   !unif=s1/2147483647.0d0                      !
      return                                    !
      end                                       !
!_______________________________________________!

!______________SUB XNOR______________________________
 !Calcula la abcisa para cualquier función           !
 ! de distribución prob                              !
      FUNCTION XNOR(prob)                            !
      implicit doubleprecision (a-h,o-z)             !
      if (prob.lt..5) then                           !
      p=prob                                         !
      else if (prob.gt..5) then                      !
      p=1.-prob                                      !
      else                                           !
      xnor=0.                                        !
      end if                                         !
      t=dsqrt (log(1./(p*p)))                        !
      x=t-(2.515517+t*(.802853+t*.010328))/&         !
          (1.+t**(1.432788+t*(.189269+t*.001308)))   !
      if (p.eq.prob) then                            !
        xnor=-x                                      !
      else                                           !
        xnor=x                                       !
      end if                                         !
      return                                         !
      END                                            !
!____________________________________________________!

!_______________X_NORMAL________________
      FUNCTION XNORMAL (x1)             !
      implicit doubleprecision (a-h,o-z)!
      xnormal=xnor(unif(x1))            !
      return                            !
      END                               !
!_______________________________________!

!_________________________________________________
      function GASDEV(x2)
!    *******************************************
!    *                                         *                          *
!    *  ESTA FUNCION GENERA NUMEROS ALEATORIOS *
!    *           DE UNA DISTRIBUCION           *
!    *     NORMAL DE MEDIA Y VARIANZA 1        *                      *
!    *                                         *                      *
!    *******************************************
      implicit doubleprecision (a-h,o-z)
	  real*8::R,unif
      DATA ISET/0/

      IF(ISET.EQ.0)THEN
    1 V1=2.0*unif(x2)-1.0
      V2=2.0*unif(x2)-1.0
      R=V1**2+V2**2
      IF(R.GE.1.)GO TO 1
      FAC=DSQRT(-2.0*DLOG(R)/R)
      GSET=V1*FAC
      gasdev=V2*FAC
      ISET=1
      ELSE
      gasdev=GSET
      ISET=0
      ENDIF
      RETURN
      END function
!___________________________________________________

      FUNCTION trun(ini,fin,mean,var,x1)
        !ojo:::: de momento no uso var porque siempre es uno para umbrales
        implicit none
        real*8 ini, fin, mean, var, trun
        real*8 p_ini, p_fin, pepe
        real*8 x1
        real*8 alnorm, ppnd, unif
        integer ifault
        logical upper
        if (ini .gt. fin)then
           trun = ( ini + fin )/2.0d0
!		 write (*,*) 'CUIDADO!! muestreando liabilities > 4 desv.tip.'
           return
        end if
        upper = .false.
!c        print *,'t-ini ', ini,' t-fin ', fin,' media ',
!c     1          mean,' varianza ',var
        !estandarizar umbrales
        p_ini = (ini - mean)/dsqrt(var)
        p_fin = (fin - mean)/dsqrt(var)
        p_ini = alnorm(p_ini, upper)
        p_fin = alnorm(p_fin, upper)
!c        print *,'probs   ',p_ini,p_fin
        trun = p_ini+(p_fin-p_ini)*unif(x1)
        pepe = ppnd(trun,ifault)
        if (ifault .eq. 0) then
            trun = mean +  pepe*dsqrt(var)
        else
            trun = ( ini + fin )/2.0d0
!		 write (*,*) 'CUIDADO!! muestreando liabilities'
!           write (*,*) 'a mas de 4 desviaciones tipicas.'
        end if
        !print *,'sampleo de la truncada ',trun
!c        pause
        return

      END FUNCTION



      FUNCTION RAN1(IDUM)
        DIMENSION R(97)
        PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=1./M1)
        PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=1./M2)
        PARAMETER (M3=243000,IA3=4561,IC3=51349)
        DATA IFF/0/
!   *******************************************************************
!   *                                                                 *
!   *   ESTA SUBRUTINA GENERA NUMERO ALEATORIOS DE UNA DISTRIBUCION   *
!   *     UNIFORME EN EL RANGO 0.0-1.0                                *
!   *                                                                 *
!   *******************************************************************
        IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
           IX1=MOD(IA1*IX1+IC1,M1)
           IX2=MOD(IA2*IX2+IC2,M2)
           R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
   11   CONTINUE
        IDUM=1
        ENDIF
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IA2*IX2+IC2,M2)
        IX3=MOD(IA3*IX3+IC3,M3)
        J=1+(97*IX3)/M3
        !IF(J.GT.97.OR.J.LT.1) PAUSE
        RAN1=R(J)
        R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
        RETURN
      END function



      SUBROUTINE gamma2(shape,rate,x2,xx)
      !**********************
      !***MUESTREA DE UNA ***
      !*** FUNCION GAMMA  ***
      !**********************
      !*  @arg rate parameter (lambda)
      !*  @arg shape parameter (r)
      !*  Algorithm to sample from a gamma distribution. Ahrens and Dieter (1974) Computing, 12:223-246.
       !*/
       real*8  :: rate ,shape
       real*8 :: xx,x2
       real*8:: unif,b,c,p,t,s,u,uu
       if (shape.le.1.d0) then !Algorithm GS
           1 u =unif(x2)
           b=(exp(1.d0)+shape)/exp(1.d0)
           p=b*u
           if (p.le.1.d0) then
                 x=p**(-shape)
                 uu=unif(x2)
                 if (uu.gt.exp(-x)) goto 1
           else
                 x=-log((b-p)/shape)
                 uu=unif(x2)
                 if (uu.gt.x**(shape -1.d0)) goto 1
           end if
       else                   !Algorithm GC
           b=shape-1
           c=shape+b
           s=dsqrt(c)
         2 u=unif(x2)
           t=s*tan(3.14159d0*(u-0.5d0))
           x=b+t
           if (x.lt.0) goto 2
           uu=unif(x2)
           if (uu.gt.exp(b*log(x/b)-t+log(1+t**2/c)) ) goto 2
       endif
       xx=(x-1.d0)*rate

       return
      END SUBROUTINE gamma2
!==================================================
FUNCTION inv_gauss(mu, lambda,x1)

! Adapted from Fortran 77 code from:
!     Michael, Schucany and Haas (1976). Generating random variates using
!         transformations with multiple roots. The American Statistician, 30:88-90


IMPLICIT NONE
REAL*8:: mu, lambda
!     Local variables
REAL*8:: x1,unif,xnormal,inv_gauss
REAL*8:: v,a,b,c,x

inv_gauss=0.d0
v=xnormal(x1)
v=v*v !sample from a chi square distribution
a=mu+0.5d0*(v*mu**2)/lambda-0.5d0*mu/lambda*sqrt(4*mu*lambda*v+(mu*v)**2)
b=mu**2/a
c=mu/(mu+a)
if (unif(x1).le.c) then
x=a
else
x=b
endif

inv_gauss = x

RETURN
END FUNCTION inv_gauss
!==================================================

! --------------------------------------
      subroutine ginv1(a,n,m,tol,irank)
! returns generalized inverse of matrix x of size n x n declared
! as m x m. tol is working zero and irank returns the rank of
! the matrix. w is a work vector of size m,
! by rohan fernando, slightly structured by i. misztal 05/05/87

      double precision:: a(m,m),w(m),re,sum,tol
      irank=n
      do 10 i=1,n
         do 20 j=1,i-1
              re=a(i,j)
              do 20 ii=i,n
20                 a(ii,i)=a(ii,i)-re*a(ii,j)
         if (a(i,i).lt.tol) then
              a(i,i)=0.0
              do 45 ii=i+1,n
45                 a(ii,i)=0.
           irank=irank-1
           else
              a(i,i)=sqrt(a(i,i))
              do 40 ii=i+1,n
40                a(ii,i)=a(ii,i)/a(i,i)
         endif
10    continue

      do 100 i=1,n
         if (a(i,i).eq.0.) then
              do 150 ii=i+1,n
150                a(ii,i)=0.
           else
              a(i,i)=1.0/ a(i,i)
              do 200 ii=i+1,n
200               w(ii)=0.0
              do 300 ii=i+1,n
                  iim1=ii-1
                  re=a(iim1,i)
                  do 400 iii=ii,n
400                   w(iii)=w(iii)-a(iii,iim1)*re
                  if (a(ii,ii).eq.0.) then
                      a(ii,i)=0.
                    else
                      a(ii,i)=w(ii)/a(ii,ii)
                  endif
300           continue
          endif
100     continue

      do 110 j=1,n
         do 110 i=j,n
              sum=0
              do 130 ii=i,n
130                sum=sum+a(ii,j)*a(ii,i)
110           a(i,j)=sum
      do 600 i=1,n
          do 600 j=i,n
600           a(i,j)=a(j,i)

      end subroutine


