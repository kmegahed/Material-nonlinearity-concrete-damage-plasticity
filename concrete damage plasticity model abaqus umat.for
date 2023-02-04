c     see reference P. Grassl, D. Xenos, U. Nyström, R. Rempling
c     , K. Gylltoft. "CDPM2: A damage-plasticity approach to 
c     modelling the failure of concrete". International Journal '
c     of Solids and Structures. Volume 50, Issue 24, pp. 3805-3816, 2013.
      SUBROUTINE UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION sig(NTENS),statev(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 stran(NTENS),dstran(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4)
      dimension cc(ntens,ntens),cin(ntens,ntens),de(ntens),rso(4),
     $eps(6),old_e(6),pj(4,4),unkn1(4),sg_p(3,3),conv_e(6),
     $dtot_e(6),e_p(6),e_rate(6),el_e(6),el_e_old(6),
     $sg_o(6),sg_eff(6),sg_old(6),sig_eff_c(6),sig_eff_t(6),
     $strain_rate(6),tot_e(6),tot_e1(6),ttot_e(6),pjo(4,4),
     $rs(4),pj_prev(4,4),sg_tr(6),dinv_dsig_pr(3,3)
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %tol=1.0d-6,tol2=10.0d0)
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,
     $iunload_flag,itypey,istrrateflg
      em=props(1);qmuo=props(2);fc=props(3);ft=props(4);
      fc0=props(5);hp=props(6);qh0=props(7);pm0=props(8);
      ah=props(9);bh=props(10);ch=props(11);dh=props(12);
      e0=props(13);ecc=props(14);gtol=props(15);itypey=props(16);
      wf=props(17); wf1=props(18);ft1=props(19); efc=props(20);
      istrrateflg=0
      em=24000.0;qmuo=0.2;fc=30.0;ft=3.0;gft=0.1;gtol=1.0d-4
      itypey=1;irate_ekffect=0;isotropic=2;imaxsubinc=20;iendflag=0
      bs=1.;df=0.85;
      ah=0.08;bh=0.003;ch=2.0; dh=1.0d-6;hp=0.01;as=15.0;
      fc0=fc**1.855/60.;qh0=fc0/fc;fb=1.50*fc**(-0.075)  
      par=ft/fb*(fb**2.-fc**2.)/(fc**2.-ft**2.)
      ecc=(1.+par)/(2.-par)
      pm0=3*(fc**2.-ft**2.)/fc/ft*ecc/(ecc+1.)
      e0=ft/em;wf=gft/0.225/ft;wf1=0.158*wf;ft1=0.3*ft;efc=0.0001
      tkp=statev(1);eps_t=statev(2);statev(28)=statev(28)
      de(1:6)=dstran;old_e(1:6)=statev(21:26);plen=CELENT;
      e_rate=de;isubinc_count=0;isubinc_flag=0;iconvrg=1
      cc=zr;cc(1:3,1:3)=p1*qmuo;
      cc(1,1)=(on-qmuo)*p1;cc(2,2)=(on-qmuo)*p1;
      cc(3,3)=(on-qmuo)*p1;Gm=Em/(on+qmuo)/tw;
      cc(4,4)=Gm;cc(5,5)=Gm;cc(6,6)=Gm;
      ddsdde=cc;
      cin=-qmuo/Em;cin(1,1)=on/Em;cin(2,2)=on/Em;
      cin(3,3)=on/Em;cin(4,4)=on/Gm;cin(5,5)=on/Gm;cin(6,6)=on/Gm;
      tot_e=old_e+de;strain_rate=de
      tot_e1=tot_e;
      ttot_e=tot_e;
      do while (iconvrg.eq.1.or.isubinc_flag.eq.1)
            el_e=ttot_e-e_p;
            sg_tr=matmul(cc,el_e)
            call khaigh(sg_tr,sv_tr,ro_tr,theta_tr,dinv_dsig_pr);
            call kff(sv_tr,ro_tr,theta_tr,tkp,yield);
            apex_sg=0.
            write(*,*) 'cin50',em
            if(yield>0.0) then
                  irtype=0
                  write(*,*) 'cin1',em
                  call kcheckvertex(sv_tr,tkp,apex_sg,irtype)
                  if (irtype.eq.1.or.irtype.eq.2) then
                        call kvertexreturn(sg_tr,apex_sg,
     $                   tkp,irtype,iconvrg,sg_tr)
                  end if
                  if(irtype.eq.0) then
                        call kregular_return(sg_tr,tkp,
     $   iendflag,gtol,sg_tr,iconvrg,tkp,pj,rs,pnorm_rs)
                  end if
            else 
                  iconvrg=0;
                  write(*,*) 'cin7',em
                  do jj=3,8
                        e_p=statev(jj)
                        write(*,*) 'cin17',em
                  end do
                  write(*,*) 'cin700',em
                  exit
            end if
            write(*,*) 'cin4',em
            if (iconvrg.eq.1) then
                  isubinc_counter=isubinc_counter+1;
                  write(*,*) 'cin5',em
                  if (isubinc_counter>imaxsubinc) then
                        write(*,*) 'erorrrr'
                        call xit
                  else if (isubinc_counter > imaxsubinc-1
     $      .and.tkp < 1.) then
                        tkp=1
                  end if
                  isubinc_flag=1;dtot_e=0.5*dtot_e;
                  ttot_e=conv_e+dtot_e;pj_prev=pj;
            else if (iconvrg.eq.0.and.isubinc_flag.eq.0) then
                  el_e=matmul(cin,sg_tr)
                  e_p=tot_e-el_e;
            else if (iconvrg.eq.0.and.isubinc_flag.eq.1) then
                  el_e=matmul(cin,sg_tr)
                  e_p=ttot_e-el_e;conv_e=ttot_e
                  dtot_e=tot_e-ttot_e;ttot_e=tot_e;
            end if
      end do
      if(itypey.eq.3) then
            wt_o=0;wc_o=0;eps_t=0;eps_c=0;pkdt=0;pkdt1=0;
            pkdt2=0
            pkdc=0;pkdc1=0;pkdc2=0;rate_fac=0;alpha=0
      else
            rate_fac=statev(17);eps_t=statev(19);eps_c=statev(20)
            pkdt=statev(9);pkdt1=statev(10);pkdt2=statev(11)
            pkdc=statev(12);pkdc1=statev(13);pkdc2=statev(14)
            wt=statev(15);wc=statev(16);alpha=statev(18)
            write(*,*) 'cin711',em
            call kc_alpha(sg_tr,sig_ekff_t,sig_ekff_c,alpha);
            write(*,*) 'cin71122',em
            do jj=3,8
                  el_e_old=statev(18+jj)-statev(jj)
            end do
            sg_old=matmul(cc,el_e_old) 
            parrr=0.
            do jj=3,8
                  parrr=parrr+(e_p(jj-2)-statev(jj))**2
            end do
            pnorm_inc_e_p=sqrt(parrr)
            write(*,*) 'cin7112',em
            call kdamage(wc_o,wt_o,strain_rate,rate_fac,
     $ alpha,eps_t,eps_c,pkdt,pkdt1,pkdt2,pkdc,pkdc1,pkdc2,sg_tr,
     $ tkp,pnorm_inc_e_p,plen,sg_old,statev(18),statev(27),
     $ wc,wt,eps_t,eps_c,pkdt,pkdt1,pkdt2,
     $ pkdc,pkdc1,pkdc2,statev(27))!statev(18) old alpha
      write(*,*) 'cin71143',em
      end if
      if (isotropic.eq.1) then
            sig=(1-wt)*sg_tr
            ddsdde=(1-wt)*ddsdde
      else 
            write(*,*) 'cin7114uu3',em
            sig=(1-wt)*sig_ekff_t+(1-wc)*sig_ekff_c
            ddsdde!!!!!!!!!!!!!!
      end if 
      write(*,*) 'cin7114vv',em
      statev(1)=tkp;statev(2)=eps_t;statev(3:8)=e_p;
      statev(9)=pkdt;statev(10)=pkdt1;statev(11)=pkdt2;
      statev(12)=pkdc;statev(13)=pkdc1;statev(14)=pkdc2;statev(15)=wt;
      statev(16)=wc;statev(17)=rate_fac;statev(18)=alpha;
      statev(19)=eps_t;statev(20)=eps_c;statev(21:26)=tot_e;
      write(*,*) 'cin88',em
      return
      end  
      subroutine kregular_return(sg,pk,iendflag,gtol,sg2,iconvrg,
     $ pkp2,pj,resd)
      INCLUDE 'ABA_PARAM.INC'
      real*8::resd(4),pnorm_res(4),pincrmt(4),pkm,
     $ qmuo,em,dl,sg(6),sv_tr1,ro_tr1,ro,sv,theta_tr1,dir(3,3),
     $ unkn1(4),tkp1,pkp,pk,gtol,sg2(6),pkp2,pj(4,4),sig_ten(3,3),
     $ dgdinv(2),ddkdldinv(2),ddgddinv(2,2),ddgdinvdk(2),sg_p(3,3),
     $ddgdinv(2,2),pi,parry,dinv_dsig_pr(3,3),
     $ dfdinv(2),ddg_ddinv(2,2),dfdk,ddg_dinvdk(2),ddk_dldinv(2)
      integer::iter,itot_iter,iendflag,iconvrg,icomputeall,i,jjj
      common/kumat/qmuo,em
      pi=4.*atan(1.d0)
      !input sg,pkp,iendflag,gtol
      iter=0;itot_iter=200;pkm=em/3./(1.-2.*qmuo);gm=em/2./(1.+qmuo)
      resd(1:4)=0.;pnorm_res(1:4)=0.;pincrmt(1:4)=0.0;dl=0.;      
      call khaigh(sg,sv_tr1,ro_tr1,theta_tr1,dinv_dsig_pr)
      call kvec_to_tens(sg,sig_ten)
      call kjacobi_eigenvalue(3,sig_ten,dir,sg_p)
      
      pkp=pk;tkp1=pk;sv=sv_tr1;ro=ro_tr1;
      unkn1(1)=sv_tr1;unkn1(2)=ro_tr1;unkn1(3)=tkp1;unkn1(4)=0.;
      call kff(sv,ro,theta_tr1,tkp1,resd(4));
      norm_res=1.;iconvrg=0;
      do while (norm_res>gtol)
            iter=iter+1
            pnorm_res(1)=resd(1)/pkm;pnorm_res(2)=resd(2)/2./gm
            pnorm_res(3)=resd(3);pnorm_res(4)=resd(4)
            norm_res=norm2(pnorm_res)
!            if((iter.eq.itot_iter).or.(isnan(norm_res))) then
            if(iter.eq.itot_iter) then
                  iconvrg=1
                  exit
            end if
            if (norm_res>gtol) then
                  icomputeall=1
                  call kderv(sv,ro,theta_tr1,tkp1,icomputeall,
     $    dgdinv,dkdl,dfdinv,ddgddinv,dfdk,ddgdinvdk,ddkdldinv)
                  pj(1,1)=1.+pkm*dl*ddgddinv(1,1);
                  pj(1,2)=pkm*dl*ddgddinv(1,2);
                  pj(1,3)=pkm*dl*ddgdinvdk(1);
                  pj(1,4)=pkm*dgdinv(1);
                  pj(2,1)=2.*gm*dl*ddgddinv(2,1);
                  pj(2,2)=1.+2.*gm*dl*ddgddinv(2,2);
                  pj(2,3)=2.*gm*dl*ddgdinvdk(2);
                  pj(2,4)=2.*gm*dgdinv(2);pj(3,1)=dl*ddkdldinv(1);
                  pj(3,2)=dl*ddkdldinv(2);pj(3,3)=dl*ddkdldk-1.;
                  pj(3,4)=dkdl;pj(4,1)=dfdinv(1);
                  pj(4,2)=dfdinv(2);pj(4,3)=dfdk;pj(4,4)=0.;
                  call Kdet44(pj,parry)
                  if (abs(parry)<1.0d-10) then
                        iconvrg=1
                        exit
                  end if
                  call kmatinv4(pj,pj)
                  pincrmt=-matmul(pj,resd)
                  unkn1=unkn1+pincrmt
                  if (unkn1(4)<0.) unkn1(4)=0.
                  if (unkn1(2)<0.) unkn1(2)=0.
                  if (unkn1(3)<pkp) unkn1(3)=pkp
                  sv=unkn1(1);ro=unkn1(2);
                  tkp1=unkn1(3);dl=unkn1(4);
                  call kderv(sv,ro,theta_tr1,tkp1,0,dgdinv,dkdl,
     $  dfdinv,ddg_ddinv,dfdk,ddg_dinvdk,ddk_dldinv)
                  resd(1)=sv-sv_tr1+km*dl*dgdinv(1)
                  resd(2)=ro-ro_tr1+2.*gm*dl*dgdinv(2)
                  resd(3)=pkp-tkp1+dl*dkdl
                  call kff(sv,ro,theta_tr1,tkp1,resd(4))
            end if 
      end do
      if(iconvrg.eq.0) then
            sg_p=0.;
            sg_p(1,1)=sv+sqrt(2./3.)*ro*cos(theta_tr1)
            sg_p(2,2)=sv+sqrt(2./3.)*ro*cos(theta_tr1-2./3.*pi)
            sg_p(3,3)=sv+sqrt(2./3.)*ro*cos(theta_tr1+2./3.*pi);
       call ktens_to_vec(matmul(transpose(dir),matmul(sig_ten,dir)),sg)
            pk=tkp1
      end if
      return
      end
      subroutine khaigh(sg,sv,ro,th,dinv_dsig_pr)
            implicit none
            real*8::sg(6),s(6),sv,ro,th,gtol,pj2,pj3,dinv_dsig_pr(3,3)
            sv=(sg(1)+sg(2)+sg(3))/3.;s=sg;s(1:3)=sg(1:3)-sv
            pj2=0.5*(s(1)**2+s(2)**2+s(3)**2)+s(4)**2+s(5)**2+s(6)**2
            if (pj2<gtol) THEN
                  th=0.;ro=0.
            else
                  ro=sqrt(2.*pj2);
                  pj3=s(1)**3.+s(2)**3.+s(3)**3.+3*s(1)*(s(4)**2.
     $             +s(6)**2.)+6*s(4)*s(5)*s(6)+3*s(2)*(s(4)**2.
     $             +s(5)**2)+3*s(3)*(s(5)**2+s(6)**2)
                  pj3=pj3/3.;th=3.*sqrt(3.)/2.*pj3/(pj2**1.5)
            end if
            if (th>1.) then
                  th=1.
            else if (th<-1.0) then
                  th=-1.0
            end if
            th= 1./3.*acos(th)
            return
            end
      subroutine kff(sv,ro,th,pkp,f)
            implicit none
            real*8::sv,ro,th,pkp,f,par1,par2,rcos,qh1,dqh1_dk,
     $ qh2,dqh2_dk,ecc,pm0,fc,e0
            common/kumat/ecc,pm0,fc,e0
            !input sv,ro,th,pkp
            !this part is repeated if you change it check others
            par1=(1.-ecc**2);par2=(2.*ecc-1.)
            rcos=(4.*par1*cos(th)**2.+par2**2.)
            rcos=rcos/(2*par1*cos(th)+par2*sqrt(4*par1*cos(th)**2.
     $       +5*ecc**2-4*ecc))
            call kcqh1(pkp,int(0),qh1,dqh1_dk)
            call kcqh2(pkp,int(0),qh2,dqh2_dk)
            par1=ro/(fc*sqrt(6.))+sv/fc;par2=(1.-qh1)*par1**2.
     $       +sqrt(1.5)*ro/fc
            f=par1**2+qh1**2*qh2*pm0*(sv/fc+ro*rcos/(sqrt(6.)*fc))
     $       -qh1**2*qh2**2
            return
            end
      subroutine kcqh1(par,icomput_deriv,qh1,dqh1_dk)
            implicit none
            real*8::hp,par,qh1,dqh1_dk
            integer::icomput_deriv
            common/kumat/hp
            if (par<=0.) then
                  qh1=1.0
            else if (par>0.0.and. par<1.0) then
                  qh1=1.0
            else
                  qh1=1.+(par-1.)*hp
            end if
            if (icomput_deriv.eq.1) then
                  if (par<=1.) then
                        dqh1_dk=0.
                  else
                        dqh1_dk=hp
                  end if
                  dqh1_dk=0.
            end if
            return
            end
      subroutine kcqh2(par,icomput_deriv,qh2,dqh2_dk)
            implicit none
            real*8::hp,qh0,par,qh2,dqh2_dk,rcos,r_top,r,par2,par1,
     $ pmQ
            integer::icomput_deriv,i
            common/kumat/hp,qh0
            if (par<=0.) then
                  qh2=qh0
            else if (par>0.0.and. par<1.0) then
                  qh2=(1.0-qh0-hp)*par**3.-3.*(1.-qh0-hp)*par**2.
     $             +(3.*(1.-qh0)-2.*hp)*par+qh0
            else
                  qh2=1.
            end if
            if (icomput_deriv.eq.1) then
                  if (par<=0.) then
                        dqh2_dk=0.
                  else if (par<=1.) then
                        dqh2_dk=3*(1.0-qh0-hp)*par**2.-6.*(1.-qh0-hp)
     $                   *par**2.+(3.*(1.-qh0)-2.*hp)
                  else
                        dqh2_dk=0.
                  end if
                  dqh2_dk=0.
            end if
            return
            end          
      subroutine kc_alpha(s,st,sc,alpha)
            implicit none
            real*8::s(6),sig_pr(3),dir(3,3),s_pt(6),s_pc(6),
     $ alpha,alpha_t,pnorm_s,st(6),sc(6),sig_ten(3,3)
            integer::i,jjj
!input s
            call kvec_to_tens(s,sig_ten)
            call kjacobi_eigenvalue(3,sig_ten,dir,sig_pr)
            pnorm_s=norm2(sig_pr);alpha_t=0
            do i=1,3
                  if (sig_pr(i)>=0.0) then
                        s_pt(i)=sig_pr(i);s_pc(i+3)=0.
                  else
                        s_pc(i)=sig_pr(i);s_pt(i+3)=0.
                  end if
                  alpha_t=alpha_t+s_pt(i)*(s_pt(i)+s_pc(i))
            end do
            if (pnorm_s>0) then
                  alpha_t=alpha_t/(pnorm_s**2.)
            end if
            alpha=1.-alpha_t
            call kvec_to_tens(s_pt,sig_ten)
       call ktens_to_vec(matmul(transpose(dir),matmul(sig_ten,dir)),st)
            call kvec_to_tens(s_pc,sig_ten)
       call ktens_to_vec(matmul(transpose(dir),matmul(sig_ten,dir)),sc)
            !check
            return
            end
      subroutine krate_fac(strain_rate,alpha,tol,rate_fac)
            implicit none
            real*8::fc,fc0,strain_rate(6),alpha,tol,rate_fac,
     $ alphas,gammas,betas,rate_t0,rate_c0,rate_t,rate_c,tmp(6),
     $ strain_pr(3),dir(3,3),pmax,pmin,rate,ratio_t,ratio_c,
     $ ratio_c0,ratio_t0,deltas,sig_ten(3,3)
            integer i
            common/kumat/fc,fc0
            !strain_rate,alpha,tol
            alphas=1./(5.+9.*fc/fc0);deltas=1./(1.+8.*fc/fc0)
            gammas=exp((6.156*alphas-2)*log(10.))!check log
            betas=exp((6*deltas-2)*log(10.))
            rate_t0=1.0d-6;rate_c0=-30.0d-6;rate_t=1.;rate_c=1.
            tmp=strain_rate/2.;tmp(1)=tmp(1)*2.;tmp(2)=tmp(2)*2.;
            tmp(3)=tmp(3)*2.
            call kvec_to_tens(tmp,sig_ten)
            call kjacobi_eigenvalue(3,sig_ten,dir,strain_pr)
            pmax=-1.0d-20;pmin=1.0d20
            do i=1,3
                  if (pmax<strain_pr(i))      pmax=strain_pr(i)
                  if (pmin>strain_pr(i))      pmin=strain_pr(i)
            end do
            if ((1-alpha)>tol) then
                  rate =pmax
            else
                  rate =pmin
            end if
            ratio_t=rate/rate_t0;ratio_c=rate/ratio_c0;
            if (rate<30.0d-6) then
                  rate_t=1.
            else if (rate>30.0d-6.and.rate<1.0) then
                  rate_t=ratio_t**deltas
            else
                  rate_t=betas*ratio_t**(0.3333)
            end if
            if (rate>-30.0d-6) then
                  rate_c=1.
            else if (rate>-30.0 .and. rate<-30.0d-6) then
                  rate_c=ratio_c**(1.026*alphas)
            else
                  rate_c=gammas*ratio_c**(0.3333)
            end if
            rate_fac=(1.-alpha)*rate_t+alpha*rate_c
           return
           end
      subroutine kratiopotential(sv,pkp,ratio)
            implicit none
            real*8::ft,fc,pm0,qmuo,df,sv,pkp,ratio,qh1,dqh1_dk,
     $qh2,dqh2_dk,o,pm,pmin_equ_e,eps_old,equ_e_minus,equ_e_new,
     $equ_e_plus,sg_minus(6),sg_plus(6),sg_old(6),dsg(6),
     $sg(6),ro,th,dgdinv(2),pmg,par1,ag,bg,par2,dmg,bg_top,
     $bg_bottom,r,b1
            common/kumat/ft,fc,pm0,qmuo,df
            !input sv,pkp
            ro=0.
            call kcqh1(pkp,int(0),qh1,dqh1_dk)
            call kcqh2(pkp,int(0),qh2,dqh2_dk)
            ag=3*ft*qh2/fc+pm0/2;
            bg_top=qh2/3.*(1.+ft/fc)
            bg_bottom=log(ag)-log(2.*df-1.)-log(3.*qh2+pm0/2.)
     $       +log(df+1.)!log vs exp
            bg=bg_top/bg_bottom
            r=(sv-qh2*ft/3.)/fc/bg
            pmg=ag*bg*fc*exp(r);par1=ro/(fc*sqrt(6.))+sv/fc!log vs exp
            dmg=ag*exp(r);
            par2=(1.-qh1)*par1**2.+sqrt(1.5)*ro/fc
            dgdinv(1)=4.*(1.-qh1)/fc*par2*b1+qh1**2.*pmg/fc
            dgdinv(2)=par2/(sqrt(6.)*fc)*(4.*(1.-qh1)*par1+6)
     $       +pm0*qh1**2/(sqrt(6.)*fc)
            ratio=dgdinv(2)/dgdinv(1)*3*(1-2.*qmuo)/(1+qmuo)
            return
            end
      subroutine kcheckunload(sg,sg_old,eps_old,gtol,pmin_equ_e,
     $       equ_e_new,iunload_flag)
            implicit none
            real*8::sg1(6),sg(6),sg_old(6),eps_old,dinv_dsig_pr(3,3)
     $,sv,sg_minus(6),sg_plus(6),pm,o,p,pn,ro,th,pmin_equ_e,
     $equ_e_plus,equ_e_minus,equ_e_new,sdg(6),equ_e1,dsg(6),
     $gtol
            integer::iunload_flag,i
            !input sg,sg_old,eps_old,gtol
            call khaigh(sg,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_new)
            dsg=sg-sg_old;sg_plus=sg_old+0.01*dsg;
            sg_minus=sg_old+0.99*dsg
            call khaigh(sg_plus,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_plus)
            call khaigh(sg_minus,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_minus)
            iunload_flag=0;pmin_equ_e=eps_old;p=equ_e_plus;
            pm=equ_e_minus;o=eps_old
            if ((p<o.and.pm<pn).and.
     $       (abs(p-o)>gtol/10.and.abs(pm-pn)>gtol/10)) then
                  iunload_flag=1
                  do i=1,100
                        sg1=sg_old+sdg*i/100
                        call khaigh(sg1,sv,ro,th,dinv_dsig_pr)
                        call kequ_e(sv,ro,th,equ_e1)
                        if (equ_e1<=pm) then
                              pm=equ_e1
                        else
                              EXIT
                        end if
                  end do
            end if
            return
            end
      subroutine kderv(sv,ro,th,tkp,icomputeall,dgdinv,dkdl,dfdinv,
     $       ddg_ddinv,dfdk,ddg_dinvdk,ddk_dldinv)
            implicit none
            real*8::sv,ro,th,tkp,dgdinv(2),dkdl,dfdinv(2),
     $       ddg_ddinv(2,2),dfdk,ddg_dinvdk(2),ddk_dldinv(2),
     $ qh1,dqh1_dk,qh2,dqh2_dk,dEquivdg_dstress_dinv(2),dductdinv(2),
     $ ddg_dinv_dk(2),pm0,equivdg_dstress,equivaplentdg_dsg,ecc,
     $ duct_m,dduct_dinv(2),dR_top_dk,dR_bottom_dk,dR_dk,dmQ_dsv,
     $ dmQ_dk,dfdqh2,dfdqh1,fc,ft,df,dEquivaplentdg_dstress_dk,
     $ddk_dldk,ddg_dsvdro,ddg_drodsv,ddg_ddsv,r,par1,par,pmQ,
     $ddg_ddro,dbl_dsv,dbl_dro,dbg_top_dk,dag_dk,dbg_bottom_dk,
     $dAl_dyieldHard,dAl_dsv,dAl_dro,Bl,bg_top,bg_bottom,
     $bg,Al,ag,pi,dg_dinv(2),rcos,R_top,r_bottom,ddg_dsv_dro,
     $dbg_dk,par2
            integer::icomputeall
            common/kumat/fc,ft,pm0,df
            !sv,ro,th,tkp,icomputeall
            pi=4.*atan(1.d0)
            call kcqh1(tkp,int(0),qh1,dqh1_dk);
            call kcqh2(tkp,int(0),qh2,dqh2_dk)
            !dgdinv 
            ag=ft*qh2*3./fc+pm0/2.;bg_top=qh2/3.*(1.+ft/fc)
            bg_bottom=log(ag)-log(2.*df-1.)-log(3.*qh2+pm0/2.)
     $       +log(df+1.)!log vs exp
            bg=bg_top/bg_bottom
            r=(sv-ft/3.*qh2)/fc/bg;pmQ=ag*exp(r);
            Bl=sv/fc+ro/(fc*sqrt(6.));
            Al=(1.-qh1)*Bl**2.+sqrt(1.5)*ro/fc;
            dgdinv(1)=4.*(1.-qh1)/fc*Al*Bl+qh1**2*pmQ/fc;
            dgdinv(2)=Al/(sqrt(6.)*fc)*(4.*(1.-qh1)*Bl+6.)+pm0
     $       *qh1**2/(sqrt(6.)*fc);
            !dkdl
            equivaplentdg_dsg=sqrt(1./3.*dgdinv(1)**2+dgdinv(2)**2);
            call kductility(sv,th,1,duct_m,dduct_dinv)
            dkdl=equivaplentdg_dsg/duct_m;
            !dfdinv
            if (icomputeall.eq.1) then
                  par1=(1.-ecc**2);par2=(2.*ecc-1.)
                  rcos=(4.*par1*cos(th)**2.+par2**2.)
                  rcos=rcos/(2.*par1*cos(th)+par2*sqrt(4.*par1
     $             *cos(th)**2.+5.*ecc**2-4.*ecc))
                  dfdinv(1)=4.*(1.-qh1)/fc*Al*Bl+qh2*qh1**2*pm0/fc
                  dfdinv(2)=Al/(sqrt(6.)*fc)*(4.*(1.-qh1)*Bl+6.)
     $             +rcos*pm0*qh2*qh1**2/(sqrt(6.)*fc);
                  !ddgddinv 
                  dmQ_dsv=ag/(bg*fc)*exp(r);dAl_dsv=2.*(1.-qh1)*Bl/fc
                  dBl_dsv=1./fc;dAl_dro=2.*(1.-qh1)*Bl/(fc*sqrt(6.))
     $             +sqrt(1.5)/fc
                  dBl_dro=1./(fc*sqrt(6.))
                  ddg_ddsv=4.*(1.-qh1)/fc*(dAl_dsv*Bl+Al*dBl_dsv)
     $             +qh1**2*dmQ_dsv/fc
                  ddg_ddro=dAl_dro/(sqrt(6.)*fc)*(4.*(1.-qh1)*Bl+6.)
     $             +Al*dBl_dro*4.*(1.-qh1)/(sqrt(6.)*fc)
                  ddg_dsvdro=4.*(1.-qh1)/fc*(dAl_dro*Bl+Al*dBl_dro)
                  ddg_drodsv=dAl_dsv/(sqrt(6.)*fc)*(4.*(1.-qh1)*Bl+6.)
     $             +Al/(sqrt(6.)*fc)*(4.*(1.-qh1)*dBl_dsv)
                  ddg_ddinv(1,1)=ddg_ddsv;ddg_ddinv(1,2)=ddg_dsv_dro
                  ddg_ddinv(2,1)=ddg_drodsv;ddg_ddinv(2,2)=ddg_ddro
                  !dfdk
                  dfdqh1=-2.*Al*(Bl**2.)+2.*qh1*qh2*pm0*(sv/fc+ro*rcos
     $             /(sqrt(6.)*fc))-2.*qh1*(qh2**2)
                  dfdqh2=(qh1**2.)*pm0*(sv/fc+ro*rcos/(sqrt(6.)*fc))-
     $             2.*qh2*(qh1**2)
                  dfdk=dqh1_dk*dfdqh1+dqh2_dk*dfdqh2
                  if(dfdk>0.)    dfdk=0.
                  !ddgdinvdk
                  dag_dk=dqh2_dk*3.*ft/fc;dbg_top_dk=dqh2_dk/3.
                  dbg_bottom_dk=-3.*dqh2_dk/(3.*qh2+pm0/2.)
                  dbg_dk=(dbg_top_dk*bg_bottom-bg_top*dbg_bottom_dk)
     $             /(bg_bottom**2.)
                  R_top=(sv-ft/3.*qh2);R_bottom=fc*bg
                  dR_top_dk=-ft/3.*dqh2_dk;dR_bottom_dk=fc*dbg_dk
                  dR_dk=(dR_top_dk*R_bottom-R_top*dR_bottom_dk)
     $             /(R_bottom**2.)
                  dmQ_dk=dag_dk*exp(r)+ag*dR_dk*exp(r);
                  dAl_dyieldHard=-(Bl**2.)
                  ddg_dinvdk(1)=(-4.*Al*Bl/fc+4.*(1-qh1)/fc
     $             *dAl_dyieldHard*Bl)*dqh1_dk
                  ddg_dinvdk(1)=ddg_dinvdk(1)+dqh1_dk*2.*qh1*pmQ
     $             /fc+qh1*dmQ_dk/fc
                  par=dAl_dyieldHard/(sqrt(6.)*fc)*(4.*(1-qh1)*Bl+6.)
                  par1=-4.*Al/(sqrt(6.)*fc)*Bl+pm0/(sqrt(6.)*fc)
                  ddg_dinvdk(2)=(par+par1)*2.*qh1*dqh1_dk;
                  !ddk_dldk
                  par1=2./3.*dg_dinv(1)/equivaplentdg_dsg*ddg_dinv_dk(1)
                  dEquivaplentdg_dstress_dk=(par1+2.*dg_dinv(2)
     $             /equivaplentdg_dsg*ddg_dinvdk(2))/2
                  ddk_dldk= dEquivaplentdg_dstress_dk/duct_m
                  !ddKappadDeltaLambdadInv
                  dEquivdg_dstress_dinv(1)=2./3.*dg_dinv(1)
     $             *ddg_ddinv(1,1)+2.*dg_dinv(2)*ddg_ddinv(2,1)
                  dEquivdg_dstress_dinv(1)=dEquivdg_dstress_dinv(1)
     $             /(2.*equivdg_dstress)
                  dEquivdg_dstress_dinv(2)=2./3.*dg_dinv(1)
     $             *ddg_ddinv(1,2)+2.*dg_dinv(2)*ddg_ddinv(2,2)
                  dEquivdg_dstress_dinv(2)=dEquivdg_dstress_dinv(2)
     $             /(2.*equivdg_dstress)
                  ddk_dldinv(1)=(dEquivdg_dstress_dinv(1)
     $             *duct_m-equivdg_dstress*dductdinv(1))/(duct_m**2)
                  ddk_dldinv(2)=(dEquivdg_dstress_dinv(2)
     $             *duct_m-equivdg_dstress*dductdinv(2))/(duct_m**2);            
            else 
                  dfdinv=0.;ddg_ddinv=0.;dfdk=0.;ddg_dinvdk=0.;
                  ddk_dldinv=0.
            end if
            return
            end
      subroutine kdamage(wc_old,wt_old,strain_rate,rate_fac,
     $       alpha,eps_t,eps_c,pkdt_old,pkdt1,pkdt2,
     $       pkdc_old,pkdc1,pkdc2,sg_ekff,tkp,pnorm_inc_e_p,
     $       plen,sg_old,alpha_old,eps_old,wc,wt,
     $       eps_t1,eps_c1,pkdt_new,pkdt1t,
     $       pkdt2t,pkdc_new,pkdc1t,pkdc2t,eps_new)
            implicit none
            real*8::wc_old,wt_old,strain_rate(6),
     $ alpha,eps_t,eps_c,pkdt_old,pkdt1,pkdt2,pkdt,dqh2_dk,
     $ pkdc_old,pkdc1,pkdc2,sg_ekff(6),tkp,pnorm_inc_e_p,
     $ plen,sg_old(6),alpha_old,eps_old,wc,wt,fac2,gtol,
     $ rate_fac,eps_t1,eps_c1,pkdt_new,pmin_equ_e,t_equ_e,
     $pkdc_new,eps_new,xs,wf1,wf,tol,theta_el,e0,fac,ft,
     $t_rate_fac,sv_el,rs1,ro_el,residualDerivative,df,
     $residual,qh2,pari,pkdt1t,pkdt2t,pkdc1t,pkdc2t,as,d_pkdc,
     $d_pkdt,dinv_dsig_pr(3,3)
            integer::istep1flag,iunload_flag,itype,istrrateflg,
     $ inewton_iter      
            common/kumat/e0,ft,as,df,iunload_flag,istrrateflg
            !input wc_old,wt_old,strain_rate,rate_fac,alpha,eps_t,eps_c,pkdt_old,pkdt1,pkdt2,
            !pkdc_old,pkdc1,pkdc2,sg_ekff,tkp,pnorm_inc_e_p,plen,sg_old,alpha_old,eps_old
            tol=gtol*10
            if (rate_fac.eq.0.) then
                  istep1flag=1;rate_fac=1.
            else 
                  istep1flag=0
            end if
            call kcheckunload(sg_ekff,sg_old,eps_old,gtol,
     $       pmin_equ_e,t_equ_e,iunload_flag)
            eps_new=t_equ_e;
            if(istrrateflg.eq.1.and.wc_old.eq.0.0.and.wt_old.eq.0.0)then
                  call krate_fac(strain_rate,alpha,tol,t_rate_fac)
            else 
                  t_rate_fac=rate_fac
            end if
            !note rate factor is calculated only once at the onset of kdamage
            eps_t1=eps_t+(t_equ_e-eps_old)/t_rate_fac
            if(iunload_flag.eq.0) then
                  eps_c1=eps_c+(t_equ_e-eps_old)*alpha/t_rate_fac
            else 
                  eps_c1=eps_c+((pmin_equ_e-eps_old)*alpha_old
     $              +(t_equ_e-pmin_equ_e)*alpha)/t_rate_fac
            end if
            if((eps_c1>e0.or.eps_t1>e0).and.
     $       (wc_old.eq.0.0.and.wt_old.eq.0.0)
     $ .and.(istrrateflg.eq.1.and.istep1flag.ne.1)) then
                  eps_t1=eps_t+(t_equ_e-eps_old)/rate_fac
                  if(iunload_flag.eq.0) then
                        eps_c1=eps_c+(t_equ_e-eps_old)*alpha/rate_fac
                  else 
                        eps_c1=eps_c+((pmin_equ_e-eps_old)*alpha_old
     $                     +(t_equ_e-pmin_equ_e)*alpha)/rate_fac
                  end if
            else
                  rate_fac=t_rate_fac
            end if
            
            call khaigh(sg_ekff,sv_el,ro_el,theta_el,dinv_dsig_pr)
            if(sv_el<0.0) then
                  rs1=-6.0**0.5*sv_el/max(ro_el,1.0d-16)
            else 
                  rs1=0.0
            end if
            xs=1.+(as-1.)*rs1;d_pkdt=(eps_t1-pkdt_old);
            d_pkdc=(eps_c1-pkdc_old);wt=0.;wc=0.;
            if(d_pkdt>=-e0*tol) then
                  if(eps_t1<e0*(1.-tol)) then
                        fac=0.
                  else if (eps_t1>e0*(1.-tol).and.pkdt_old
     $                  <e0*(1.-tol)) then
                        fac=(1.-(e0-pkdt_old)/(eps_t1-pkdt_old))
                  else 
                        fac=1.
                  end if
                  pkdt1t=pkdt1+pnorm_inc_e_p*fac/xs/rate_fac;
                  pkdt2t=pkdt2+d_pkdt/xs
                  pkdt_new=eps_t1
                  call kdamaget(pkdt_new,pkdt1t,pkdt2t,plen,wt_old,wt)
            else if(d_pkdc>=-e0*tol) then
                  if (eps_c1<e0*(1.-tol)) then
                        fac=0.
                  elseif (eps_c1>e0*(1.-tol).and.pkdc_old
     $                        <e0*(1.-tol)) then
                        fac=(1.-(e0-pkdc_old)/(eps_c1-pkdc_old))
                  else 
                        fac=1.
                  end if
                  call kcqh2(tkp,int(0),qh2,dqh2_dk)
                  fac2=ft*qh2*sqrt(0.66666)
     $             /max(ro_el,1.0d-16)/sqrt(1.+2.*(df**2.));
                  pkdc1t=pkdc1+pnorm_inc_e_p*fac*fac2*alpha/xs/rate_fac;
                  pkdc2t=pkdc2+d_pkdc/xs;
                  pkdc_new=eps_c1
                  call kdamagec(pkdc_new,pkdc1t,pkdc2t,wc_old,wc)
            end if
            return
            end
      subroutine kdamaget(pkdt,pkdt1,pkdt2,plen,wt_old,wt)
            implicit none
            real*8::pkdt,pkdt1,pkdt2,plen,wt_old,wt,Em,wf,ft,ft1,
     $ wf1,pari,residual,residualDerivative,ytol,gtol,e0
            integer::itype,iter,inewton_iter
            common/kumat/Em,wf,ft,ft1,wf1,e0
            !pkdt,pkdt1,pkdt2,plen,wt_old output wt
            inewton_iter=100;ytol=gtol*10.;
            if (pkdt>e0*(1.-ytol)) then
                  if (itype.eq.0) then
                        wt=(Em*pkdt*wf-ft*wf+ft*pkdt1*plen)
     $                        /(Em*pkdt*wf-ft*plen*pkdt2)
                  else if (itype.eq.1) then
                        wt=(Em*pkdt*wf1-ft*wf1-(ft1-ft)*pkdt1*plen)
     $                        /(Em*pkdt*wf1+(ft1-ft)*plen*pkdt2)
                        pari=plen*pkdt1+plen*wt*pkdt2
                        if (pari>wf1.and.pari<wf) then
                              wt=(Em*pkdt*(wf-wf1)-ft1*(wf-wf1)+ft1
     $                              *pkdt1*plen-ft1*wf1)
                              wt=wt/(Em*pkdt*(wf-wf1)-ft1*plen*pkdt2)
                              pari=plen*pkdt1+plen*wt*pkdt2
                        else if (pari>wf) then
                              wt=1.
                        end if
                  else if (itype.eq.2) then
                        !Exponential: Iterative solution with N-R procedure
                        wt=1.;residual=0.;residualDerivative=0.;iter=0;
                        pari=1.;
                        do while (pari.eq.1.)
                              iter=iter+1;residual=(1.-wt)*Em*pkdt-ft
     $                              *exp(-plen*(wt*pkdt2+pkdt1)/wf)
                              residualDerivative=-Em*pkdt+ft*plen*pkdt2/wf
     $                              *exp(-plen*(wt*pkdt2+pkdt1)/wf)
                              wt=wt-residual/residualDerivative
                              !if(iter>inewton_iter);  
                              !disp('* Algorithm for tensile kdamage-No convergence reached after 100 iterations *')
                              !error('erorrre');end
                              if(abs(residual/ft)<1.0d-8) pari=0.0
                        end do
                  else
                        wt=0.
                  end if
                  if(wt>1.) wt=1.
                  if(wt<0.0.or.wt<wt_old) wt=wt_old
            else 
                  wt=0.
            end if
            return
            end
      subroutine kdamagec(pkdc,pkdc1,pkdc2,wc_old,wc)
            implicit none
            real*8::pkdc,pkdc1,pkdc2,wc_old,wc,pari,gtol,ytol,tol,
     $ residual,residualDerivative,dResdw,errorOld,efc,ft,e0,em
            integer::itype,iter,nite,inewton_iter,isotropic
            common/kumat/em,ft,efc,e0,isotropic
            !pkdc,pkdc1,pkdc2,wc_old
            if (isotropic.eq.1) then
                  wc=0.;pkdc1=0.;pkdc2=0.;pkdc=0.
            else
                  inewton_iter=200;tol=gtol;ytol=gtol*10.;nite=0;
                  residual=0.;dResdw=0.;
                  if (pkdc>e0*(1.-ytol)) then
                        do while(nite<inewton_iter)
                              nite=nite+1;residual=(1.-wc)*em*pkdc-ft
     $                              *exp(-(pkdc1+wc*pkdc2)/efc)
                              dResdw=-em*pkdc+ft*pkdc2/efc*exp(-(pkdc1+
     $                              wc*pkdc2)/efc);
                              wc=wc-residual/dResdw
                              errorOld = residual/ft
                              if(wc<0.0) then
                                    wc=0.
                                    exit
                              end if
                              if(nite.eq.inewton_iter) then
                                    if(residual<0.) then
                                          wc=wc_old
                                          exit
                                    else !disp('*** Algorithm for compressive kdamage-No convergence reached.');error('errerrr')
                                    end if
                              end if
                              if(abs(residual/ft)<tol) then
                                    exit
                              end if
                        end do
                  else 
                        wc=0.
                  end if
            end if
            if(wc>1.) wc=1.
            if(wc<0.0.or.wc<wc_old) wc=wc_old
            return
            end
      subroutine kcheckvertex(sv_tr,tkp,apex_sg,irtype)
            implicit none
            real*8::sv_tr,sv,tkp,apex_sg(6),qh22,dqh2_dk,fc,pm0,qh2
            integer::irtype
            common/kumat/fc,pm0
            !sv_tr,tkp
            if(sv_tr>0.) then
                  irtype=1
                  if(tkp<1) then
                        apex_sg=0.
                  else 
                        call kcqh2(tkp,int(0),qh2,dqh2_dk)
                        apex_sg=qh2*fc/pm0
                  end if
            else if (sv<0. .and. tkp<1.) then
                  irtype=2;apex_sg=0.
            else
                  irtype=0;apex_sg=0.
            end if
            return
            end
      subroutine kvertexreturn(sg,apex_sg,tkp,irtype,iconvrg,sg_ekff)
            implicit none
            real*8::pkp,sv,sv2,yvalue_mid,gtol,pkp0,tkp,sg(6),apex_sg,
     $sg_ekff(6),dsv,pari,ratioPotent1, ratiotrial,ro,sgAnswer(3),
     $stress(6),sv_mid,svAnswer,theta,tk,ytol,yvalue,dinv_dsig_pr(3,3)
            integer::irtype,iconvrg,j,maxiter
            !sg,apex_sg,tkp,irtype,iconvrg
            ytol=gtol;yvalue=0.;yvalue_mid=0.;sv2=0.;pkp0=tkp;pkp=tkp
            maxiter=250
            call khaigh(sg,sv,ro,theta,dinv_dsig_pr)
            sv2 = apex_sg;
            call kpp(pkp0,sv,ro,sv,tkp);
            call kff(sv,0.0d0,0.0d0,tkp,yvalue)
            call kpp(pkp0,sv,ro,sv2,tkp);
            call kff(sv2,0.0d0,0.0d0,tk,yvalue_mid)
            pari=0.
            if(yvalue*yvalue_mid>=0) then
                  iconvrg=1;irtype=0;
            else
                  if(yvalue<0.) then
                        dsv=sv2-sv;svAnswer=sv2
                  else 
                        dsv=sv-sv2;svAnswer=sv2
                  end if
                  do j=1,maxiter
                        dsv=0.5*dsv
                        sv_mid=svAnswer+dsv
                        call kpp(pkp0,sv,ro,sv_mid,tkp)
                        call kff(sv_mid,0.0d0,0.0d0,tkp,yvalue_mid)
                        if(yvalue_mid<=.0d0) then
                              svAnswer=sv_mid
                        end if
                        if (abs(yvalue_mid)<ytol.and.
     $                   yvalue_mid<=.0d0) then 
                              call kratiopotential(svAnswer,tkp,
     $                              ratioPotent1)
                              ratiotrial=ro/(sv-svAnswer)
                              if((ratioPotent1>=ratiotrial.and.irtype
     $                           .eq.1).or.(ratioPotent1<=ratiotrial
     $                           .and.irtype.eq.2)) then
                                    exit
                              else
                                    iconvrg=1.;irtype=0.;pari=1.
                             end if
                        end if
                  end do
                  if (pari.eq.0.) then
                        stress(1:3)=sgAnswer;stress(4:6)=0.;pkp=tkp
                        iconvrg=0.
                  end if
            end if
            return
            end
      subroutine kpp(pkp_old,sv1,dro,sv2,pkp)
            implicit none
            real*8::pkp_old,pkm3,em,qmuo,dro,sv1,sv2,
     $ equ_d_e_p,duct_m,dduct_dinv(2),pkp,pi,gm2
            common/kumat/em,qmuo
            pkm3=em/(1.-2.*qmuo);gm2=em/(1.+qmuo)
            equ_d_e_p=sqrt((((sv1-sv2)/pkm3)**2+(dro/gm2)**2))
            pi=4.*atan(1.d0)
            call kductility(sv2,pi/3.,int(0),duct_m,dduct_dinv)
            pkp=pkp_old+equ_d_e_p/duct_m
            return
            end
      subroutine kductility(sv,th,icomput_deriv,duct_m,dduct_dinv)
            implicit none
            real*8::par1,sv,th,duct_m,dduct_dinv(2),fc,ah,bh,ch,
     $dh,eh,x,fh
            common/kumat/fc,ah,bh,ch,dh
            integer::icomput_deriv
            !input sv,th,icomput_deriv (integer important)
            par1=(2*cos(th))**2.;x=-1.*(sv+fc/3.)/fc;
            if (x<0.) then
                  eh=bh-dh;fh=eh*ch/(ah-bh)
                  duct_m=(eh*exp(x/fh)+dh)/par1
                  if (icomput_deriv.eq.1) then
                        dduct_dinv(1)=eh/fh*exp(x/fh)/par1*(-1.)/fc
                  else
                        dduct_dinv(1)=0.
                  end if
            else
                  duct_m=(ah-(ah-bh)*exp(-x/fh))/par1
                  if (icomput_deriv.eq.1) then
                        dduct_dinv(1)=(bh-ah)/ch*exp(-x/ch)/par1/fc
                  else
                        dduct_dinv(1)=0.
                  end if
            end if
            return
            end  
      subroutine kequ_e(sv,ro,th,equ_e)
            implicit none
            real*8::sv,ro,th,equ_e,ecc,pm0,fc,e0,par1,par2,
     $ rcos,par_p,par_q    
            common/kumat/ecc,pm0,fc,e0
            par1=(1.-ecc**2);par2=(2.*ecc-1.)
            rcos=(4.*par1*cos(th)**2.+par2**2.)
            rcos=rcos/(2.*par1*cos(th)+par2*sqrt(4.
     $       *par1*cos(th)**2.+5.*ecc**2-4.*ecc))
            par_p=-pm0*(ro*rcos/sqrt(6.)/fc+sv/fc);
            par_q=-1.5*ro**2./fc**2.
            equ_e=e0*(-0.5*par_p+sqrt(par_p**2/4-par_q))
            if (equ_e<0.) equ_e=0
            return
            end
      subroutine kvec_to_tens(v,t)
            implicit none
            real*8,intent(in)      ::v(6)
            real*8            ::  t(3,3)
            t(1,1)=v(1)
            t(2,2)=v(2)
            t(3,3)=v(3)
            t(1,2)=v(4)
            t(1,3)=v(5)
            t(3,2)=v(6)
            t(2,1)=v(4)
            t(3,1)=v(5)
            t(2,3)=v(6)
            end 
      subroutine ktens_to_vec(t,v)
            implicit none
            real*8,intent(in) ::t(3,3)
            real*8            ::v(6)
            v(1)=t(1,1)
            v(2)=t(2,2)
            v(3)=t(3,3)
            v(4)=t(1,2)
            v(5)=t(1,3)
            v(6)=t(2,3)
            end
!
      subroutine kjacobi_eigenvalue (n, a,v,d)

            !*****************************************************************************80
            !
            !! kjacobi_eigenvalue carries out the Jacobi eigenvalue iteration.
            !
            !  Discussion:
            !
            !    This function computes the eigenvalues and eigenvectors of a
            !    real symmetric matrix, using Rutishauser's modfications of the classical
            !    Jacobi rotation method with threshold pivoting. 
            !
            !  Licensing:
            !
            !    This code is distributed under the GNU LGPL license.
            !
            !  Modified:
            !
            !    17 September 2013
            !
            !  Author:
            !
            !    FORTRAN90 version by John Burkardt
            !
            !  Parameters:
            !
            !    Input, integer ( kind = 4 ) N, the order of the matrix.
            !
            !    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
            !    and symmetric.
            !
            !    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
            !
            !    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
            !
            !    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
            !
            !    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
            !
            !    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
            !
      implicit none
            
      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n,n)
      real ( kind = 8 ) bw(n)
      real ( kind = 8 ) c
      real ( kind = 8 ) d(n)
      real ( kind = 8 ) g
      real ( kind = 8 ) gapq
      real ( kind = 8 ) h
      integer ( kind = 4 ) i
      integer ( kind = 4 ) it_max
      integer ( kind = 4 ) it_num
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) m
      integer ( kind = 4 ) p
      integer ( kind = 4 ) q
      integer ( kind = 4 ) rot_num
      real ( kind = 8 ) s
      real ( kind = 8 ) t
      real ( kind = 8 ) tau
      real ( kind = 8 ) term
      real ( kind = 8 ) termp
      real ( kind = 8 ) termq
      real ( kind = 8 ) theta
      real ( kind = 8 ) thresh
      real ( kind = 8 ) v(n,n)
      real ( kind = 8 ) w(n)
      real ( kind = 8 ) zw(n)

      do j = 1, n
            do i = 1, n
                  v(i,j) = 0.0D+00
            end do
            v(j,j) = 1.0D+00
      end do

      do i = 1, n
            d(i) = a(i,i)
      end do

      bw(1:n) = d(1:n)
      zw(1:n) = 0.0D+00
      it_num = 0
      rot_num = 0

      do while ( it_num < it_max )
            it_num = it_num + 1
      !
      !  The convergence threshold is based on the size of the elements in
      !  the strict upper triangle of the matrix.
      !
            thresh = 0.0D+00
            do j = 1, n
                  do i = 1, j - 1
                        thresh = thresh + a(i,j) ** 2
                  end do
            end do

            thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

            if ( thresh == 0.0D+00 ) then
                  exit 
            end if

            do p = 1, n
                  do q = p + 1, n

                        gapq = 10.0D+00 * abs ( a(p,q) )
                        termp = gapq + abs ( d(p) )
                        termq = gapq + abs ( d(q) )
      !
      !  Annihilate tiny offdiagonal elements.
      !
                        if ( 4 < it_num .and.
     $                   termp .eq. abs ( d(p) ) .and.
     $                   termq .eq. abs ( d(q) ) ) then

                              a(p,q) = 0.0D+00
      !
      !  Otherwise, apply a rotation.
      !
                        else if ( thresh <= abs ( a(p,q) ) ) then

                              h = d(q) - d(p)
                              term = abs ( h ) + gapq

                              if ( term == abs ( h ) ) then
                                    t = a(p,q) / h
                              else
                                    theta = 0.5D+00 * h / a(p,q)
                                    t = 1.0D+00 / ( abs ( theta) 
     $                         + sqrt ( 1.0D+00 + theta * theta ) )
                                    if ( theta < 0.0D+00 ) then 
                                          t = - t
                                    end if
                              end if
            
                              c = 1.0D+00 / sqrt ( 1.0D+00 +t*t)
                              s = t * c
                              tau = s / ( 1.0D+00 + c )
                              h = t * a(p,q)
            !
            !  Accumulate corrections to diagonal elements.
            !
                              zw(p) = zw(p) - h                  
                              zw(q) = zw(q) + h
                              d(p) = d(p) - h
                              d(q) = d(q) + h
            
                              a(p,q) = 0.0D+00
            !
            !  Rotate, using information from the upper triangle of A only.
            !
                              do j = 1, p - 1
                                    g = a(j,p)
                                    h = a(j,q)
                                    a(j,p) = g -s*( h + g * tau )
                                    a(j,q) = h +s*( g - h * tau )
                              end do
            
                              do j = p + 1, q - 1
                                    g = a(p,j)
                                    h = a(j,q)
                                    a(p,j) = g - s*(h+ g * tau )
                                    a(j,q) = h + s*(g- h * tau )
                              end do
            
                              do j = q + 1, n
                                    g = a(p,j)
                                    h = a(q,j)
                                    a(p,j) = g - s*(h+ g * tau )
                                    a(q,j) = h + s*(g- h * tau )
                              end do
            !
            !  Accumulate information in the eigenvector matrix.
            !
                              do j = 1, n
                                    g = v(j,p)
                                    h = v(j,q)
                                    v(j,p) = g - s * (h+g*tau )
                                    v(j,q) = h + s * (g-h*tau )
                              end do
            
                              rot_num = rot_num + 1
            
                        end if
            
                  end do
            end do

            bw(1:n) = bw(1:n) + zw(1:n)
            d(1:n) = bw(1:n)
            zw(1:n) = 0.0D+00

      end do
      !
      !  Restore upper triangle of input matrix.
      !
      do j = 1, n
            do i = 1, j - 1
                  a(i,j) = a(j,i)
            end do
      end do
      !
      !  Ascending sort the eigenvalues and eigenvectors.
      !
      do k = 1, n - 1

            m = k

            do l = k + 1, n
                  if ( d(l) < d(m) ) then
                        m = l
                  end if
            end do

            if ( m /= k ) then

                  t    = d(m)
                  d(m) = d(k)
                  d(k) = t
            
                  w(1:n)   = v(1:n,m)
                  v(1:n,m) = v(1:n,k)
                  v(1:n,k) = w(1:n)

            end if

      end do
      return
      end
      subroutine Kdet33 (A,det3)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
            real*8::det3
            det3 =A(1,1)*A(2,2)*A(3,3)- A(1,1)*A(2,3)*A(3,2)
            det3=det3-A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
            det3=det3+A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
            RETURN
            end
      subroutine Kdet44 (A,det4)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
            real*8::det4
            det4=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))
     $+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)
     $-A(3,3)*A(4,2)))- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))
     $+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)
     $-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))
     $+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)
     $-A(3,2)*A(4,1)))- A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
     $+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)
     $-A(3,2)*A(4,1)))
            RETURN
            end
      subroutine kmatinv4(A,B)
            !! Performs a direct calculation of the kinverse of a 4×4 matrix.
            real*8, intent(in) :: A(4,4)   !! Matrix
            real*8             :: B(4,4)   !! kinverse matrix
            real*8             :: detinv

            ! Calculate the kinverse determinant of the matrix
            detinv =1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))
     $+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)
     $-A(3,3)*A(4,2)))- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))
     $+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)
     $-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))
     $+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)
     $-A(3,2)*A(4,1)))- A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
     $+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)
     $-A(3,2)*A(4,1))))

            ! Calculate the kinverse of the matrix
            B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)
     $       -A(3,4)*A(4,3))+A(2,3)
     $*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)
     $*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
            B(2,1) = detinv*(A(2,1)*(A(3,4)
     $       *A(4,3)-A(3,3)*A(4,4))+A(2,3)
     $*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)
     $*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))       
            B(3,1) = detinv*(A(2,1)*(A(3,2)
     $       *A(4,4)-A(3,4)*A(4,2))+A(2,2)
     $*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)
     $*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))       
            B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)
     $*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)
     $*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
            B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)
     $*(A(3,2)*A(4,4)-A(3,4)*A(4,2))
     $+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))       
            B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)
     $*(A(3,4)*A(4,1)-A(3,1)*A(4,4))
     $+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))       
            B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)
     $*(A(3,1)*A(4,4)-A(3,4)*A(4,1))
     $+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))       
            B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)
     $*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
     $+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))       
            B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)
     $*(A(2,4)*A(4,2)-A(2,2)*A(4,4))
     $+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))       
            B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)
     $*(A(2,1)*A(4,4)-A(2,4)*A(4,1))
     $+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))       
            B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)
     $*(A(2,4)*A(4,1)-A(2,1)*A(4,4))
     $+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))       
            B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)
     $*(A(2,1)*A(4,3)-A(2,3)*A(4,1))
     $+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))       
            B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)
     $*(A(2,2)*A(3,4)-A(2,4)*A(3,2))
     $+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))       
            B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)
     $*(A(2,4)*A(3,1)-A(2,1)*A(3,4))
     $+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))       
            B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)
     $*(A(2,1)*A(3,4)-A(2,4)*A(3,1))
     $+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))       
            B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)
     $*(A(2,3)*A(3,1)-A(2,1)*A(3,3))
     $+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))       
            end            
      !integer:: mu,tr,kvec_to_tens,ktens_to_vec
      subroutine kinverse(a,c,n)
      !============================================================
      ! kinverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - kinverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

      ! step 1: forward elimination
      do k=1, n-1
            do i=k+1,n
                  coeff=a(i,k)/a(k,k)
                  L(i,k) = coeff
                  do j=k+1,n
                        a(i,j) = a(i,j)-coeff*a(k,j)
                  end do
            end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
            L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
            do i=1,j
                  U(i,j) = a(i,j)
            end do
      end do

      ! Step 3: compute columns of the kinverse matrix C
      do k=1,n
            b(k)=1.0
            d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
            do i=2,n
                  d(i)=b(i)
                  do j=1,i-1
                        d(i) = d(i) - L(i,j)*d(j)
                  end do
            end do
      ! Step 3b: Solve Ux=d using the back substitution
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                  x(i) = d(i)
                  do j=n,i+1,-1
                        x(i)=x(i)-U(i,j)*x(j)
                  end do
                  x(i) = x(i)/u(i,i)
            end do
      ! Step 3c: fill the solutions x(n) into column k of C
            do i=1,n
                  c(i,k) = x(i)
            end do
            b(k)=0.0
      end do
      return
      end
      subroutine gauss_2(a,b,x,n)
            !===========================================================
            ! Solutions to a system of linear equations A*x=b
            ! Method: Gauss elimination (with scaling and pivoting)
            ! Alex G. (November 2009)
            !-----------------------------------------------------------
            ! input ...
            ! a(n,n) - array of coefficients for matrix A
            ! b(n) - array of the right hand coefficients b
            ! n - number of equations (size of matrix A)
            ! output ...
            ! x(n) - solutions
            ! coments ...
            ! the original arrays a(n,n) and b(n) will be destroyed
            ! during the calculation
            !===========================================================
      implicit none
      integer n
      double precision a(n,n), b(n), x(n)
      double precision s(n)
      double precision c, pivot, store
      integer i, j, k, l
      ! step 1: begin forward elimination
      do k=1, n-1
      ! step 2: "scaling"
      ! s(i) will have the largest element from row i
            do i=k,n ! loop over rows
                  s(i) = 0.0
                  do j=k,n ! loop over elements of row i
                        s(i) = max(s(i),abs(a(i,j)))
                  end do
            end do
      ! step 3: "pivoting 1"
      ! find a row with the largest pivoting element
            pivot = abs(a(k,k)/s(k))
            l=k
            do j=k+1,n
                  if(abs(a(j,k)/s(j)) > pivot) then
                        pivot = abs(a(j,k)/s(j))
                        l=j
                  end if
            end do
      ! Check if the system has a sigular matrix
            if(pivot == 0.0) then
                  write(*,*) 'The matrix is sigular '
                  return
            end if
      ! step 4: "pivoting 2" interchange rows k and l (if needed)
            if (l /= k) then
                  do j=k,n
                        store = a(k,j)
                        a(k,j) = a(l,j)
                        a(l,j) = store
                  end do
                  store = b(k)
                  b(k) = b(l)
                  b(l) = store
            end if
      ! step 5: the elimination (after scaling and pivoting)
            do i=k+1,n
                  c=a(i,k)/a(k,k)
                  a(i,k) = 0.0
                  b(i)=b(i)- c*b(k)
                  do j=k+1,n
                        a(i,j) = a(i,j)-c*a(k,j)
                  end do
            end do
      end do
      ! step 6: back substiturion
      x(n) = b(n)/a(n,n)
      do i=n-1,1,-1
            c=0.0
            do j=i+1,n
                  c= c + a(i,j)*x(j)
            end do
            x(i) = (b(i)- c)/a(i,i)
      end do
      return
      end
      
