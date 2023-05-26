!  codey.f90 
      PROGRAM HELLO
      implicit double precision(a-h,o-z)
c      implicit real(a-h,o-z)
      parameter (nprecd=2)
      CHARACTER*80 CMNAME
      DIMENSION sig(6),statev(26),PROPS(19),
     1 DDSDDE(6,6),
     2 DDSDDT(6),DRPLDE(6),
     3 stran(6),dstran(6),TIME(2),PREDEF(1),DPRED(1),
     4 COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4),sig_ten(3,3),sg(6),dir(3,3),sg_p(3),
     6 sig1(1),sig2(1),sig3(1),sig4(1),sig5(1),sig6(1),
     7 d1(1),d2(1),d3(1),d4(1),d5(1),d6(1),hsvs(1,27),cm(24),
     8 idelev(1),elsizv(1),qmat(1,3,3),crv(1,2,1),nnpcrv(1),cma(1),
     9 dt1siz(1),temps(1),epsps(1)
      logical failels(1),reject
      character*5 etype
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %tol=1.0d-6,tol2=10.0d0,sx=6.0d0,half=0.5d0)
      stran(1)=0.d0;dstran=0.d0;dstran(1)=-5.0d-3;dstran(2)=1.5d-3
c      stran(1)=0.d0;dstran=0.d0;dstran(1)=-2.0d-3;dstran(2:3)=4.0d-4
      dstran(3)=2.5d-3;dstran(4)=7.0d-4;
      dstran=dstran*1.d0
c      stran(1)=0.d0;dstran=0.d0;dstran(1)=-1.0d-3;dstran(2:3)=2.0d-4
c      stran(1)=0.;dstran=0.;dstran(1)=5.0d-3;dstran(2:3)=-2.0d-3

      
      WRITE(*,*) 'Hello, world!'
      sg(1)=3.0d0/2.0d0;sg(2)=11.0d0/4.0d0;sg(3)=11.0d0/4.0d0;
      sg(4)=-1.0d0/2.d0/sqrt(2.0d0);sg(5)=-1.0d0/2.d0/sqrt(2.0d0);
      sg(6)=-5.0d0/4.0d0;
c      call kvec_to_tens(sg,sig_ten)
c      call kjacobi_eigenvalue(3,sig_ten,dir,sg_p)
c                  sig_ten=0.d0
c            sig_ten(1,1)=sg_p(1);sig_ten(2,2)=sg_p(2);
c            sig_ten(3,3)=sg_p(3);
c      call ktens_to_vec(matmul(dir,matmul(sig_ten,transpose(dir))),sg)
      em=24000.0d0;qmuo=0.2d0;fc=30.0d0;ft=3.0d0;
      fc0=fc**1.855d0/60.d0;itypey=1
      hp=0.01d0;qh0=0.3d0;ah=0.08d0;bh=0.003d0;ch=2.0d0;dh=1.0d-6;
      as=1.5d0;gft=0.1d0;wf=gft/0.225d0/ft;wf1=0.158d0*wf;ft1=0.3d0*ft;
      
      itypey=1;
      
      gtol=1.d-4
      fb=1.50d0*fc**(-0.075d0)*fc  
      par=ft/fb*(fb**tw-fc**tw)/(fc**tw-ft**tw)
      ecc=(on+par)/(tw-par)
      pm0=thr*(fc**tw-ft**tw)/fc/ft*ecc/(ecc+on)

      props(1)=em;props(2)=qmuo;props(3)=fc;props(4)=ft;
      props(5)=fc0;props(6)=hp;props(7)=qh0;props(8)=pm0;
      props(9)=ah;props(10)=bh;props(11)=ch;props(12)=dh;
      props(13)=as;props(14)=ecc;props(15)=gtol;props(16)=itypey;
      props(17)=wf; props(18)=wf1;props(19)=ft1;
      statev=0.;SSE=0.;SPD=0.;SCD=0.;RPL=0.;DDSDDT=0.
      DRPLDE=0.;DRPLDT=0.;time=0.;DTIME=0.;TEMP=0.;DTEMP=0.;PREDEF=0.
      DPRED=0.;NDI=3;NSHR=3;NTENS=6;NSTATV=26;NPROPS=19;
      CMNAME='( ''temperature = '')';
      COORDS=0.;DROT=0.;PNEWDT=0.;CELENT=10.0d0;DFGRD0=0.;DFGRD1=0.;
      NOEL=1;NPT=1;LAYER=1;KSPT=1;JSTEP=1;KINC=1;
      
      efc=wf/celent
      
      bs=1.d0;df=0.85d0;
      
      
      
      

      
      sig1(1)=0.;sig2(1)=0.;sig3(1)=0.;sig4(1)=0.;sig5(1)=0.;sig6(1)=0.
      d1(1)=dstran(1);d2(1)=dstran(2);d3(1)=dstran(3);
      d4(1)=dstran(4);d5(1)=dstran(5);d6(1)=dstran(6);
      dt1siz=1;lft=1;llt=1;
      epsps(1)=0.;hsvs(1,1:26)=statev(1:26);hsvs(1,27)=0.0;capa=0.;
      cm(1)=em;cm(2)=qmuo;cm(3)=ecc;cm(4)=qh0;cm(5)=ft;cm(6)=fc
      cm(7)=hp;cm(8)=ah;cm(9)=bh;cm(10)=ch;cm(11)=dh;cm(12)=as
      cm(13)=df;cm(14)=0;cm(15)=1;cm(16)=bs
      cm(17)=wf;cm(18)=wf1;cm(19)=ft1;cm(20)=1
      cm(21)=0;cm(22)=efc;cm(23)=0;cm(24)=0;nlqa=1;crv=0.;nnpcrv(1)=1
      qmat(1,1:3,1:3)=0.;elsizv(1)=CELENT;idelev(1)=1;
      cma=0.;lq1=1;
      etype='( ''temperature = '')';
c            call umat50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
c     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
c     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
c     . reject,nlq,lq1)
      

      WRITE(*,*) 'Hello, world!'
      dstran=0.d0;dstran(1:2)=3.0d-4;dstran(3)=3.2d-4
      dstran(4)=1.5d-8;
      call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      dstran=0.d0;dstran(1)=-5.0d-5;dstran(2)=1.5d-5
      dstran(3)=3.5d-5;dstran(4)=7.0d-6;
      dstran=dstran*1.d0
      call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      stran(1)=0.d0;dstran=0.d0;dstran(1)=-5.0d-3;dstran(2)=1.5d-3
      dstran(3)=2.5d-3;dstran(4)=7.0d-4;
      call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
kcheckvertex      dstran=0.d0;dstran(1)=-5.0d-3;dstran(2:3)=1.5d-3;dstran(4)=7.0d-4;
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      dstran=0.d0;dstran(1)=1.0d-3;dstran(2:3)=-3.5d-3;dstran(4)=7.0d-4;
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)            
      dstran=0.d0;dstran(1)=-5.0d-3;dstran(2:3)=1.5d-3;dstran(4)=7.0d-4;
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)  
      dstran=0.d0;dstran(1)=10.d-3;dstran(2:3)=-3.5d-3;dstran(4)=-5.d-4
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)  
      dstran=0.d0;dstran(1)=2.d-3;dstran(2)=-1.d-3;dstran(4)=-5.d-3
      dstran(3)=-2.d-3;
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)  
      dstran=0.d0;dstran(1)=-5.d-3;dstran(2)=2.d-3;dstran(4)=-5.d-3
      dstran(3)=-2.d-3;
            call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)  
      dstran=0.d0;dstran(1)=-5.d-3;dstran(2)=4.d-3;dstran(4)=-5.d-3
      dstran(3)=-2.d-3;
      call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)  
      dstran=0.d0;dstran(1)=-5.d-3;dstran(2)=4.d-3;dstran(4)=-5.d-3
      dstran(3)=-2.d-3;
      call UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      write(*,*) 'cin50',sig,DDSDDE
      END
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
      implicit double precision(a-h,o-z)
      parameter (nprecd=2)
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
     $rs(4),pj_prev(4,4),sg_tr(6),dinv_dsig_pr(3,3),pp(3,3),pxx(3,3),
     $pnorm_rs(4),sig_ekff_t(6),sig_ekff_c(6),sigi(6),sigii(6)
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %tol=1.0d-6,tol2=10.0d0,sx=6.0d0,half=0.5d0)
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      em=props(1);qmuo=props(2);fc=props(3);ft=props(4);
      fc0=props(5);hp=props(6);qh0=props(7);pm0=props(8);
      ah=props(9);bh=props(10);ch=props(11);dh=props(12);
      as=props(13);ecc=props(14);gtol=props(15);itypey=props(16);
      wf=props(17); wf1=props(18);ft1=props(19); efc=wf/CELENT;
      gtol=1.0d-4;irate_ekffect=0;isotropic=2;imaxsubinc=20;iendflag=0
      bs=1.d0;df=0.85d0;
      fb=1.50d0*fc**(-0.075d0)*fc  
      par=ft/fb*(fb**tw-fc**tw)/(fc**tw-ft**tw)
      ecc=(on+par)/(tw-par)
      pm0=thr*(fc**tw-ft**tw)/fc/ft*ecc/(ecc+on)
      e0=ft/em;
      istrrateflg=0
      tkp=statev(1);
      de(1:6)=dstran;old_e(1:6)=statev(21:26);plen=CELENT;
      e_rate=de;isubinc_count=0;isubinc_flag=0;iconvrg=1
      cc=zr;p1=em/(on+qmuo)/(on-tw*qmuo);cc(1:3,1:3)=p1*qmuo;
      cc(1,1)=(on-qmuo)*p1;cc(2,2)=(on-qmuo)*p1;
      cc(3,3)=(on-qmuo)*p1;Gm=Em/(on+qmuo)/tw;
      cc(4,4)=Gm;cc(5,5)=Gm;cc(6,6)=Gm;
      ddsdde=cc;
      cin=zr;cin(1:3,1:3)=-on*qmuo/Em;cin(1,1)=on/Em;cin(2,2)=on/Em;
      cin(3,3)=on/Em;cin(4,4)=on/Gm;cin(5,5)=on/Gm;cin(6,6)=on/Gm;
      tot_e=old_e+de;strain_rate=de;e_p=statev(3:8)
      tot_e1=tot_e;
      ttot_e=tot_e;
      do while (iconvrg.eq.1.or.isubinc_flag.eq.1)
            el_e=ttot_e-e_p;
            sg_tr=matmul(cc,el_e)
            call khaigh(sg_tr,sv_tr,ro_tr,theta_tr,dinv_dsig_pr);
            call kff(sv_tr,ro_tr,theta_tr,tkp,yield);
            apex_sg=zr
            write(*,*) 'cin50',em
            if(yield>zr) then
                  irtype=0
                  write(*,*) 'cin1',em
                  call kcheckvertex(sv_tr,tkp,apex_sg,irtype)
                  if (irtype.eq.1.or.irtype.eq.2) then
                        call kvertexreturn(sg_tr,apex_sg,
     $                   tkp,irtype,iconvrg,sg_tr)
                  end if
                  if(irtype.eq.0) then
                        call kregular_return(sg_tr,tkp,
     $   iendflag,sg_tr,iconvrg,tkp,pj,rs,pnorm_rs,pp,pxx)
                  end if
            else 
                  iconvrg=0;
                  write(*,*) 'cin7',em
                  do jj=1,6
                        e_p(jj)=statev(jj+2)
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
c                        call xit
                  else if (isubinc_counter > imaxsubinc-1
     $      .and.tkp < on) then
                        tkp=on
                  end if
                  isubinc_flag=1;dtot_e=half*dtot_e;
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
            wt_o=0.d0;wc_o=0.d0;eps_t=0.d0;eps_c=0.d0;pkdt=0.d0;
            pkdt1=0.d0;pkdt2=0.d0
            pkdc=0.d0;pkdc1=0.d0;pkdc2=0.d0;rate_fac=0.d0;alpha=0.d0
      else
            rate_fac=statev(17);eps_t=statev(19);eps_c=statev(20)
            pkdt=statev(9);pkdt1=statev(10);pkdt2=statev(11)
            pkdc=statev(12);pkdc1=statev(13);pkdc2=statev(14)
            wt_o=statev(15);wc_o=statev(16);alpha=statev(18)
            write(*,*) 'cin711',em
            call kc_alpha(sg_tr,sig_ekff_t,sig_ekff_c,alpha);
            write(*,*) 'cin71122',em
            do jj=1,6
                  el_e_old(jj)=statev(20+jj)-statev(jj+2)
            end do
            sg_old=matmul(cc,el_e_old) 
            parrr=0.d0
            do jj=3,8
                  parrr=parrr+(e_p(jj-2)-statev(jj))**tw
            end do
            pnorm_inc_e_p=sqrt(parrr)
            write(*,*) 'cin7112',em
            call kdamage(wc_o,wt_o,strain_rate,rate_fac,
     $ alpha,eps_t,eps_c,pkdt,pkdt1,pkdt2,pkdc,pkdc1,pkdc2,sg_tr,
     $ tkp,pnorm_inc_e_p,plen,sg_old,statev(18),statev(2),
     $ wc,wt,eps_t,eps_c,pkdt,pkdt1,pkdt2,
     $ pkdc,pkdc1,pkdc2,eps_new)!statev(18) old alpha
      write(*,*) 'cin71143',em
      end if
      if (isotropic.eq.0) then
            sig=(on-wt)*sig_ekff_t+(on-wc)*sig_ekff_c
      else if (isotropic.eq.1) then
          sig=(on-wt)*sg_tr
      else
            sig=(on-wt*(on-alpha))*(on-wc*alpha)*sg_tr
c     modified from original code
      end if 
      write(*,*) 'cin7114vv',em
      statev(1)=tkp;statev(3:8)=e_p;statev(2)=eps_new
      statev(9)=pkdt;statev(10)=pkdt1;statev(11)=pkdt2;
      statev(12)=pkdc;statev(13)=pkdc1;statev(14)=pkdc2;statev(15)=wt;
      statev(16)=wc;statev(17)=rate_fac;statev(18)=alpha;
      statev(19)=eps_t;statev(20)=eps_c;statev(21:26)=tot_e;     
      
      
      
c      stiffness
c      d_strain=gtol/10.d0;
c      do iij=1,6
c      de=zr;sigi=sig;
c      if (dstran(iij)>zr) then
c          de(iij)=d_strain
c      else
c          de(iij)=-d_strain
c      end if
c      tkp=statev(1);e_p=statev(3:8)
c      old_e(1:6)=statev(21:26);plen=CELENT;
c      e_rate=de;isubinc_count=0;isubinc_flag=0;iconvrg=1
c      tot_e=old_e+de;strain_rate=de
c      tot_e1=tot_e;
c      ttot_e=tot_e;
c      do while (iconvrg.eq.1.or.isubinc_flag.eq.1)
c            el_e=ttot_e-e_p;
c            sg_tr=matmul(cc,el_e)
c            call khaigh(sg_tr,sv_tr,ro_tr,theta_tr,dinv_dsig_pr);
c            call kff(sv_tr,ro_tr,theta_tr,tkp,yield);
c            apex_sg=zr
c            write(*,*) 'cin50',em
c            if(yield>zr) then
c                  irtype=0
c                  write(*,*) 'cin1',em
c                  call kcheckvertex(sv_tr,tkp,apex_sg,irtype)
c                  if (irtype.eq.1.or.irtype.eq.2) then
c                        call kvertexreturn(sg_tr,apex_sg,
c     $                   tkp,irtype,iconvrg,sg_tr)
c                  end if
c                  if(irtype.eq.0) then
c                        call kregular_return(sg_tr,tkp,
c     $   iendflag,sg_tr,iconvrg,tkp,pj,rs,pnorm_rs)
c                  end if
c            else 
c                  iconvrg=0;
c                  write(*,*) 'cin7',em
c                  do jj=1,6
c                        e_p(jj)=statev(jj+2)
c                        write(*,*) 'cin17',em
c                  end do
c                  write(*,*) 'cin700',em
c                  exit
c            end if
c            write(*,*) 'cin4',em
c            if (iconvrg.eq.1) then
c                  isubinc_counter=isubinc_counter+1;
c                  write(*,*) 'cin5',em
c                  if (isubinc_counter>imaxsubinc) then
c                       write(*,*) 'erorrrr'
c                        call xit
c                  else if (isubinc_counter > imaxsubinc-1
c     $      .and.tkp < on) then
c                        tkp=on
c                  end if
c                  isubinc_flag=1;dtot_e=half*dtot_e;
c                  ttot_e=conv_e+dtot_e;pj_prev=pj;
c            else if (iconvrg.eq.0.and.isubinc_flag.eq.0) then
c                  el_e=matmul(cin,sg_tr)
c                  e_p=tot_e-el_e;
c            else if (iconvrg.eq.0.and.isubinc_flag.eq.1) then
c                  el_e=matmul(cin,sg_tr)
c                  e_p=ttot_e-el_e;conv_e=ttot_e
c                 dtot_e=tot_e-ttot_e;ttot_e=tot_e;
c            end if
c     end do
c    if(itypey.eq.3) then
c            wt_o=0.d0;wc_o=0.d0;eps_t=0.d0;eps_c=0.d0;pkdt=0.d0;
c            pkdt1=0.d0;pkdt2=0.d0;pkdc=0.d0;pkdc1=0.d0;
c            pkdc2=0.d0;rate_fac=0.d0;alpha=0.d0
c      else
c            rate_fac=statev(17);eps_t=statev(19);eps_c=statev(20)
c            pkdt=statev(9);pkdt1=statev(10);pkdt2=statev(11)
c            pkdc=statev(12);pkdc1=statev(13);pkdc2=statev(14)
c            wt_o=statev(15);wc_o=statev(16);alpha=statev(18)
c            write(*,*) 'cin711',em
c            call kc_alpha(sg_tr,sig_ekff_t,sig_ekff_c,alpha);
c            write(*,*) 'cin71122',em
c            do jj=1,6
c                  el_e_old(jj)=statev(20+jj)-statev(jj+2)
c            end do
c            sg_old=matmul(cc,el_e_old) 
c            parrr=zr
c            do jj=3,8
c                  parrr=parrr+(e_p(jj-2)-statev(jj))**tw
c            end do
c            pnorm_inc_e_p=sqrt(parrr)
c            write(*,*) 'cin7112',em
c            call kdamage(wc_o,wt_o,strain_rate,rate_fac,
c     $ alpha,eps_t,eps_c,pkdt,pkdt1,pkdt2,pkdc,pkdc1,pkdc2,sg_tr,
c     $ tkp,pnorm_inc_e_p,plen,sg_old,statev(18),statev(2),
c     $ wc,wt,eps_t,eps_c,pkdt,pkdt1,pkdt2,
c     $ pkdc,pkdc1,pkdc2,eps_new)!statev(18) old alpha
c            
c      write(*,*) 'cin71143',em
c      end if
c     if (isotropic.eq.1) then
c            sigii=(on-wt)*sg_tr
cc            ddsdde=(1-wt)*ddsdde
c      else 
c            write(*,*) 'cin7114uu3',em
c            sigii=(on-wt)*sig_ekff_t+(on-wc)*sig_ekff_c
Cc            ddsdde!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      end if
c      ddsdde(1:6,iij)=(sigii-sigi)/de(iij)
c      end do
c          
c      write(*,*) 'cin88',em
      return
      end  
      subroutine kregular_return(sg,pk,iendflag,sg2,iconvrg,
     $ pkp2,pj,resd,pnorm_res,pp,pxx)
      implicit double precision(a-h,o-z)
      parameter (nprecd=2)
      real*8::resd(4),pnorm_res(4),pincrmt(4),pkm,temp_sg,dth(3),
     $ dl,sg(6),sv_tr1,ro_tr1,ro,sv,theta_tr1,dir(3,3),dir_tem(3),
     $ unkn1(4),tkp1,pkp,pk,sg2(6),pkp2,pj(4,4),sig_ten(3,3),pjr(3,3),
     $ dgdinv(2),ddkdldinv(2),ddgddinv(2,2),ddgdinvdk(2),sg_p(3),
     $ddgdinv(2,2),pi,parry,dinv_dsig_pr(3,3),ppnorm_res,pj2(4,4),
     $ dfdinv(2),ddg_ddinv(2,2),dfdk,ddg_dinvdk(2),ddk_dldinv(2),sg3(6),
     $ dthi(3),sg_p_tr(3),pp(3,3),pt(4,3),pr(3,3),pjr1(3),pxx(3,3)
      integer::iter,itot_iter,iendflag,iconvrg,icomputeall,i,jjj
      real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     $ sx=6.0d0,half=0.5d0,fv=5.0d0)
      pi=fr*datan(on)
      !input sg,pkp,iendflag,gtol
      iter=0;itot_iter=200;pkm=em/thr/(on-tw*qmuo);gm=em/tw/(on+qmuo)
      resd(1:4)=zr;pnorm_res(1:4)=zr;pincrmt(1:4)=zr;dl=zr;      
      call khaigh(sg,sv_tr1,ro_tr1,theta_tr1,dinv_dsig_pr)
      call kvec_to_tens(sg,sig_ten)
      call kjacobi_eigenvalue(3,sig_ten,dir,sg_p)
c      call cdpm2u_computePrincValues(sg,sg_p,0,dir)
c      do j = 1, 3
c          do i = j+1, 3
c              if (sg_p(i) > sg_p(j)) then
c               temp_sg = sg_p(i);dir_tem = dir(1:3,i)
c               sg_p(i) = sg_p(j);dir(1:3,i) = dir(1:3,j)
c               sg_p(j) = temp_sg;dir(1:3,j) = dir_tem
c            end if
c          end do
c      end do
    
      sg_p_tr=sg_p
      pkp=pk;tkp1=pk;sv=sv_tr1;ro=ro_tr1;
      unkn1(1)=sv_tr1;unkn1(2)=ro_tr1;unkn1(3)=tkp1;unkn1(4)=zr;
      call kff(sv,ro,theta_tr1,tkp1,resd(4));
      ppnorm_res=on;ppnorm_res1=on;iconvrg=0;ggtol=gtol*1.d-2
      do while (ppnorm_res>ggtol)
            iter=iter+1
            pnorm_res(1)=resd(1)/pkm;pnorm_res(2)=resd(2)/tw/gm
            pnorm_res(3)=resd(3);pnorm_res(4)=resd(4)
            ppnorm_res1=(pnorm_res(1)**tw+pnorm_res(2)**tw)**half
            ppnorm_res=(pnorm_res(1)**tw+pnorm_res(2)**tw+
     $                pnorm_res(3)**tw+pnorm_res(4)**tw)**half
            if (iter.gt.1) then
                if (ppnorm_res1<gtol*gtol*10.d0) then
                    exit
                end if
            end if            
!            if((iter.eq.itot_iter).or.(isnan(norm_res))) then
            if(iter.eq.itot_iter) then
                  if (ppnorm_res<gtol*1.d-2) then
                      exit
                  end if
                  iconvrg=1
                  exit
            end if
            if (ppnorm_res>ggtol) then
                  icomputeall=1
                  call kderv(sv,ro,theta_tr1,tkp1,icomputeall,
     $    dgdinv,dkdl,dfdinv,ddgddinv,dfdk,ddgdinvdk,ddkdldinv,ddk_dldk)
                  pj(1,1)=on+pkm*dl*ddgddinv(1,1);
                  pj(1,2)=pkm*dl*ddgddinv(1,2);
                  pj(1,3)=pkm*dl*ddgdinvdk(1);
                  pj(1,4)=pkm*dgdinv(1);
                  pj(2,1)=tw*gm*dl*ddgddinv(2,1);
                  pj(2,2)=on+tw*gm*dl*ddgddinv(2,2);
                  pj(2,3)=tw*gm*dl*ddgdinvdk(2);
                  pj(2,4)=tw*gm*dgdinv(2);pj(3,1)=dl*ddkdldinv(1);
                  pj(3,2)=dl*ddkdldinv(2);pj(3,3)=dl*ddk_dldk-on;
                  pj(3,4)=dkdl;pj(4,1)=dfdinv(1);
                  pj(4,2)=dfdinv(2);pj(4,3)=dfdk;pj(4,4)=zr;
                  call Kdet44(pj,parry)
                  if (abs(parry)<1.0d-10) then
                        iconvrg=1
                        exit
                  end if
c                  call computeInverseJac(pj2,pj,error)
c                  if (error.eq.-1) then
c                      converged=1;exit
c                  end if
                  call kmatinv4(pj,pj2)
                  pincrmt=-matmul(pj2,resd)
                  unkn1=unkn1+pincrmt
                  if (unkn1(4)<zr) unkn1(4)=zr
                  if (unkn1(2)<zr) unkn1(2)=zr
                  if (unkn1(3)<pkp) unkn1(3)=pkp
                  sv=unkn1(1);ro=unkn1(2);
                  tkp1=unkn1(3);dl=unkn1(4);
                  call kderv(sv,ro,theta_tr1,tkp1,0,dgdinv,dkdl,
     $  dfdinv,ddg_ddinv,dfdk,ddg_dinvdk,ddk_dldinv,ddk_dldk)
                  resd(1)=sv-sv_tr1+pkm*dl*dgdinv(1)
                  resd(2)=ro-ro_tr1+tw*gm*dl*dgdinv(2)
                  resd(3)=pkp-tkp1+dl*dkdl
                  call kff(sv,ro,theta_tr1,tkp1,resd(4))
            end if 
      end do
      if(iconvrg.eq.0) then
            sg_p=zr;
            sg_p(1)=sv+sqrt(tw/thr)*ro*dcos(theta_tr1)
            sg_p(2)=sv+sqrt(tw/thr)*ro*dcos(theta_tr1-tw/thr*pi)
            sg_p(3)=sv+sqrt(tw/thr)*ro*dcos(theta_tr1+tw/thr*pi);
            sig_ten=zr;
            sig_ten(1,1)=sg_p(1);sig_ten(2,2)=sg_p(2);
            sig_ten(3,3)=sg_p(3);
       call ktens_to_vec(matmul(dir,matmul(sig_ten,transpose(dir))),sg)
            pk=tkp1           
      end if
      return
      end
      subroutine khaigh(sg,sv,ro,thi,dinv_dsig_pr)
            implicit none
            real*8::sg(6),s(6),sv,ro,thi,gtol,pj2,pj3,dinv_dsig_pr(3,3),
     $ zr,on,tw,thr,sx,half,em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ ecc,wf,wf1,ft1,efc,bs,df,fb,gft
          common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol
            PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,sx=6.0d0,
     $       half=0.5d0)
            sv=(sg(1)+sg(2)+sg(3))/thr;s=sg;s(1:3)=sg(1:3)-sv
            pj2=half*(s(1)**tw+s(2)**tw+s(3)**tw)
            pj2=pj2+s(4)**tw+s(5)**tw+s(6)**tw
            if (pj2<=gtol) THEN
                  thi=zr;ro=zr
            else
                  ro=sqrt(tw*pj2);
                  pj3=s(1)**thr+s(2)**thr+s(3)**thr+thr*s(1)*(s(4)**tw
     $             +s(6)**tw)+sx*s(4)*s(5)*s(6)+thr*s(2)*(s(4)**tw
     $             +s(5)**tw)+thr*s(3)*(s(5)**tw+s(6)**tw)
                  pj3=pj3/thr;
                  thi=thr*sqrt(thr)/tw*pj3/(pj2**(1.5d0));
            end if
            if (thi>on) then
                  thi=on
            else if (thi<-on) then
                  thi=-on
            end if
            thi= on/thr*dacos(thi)
            return
            end
      subroutine kff(sv,ro,th,pkp,f)
            implicit none
            real*8::sv,ro,th,pkp,f,par1,par2,rcos,qh1,dqh1_dk,
     $ qh2,dqh2_dk,ecc,pm0,fc,e0,qh0
            real*8::em,qmuo,ft,fc0,hp,ah,bh,ch,dh,as,
     $ gtol,wf,wf1,ft1,efc,bs,df,fb,gft,on,tw,thr,fr,sx,half
            integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
            PARAMETER(on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %sx=6.0d0,half=0.5d0)
            !input sv,ro,th,pkp
            !this part is repeated if you change it check others
            par1=(on-ecc**tw);par2=(tw*ecc-on)
            rcos=(fr*par1*dcos(th)**tw+par2**tw)
            rcos=rcos/(tw*par1*dcos(th)+par2*sqrt(fr*par1*dcos(th)**tw
     $       +5.d0*ecc**tw-fr*ecc))
            call kcqh1(pkp,int(0),qh1,dqh1_dk)
            call kcqh2(pkp,int(0),qh2,dqh2_dk)
            par1=ro/(fc*sqrt(sx))+sv/fc;par1=(on-qh1)*par1**tw
     $       +sqrt(1.5d0)*ro/fc
            f=par1**tw+qh1**tw*qh2*pm0*(sv/fc+ro*rcos/(sqrt(tw*thr)*fc))
     $       -qh1**tw*qh2**tw
            return
            end
      subroutine kcqh2(par,icomput_deriv,qh1,dqh1_dk)
            implicit none
            real*8::hp,par,qh1,dqh1_dk
            integer::icomput_deriv
            real*8::em,qmuo,fc,ft,fc0,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,on,zr
            integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
              PARAMETER(zr=0.0D0,on=1.0D0)
            if (par<on) then
                  qh1=on
            else
                  qh1=on+(par-on)*hp
            end if
            dqh1_dk=zr
            if (icomput_deriv.eq.1) then
                  if (par<=on) then
                        dqh1_dk=zr
                  else
                        dqh1_dk=hp
                  end if                 
            end if
            return
            end
      subroutine kcqh1(par,icomput_deriv,qh2,dqh2_dk)
            implicit none
            real*8::hp,qh0,par,qh2,dqh2_dk,rcos,r_top,r,par2,par1,
     $ pmQ
            integer::icomput_deriv,i
            real*8::em,qmuo,fc,ft,fc0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr,r1
            integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
            PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0)
            if (par<=zr) then
                  qh2=qh0
            else if (par>zr.and. par<on) then
                  qh2=(on-qh0-hp)*par**thr-thr*(on-qh0-hp)*par**tw
     $             +(thr*(on-qh0)-tw*hp)*par+qh0
            else
                  qh2=on
            end if
            dqh2_dk=zr
            if (icomput_deriv.eq.1) then
                  if (par<=on) then
                    dqh2_dk=thr*(on-qh0-hp)*par**tw-tw*thr*(on-qh0-hp)
     $                   *par+(thr*(on-qh0)-tw*hp)
                  else
                        dqh2_dk=zr
                  end if
            end if
          
            
c            r1=1.1d0
c            if (par<=gtol) then
c                  qh2=qh0
c            else if (par>gtol.and. par<on) then
c                  qh2=qh0+(on-qh0)*par*r1/(r1-on+par**r1)
c     $             +hp*(par**tw-par)
c            else
c                  qh2=on
c            end if
c            dqh2_dk=zr
c            if (icomput_deriv.eq.1) then
c                  if (par<=on) then
c                      if (par<=gtol) then 
c                          par=gtol
c                      end if
c           dqh2_dk=(qh0-on)*r1*(par**r1-on)*(r1-on)/
c     $     (r1-on+par**r1)**2+hp*(tw*par-on)
c                  else
c                        dqh2_dk=zr
c                  end if
c            end if
            return
            end          
      subroutine kc_alpha(s,st,sc,alpha)
            implicit none
            real*8::s(6),sig_pr(3),dir(3,3),s_pt(6),s_pc(6),
     $ alpha,alpha_t,pnorm_s,st(6),sc(6),sig_ten(3,3),zr,on,tw,half
            integer::i,jjj
            PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,half=0.5d0)
            !input s
            call kvec_to_tens(s,sig_ten)
            call kjacobi_eigenvalue(3,sig_ten,dir,sig_pr)
c            call cdpm2u_computePrincValues(sig_ten,sig_pr,0,dir)
            pnorm_s=(sig_pr(1)**tw+sig_pr(2)**tw+sig_pr(3)**tw)
            alpha_t=zr;s_pt=zr;s_pc=zr;
            if (pnorm_s>zr) then
               do i=1,3
                  if (sig_pr(i)>=zr) then
                        s_pt(i)=sig_pr(i);s_pc(i)=zr
                  else
                        s_pc(i)=sig_pr(i);s_pt(i)=zr
                  end if
                  alpha_t=alpha_t+s_pt(i)*(s_pt(i)+s_pc(i))
               end do
               alpha_t=alpha_t/pnorm_s
            end if
            alpha=on-alpha_t

            call kvec_to_tens(s_pt,sig_ten)
c       call ktens_to_vec(matmul(transpose(dir),matmul(sig_ten,dir)),st)
        call ktens_to_vec(matmul(dir,matmul(sig_ten,transpose(dir))),st)
            call kvec_to_tens(s_pc,sig_ten)
c       call ktens_to_vec(matmul(transpose(dir),matmul(sig_ten,dir)),sc)
        call ktens_to_vec(matmul(dir,matmul(sig_ten,transpose(dir))),sc)
            !check
            return
            end
      subroutine krate_fac(strain_rate,alpha,tol,rate_fac)
            implicit none
            real*8::strain_rate(6),alpha,tol,rate_fac,
     $ alphas,gammas,betas,rate_t0,rate_c0,rate_t,rate_c,tmp(6),
     $ strain_pr(3),dir(3,3),pmax,pmin,rate,ratio_t,ratio_c,
     $ ratio_c0,ratio_t0,deltas,sig_ten(3,3)
            integer i
            real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,on,tw,sx
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
       PARAMETER(on=1.0D0,tw=2.0D0,sx=6.0d0)
            !strain_rate,alpha,tol
            alphas=on/(5.d0+9.d0*fc/fc0);deltas=on/(on+8.d0*fc/fc0)
            gammas=exp((6.156d0*alphas-tw)*log(10.d0))!check log
            betas=exp((sx*deltas-tw)*log(10.d0))
            rate_t0=1.0d-6;rate_c0=-30.0d-6;rate_t=on;rate_c=on
            tmp=strain_rate/tw;tmp(1)=tmp(1)*tw;tmp(2)=tmp(2)*tw;
            tmp(3)=tmp(3)*tw
            call kvec_to_tens(tmp,sig_ten)
            call kjacobi_eigenvalue(3,sig_ten,dir,strain_pr)
c            call cdpm2u_computePrincValues(sig_ten,strain_pr,0,dir)
            pmax=-1.0d-20;pmin=1.0d20
            do i=1,3
                  if (pmax<strain_pr(i))      pmax=strain_pr(i)
                  if (pmin>strain_pr(i))      pmin=strain_pr(i)
            end do
            if ((on-alpha)>tol) then
                  rate =pmax
            else
                  rate =pmin
            end if
            ratio_t=rate/rate_t0;ratio_c=rate/ratio_c0;
            if (rate<30.0d-6) then
                  rate_t=on
            else if (rate>30.0d-6.and.rate<on) then
                  rate_t=ratio_t**deltas
            else
                  rate_t=betas*ratio_t**(0.3333d0)
            end if
            if (rate>-30.0d-6) then
                  rate_c=on
            else if (rate>-30.0d0 .and. rate<-30.0d-6) then
                  rate_c=ratio_c**(1.026d0*alphas)
            else
                  rate_c=gammas*ratio_c**(0.3333d0)
            end if
            rate_fac=(on-alpha)*rate_t+alpha*rate_c
           return
           end
      subroutine kratiopotential(sv,pkp,ratio)
            implicit none
            real*8::sv,pkp,ratio,qh1,dqh1_dk,
     $qh2,dqh2_dk,o,pm,pmin_equ_e,eps_old,equ_e_minus,equ_e_new,
     $equ_e_plus,sg_minus(6),sg_plus(6),sg_old(6),dsg(6),
     $sg(6),ro,th,dgdinv(2),pmg,par1,ag,bg,par2,dmg,bg_top,
     $bg_bottom,r,b1
             real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr,fr,sx
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,sx=6.0d0)
            !input sv,pkp
            ro=zr
            call kcqh1(pkp,int(0),qh1,dqh1_dk)
            call kcqh2(pkp,int(0),qh2,dqh2_dk)
            ag=thr*ft*qh2/fc+pm0/tw;
            bg_top=qh2/thr*(on+ft/fc)
            bg_bottom=log(ag)-log(tw*df-on)-log(thr*qh2+pm0/tw)
     $       +log(df+on)!log vs exp
            bg=bg_top/bg_bottom
            r=(sv-qh2*ft/thr)/fc/bg
            pmg=ag*bg*fc*exp(r);par1=ro/(fc*sqrt(sx))+sv/fc!log vs exp
            dmg=ag*exp(r);
            par2=(on-qh1)*par1**tw+sqrt(1.5d0)*ro/fc
            dgdinv(1)=fr*(on-qh1)/fc*par2*b1+qh1**tw*pmg/fc
            dgdinv(2)=par2/(sqrt(sx)*fc)*(fr*(on-qh1)*par1+6)
     $       +pm0*qh1**tw/(sqrt(sx)*fc)
            ratio=dgdinv(2)/dgdinv(1)*thr*(on-tw*qmuo)/(on+qmuo)
            return
            end
      subroutine kcheckunload(sg,sg_old,eps_old,gtol,pmin_equ_e,
     $       equ_e_new,iunload_flag)
            implicit none
            real*8::sg1(6),sg(6),sg_old(6),eps_old,dinv_dsig_pr(3,3)
     $,sv,sg_minus(6),sg_plus(6),pm,o,p,pn,ro,th,pmin_equ_e,
     $equ_e_plus,equ_e_minus,equ_e_new,sdg(6),equ_e1,dsg(6),
     $gtol,grgtol
            integer::iunload_flag,i
            !input sg,sg_old,eps_old,gtol
            call khaigh(sg,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_new)
            dsg=sg-sg_old;sg_plus=sg_old+0.01d0*dsg;
            sg_minus=sg_old+0.99d0*dsg
            call khaigh(sg_plus,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_plus)
            call khaigh(sg_minus,sv,ro,th,dinv_dsig_pr);
            call kequ_e(sv,ro,th,equ_e_minus)
            iunload_flag=0;pmin_equ_e=eps_old;p=equ_e_plus;
            pm=equ_e_minus;o=eps_old;pn=equ_e_new;grgtol=gtol*1.d-3
            if ((p<o.and.pm<pn).and.
     $       (abs(p-o)>grgtol.and.abs(pm-pn)>grgtol)) then
                  iunload_flag=1
                  do i=1,100
                        sg1=sg_old+dsg*i/100.0d0
                        call khaigh(sg1,sv,ro,th,dinv_dsig_pr)
                        call kequ_e(sv,ro,th,equ_e1)
                        if (equ_e1<=pmin_equ_e) then
                              pmin_equ_e=equ_e1
                        else
                              EXIT
                        end if
                  end do
            end if
            return
            end
      subroutine kderv(sv,ro,th,tkp,icomputeall,dg_dinv,dkdl,dfdinv,
     $       ddg_ddinv,dfdk,ddg_dinvdk,ddk_dldinv,ddk_dldk)
            implicit none
            real*8::sv,ro,th,tkp,dg_dinv(2),dkdl,dfdinv(2),
     $       ddg_ddinv(2,2),dfdk,ddg_dinvdk(2),ddk_dldinv(2),
     $ qh1,dqh1_dk,qh2,dqh2_dk,dEquivdg_dstress_dinv(2),dduct_dinv(2),
     $ equivdg_dstress,equivaplentdg_dsg,ecc,
     $ duct_m,dR_top_dk,dR_bottom_dk,dR_dk,dmQ_dsv,
     $ dmQ_dk,dfdqh2,dfdqh1,dEquivaplentdg_dstress_dk,
     $ddk_dldk,ddg_dsvdro,ddg_drodsv,ddg_ddsv,r,par1,par,pmQ,
     $ddg_ddro,dbl_dsv,dbl_dro,dbg_top_dk,dag_dk,dbg_bottom_dk,
     $dAl_dyieldHard,dAl_dsv,dAl_dro,Bl,bg_top,bg_bottom,
     $bg,Al,ag,pi,rcos,R_top,r_bottom,ddg_dsv_dro,
     $dbg_dk,par2
            integer::icomputeall
                  real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr,fr,sx,half
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %sx=6.0d0,half=0.5d0)
            !sv,ro,th,tkp,icomputeall
            pi=fr*datan(on)
            call kcqh1(tkp,int(1),qh1,dqh1_dk);
            call kcqh2(tkp,int(1),qh2,dqh2_dk)
            !dgdinv 
            ag=ft*qh2*thr/fc+pm0/tw;bg_top=qh2/thr*(on+ft/fc)
            bg_bottom=log(ag)-log(tw*df-on)-log(thr*qh2+pm0/tw)
     $       +log(df+on)!log vs exp
            bg=bg_top/bg_bottom
            r=(sv-ft/thr*qh2)/fc/bg;pmQ=ag*exp(r);
            Bl=sv/fc+ro/(fc*sqrt(sx));
            Al=(on-qh1)*Bl**tw+sqrt(1.5d0)*ro/fc;
            dg_dinv(1)=fr*(on-qh1)/fc*Al*Bl+qh1**tw*pmQ/fc;
            dg_dinv(2)=Al/(sqrt(sx)*fc)*(fr*(on-qh1)*Bl+sx)+pm0
     $       *qh1**tw/(sqrt(sx)*fc);
            !dkdl
            equivaplentdg_dsg=sqrt(on/thr*dg_dinv(1)**tw+dg_dinv(2)**tw)
            call kductility(sv,th,1,duct_m,dduct_dinv)
            dkdl=equivaplentdg_dsg/duct_m;
            !dfdinv
            if (icomputeall.eq.1) then
                  par1=(on-ecc**tw);par2=(tw*ecc-on)
                  rcos=(fr*par1*dcos(th)**tw+par2**tw)
                  rcos=rcos/(tw*par1*dcos(th)+par2*sqrt(fr*par1
     $             *dcos(th)**tw+5.d0*ecc**tw-fr*ecc))
                  dfdinv(1)=fr*(on-qh1)/fc*Al*Bl+qh2*qh1**tw*pm0/fc
                  dfdinv(2)=Al/(sqrt(sx)*fc)*(fr*(on-qh1)*Bl+sx)
     $             +rcos*pm0*qh2*qh1**tw/(sqrt(sx)*fc);
                  !ddgddinv 
                  dmQ_dsv=ag/(bg*fc)*exp(r);dAl_dsv=tw*(on-qh1)*Bl/fc
                  dBl_dsv=on/fc;dAl_dro=tw*(on-qh1)*Bl/(fc*sqrt(sx))
     $             +sqrt(1.5d0)/fc
                  dBl_dro=on/(fc*sqrt(sx))
                  ddg_ddsv=fr*(on-qh1)/fc*(dAl_dsv*Bl+Al*dBl_dsv)
     $             +qh1**tw*dmQ_dsv/fc
                  ddg_ddro=dAl_dro/(sqrt(sx)*fc)*(fr*(on-qh1)*Bl+sx)
     $             +Al*dBl_dro*fr*(on-qh1)/(sqrt(sx)*fc)
                  ddg_dsvdro=fr*(on-qh1)/fc*(dAl_dro*Bl+Al*dBl_dro)
                  ddg_drodsv=dAl_dsv/(sqrt(sx)*fc)*(fr*(on-qh1)*Bl+sx)
     $             +Al/(sqrt(sx)*fc)*(fr*(on-qh1)*dBl_dsv)
                  ddg_ddinv(1,1)=ddg_ddsv;ddg_ddinv(1,2)=ddg_drodsv
                  ddg_ddinv(2,1)=ddg_drodsv;ddg_ddinv(2,2)=ddg_ddro
                  !dfdk
                  dfdqh1=-tw*Al*(Bl**tw)+tw*qh1*qh2*pm0*(sv/fc+ro*rcos
     $             /(sqrt(sx)*fc))-tw*qh1*(qh2**tw)
                  dfdqh2=(qh1**tw)*pm0*(sv/fc+ro*rcos/(sqrt(sx)*fc))-
     $             tw*qh2*(qh1**tw)
                  dfdk=dqh1_dk*dfdqh1+dqh2_dk*dfdqh2
                  if(dfdk>zr)    dfdk=zr
                  !ddgdinvdk
                  dag_dk=dqh2_dk*thr*ft/fc;dbg_top_dk=dqh2_dk/thr
                  dbg_bottom_dk=-thr*dqh2_dk/(thr*qh2+pm0/tw)
                  dbg_dk=(dbg_top_dk*bg_bottom-bg_top*dbg_bottom_dk)
     $             /(bg_bottom**tw)
                  R_top=(sv-ft/thr*qh2);R_bottom=fc*bg
                  dR_top_dk=-ft/thr*dqh2_dk;dR_bottom_dk=fc*dbg_dk
                  dR_dk=(dR_top_dk*R_bottom-R_top*dR_bottom_dk)
     $             /(R_bottom**tw)
                  dmQ_dk=dag_dk*exp(r)+ag*dR_dk*exp(r);
                  dAl_dyieldHard=-(Bl**tw)
                  ddg_dinvdk(1)=(-fr*Al*Bl/fc+fr*(on-qh1)/fc
     $             *dAl_dyieldHard*Bl)*dqh1_dk
                  ddg_dinvdk(1)=ddg_dinvdk(1)+dqh1_dk*tw*qh1*pmQ
     $             /fc+qh1*dmQ_dk/fc
                  par=dAl_dyieldHard/(sqrt(sx)*fc)*(fr*(on-qh1)*Bl+sx)
                  par1=-fr*Al/(sqrt(sx)*fc)*Bl+pm0/(sqrt(sx)*fc)
                  ddg_dinvdk(2)=(par+par1)*tw*qh1*dqh1_dk;
                  !ddk_dldk
                  par1=dg_dinv(1)/equivaplentdg_dsg*ddg_dinvdk(1)
                  par1=tw/thr*par1
                  dEquivaplentdg_dstress_dk=(par1+tw*dg_dinv(2)
     $             /equivaplentdg_dsg*ddg_dinvdk(2))/tw
                  ddk_dldk= dEquivaplentdg_dstress_dk/duct_m
                  !ddKappadDeltaLambdadInv
                  dEquivdg_dstress_dinv(1)=tw/thr*dg_dinv(1)
     $             *ddg_ddinv(1,1)+tw*dg_dinv(2)*ddg_ddinv(2,1)
                  dEquivdg_dstress_dinv(1)=dEquivdg_dstress_dinv(1)
     $             /(tw*equivaplentdg_dsg)
                  dEquivdg_dstress_dinv(2)=tw/thr*dg_dinv(1)
     $             *ddg_ddinv(1,2)+tw*dg_dinv(2)*ddg_ddinv(2,2)
                  dEquivdg_dstress_dinv(2)=dEquivdg_dstress_dinv(2)
     $             /(tw*equivaplentdg_dsg)
                  ddk_dldinv(1)=(dEquivdg_dstress_dinv(1)
     $             *duct_m-equivaplentdg_dsg*dduct_dinv(1))/(duct_m**2)
                  ddk_dldinv(2)=(dEquivdg_dstress_dinv(2)
     $             *duct_m-equivaplentdg_dsg*dduct_dinv(2))/(duct_m**2);
            else 
                  dfdinv=zr;ddg_ddinv=zr;dfdk=zr;ddg_dinvdk=zr;
                  ddk_dldinv=zr;ddk_dldk=zr
            end if
            return
            end
      subroutine kdamage(wc_old,wt_old,strain_rate,rate_fc,
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
     $ rate_fc,eps_t1,eps_c1,pkdt_new,pmin_equ_e,t_equ_e,
     $pkdc_new,eps_new,xs,wf1,wf,tol,theta_el,e0,fac,ft,
     $t_rate_fac,sv_el,rs1,ro_el,residualDerivative,df,
     $residual,qh2,pari,pkdt1t,pkdt2t,pkdc1t,pkdc2t,as,d_pkdc,
     $d_pkdt,dinv_dsig_pr(3,3),dc,dt
            integer::istep1flag,iunload_flag,itypey,istrrateflg,
     $ inewton_iter      
                  real*8::em,qmuo,fc,fc0,hp,qh0,pm0,ah,bh,ch,dh,
     $ecc,ft1,efc,bs,fb,gft,zr,on,tw,sx,half,thr
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0d0,sx=6.d0,half=0.5d0)
            !input wc_old,wt_old,strain_rate,rate_fc,alpha,eps_t,eps_c,pkdt_old,pkdt1,pkdt2,
            !pkdc_old,pkdc1,pkdc2,sg_ekff,tkp,pnorm_inc_e_p,plen,sg_old,alpha_old,eps_old
            tol=gtol*10
            if (rate_fc.eq.zr) then
                  istep1flag=1;rate_fc=on
            else 
                  istep1flag=0
            end if
            call kcheckunload(sg_ekff,sg_old,eps_old,gtol,
     $       pmin_equ_e,t_equ_e,iunload_flag)
            eps_new=t_equ_e;
            if(iunload_flag.eq.0) then
              pmin_equ_e=eps_old
            end if
            dt=t_equ_e-eps_old;
            dc=(pmin_equ_e-eps_old)*alpha_old+(t_equ_e-pmin_equ_e)*alpha
            if(istrrateflg.eq.1.and.wc_old.eq.zr.and.wt_old.eq.zr)then
                  call krate_fac(strain_rate,alpha,tol,t_rate_fac)
                  eps_t1=eps_t+dt/t_rate_fac;eps_c1=eps_c+dc/t_rate_fac
                  if((eps_c1>e0.or.eps_t1>e0).and.istep1flag.ne.1) then
                        eps_t1=eps_t+dt/rate_fc;eps_c1=eps_c+dc/rate_fc
                  else
                        rate_fc=t_rate_fac
                  end if
            else
                  eps_t1=eps_t+dt/rate_fc;eps_c1=eps_c+dc/rate_fc
            end if
            
            call khaigh(sg_ekff,sv_el,ro_el,theta_el,dinv_dsig_pr)
            if(sv_el<zr) then
                  rs1=-(sx)**(half)*sv_el/max(ro_el,1.0d-16)
            else 
                  rs1=zr
            end if
            xs=on+(as-on)*rs1;d_pkdt=(eps_t1-pkdt_old);
            d_pkdc=(eps_c1-pkdc_old);wt=zr;wc=zr;
            if(d_pkdt>=-e0*tol) then
                  if(eps_t1<e0*(on-tol)) then
                        fac=zr
                  else if (eps_t1>e0*(on-tol).and.pkdt_old
     $                  <e0*(on-tol)) then
                        fac=(on-(e0-pkdt_old)/(eps_t1-pkdt_old))
                  else 
                        fac=on
                  end if
                  pkdt1t=pkdt1+pnorm_inc_e_p*fac/xs/rate_fc;
                  pkdt2t=pkdt2+d_pkdt/xs
                  pkdt_new=eps_t1
                  call kdamaget(pkdt_new,pkdt1t,pkdt2t,plen,wt_old,wt)
            end if
            if(d_pkdc>=-e0*tol) then
                  if (eps_c1<e0*(on-tol)) then
                        fac=zr
                  elseif (eps_c1>e0*(on-tol).and.pkdc_old
     $                        <e0*(on-tol)) then
                        fac=(on-(e0-pkdc_old)/(eps_c1-pkdc_old))
                  else 
                        fac=on
                  end if
                  call kcqh2(tkp,int(0),qh2,dqh2_dk)
                  fac2=ft*qh2*sqrt(tw/thr)
     $             /max(ro_el,1.0d-16)/sqrt(on+tw*(df**tw));
                  pkdc1t=pkdc1+pnorm_inc_e_p*fac*fac2*alpha/xs/rate_fc;
                  pkdc2t=pkdc2+d_pkdc/xs;
                  pkdc_new=eps_c1
                  call kdamagec(pkdc_new,pkdc1t,pkdc2t,wc_old,wc)
            end if
            return
            end
      subroutine kdamaget(pkdt,pkdt1,pkdt2,plen,wt_old,wt)
            implicit none
            real*8::pkdt,pkdt1,pkdt2,plen,wt_old,wt,
     $ pari,residual,residualDerivative,ytol
            integer::iter,inewton_iter
                  real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr,fr,sx,half
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %sx=6.0d0,half=0.5d0)
            !pkdt,pkdt1,pkdt2,plen,wt_old output wt
            inewton_iter=100;ytol=gtol*10.d0;
            if (pkdt>e0*(on-ytol)) then
                  if (itypey.eq.0) then
                        wt=(Em*pkdt*wf-ft*wf+ft*pkdt1*plen)
     $                        /(Em*pkdt*wf-ft*plen*pkdt2)
                  else if (itypey.eq.1) then
                        wt=(Em*pkdt*wf1-ft*wf1-(ft1-ft)*pkdt1*plen)
     $                        /(Em*pkdt*wf1+(ft1-ft)*plen*pkdt2)
                        pari=plen*pkdt1+plen*wt*pkdt2
                        if (pari>wf1.and.pari<wf) then
                              wt=(Em*pkdt*(wf-wf1)-ft1*(wf-wf1)+ft1
     $                              *pkdt1*plen-ft1*wf1)
                              wt=wt/(Em*pkdt*(wf-wf1)-ft1*plen*pkdt2)
                              pari=plen*pkdt1+plen*wt*pkdt2
                        else if (pari>wf) then
                              wt=on
                        end if
                  else if (itypey.eq.2) then
                        !Exponential: Iterative solution with N-R procedure
                        wt=on;residual=zr;residualDerivative=zr;iter=0;
                        pari=on;
                        do while (pari.eq.on)
                              iter=iter+1;residual=(on-wt)*Em*pkdt-ft
     $                              *EXP(-plen*(wt*pkdt2+pkdt1)/wf)
                              residualDerivative=-Em*pkdt+ft*plen*pkdt2
     $                              *EXP(-plen*(wt*pkdt2+pkdt1)/wf)/wf
                              wt=wt-residual/residualDerivative
                              !if(iter>inewton_iter);  
                              !disp('* Algorithm for tensile kdamage-No convergence reached after 100 iterations *')
                              !error('erorrre');end
                              if(abs(residual/ft)<1.0d-8) pari=zr
                        end do
                  else
                        wt=zr
                  end if
                  if(wt>on) wt=on
                  if(wt<zr.or.wt<wt_old) wt=wt_old
            else 
                  wt=zr
            end if
            return
            end
      subroutine kdamagec(pkdc,pkdc1,pkdc2,wc_old,wc)
            implicit none
            real*8::pkdc,pkdc1,pkdc2,wc_old,wc,pari,gtol,ytol,tol,
     $ residual,residualDerivative,dResdw,errorOld,efc,ft,e0,em
            integer::iter,nite,inewton_iter,isotropic
            real*8::qmuo,fc,fc0,hp,qh0,pm0,ah,bh,ch,dh,
     $ecc,wf,wf1,ft1,bs,df,fb,gft,as,zr,on,tw,thr,fr,sx,half
      integer::iunload_flag,itypey,istrrateflg
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %sx=6.0d0,half=0.5d0)
            !pkdc,pkdc1,pkdc2,wc_old
            if (isotropic.eq.1) then
                  wc=zr;pkdc1=zr;pkdc2=zr;pkdc=zr
            else
                  inewton_iter=200;tol=gtol;ytol=gtol*10.d0;nite=0;
                  residual=zr;dResdw=zr;
                  if (pkdc>e0*(on-ytol)) then
                        do while(nite<inewton_iter)
                              nite=nite+1;residual=(on-wc)*em*pkdc-ft
     $                              *exp(-(pkdc1+wc*pkdc2)/efc)
                              dResdw=-em*pkdc+ft*pkdc2/efc*exp(-(pkdc1+
     $                              wc*pkdc2)/efc);
                              wc=wc-residual/dResdw
                              errorOld = residual/ft
                              if(wc<zr) then
                                    wc=zr
                                    exit
                              end if
                              if(nite.eq.inewton_iter) then
                                    if(residual<zr) then
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
                        wc=zr
                  end if
            end if
            if(wc>on) wc=on
            if(wc<zr.or.wc<wc_old) wc=wc_old
            return
            end
      subroutine kcheckvertex(sv_tr,tkp,apex_sg,irtype)
            implicit none
            real*8::sv_tr,sv,tkp,apex_sg,qh22,dqh2_dk,qh2
            integer::irtype
            real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on
      integer::iunload_flag,itypey,istrrateflg,isotropic
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0)
            !sv_tr,tkp
            if(sv_tr>zr) then
                  irtype=1
                  if(tkp<on) then
                        apex_sg=zr
                  else 
                        call kcqh2(tkp,int(0),qh2,dqh2_dk)
                        apex_sg=qh2*fc/pm0
                  end if
            else if (sv_tr<zr .and. tkp<on) then
                  irtype=2;apex_sg=zr
            else
                  irtype=0;apex_sg=zr
            end if
            return
            end
      subroutine kvertexreturn(sg,apex_sg,tkp,irtype,iconvrg,sg_ekff)
            implicit none
            real*8::pkp,sv,sv2,yvalue_mid,gtol,pkp0,tkp,sg(6),apex_sg,
     $sg_ekff(6),dsv,pari,ratioPotent1, ratiotrial,ro,tkpi,
     $stress(6),sv_mid,svAnswer,theta,tk,ytol,yvalue,dinv_dsig_pr(3,3),
     $zr,on,tw,thr,fr,sx,half,em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,
     $e0,ecc,wf,wf1,ft1,efc,bs,dffb,gft,as,df,fb
            integer::irtype,iconvrg,j,maxiter,iunload_flag,itypey,
     $istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %sx=6.0d0,half=0.5d0)
            common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
            !sg,apex_sg,tkp,irtype,iconvrg
            ytol=gtol*1.d-2;yvalue=zr;yvalue_mid=zr;sv2=zr;pkp0=tkp;
            maxiter=250
            call khaigh(sg,sv,ro,theta,dinv_dsig_pr)
            sv2 = apex_sg;
            call kpp(pkp0,sv,ro,sv,tkpi);
            call kff(sv,zr,zr,tkpi,yvalue)
            call kpp(pkp0,sv,ro,sv2,tkpi);
            call kff(sv2,zr,zr,tkpi,yvalue_mid)
            pari=zr
            if(yvalue*yvalue_mid>=zr) then
                  iconvrg=1;irtype=0;
            else
                  if(yvalue<zr) then
                        dsv=sv2-sv;svAnswer=sv2
                  else 
                        dsv=sv-sv2;svAnswer=sv2
                  end if
                  do j=1,maxiter
                        dsv=half*dsv
                        sv_mid=svAnswer+dsv
                        call kpp(pkp0,sv,ro,sv_mid,tkpi)
                        call kff(sv_mid,zr,zr,tkpi,yvalue_mid)
                        if(yvalue_mid<=zr) then
                              svAnswer=sv_mid
                        end if
                        if (abs(yvalue_mid)<ytol.and.
     $                   yvalue_mid<=zr) then 
                              call kratiopotential(svAnswer,tkpi,
     $                              ratioPotent1)
                              ratiotrial=ro/(sv-svAnswer)
                              if((ratioPotent1>=ratiotrial.and.irtype
     $                           .eq.1).or.(ratioPotent1<=ratiotrial
     $                           .and.irtype.eq.2)) then
                                    exit
                              else
                                    iconvrg=1;irtype=0;pari=on
                                    exit
                             end if
                        end if
                  end do
                  if (pari.eq.zr) then
                        sg_ekff(1:3)=svAnswer;sg_ekff(4:6)=zr;tkp=tkpi
                        iconvrg=zr
                  end if
            end if
            return
            end
      subroutine kpp(pkp_old,sv1,dro,sv2,pkp)
            implicit none
            real*8::pkp_old,pkm3,em,qmuo,dro,sv1,sv2,
     $ equ_d_e_p,duct_m,dduct_dinv(2),pkp,pi,gm2
            real*8::fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr,fr
      integer::iunload_flag,itypey,istrrateflg,isotropic
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0d0)
            pkm3=em/(on-tw*qmuo);gm2=em/(on+qmuo)
            equ_d_e_p=sqrt((((sv1-sv2)/pkm3)**tw+(dro/gm2)**tw))
            pi=fr*datan(on)
            call kductility(sv2,pi/thr,int(0),duct_m,dduct_dinv)
            pkp=pkp_old+equ_d_e_p/duct_m
            return
            end
      subroutine kductility(sv,th,icomput_deriv,duct_m,dduct_dinv)
            implicit none
            real*8::par1,sv,th,duct_m,dduct_dinv(2),eh,x,fh
            real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,thr
      integer::iunload_flag,itypey,istrrateflg,isotropic,icomput_deriv
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0)
            !input sv,th,icomput_deriv (integer important)
            par1=(tw*dcos(th))**tw;x=-on*(sv+fc/thr)/fc;
            if (x<zr) then
                  eh=bh-dh;fh=eh*ch/(ah-bh)
                  duct_m=(eh*exp(x/fh)+dh)/par1
                  if (icomput_deriv.eq.1) then
                        dduct_dinv(1)=eh/fh*exp(x/fh)/par1*(-on)/fc
                  else
                        dduct_dinv(1)=zr
                  end if
            else
                  duct_m=(ah-(ah-bh)*exp(-x/ch))/par1
                  if (icomput_deriv.eq.1) then
                        dduct_dinv(1)=(bh-ah)/ch*exp(-x/ch)/par1/fc
                  else
                        dduct_dinv(1)=zr
                  end if
            end if
            return
            end  
      subroutine kequ_e(sv,ro,th,equ_e)
            implicit none
            real*8::sv,ro,th,equ_e,par1,par2,
     $ rcos,par_p,par_q    
            real*8::em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,gtol,wf,wf1,ft1,efc,bs,df,fb,gft,as,zr,on,tw,fr,half,sx
      integer::iunload_flag,itypey,istrrateflg,isotropic
      common/kumat/em,qmuo,fc,ft,fc0,hp,qh0,pm0,ah,bh,ch,dh,e0,
     $ecc,wf,wf1,ft1,efc,bs,df,fb,gft,gtol,as,
     $iunload_flag,itypey,istrrateflg,isotropic
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,fr=4.0D0,sx=6.0d0,half=0.5d0)
            par1=(on-ecc**tw);par2=(tw*ecc-on)
            rcos=(fr*par1*dcos(th)**tw+par2**tw)
            rcos=rcos/(tw*par1*dcos(th)+par2*sqrt(fr
     $       *par1*dcos(th)**tw+5.d0*ecc**tw-fr*ecc))
            par_p=-pm0*(ro*rcos/sqrt(sx)/fc+sv/fc);
            par_q=-1.5d0*ro**tw/fc**tw
            equ_e=e0*(-half*par_p+sqrt(par_p**tw/fr-par_q))
            if (equ_e.le.zr) equ_e=zr
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
      it_num = 0;it_max=100;
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
                  if ( d(l) > d(m) ) then
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
      subroutine computeInverseJac(jacinv,jac,error)
c     Subroutine to calculate the inverse of a 4x4 matrix
      real jacinv(4,4),jac(4,4),tmp(4,4)
      real  piv, linkomb,dtol
      integer i,j,k,error
c     jac(4,4)    ---------  matrix whose inverse will be calculated <Input>
c     jacinv(4,4) ---------  inverse matrix  <Output>
c     error       --------- integer showing whether solution has converged <Output>
c                            = 0 solution has converged
c                            =-1 solution has not converged
c     tmp(4,4)    --------- temporary matrix
c     piv,lincomb --------- variables used in gaussian elimination
c     i,j,k       --------- counters
c     dtol        --------- tolerance of the algorithm

      dtol=1.e-20           
      do i=1,4
         do j=1,4
            tmp(i,j)=jac(i,j)
            jacinv(i,j)=0.
         enddo
         jacinv(i,i)=1.
      enddo
      
      do i=1,3
         piv = tmp(i, i)
         if (abs(piv) .lt. dtol) then
            error=-1
            goto 452
         endif

         do j=i+1,4
            linkomb = tmp(j, i) / tmp(i, i)
            do k=i,4
               tmp(j, k) = tmp(j, k) - tmp(i, k) * linkomb
            enddo
            do k=1,4
               jacinv(j, k) =jacinv(j, k)-jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=4,2,-1
         piv = tmp(i, i)
         do j=i-1,1,-1
            linkomb = tmp(j, i) / piv
            do k=i,1,-1
               tmp(j, k) =  tmp(j, k) -tmp(i, k) * linkomb
            enddo
            do k=4,1,-1
               jacinv(j, k) =jacinv(j, k) - jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=1,4
         do j=1,4
            jacinv(i, j) = jacinv(i, j) / tmp(i, i)
         enddo
      enddo 
      error=0
 452  continue
      return
      end

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

            subroutine umat50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject,nlq,lq1)
c      include 'nlqparm'
      implicit real(a-h,o-z)
      dimension d1(1),d2(1),d3(1),d4(1),d5(1),d6(1)
      dimension sig1(1),sig2(1),sig3(1),sig4(1),sig5(1),sig6(1)
      dimension cm(24),epsps(1),hsvs(1,27),dt1siz(1)
      dimension temps(1),crv(lq1,2,1),cma(1),qmat(nlq,3,3),elsizv(1)
      integer nnpcrv(1)
      integer idelev(1)
      logical failels(1),reject
      character*5 etype
      real eps(6);
      integer maxnip,nnm1,lft,llt

      integer mx,i,j
c
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      real sig(6),sigEff(6),oldStrain(6),convStrain(6),
     $     deltaTotStrain(6),
     $     tempTotStrain(6),elStrain(6),princStress(3),
     $     princDir(3,3),totStrain(6),plastStrain(6),
     $     sum,tempTheta,
     $     strain(6), sigVTrial,rhoTrial,thetaTrial,apexStress,
     $     tempkappaP,yieldval
      integer subincCounter,rtype,subincFlag,converged,l,k

      real alpha,omegaT,omegaC,effStressT(6),effStressC(6),
     $     rateFactor,cdpm2u_computeRateFactor,strainrate(6),epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     length,deltaElStrain(6),stressOld(6), damagedGP
c$omp threadprivate (/cdpmc/)
      ym=cm(1)
      pr=cm(2)
      ecc=cm(3)
      qh0=cm(4)
      ft=cm(5)
      fc=cm(6)
      hp=cm(7)
      ah=cm(8)
      bh=cm(9)
      ch=cm(10)
      dh=cm(11)
      as=cm(12)
      df=cm(13)
      eratetype=cm(14)
      type=cm(15)
      bs=cm(16)
      wf=cm(17)
      wf1=cm(18)      
      ft1=cm(19)
      sratetype=cm(20)
      failflg=cm(21)
      efc=cm(22)
      damageflag=cm(23)
      printflag=cm(24)
c This is the global tolerance. All other tolerances are made relative to this global tolerance
      gTol=1.e-6
      damagedGP=0.

c This function needs allows us to run the code without having to enter all 
c default parameters
      call cdpm2u_giveDefaultValuesIfNotGiven()
      
      cm(24) = printflag;
      
      m0 = 3. * ( fc** 2. - ft**2. ) / ( fc * ft ) * ecc / ( ecc + 1. );
      do i=lft,llt
         tempkappaP=hsvs(i,1)
c Write strain increment vector
         eps(1) = d1(i)
         eps(2) = d2(i)
         eps(3) = d3(i)
         eps(4) = d4(i)
         eps(5) = d5(i)
         eps(6) = d6(i)

         do l=1,6
            totStrain(l)   = eps(l) + hsvs(i,l+20)
            strainrate(l)  = eps(l)/dt1siz(i)
            plastStrain(l) =  hsvs(i,l+2)
            oldStrain(l)  = hsvs(i,l+20)
            convStrain(l) = oldStrain(l)
            tempTotStrain(l) = totStrain(l)
            deltaTotStrain(l) = eps(l)
         enddo

         subincCounter=0
         subincFlag=0
         
         converged=1
                 
         do while ( converged .eq. 1 .or. subincFlag .eq. 1 ) 
            
            do l=1,6
               elStrain(l) = tempTotstrain(l) - plastStrain(l)    
            enddo   

            call cdpm2u_computeStressesfromStrains(sigEff,elStrain,
     $           ym,pr)
            
            call cdpm2u_computeTrialCoordinates(sigEff,sigVTrial,
     $           rhoTrial,tempTheta)
            thetaTrial=tempTheta
            call cdpm2u_computeYieldValue(yieldval,sigVTrial,rhoTrial,
     $           thetaTrial,tempKappaP)

            apexStress=0.

            if (yieldval .gt. 0.) then
               call cdpm2u_checkForVertexCase(apexStress,sigVTrial,
     $              tempKappaP,rtype)
               if (rtype.eq.1 .or. rtype .eq. 2) then
                  call cdpm2u_performVertexReturn(sigEff,
     $                 apexStress,tempKappaP,rtype,converged)
               end if
               
               if (rtype.eq.0) then
                  call cdpm2u_performRegularReturn(sigEff, 
     $                 tempKappaP,converged,ym,pr,gTol)
               end if
            else
               converged=0
               do l=1,6        
                  plastStrain(l)=hsvs(i,l+2)
               enddo
               goto 925
            end if                     
            
            if ( converged .eq. 1 ) then
               subincCounter=subincCounter+1
               if ( subincCounter .gt. 10 ) then
                  write(*,*) '*** Perform Plasticity return with' 
                  write(*,*) 'subincrementation methodology'
                  write(*,*) 'No convergence reached !***'
                  stop
               else if (subincCounter .gt. 9 .and. 
     $                 tempKappaP .lt. 1.0 ) then
                  tempKappaP=1.
               end if
               subIncFlag = 1
               do l=1,6
                  deltaTotStrain(l)=deltaTotStrain(l)*0.5
                  tempTotStrain(l)=convStrain(l)+deltaTotStrain(l)
               enddo
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 0) then
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $          ym,pr)
               do l=1,6        
                  plastStrain(l)=totStrain(l)-elStrain(l)
               enddo
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 1) then
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $              ym,pr)
               do l=1,6
                  plastStrain(l)=tempTotStrain(l)-elStrain(l)
                  convStrain(l)=tempTotStrain(l)
                  deltaTotStrain(l)=totStrain(l)-convStrain(l)
                  tempTotStrain(l)=totStrain(l)
               enddo
               subincCounter = 0
               subincFlag=0
               converged=1            
            end if
         end do
         
                     
 925     continue
         if (damageflag .eq. 3.0) then
            omegaT=0.0
            omegaC=0.0            
            epsilonT=0.0
            kappaDT=0.0
            kappaDT1=0.0
            kappaDT2=0.0
            kappaDC=0.0
            kappaDC1=0.0
            kappaDC2=0.0
            omegaT=0.0
            omegaC=0.0
            rateFactor=0.0
            alpha=0.0
            epsilonT=0.0
            epsilonC=0.0
            do l=1,6
c               sig(i,l)=sigEff(l)
               sig(l)=sigEff(l)
            enddo            
            goto 152
         end if

c Initialize parameters used in the damage algorithm              
         rateFactor=hsvs(i,17)

         epsilonT=hsvs(i,19)
         epsilonC=hsvs(i,20)
         kappaDT=hsvs(i,9)
         kappaDT1=hsvs(i,10)
         kappaDT2=hsvs(i,11)
         kappaDC=hsvs(i,12)
         kappaDC1=hsvs(i,13)
         kappaDC2=hsvs(i,14)
         omegaT=hsvs(i,15)
         omegaC=hsvs(i,16)
         alpha=hsvs(i,18)
         call cdpm2u_computeAlpha(effStressT,effStressC,sigEff,alpha)
         length = elsizv(i)
         sum=0.
         do l=1,6
c     Compute norm of increment of plastic strains       
            sum=sum+(plastStrain(l)-hsvs(i,l+2))**2
            deltaElStrain(l)=(hsvs(i,l+20)-hsvs(i,l+2))
         enddo
         call cdpm2u_computeStressesfromStrains(stressOld,deltaElStrain,
     $        ym,pr)
         sum=sqrt(sum)
         rateFactor=hsvs(i,17)

         call cdpm2u_computeDamage(omegaC,omegaT,strainrate,
     $        rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,
     $        kappaDT2,kappaDC,kappaDC1,kappaDC2,sigEff,sum,
     $        tempKappaP,length,stressOld,hsvs(i,18),hsvs(i,27))
         
         do l=1,6
            if (damageflag .eq. 0.0) then
               sig(l)=(1.-omegaT)*effStressT(l)+
     $              (1.-omegaC)*effStressC(l)               
            else if (damageflag .eq. 1.0) then 
                sig(l)=(1.-omegaT)*sigEff(l)
           else if (damageflag .eq. 2.0) then
               sig(l)=1.-(1.-omegaT*(1-alpha))*
     $              (1-omegaC*alpha)*sigEff(l)
            end if
         enddo

 152     continue
         
c     Write the history variable at the end of the routine
         epsps(i)= tempkappaP 
         hsvs(i,1)= tempkappaP 
         hsvs(i,2)= epsilonT

         hsvs(i,3)= plastStrain(1)
         hsvs(i,4)= plastStrain(2)
         hsvs(i,5)= plastStrain(3)
         hsvs(i,6)= plastStrain(4)
         hsvs(i,7)= plastStrain(5)
         hsvs(i,8)= plastStrain(6)
         
         hsvs(i,9)=  kappaDT
         hsvs(i,10)= kappaDT1
         hsvs(i,11)= kappaDT2
         hsvs(i,12)= kappaDC
         hsvs(i,13)= kappaDC1
         hsvs(i,14)= kappaDC2
         hsvs(i,15)= omegaT
         hsvs(i,16)= omegaC
         hsvs(i,17)=rateFactor
         hsvs(i,18)=alpha
         hsvs(i,19)=epsilonT
         hsvs(i,20)=epsilonC

         hsvs(i,21)=totStrain(1)
         hsvs(i,22)=totStrain(2)
         hsvs(i,23)=totStrain(3)
         hsvs(i,24)=totStrain(4)
         hsvs(i,25)=totStrain(5)
         hsvs(i,26)=totStrain(6)

         if ( hsvs(i,15).gt. 0.9995 .and. hsvs(i,16).gt. 0.9995 ) then
            damagedGP=damagedGP+1.0
         end if

         if (failflg .gt. 0.0) then
            damagedGP=damagedGP/FLOAT(maxnip)
            if (damagedGP .ge. failflg) then
               failur=1
                  failels(i)=1
                  do l=1,6
                     sig(l)=0.0
                  enddo
            end if
         end if
         
c     Write sig components
         sig1(i) = sig(1)
         sig2(i) = sig(2)
         sig3(i) = sig(3)
         sig4(i) = sig(4)
         sig5(i) = sig(5)
         sig6(i) = sig(6)
         enddo
c       sig(nlq,6)   --- nominal stress in voigt notation xx,yy,zz,xy,yz,xz
c       epsps(nlq)   --- cummulative plastic strain (kappa)
c       hsvs(nlq,35) - history variables
c       hsvs(nlq,1) --- kappa
c       hsvs(nlq,2) --- equivalent strain
c       hsvs(nlq,3) --- plastic strain xx
c       hsvs(nlq,4) --- plastic strain yy
c       hsvs(nlq,5) --- plastic strain zz
c       hsvs(nlq,6) --- plastic strain xy
c       hsvs(nlq,7) --- plastic strain zy
c       hsvs(nlq,8) --- plastic strain xz
c       hsvs(nlq,9) --- kappa tension kdt
c       hsvs(nlq,10) -- kappa tension 1 kdt1
c       hsvs(nlq,11) -- kappa tension 2 kdt2
c       hsvs(nlq,12) -- kappa compression kdc
c       hsvs(nlq,13) -- kappa compression 1 kdc1
c       hsvs(nlq,14) -- kappa compression 2 kdc2
c       hsvs(nlq,15) -- damage variable tension omegaT (in LS-DYNA omegaT is written as wt)
c       hsvs(nlq,16) -- damage variable tension omegaC (in LS-DYNA omegaC is written as wc)
c       hsvs(nlq,17) -- strain rate factor used for the dynamic formulation based on quasi-static analysis  
c       hsvs(nlq,18) -- alphac is the compression factor given in paper in IJSS by Grassl et al. in equation 46
c       hsvs(nlq,19) -- equivalent strain tension eqstrT
c       hsvs(nlq,20) -- equivalent strain compression eqstrC
c       hsvs(nlq,21) -- total strain along xx
c       hsvs(nlq,22) -- total strain along yy
c       hsvs(nlq,23) -- total strain along zz
c       hsvs(nlq,24) -- total strain along xy
c       hsvs(nlq,25) -- total strain along yz
c       hsvs(nlq,26) -- total strain along xz
c       hsvs(nlq,27) -- equivalent strain (without rate factor influence)

c We assume that elen contains the length associated with one 
c integration point. This needs to be improved for triangular elements.
         
c       cm(1)   YM (Youngs modulus)
c	cm(2)	PR (Poissons ratio)
c	cm(3)	ECC (Eccentricity)
c	cm(4)	QH0 (Initial hardening)
c	cm(5)	FT (Uniaxial tension strength)
c	cm(6)	FC (Uniaxial compression strength)
c	cm(7)	HP (Hardening parameter)
c	cm(8)	AH (Hardening ductility measure)
c	cm(9)	BH (Hardening ductility measure)
c	cm(10)	CH (Hardening ductility measure)
c	cm(11)	DH (Hardening ductility measure)
c	cm(12)	AS (Damage ductility measure)
c	cm(13)	DF (Dilation constant)
c       cm(14)  ERT (Energy strain rate type)
c		     = 0.0: Rate does not affect fracture energy        
c                    = 1.0: Rate effect on fracture energy 
c                    = 2.0: Square of rate effect on fracture energy
c       cm(15)	TYPE (tensile damage type)
c		     = 0.0: Linear softening     
c                    = 1.0: Bi-Linear softening
c                    = 2.0: Exponential softening
c	cm(16)	BS (Damage: ductility parameter) 
c	cm(17)	WF (Damage: disp threshold 0)
c       cm(18)	WF1 (Damage: disp threshold 1)
c       cm(19)	FT1 (Damage: stress threshold 1)
c       cm(20)	SRT (Strength strain rate type)
c		     = 0.0: No rate effects        
c                    = 1.0: Model code 2010 first branch only 
c                    = 2.0: Model code 2010 first and second branch
c       cm(21)	FAILFLG 
c		     = 0.0: if not ALL  gausspoints of the element are damaged        
c                    = 1.0: if ALL gausspoints of the element are damaged
c     cm(22)	efc (Damage Compression: Strain/displacement threshold )
c     cm(23)	DAMFLG (damage flag)
c     cm(24)	pflag (print flag)
c     bqs(nlq)  - pressure at each gauss point
c     maxnip    ---- variable denoting max number of integration points
c     ym   --------------- Young's modulus
c     pr   --------------- Poisson's ratio
c     ecc  --------------- eccentricity parameter used in the plasticity law. Its calibration is described in Jirasek & Bazant (2002)
c     qh0  --------------- value of the derivative of the 1st hardening function at kappa=0 as described in eq. (30) of the IJSS paper by Grassl et al.
c     ft   --------------- tensile strength
c     fc   --------------- compressive strength
c     hp   --------------- value of the derivative of the 2nd hardening function as described in eq. (31) of the IJSS paper by Grassl et al.
c     ah,bh,ch,dh -------- hardening parameters used in eq. (33) of the IJSS paper by Grassl et al.
c     as,bs,df ----------- softening  parameters used in eqs. (56) and (50) of the IJSS paper by Grassl et al.
c     eratetype ---------  parameter denoting whether and how strain rate dependence is taken into account for fracture energy
c     type --------------- softening type used in the formulation of the tensile damage variable
c     wf   --------------- max crack opening displacement used in the tensile damage law
c     wf1  --------------- max crack opening displacement used in the bilinear tensile damage law
c     ft1  --------------- tensile stress threshold used in the billinear tensile damage law
c     sratetype ---------  parameter denoting whether and how strain rate dependence is taken into account for strength
c     failflg ------------ flag denoting when an element should be deleted
c     efc  --------------- compressive strain/displacement threshold used as a parameter in the compressive damage law
c     m0   --------------- parameter used in the plasticity law calculated in eq.(20) of the IJSS paper by Grassl et al.
c     damageflag ------------ flag to denote whether isotropic or anisotropic damage law is used
c     printflag ---------- flag to print input only once
c
c Variables used for the evaluation of the plasticity algorithm
c       totstrain      -----------  total strain vector equal to the sum of  plastic and elastic strains
c       sigVTrial      ----------- trial volumetric stress 
c       rhoTrial       ----------- trial deviatoric stress
c       thetaTrial     ----------- trial Lode angle
c       apexStress     ----------- apexstress. Used only in the vertex case of the plasticity algorithm
c       yieldval   - value of the yield function
c       subincCounter   - counter of subincrementations performed in the plasticity algorithm
c       subincFlag   - flag denoting whether or not subincrementation is taking place
c                      =0 no subincrementation is taking place
c                      =1 subincrementation is taking place
c       rtype      - return type of the yield surface
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged       - integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     sigEff(6)        ---------- effective stress (no damage included)
c     oldStrain(6)     ---------- old elastic strain 
c     convStrain(6)    ---------- last elastic strain vector for which the plasticity algorithm has converged
c     deltaTotStrain(6)   ---------- difference between last and new strain vector
c     tempTotStrain(6)  ---------- temporary elastic strain
c     elStrain(6)      ----------  elastic strain
c     princStress(3)   ----------  array containing principal stresses
c     princDir(3,3)    ----------  matrix containing principal eigenvectors stored columnwise
c     plastStrain(6)   ----------  plastic strain vector
c     sum              ----------  variable used for summation
c     tempTheta        ---------- temporary Lode angle
c     strain(6)        ---------- strain array
c     l,k              ---------- counters used in various iterations
c
c Variables used in the damage algorithm
c     alpha             --------- variable used to identify amount of contribution of compressive stresses and is defined in eq. (46) of the IJSS paper by Grassl et al.
c     omegaT            --------- tensile damage variable
c     omegaC            --------- compressive damage variable
c     effStressT        --------- effective tensile part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     effStressC        --------- effective compressive part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     rateFactor       ---------- variable used to incorporate impact effects on the constitutive law
c     cdpm2u_computeRateFactor --------- function used to calculate the rate factor
c     strainrate(6)    ---------- rate of the strain tensor used for the calculation of the rateFacto
c     epsilonT         ---------- tensile equivalent strain
c     epsilonC         ---------- compressive equivalent strain
c     kappaDT          ---------- history parameter kappaDT        
c     kappaDT1         ---------- history parameter kappaDT1       
c     kappaDT2         ---------- history parameter kappaDT2       
c     kappaDC          ---------- history parameter kappaDC        
c     kappaDC1         ---------- history parameter kappaDC1       
c     kappaDC2         ---------- history parameter kappaDC2       
c     length           ---------- characteristic length used for the formulation of the tensile damage function based on the crack-band approach
c     deltaElStrain(6) ---------- elastic strains of the previous step
c     damagedGP        ---------- number of damaged gausspoints within analysed element 
c
c ----------------------------------- Variable initialisation -----------------------------------------------
c
c      mx=48*(mxt(lft)-1)
c
c     Material constants
c
      return
      end
      
      subroutine cdpm2u_giveDefaultValuesIfNotGiven()
c     Subroutine to check if all model parameters have values. If a parameter does not have any values a default value is provided.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c
c$omp threadprivate (/cdpmc/)

      real epsilon
      if( pr .lt. 0) then
         pr = 0.2
      end if
      if (wf1 .le. 0.) then
         wf1= 0.15*wf
      end if
      if (ft1 .le. 0.) then
         ft1= 0.3*ft
      end if
      if (qh0 .le. 0.) then
         qh0=0.3
      end if
      if (hp .lt. 0.) then
         hp=0.01
      end if
      if (ah .le. 0.) then
         ah=8.e-2
      end if
      if (bh .le. 0.) then
         bh=3.e-3
      end if
      if (ch .le. 0.) then
         ch=2.
      end if
      if (dh .le. 0.) then
         dh=1.e-6
      end if
      if (df .le. 0.) then
         df=0.85
      end if
      if (type .le. 0.) then
         type=0.
      end if
      if (eratetype .le. 0.) then
         eratetype=0.
      end if
      if (efc .le. 0. ) then
         efc=1.e-4
      end if
      if (as .le. 0.) then
         as=15.
      end if
      if (bs .le. 0.) then
         bs=1.
      end if
      if (sratetype .le. 0.) then
         sratetype=0.
      end if
      if (damageflag .lt. 0.) then
         damageflag=1.
      end if
      if (gTol .le. 0.) then
         gTol=1.e-6
      end if
      if (ecc .le. 0.) then
         epsilon=ft*((1.16*fc)**2.-fc**2.)/(1.16*fc*(fc**2.-ft**2.))
         ecc=(1.+epsilon)/(2.-epsilon)
      end if
      if (printflag .lt. 0.) then
         printflag = 0
      end if

c Write output
c     Print all input variables
c
      if (printflag .lt. 0.5) then
      
      printflag=1.0
      end if
c
      return
      end    


      subroutine cdpm2u_computeStrainsfromStresses(stress,strain,ym,pr)
c     Subroutine to calculate elastic strains from the stress tensor. Performs operation epsilon = D : sigma
      real stress(6),strain(6),ym,pr
c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      strain(1)=(stress(1) - pr * stress(2) - pr * stress(3))/ym
      strain(2)=(-pr*stress(1) + stress(2) - pr * stress(3))/ym
      strain(3)=(-pr*stress(1) - pr * stress(2) + stress(3))/ym
      strain(4)=(2. * (1+pr) * stress(4) )/ym
      strain(5)=(2. * (1+pr) * stress(5) )/ym
      strain(6)=(2. * (1+pr) * stress(6) )/ym
      return
      end
      
      subroutine cdpm2u_computeStressesfromStrains(stress,strain,ym,pr)
c     Subroutine to calculate strains from the elastic strain tensor. Performs operation sigma = C : epsilon

c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      real factor,  stress(6),strain(6),ym,pr
      factor = ym/((1.+pr)*(1.-2.*pr))
      stress(1)=factor*((1.-pr)*strain(1) + pr * strain(2) + 
     $     pr * strain(3))
      stress(2)=factor*(pr*strain(1) + (1.-pr) * strain(2) + 
     $     pr * strain(3))
      stress(3)=factor*(pr*strain(1) + pr * strain(2) + 
     $     (1.-pr)*strain(3))
      stress(4)=factor*(((1.-2.*pr)/2.) * strain(4) )
      stress(5)=factor*(((1.-2.*pr)/2.) * strain(5) )
      stress(6)=factor*(((1.-2.*pr)/2.) * strain(6) )
      return
      end
c ---------------------------------------------------------VERTEX RETURN FUNCTIONS ---------------------------------------------------

      subroutine cdpm2u_checkForVertexCase(apexStress,sigV,tempkappa,
     $     rtype)
c     Subroutine that check whether the current stress state requires plasticity return to the vertex, at which
c     derivative of plastic potential and yield surface are discontinuous.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c
c$omp threadprivate (/cdpmc/)
      real apexStress,sigV,tempkappa,qh2,cdpm2u_qh2fun 
      integer rtype
c     sigV              -----------  volumetric stress <Input>
c     tempKappa         -----------  cummulative plastic strain  <Input>
c     apexStress        -----------  sigmaV of the yield surface for the current 
c                                    tempkappa and rho=theta=0 <Output>
c       rtype           -----------  return type of the yield surface  <Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     qh2               -----------  variables containing the results of the hardening functions 
c     cdpm2u_qh2fun     -----------  function to calculate the hardening function  given in eq. (31) of IJSS paper by P. Grassl et al.

      if ( sigV .gt. 0. ) then
         rtype = 1
         if (tempKappa .lt. 1.) then
            apexStress = 0.
         else
            qh2=cdpm2u_qh2fun(tempKappa,hp)
            apexStress=qh2*fc/m0
         end if        
      else if ( sigV .lt. 0. .and. tempKappa .lt. 1.) then
         rtype = 2
         apexStress = 0.
      else 
         rtype = 0
         apexStress=0.
      end if
      
      return 
      end
      
      
      subroutine cdpm2u_performVertexReturn(stress,apexStress,
     $     kappa,rtype,converged)
c     Subroutine that performs plasticity return close whenever the stress state needs
c     to be returned to the apex. If the stress state is not an actual vertex case
c     rtype=0 is returned.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real sigV,rho,apexStress,kappa,stress(6),theta,
     $     yieldValue,yieldValueMid,sig2,dSig,sigMid,sigAnswer,
     $     ratioPotential,kappa0,tempKappaP, cdpm2u_computeTempKappa,
     $     cdpm2u_computeRatioPotential,ratioTrial,yldTol
      integer i,j,k,rtype,maxiter,converged
c     apexStress     ------------------ variable containing the apex to return <Input>
c     kappa          ------------------ cummulative plastic strain  <Input>
c     rtype          ------------------  return type of the yield surface  <Input/Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged      ------------------ integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     kappa0,tempKappa---------------- variables for cummulative plastic strains
c     yieldValue,yieldValueMid ------- variables with yield values 
c     sigV           ----------------- volumetric stress 
c     rho            ----------------- deviatoric stress 
c     sig2,sigMid,sigAnswer,dSig ----- variables containing volumetric stress units
c     ratioPotential ----------------- variable containing the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     ratioTrial     ----------------- ratio of rho/sigmaV
c     yldTol         ----------------- tolerance used in bisection method solver
c     cdpm2u_computeTempKappa--------- function to calculate tempKappa according to eq.(32) of the IJSS paper by Grassl et al.
c     cdpm2u_computeRatioPotential --- function to calculate the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     maxiter         ---------------- parameter denoting max number of iterations performed using the bisection method
      yldTol=gTol
      yieldValue = 0.
      yieldValueMid = 0.
      sig2 = 0.
      kappa0=kappa
      tempKappaP=kappa
      maxiter=250

      call cdpm2u_computeTrialCoordinates(stress,sigV,rho,theta)

      sig2 = apexStress
      
      tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigV,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValue,sigV, 0., 0., tempKappaP)
      
      tempKappaP =
     $     cdpm2u_computeTempKappa(kappa0, sigV, rho, sig2,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValueMid,sig2, 0.,0.,
     $     tempKappaP)
      
      if ( yieldValue * yieldValueMid .ge. 0. )  then
         converged=1
         rtype = 0  
         goto 501
      end if
      
      if ( yieldValue .lt. 0.0 ) then
         dSig = sig2 - sigV
         sigAnswer = sig2
      else 
         dSig = sigV - sig2
         sigAnswer = sig2
      end if
      
      do  j = 1, maxiter
         dSig = 0.5 * dSig
         sigMid = sigAnswer + dSig
         tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigMid,
     $        ym,pr)
         
         call cdpm2u_computeYieldValue(yieldValueMid,sigMid, 0., 0.,
     $        tempKappaP)
         
        if ( yieldValueMid .le. 0. ) then
            sigAnswer = sigMid
         end if
         if (abs(yieldValueMid) .lt. yldTol .and. 
     $        yieldValueMid .le. 0.) then

            ratioPotential =
     $           cdpm2u_computeRatioPotential(sigAnswer, tempKappaP)
            
            ratioTrial = rho / ( sigV - sigAnswer );
            
            if ( ( ( ( ratioPotential .ge. ratioTrial ) .and. 
     $           rtype .eq. 1 ) ) .or.
     $           ( ( ratioPotential .le. ratioTrial ) .and. 
     $           rtype .eq. 2  ) ) then
               goto 500
            else    
               converged=1
               rtype = 0           
               goto 501
            endif
         endif
      enddo
 500  do k = 1, 3
         stress(k) = sigAnswer
         stress(k+3) = 0.
      enddo
      kappa=tempKappaP                          
      converged=0
 501  continue
      return
      end


      real function cdpm2u_computeTempKappa(kappaInitial,sigV1,rho,
     $     sigV2,ym,pr)
c     Function to calculate the tempKappa whenever requested from the performVertexReturn function.
c     TempKappa is calculated according to eq. (32) of the IJSS paper by P. Grassl et al.

      real sigV1,rho,sigV2,kappaInitial,ym,pr,
     $     equivalentDeltaPlasticStrain,kM,gM,ducMeas, 
     $     cdpm2u_computeDucMeas
c     sigV1 -------------- volumetric stress in the previous stress state <Input>
c     sigV2 -------------- volumetric stress in the current stress state <Input>
c     rho   -------------- deviatoric stress  <Input>
c     kappaInitial ------- previous kappaP (cummulative plastic strain) <Input>
c     ym    -------------- Young's modulus  <Input>
c     pr    -------------- Poisson's ratio <Input>
c     kM,gM -------------- bulk and shear moduli
c     equivalentDeltaPlasticStrain  -----  Increase of the plastic strains
c     ducMeas------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas------ function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM = ym / ( 2. * ( 1. + pr ) )

      equivalentDeltaPlasticStrain = sqrt( 1. / 9. *  (( sigV1 - sigV2 ) 
     $     /  kM )** 2.  + (rho / ( 2. * gM ))** 2. )
                                         
      ducMeas = cdpm2u_computeDucMeas(sigV2, 0., 3.141592653589793/3.)

      cdpm2u_computeTempKappa=kappaInitial+equivalentDeltaPlasticStrain/
     $     ducMeas
      return
      end

      real function cdpm2u_computeRatioPotential(sig ,kappa)
c     Function to calculate the ratio of the derivatives of the plastic potential, given in eq.(22) of the IJSS paper by P. Grassl et al. with respect to the deviatoric and volumetric stress respectively.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real AGParam,BGParam,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,R,mQ,kappa,
     $     sig,rho,
     $     dgdsig,dgdrho,Al,Bl
      integer j
c     sig          ---------------- volumetric stress <Input>
c     rho          ---------------- deviatoric stress <Input>
c     kappa        ---------------- cummulative plastic strain kappaP <Input>
c     dgdsig,dgdrho --------------- derivatives of the plastic potential
c     AGParam,BGParam ------------- components of the plastic potential function
c     qh1,qh2          ------------ variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun - functions to calculate the hardening functions  given in eqs. (30) in IJSS paper by Grassl et al.
c     Al,Bl        ---------------- components of the function given in eq.(23) of the IJSS paper by P. Grassl et al.
c     R,mQ         ---------------- variables 
      rho=0.
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =    qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2.*df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(1.5) *rho/fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * (1. - qh1 ) * Bl + 6. ) +
     $     m0 * (qh1**2.) / ( sqrt(6.) * fc )
      
      cdpm2u_computeRatioPotential= dgdrho/dgdsig*3.*(1.-2.*pr ) / 
     $     ( 1. + pr )
      return 
      end

c ------------------------------------ Plasticity Algorithm functions (regular return functions) ---------------------------------
      subroutine cdpm2u_performRegularReturn(stress,kappa,converged,ym,
     $     pr,yieldTol)
c     Subroutine to perform regular plasticity return
      
      real stress(6),tempKappa,ym,pr, resid(4), normalisedResid(4),PI,
     $     jacobian(4,4), inverseJac(4,4),increment(4),
     $     deltaLambda,normOfResiduals, unknowns(4),
     $     trialSig,trialRho,trialTheta,kappa,kappaP,tempKappaP,
     $     yieldTol,sig,rho,ddkappadDeltaLambdadInv(2),
     $     dgdInv(2),ddgddInv(2,2), dkappadDeltaLambda,
     $     dfdkappa, ddgdInvDkappa(2), ddkappadDeltaLambdadKappa,
     $     ddkappaddDeltaLambdadInv(2),kM,gM,dfdInv(2),sum,
     $     stressPrincipal(3),princDir(3,3)
      integer i,j,iterations,totIter,converged,error
      
c     stress(6)      ------------------  effective stress components are in order xx,yy,zz,xy,yz,xz  <Input/Output>
c     kappa          ------------------  cummulative plastic strain (kappa) <Input>
c     ym             ------------------  Young's modulus  <Input>
c     pr             ------------------  Poisson's ratio <Input>
c     yieldTol       ------------------ tolerance of the N-R solver <Input>
c     converged      ------------------ integer showing whether solution has converged <Output>
c                                        = 0 solution has converged
c                                        = 1 solution has not converged
c     resid(4)       ------------------ array containing the residuals of the N-R solver
c     normalisedResid(4) -------------- array containing the normalised residuals of the N-R solver
c     PI             ------------------ parameter defined as the pi=3.14....
c     jacobian(4,4)  ------------------ matrix containing the jacobian of the problem in order to calculate the solution
c     inversJac(4,4) ------------------ matrix containing the inverse matrix of the jacobian of the problem in order to calculate the solution
c     increment(4)   ------------------ array containing the increments of the unknowns at each N-R iteration
c     deltaLambda    ------------------ plastic multiplier
c     normOfResiduals ----------------- norm the array resid(4)
c     unknowns(4)    ------------------ array with the unknowns in order sigV,rho,kappa,deltaLambda
c     trialSig       ------------------ trial volumetric stress (initial guess)
c     trialRho       ------------------ trial deviatoric stress (initial guess)
c     trialTheta     ------------------ trial Lode angle (initial guess)
c     kappaP          ------------------ cummulative plastic strain kappa (initial guess)
c     tempKappa      ------------------ temporary cummulative plastic strain kappa (each iteration)
c     sig            ------------------ temporary volumetric stress (each iteration)
c     rho            ------------------ temporary deviatoric stress (each iteration)
c     ddkappadDeltaLambdadInv(2) ------ derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric stress
c     dfdInv(2)      ------------------ derivative of the yield function with respect to the volumetric and deviatoric stress
c     dgdInv(2)      ------------------ derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     ddgddInv(2,2)  ------------------ second derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     dkappadDeltaLambda -------------- derivative of kappa with respect to the plastic multiplier
c     dfdkappa       ------------------ derivative of the yield function with respect to kappa
c     ddgdInvDkappa(2) ---------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and to kappa
c     ddkappadDeltaLambdadKappa ------- derivative of kappa with respect to the plastic multiplier and kappa
c     ddkappaddDeltaLambdadInv(2) ----- derivative of kappa with respect to the plastic multiplier and deviatoric and volumetric stress
c     kM             ------------------ bulk modulus
c     gM             ------------------ shear modulus
c     sum            ------------------ variable used for summation
c     stressPrincipal(3)  ------------- array containing the principal stresses 
c     princDir(3,3)  ------------------ matrix containing eigenvectors of the effective stress tensor stored columnwise 
c     i,j,iterations ------------------ integers used as counters
c     totIter        ------------------ maximum number of iterations of the N-R algorithm
c     error          ------------------ integer indicating whether the inversion of the jacobian matrix was successful
      iterations=0
      totIter=100

      PI=3.1415926535897932384626433832795029
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM =  ym / ( 2. * ( 1. + pr ) )
      normOfResiduals=1.
      do i=1,4
         resid(i)=0.
         normalisedResid(i)=0.
         unknowns(i)=0.
         increment(i)=0.
      enddo
      deltaLambda=0.
      call cdpm2u_computePrincValues(stress,stressPrincipal,0,princDir)
      call cdpm2u_computeTrialCoordinates(stress,trialSig,trialRho,
     $     trialTheta)     
      kappaP=kappa
      tempKappaP=kappa
      sig=trialSig
      rho=trialRho
      unknowns(1)=trialSig
      unknowns(2)=trialRho
      unknowns(3)=tempKappaP
      unknowns(4)=0.
      
      call cdpm2u_computeYieldValue(resid(4),sig, rho,trialTheta,
     $     tempKappaP)
      normOfResiduals=1.
      do while (normOfResiduals .gt. yieldTol)
c      write(*,*) 'normOfResiduals = ',normOfResiduals
         iterations=iterations+1
         if (iterations .eq. totIter) then
            converged=1
            goto 600
         end if 
         normalisedResid(1)=resid(1)/kM
         normalisedResid(2)=resid(2)/2./gM
         normalisedResid(3)=resid(3)
         normalisedResid(4)=resid(4)
         normOfResiduals=sqrt(normalisedResid(1)**2.+
     $        normalisedResid(2)**2.+normalisedResid(3)**2. +
     $        normalisedResid(4)**2.)
         
         if (isnan(normOfResiduals)) then
            converged=1
            goto 600
         end if

         if (normOfResiduals .gt. yieldTol) then
c     ----------------- compute jacobian ---------------------------------
           call cdpm2u_computedfdInv(dfdInv,sig,rho,trialTheta,
     $           tempKappaP) 
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgddInv(ddgddInv,sig,rho,trialTheta,
     $          tempKappaP)
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP)
           call cdpm2u_computedfdKappa(dfdkappa,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $          trialTheta,tempKappaP)
           call cdpm2u_computeddKappadDeltaLambdadKappa(
     $          ddkappadDeltaLambdadKappa,sig,rho,tempKappaP,trialTheta)

           call cdpm2u_computeddKappadDeltaLambdadInv(
     $          ddKappaddDeltaLambdadInv,sig,rho,tempKappaP,trialTheta)

           jacobian(1,1) = 1. + kM * deltaLambda *   ddgddInv(1, 1)
           jacobian(1, 2) = kM * deltaLambda * ddgddInv(1, 2)
           jacobian(1, 3) = kM * deltaLambda * ddgdInvdKappa(1)
           jacobian(1, 4) = kM * dgdInv(1)
           
           jacobian(2, 1) = 2. *gM *deltaLambda *ddgddInv(2, 1)
           jacobian(2, 2) = 1. + 2. *gM *deltaLambda *  ddgddInv(2, 2)
           jacobian(2, 3) = 2. *gM *deltaLambda * ddgdInvdKappa(2)
           jacobian(2, 4) = 2. *gM *dgdInv(2)
           
           jacobian(3, 1) = deltaLambda * ddKappaddDeltaLambdadInv(1)
           jacobian(3, 2) = deltaLambda * ddKappaddDeltaLambdadInv(2)
           jacobian(3, 3) = deltaLambda * ddkappadDeltaLambdadKappa - 1.
           jacobian(3, 4) = dkappadDeltaLambda
           
           jacobian(4, 1) = dfdInv(1)
           jacobian(4, 2) = dfdInv(2)
           jacobian(4, 3) = dfdKappa
           jacobian(4, 4) = 0.

           call cdpm2u_computeInverseJac(inverseJac,jacobian,error)
           if (error.eq.-1) then
              converged=1
              goto 600
           end if
           do i=1,4
              sum = 0.
              do j = 1,4
                 sum =sum+ inverseJac(i, j) * resid(j)
              enddo
              increment(i) =0.-sum
              unknowns(i)=unknowns(i)+increment(i)
           enddo
           if (unknowns(4) .le. 0.) then
              unknowns(4)=0.
           end if
           if (unknowns(2) .le. 0.) then
              unknowns(2)=0.
           end if
           if (unknowns(3)-kappaP .le. 0.) then
              unknowns(3)=kappaP
           end if
           sig = unknowns(1)
           rho = unknowns(2)
           tempKappaP=unknowns(3)
           deltaLambda=unknowns(4)
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP)
           resid(1) = sig - trialSig + kM *deltaLambda * dgdInv(1)
           resid(2) = rho - trialRho +  2.* gM *deltaLambda * dgdInv(2)
           resid(3) = -tempKappaP +kappaP+deltaLambda*dkappadDeltaLambda
           call cdpm2u_computeYieldValue(resid(4),sig,rho,trialTheta,
     $          tempKappaP)
        end if
      end do
      converged=0
      
      
      stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(trialTheta)
      stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta - 2. * PI/ 3.)
      stressPrincipal(3) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta + 2. * PI / 3.)

      call cdpm2u_transformStressVectorTo(stress,princDir,
     $     stressPrincipal)
      
      kappa=tempKappaP
      write(*,*) 'cin50',stressPrincipal
      write(*,*) 'cin50',kappa
 600  continue
      return
      end
      
c ------------------------------------ Plasticity Algorithm functions (general functions) ----------------------------------------
      subroutine cdpm2u_computedfdInv(dfdInv,sig,rho,theta,kappa) 
c     Subroutine to calculate the derivative of the yield function with respect to the volumetric and deviatoric stresses respectively 
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)      
      real rFunction,qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     dfdsig,dfdrho,sig,dfdInv(2),Bl
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdInv(2)        -------------  derivatives of yield function with respect to volumeric and deviatoric stress <Output>
c     dfdsig           -------------  derivative of the yield function with respect to volumetric stress
c     dfdrho           -------------  derivative of the yield function with respect to deviatoric stress
c     Al,Bl            -------------  variables corresponding to components of the yield function
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.       
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      
      
      Al =( 1. - qh1 ) * ( sig / fc + rho / ( sqrt(6.) * fc ) )** 2. +
     $     sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      
      dfdsig= 4. * ( 1. - qh1 ) / fc * Al * Bl + qh2* qh1**2. * m0 / fc
      dfdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6. )+ 
     $     rFunction * m0 * qh2 * qh1**2./ ( sqrt(6.) * fc )
      
      dfdInv(1) = dfdsig
      dfdInv(2) = dfdrho
      return
      end
      
      subroutine cdpm2u_computedgdInv(dgdInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the plastic potential function with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)      
      real qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,cdpm2u_qh2fun,kappa,
     $     Bl,AGParam,BGParam,R,mQ,dgdsig,dgdrho,sig,dgdInv(2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2)        -------------  derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     dgdsig           -------------  derivative of the plastic potential with respect to volumetric stress
c     dgdrho           -------------  derivative of the plastic potential with respect to deviatoric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.      
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) /  ( log(AGParam)+
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6.) + 
     $     m0 * qh1**2. / ( sqrt(6.) * fc )
      
      dgdInv(1) = dgdsig
      dgdInv(2) = dgdrho
      return 
      end 
      
      subroutine cdpm2u_computeddgddInv(ddgddInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the derivatives of the plastic potential with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)      
      real qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa, Bl,AGParam,
     $     BGParam,R,mQ,sig, dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho,
     $     ddgddSig,ddgddRho,ddgdSigdRho,ddgdRhodSig,ddgddInv(2,2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2,2)      -------------  derivatives of derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     ddgddsig         -------------  second derivative of the plastic potential with respect to volumetric stress
c     ddgddrho         -------------  second derivative of the plastic potential with respect to deviatoric stress
c     ddgdsigdrho      -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     ddgdrhodsig      -------------  derivative of the plastic potential with respect to deviatoric and volumetric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho, ---  variables corresponding to derivatives of components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.            
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) / ( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      dMQDSig = AGParam / ( BGParam * fc ) * exp(R)
      dAlDSig = 2. * ( 1. - qh1 ) * Bl / fc
      dBlDSig = 1. / fc
      dAlDRho = 2. * ( 1. - qh1 ) * Bl / ( fc * sqrt(6.) ) + 
     $     sqrt(3. / 2.) / fc;
      dBlDRho = 1. / ( fc * sqrt(6.) )
      
      ddgddSig = 4. * ( 1. - qh1 ) / fc * ( dAlDSig * Bl + Al * 
     $     dBlDSig ) +qh1**2. * dMQDSig / fc
      ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * Bl + 6. ) +Al * dBlDRho * 4. *
     $     ( 1. - qh1 ) / ( sqrt(6.) * fc )
      ddgdSigdRho = 4. * (1. - qh1 )/fc *( dAlDRho * Bl + Al * dBlDRho )
      ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - 
     $     qh1 ) * Bl + 6. ) + Al / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * dBlDSig )
      
      ddgddInv(1, 1) = ddgddSig
      ddgddInv(1, 2) = ddgdSigdRho
      ddgddInv(2, 1) = ddgdRhodSig
      ddgddInv(2, 2) = ddgddRho
      
      return
      end
      
      subroutine cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,
     $     sig,rho,theta,kappa)
c     Subroutine to calculate the derivative of the hardening variable with respect to plastic multiplier
      real rho,sig,theta,kappa,dkappadDeltaLambda,equivalentDGDStress,
     $     ductilityMeasure,dgdInv(2),cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     dkappadDeltaLambda -----------  derivative of the hardening variable with respect to plastic multiplier <Output>
c     dgdInv(2)        -------------  derivative of plastic potential with respect to volumeric and deviatoric stress
c     ductilityMeasure -------------  ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- norm of the derivative of the plastic potential with respect to volumetric and deviatoric stress
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      equivalentDGDStress = sqrt( 1. / 3.*dgDInv(1)** 2.+dgdInv(2)** 2.)
      ductilityMeasure = cdpm2u_computeDucMeas(sig, rho,theta)
      dkappadDeltaLambda = equivalentDGDStress / ductilityMeasure
      return
      end
      
      subroutine cdpm2u_computedfdKappa(dfdkappa,sig,rho,theta,kappa)
c     Subroutine to calculate the derivative of the yield function with respect to the cummulative plastic strain(kappa)
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real dfdkappa,sig,rho,theta,kappa,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,dfdqh1,
     $     dfdqh2,dqh1dkappa,dqh2dkappa,cdpm2u_dqh1dkappaFun,
     $     cdpm2u_dqh2dkappaFun,Al,
     $     Bl,rFunction
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdkappa         -------------  derivative of the yield function with respect to cummulative plastic strain (kappa) <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     Al,Bl            -------------  variables corresponding to components of the yield function
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      Al = ( 1. - qh1 ) * ( ( sig / fc + rho / ( sqrt(6.) * 
     $     fc ) )) **2.  + sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      dfdqh1 = -2. *Al *(Bl** 2.) + 2. * qh1 * qh2 *   m0 * ( sig / fc +
     $     rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *qh1*(qh2**2.)

      dfdqh2 = (qh1**2.) * m0 * ( sig / fc + rho * rFunction /
     $     (sqrt(6.) * fc)) -  2. *qh2 *(qh1** 2.)
      dfdkappa =  dqh1dkappa * dfdqh1 + dqh2dkappa * dfdqh2
      
      if ( dfdkappa .gt. 0. ) then
         dfdkappa = 0.
      end if
      
      return
      end      
      
      subroutine cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $     theta,kappa)
c     Subroutine to calculate the derivative of the plastic potential function with respect to the volumetric and deviatoric stresses and the cummulative plastic strain (kappa)      
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      
      real qh1,cdpm2u_qh1fun,qh2,cdpm2u_qh2fun,dqh1dkappa,
     $     cdpm2u_dqh1dkappaFun,dqh2dkappa,
     $     cdpm2u_dqh2dkappaFun, AGParam,BGParam,R,mQ,dAGParamdKappa,
     $     BGParamTop,BGParamBottom,dBGParamTopDKappa,
     $     dBGParamBottomDKappa,dBGParamDKappa,RTop,RBottom,dRTopDKappa,
     $     dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard,
     $     dDGDSigDKappa,ddgdInvdKappa(2),kappa,sig,rho,theta,
     $     dDGDRhoDKappa
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     ddgdInvdKappa(2) -------------  derivative of the plastic potential with respect to deviatoric and volumetric strains and the cummulative plastic strain (kappa) <Output>
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     dDGDSigDKappa    -------------  derivative of the plastic potential with respect to volumetric stress and cummulative plastic strain (kappa)
c     dDGDRhoDKappa    -------------  derivative of the plastic potential with respect to deviatoric stress and cummulative plastic strain (kappa)
c     Al,Bl,AGParam,BGParam,R,mQ,dAGParamdKappa,BGParamTop,BGParamBottom,dBGParamTopDKappa,  dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard           -------------  variables corresponding to components and their derivatives of the plastic potential
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)
      
      AGParam = ft * qh2 * 3 / fc + m0 / 2
      BGParam =  qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2+ m0 / 2) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      dAGParamDKappa = dqh2dkappa * 3. * ft / fc
      BGParamTop = qh2 / 3. * ( 1. + ft / fc );
      BGParamBottom = ( log(AGParam) + log(df + 1.) - 
     $     log(2 * df - 1.) - log(3. * qh2 + m0 / 2) )
      dBGParamTopDKappa = dqh2dkappa / 3.
      dBGParamBottomDKappa = -3. * dqh2dkappa / ( 3 * qh2 + m0 / 2. )
      dBGParamDKappa =( dBGParamTopDKappa * BGParamBottom - BGParamTop * 
     $     dBGParamBottomDKappa ) / (BGParamBottom**2.)
      RTop = ( sig - ft / 3. * qh2 )
      RBottom = fc * BGParam
      dRTopDKappa = -ft / 3. * dqh2dkappa
      dRBottomDKappa = fc * dBGParamDKappa
      dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / 
     $     (RBottom** 2.)
      dMQDKappa = dAGParamDKappa * exp(R) + AGParam *dRDKappa *exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * (Bl** 2.) + sqrt(3. / 2.) * rho / fc
      dAlDYieldHard = -(Bl** 2.)
      
      dDGDSigDKappa =  ( -4. * Al * Bl / fc + 4. * ( 1 - qh1 ) / fc * 
     $     dAlDYieldHard * Bl ) * dqh1dkappa +
     $     dqh1dkappa * 2 * qh1 * mQ / fc + qh1 * dMQDKappa / fc
      dDGDRhoDKappa =
     $     ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 )* 
     $     Bl + 6. ) - 4. * Al / ( sqrt(6.) * fc ) * Bl + m0/( sqrt(6.)* 
     $     fc ) ) * 2 * qh1 * dqh1dkappa
      
      ddgdInvdKappa(1) = dDGDSigDKappa
      ddgdInvdKappa(2) = dDGDRhoDKappa
      
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadKappa(
     $     ddkappadDeltaLambdadKappa,sig,rho,kappa,theta)
c     Subroutine to compute the derivative of the cummulative plastic strain(kappa) with respect to the plastic multiplier and the hardening variable kappa
      real equivalentDGDStress,dEquivalentDGDStressDKappa,ducMeas,
     $     ddkappadDeltaLambdadKappa,dgdInv(2),ddgdInvdKappa(2),
     $     sig,rho,kappa,theta,cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddkappadDeltaLambdadKappa ----  derivative of the kappa with respect to plastic multiplier and kappa <Output> 
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     ddgdInvdKappa(2) ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and cummulative plastic strain (kappa)
c     dgdInv(2)        ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDKappa ---  scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress and kappa which is the norm of the derivative of the increment of the plastic strains with respect to kappa 
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa, sig, rho,theta,
     $     kappa)      
      equivalentDGDStress =sqrt( 1./ 3.*( dgdInv(1)**2.)+ dGDInv(2)**2.)
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivalentDGDStressDKappa = (2./3.*dgdInv(1) * ddgdInvdKappa(1) +
     $     2. * dgdInv(2) *ddgdInvdKappa(2) ) / 2./equivalentDGDStress
      ddkappadDeltaLambdadKappa= dEquivalentDGDStressDKappa/ducMeas
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadInv(
     $     ddKappaddDeltaLambdadInv, sig,rho,kappa,theta)
c     Subroutine to compute the derivative of the cummulative plastic strain (kappa) with respect to the plastic multiplier and the volumetric and deviatoric stress
      real equivDGDStress, dgdInv(2),ddgddInv(2, 2),
     $     dEquivDGDStressDInv(2),ducMeas,
     $     ddKappaddDeltaLambdadInv(2),sig,rho,kappa,theta,
     $     dDuctilityMeasureDInv(2),cdpm2u_computeDucMeas
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddKappaddDeltaLambdadInv -----  derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric strain <Output> 
c     ddgdd(2,2)       -------------  second derivative of the plastic potential with respect to volumetric and deviatoric stress 
c     dgdInv(2)        -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     dDuctilityMeasureDInv(2) -----  derivative of the ductility measure with respect to volumetric and deviatoric stress
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDInv(2) -- scalar function of the second derivative of the plastic potential with respect to volumetric and deviatoric stress which is the derivative of the norm of the increment of the plastic strains with respect to deviatoric and volumetric plastic strains 
      
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgddInv(ddgddInv, sig, rho,theta, kappa)
      
      equivDGDStress = sqrt( 1. / 3. * (dGDInv(1)** 2.) +
     $     dGDInv(2)**2. )
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivDGDStressDInv(1) =  ( 2. / 3.*dgdInv(1)*ddgddInv(1,1) + 
     $     2. * dgdInv(2)*ddgddInv(2, 1) ) / ( 2. * equivDGDStress)
      dEquivDGDStressDInv(2) = ( 2. / 3.*dgdInv(1)*ddgddInv(1,2) +
     $     2. * dgdInv(2)*ddgddInv(2, 2) ) /( 2. * equivDGDStress )
      call cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, rho,
     $     theta, kappa)
      
      ddKappaddDeltaLambdadInv(1) = ( dEquivDGDStressDInv(1) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(1) ) / (ducMeas** 2.)
      ddKappaddDeltaLambdadInv(2) = ( dEquivDGDStressDInv(2) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(2) ) / (ducMeas** 2.)
      
      return
      end
      
      
      subroutine cdpm2u_computeTrialCoordinates(stress, sigV,rho, theta)
c     Subroutine which returns volumetric and deviatoric stress and the Lode angle
c     based on the given stress tensor
      real stress(6),rho,sigV,theta,tempdevSig(6),j2,j3
      integer i,j
      
c     stress(6)     -------------  Stress array. stress={sigXX,sigYY,sigZZ,sigXY,sigYZ,sigXZ} <Input>
c     sigV          -------------  Volumetric stress. <Output>
c     rho           -------------  Norm of deviatoric stress. <Output>
c     theta         -------------  Lode angle.  <Output>
c     tempdevSig(6) -------------  Array containing  deviatoric stress tensor
c     j2            -------------  Second invariant of the deviatoric stress tensor
c     j3            -------------  Third invariant of the deviatoric stress tensor
c     i,j           -------------  counters used in iterations
      
      sigV= (stress(1)+stress(2)+stress(3))/3.
      do j=1,3
         tempdevSig(j)=stress(j)-sigV
         tempdevSig(j+3)=stress(j+3)
      enddo
      j2=0.5*(tempdevSig(1)*tempdevSig(1)+tempdevSig(2)*tempdevSig(2)+ 
     $     tempdevSig(3)*tempdevSig(3))+ tempdevSig(4)*tempdevSig(4)+ 
     $     tempdevSig(5)*tempdevSig(5)+tempdevSig(6)*tempdevSig(6)
      
      if(j2 .eq. 0.) then
         theta=0.
         rho=0.
      else   
         rho=sqrt(2.*j2)
         j3= (1./3.) * ( tempdevSig(1)*tempdevSig(1)*tempdevSig(1) + 
     $        3.*tempdevSig(1) *tempdevSig(4)*tempdevSig(4)+
     $        3.*tempdevSig(1) * tempdevSig(6) * tempdevSig(6) + 
     $        6. * tempdevSig(5) * tempdevSig(4)  *  tempdevSig(6) +
     $        3. * tempdevSig(2) * tempdevSig(4)**2.+
     $        3 * tempdevSig(3) *tempdevSig(6)*tempdevSig(6)+
     $        tempdevSig(2)*tempdevSig(2)*tempdevSig(2)+ 
     $        3. * tempdevSig(2) * tempdevSig(5)*tempdevSig(5)+
     $        3. * tempdevSig(3)*tempdevSig(5)*tempdevSig(5)+ 
     $        tempdevSig(3)*tempdevSig(3)*tempdevSig(3))
         theta=(3.*sqrt(3.)/2.)*j3/(j2**(3./2.))
      end if
      if (theta .gt. 1.) then
         theta=1.
      else  if (theta .lt. -1.) then
         theta=-1.
      end if
      theta=1./3.*acos(theta)
      return
      end
      
      subroutine cdpm2u_computeYieldValue(answer,sigV, rho,theta,kappa)
c     Function to evaluate the yield function f based on the given stress High -Westergaard coordinates and kappa.
c     The equation is given in eq. (18) of IJSS paper by P. Grassl et al. Returns a real number 
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      
      real rFunction,qh1,qh2, Al,sigV,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     answer

c     sigV             -------------  volumetric stress <Input>
c     rho              -------------  deviatoric stress <Input>
c     theta            -------------  Lode angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     answer           -------------  variable containing the result of the yield function <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.
c     Al               -------------  variable used to siplify and facilitate the calculation of the yield surface
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     1     ( 2. * ecc - 1. )**2.  ) /
     2     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     3     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     4     cos(theta)**2.+ 5. * ecc**2. - 4. * ecc) )                                                              

      Al = ( 1. - qh1 ) * (sigV / fc + rho / (sqrt(6.) * fc) ) ** 2. +
     &     sqrt(1.5) * rho / fc
      
      answer=Al**2. + qh1**2. * qh2* m0 * ( sigV / fc + rho *
     &     rFunction / ( sqrt(6.) * fc ) ) - qh1**2. * qh2**2.      
      return
      end
      
      
      real function cdpm2u_qh1fun(kappa,qh0,hp)
c     Function to calculate the first hardening function given in eq. (30) of the IJSS paper 
c     by Grassl et al.

      real hp,qh0, kappa,answer

c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     qh0    ---------- Initial hardening parameter <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer=qh0
      else if ((kappa .gt. 0.0) .and. (kappa .lt. 1.0)) then
         answer=
     $        (1.0 - qh0 - hp) *(kappa**3.0)-
     $        ( 3.0 * ( 1.0 - qh0 ) - 3.0 * hp ) * (kappa**2.0) +
     $        ( 3.0 * ( 1.0 - qh0 ) - 2.0 * hp ) * kappa + qh0
      else 
         answer=1.
      endif      
      cdpm2u_qh1fun=answer      
      return
      end
      
      real function cdpm2u_qh2fun(kappa,hp)
c     Function to calculate the first hardening function given in eq. (31) of the IJSS paper 
c     by Grassl et al.
      
      real kappa,answer,hp
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      
      if ( kappa .le. 0. ) then
         answer=1.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=1.
      else 
         answer= 1.+(kappa-1.)*hp
      endif      
      cdpm2u_qh2fun=answer     
      return 
      end
      
      real function cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
c     Function to calculate the derivative of the first hardening function, given in eq. (30) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real kappa,answer,hp,qh0
c     kappa  ---------- Cummulative plastic strain <Input> 
c     qh0    ---------- Initial hardening parameter <Input>
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer= 3. * ( 1 - qh0 ) - 2. * hp
      else if ( kappa .ge. 0. .and. kappa .lt. 1. ) then
         answer=   3. * ( 1. - qh0 - hp ) * (kappa**2.)
     $        - 2. * ( 3. * ( 1. - qh0 ) - 3. * hp ) * kappa
     $        + ( 3. * ( 1. - qh0 ) - 2. * hp )
      else 
         answer=  0.
      endif
      cdpm2u_dqh1dkappaFun=answer
      return
      end
      
      real function cdpm2u_dqh2dkappaFun(kappa,hp)
c     Function to calculate the derivative of the second hardening function, given in eq. (31) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real kappa,hp,answer
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      if ( kappa .le. 0. ) then
         answer=0.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=0.
      else 
         answer=hp
      endif
      cdpm2u_dqh2dkappaFun=answer
      return 
      end  


      real function cdpm2u_computeDucMeas( sigV,rho,theta)
c     Function to calculate ductility measure according to eq. (33) of IJSS paper by P. Grassl et al.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)

      real thetaConst,x,  sigV,rho,theta ,eh,fh,answer
c     sigV ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta---------------- Lode angle <Input>
c     eh,fh---------------- hardening parameters
c     thetaConst,x -------- variables
c     answer--------------- variable containing the answer of the function

      thetaConst = (2. * cos(theta))**2.
      x = -( sigV + fc / 3 ) / fc
      if ( x .lt. 0. ) then
         eh = bh - dh
         fh = ( bh - dh ) * ch / ( ah - bh )
         answer = ( eh * exp(x / fh) + dh ) / thetaConst
      else 
         answer = ( ah + ( bh - ah ) * exp( -x / ( ch ) ) ) / thetaConst
      endif
      cdpm2u_computeDucMeas=answer
      return
      end
      
      
      
      subroutine cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, 
     $     rho,theta, kappa)
c     Subroutine to compute the derivative of the ductility measure with respect to volumetric and deviatoric stress
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)      
      real dDuctilityMeasureDInv(2), sig,rho, kappa,theta,x,dXDSig,
     $     EHard,FHard,dDuctilityMeasureDX,theta1
c     sig  ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta --------------- Lode angle <Input>
c     dDuctilityMeasureDInv(2) --  derivative of the ductility measure with respect to volumetric and deviatoric stress <Output>
c     eh,fh --------------- hardening parameters
c     theta1,x ------------ variables
c     EHard,FHard --------- hardening parameters
c     dXDSig -------------- derivative of the function R, given in eq. (34) in IJSS paper by P.Grassl et al., with respect to volumetric stress
c     dDuctilityMeasureDX - derivative of the ductility measure with respect to the function R

      theta1 = (2. * cos(theta))** 2.
      x = ( -( sig + fc / 3. ) ) / fc
      
      if ( x .lt. 0. ) then
         dXDSig = -1. / fc
         EHard = bh - dh
         FHard = ( bh - dh ) * ch / ( ah - bh )
         
         dDuctilityMeasureDX = EHard  / FHard *exp(x / FHard) / theta1
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      else 
         dXDSig = -1. / fc
         dDuctilityMeasureDX = -( bh - ah ) / ( ch ) / theta1 *
     $        exp( -x / ( ch ) )
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      endif
      return
      end
c     --------------------------------------------- Functions used in the damage algorithm ---------
      subroutine cdpm2u_computeAlpha(stressT,stressC,stress,alpha)
      
c     Subroutine calculating the alpha according to eq.46 of the IJSS paper by Grassl et al. and splitting stress tensor in tensile and compressive part in the principal effective stress coordinate system. Finally principal effective tensile and compressive stress tensors are rotated back to the original system

      real stress(6),stressT(6),stressC(6)
      real princDir(3, 3),princStress(3),princStressT(3),alpha,
     $     princStressC(3), squareNormOfPrincipalStress,alphaTension
      integer i
c     stress(6)       -------------  effective stress tensor (no damage included) <Input>
c     stressT(6)      -------------  effective tensile stress tensor (no damage included) <Output>
c     alpha           ------------- alphaC used in the material model <Output>
c     stressC(6)      -------------  effective compressive stress tensor (no damage included) 
c     alphaTension    ------------- 1-alphaC <Output>
c     princDir(3,3)   -------------  matrix containing the principal directions of original stress tensor stored columnwise
c     princStress(3)  ------------- array containing principal stresses {sigma1,sigma2,sigma3} 
c     princStressT(3) ------------- array containing principal tensile stresses {sigmaT1,sigmaT2,sigmaT3} 
c     princStressC(3) ------------- array containing principal compressive stresses {sigmaC1,sigmaC2,sigmaC3} 
c     squareNormOfPrincipalStress - norm of princStress 

      call cdpm2u_computePrincValues(stress,princStress,0,princDir)
      
c     Split the principal values in a tension and a compression part
      do i = 1,3
         if ( princStress(i) .ge. 0. ) then
            princStressT(i) = princStress(i)
            princStressC(i) = 0.
         else 
            princStressC(i) = princStress(i)
            princStressT(i) = 0.
         endif
      enddo
      
c     Transform the tension and compression principal stresses back to the original coordinate system

      call cdpm2u_transformStressVectorTo(stressT, princDir,
     $     princStressT)
      call cdpm2u_transformStressVectorTo(stressC, princDir,
     $     princStressC)
      
c     Determine the two factors from the stress
      squareNormOfPrincipalStress = 0.

         squareNormOfPrincipalStress = princStress(1)**2.+
     $     princStress(2)**2.+princStress(3)**2.

      
      alphaTension = 0.
      
      if ( squareNormOfPrincipalStress .gt. 0. ) then
         do i = 1,3
            alphaTension = alphaTension+ princStressT(i) *
     $           ( princStressT(i) + princStressC(i) ) /
     $           squareNormOfPrincipalStress
         enddo
      endif

      alpha= 1. - alphaTension

      return
      end

      real function cdpm2u_computeRateFactor(alpha,strainrate,fc,
     $     tol,sratetype)
c     Function to incorporate impact effects in the constitutive law and calculate rateFactor. 
c     All functions used are based on Model Code 2010. 
      
      real strainrate(6),princDir(3,3),princStrainRate(3),max,min,tol,
     $     alphaS,gammaS, deltaS,betaS,strainRateTension0,
     $     strainRateCompression0,rate,ratioT,ratioC,rateFactorTension,
     $     rateFactorCompression,alpha,rateFactor,
     $     strainRateRatioCompression,sratetype
      integer k
c     alpha                     --------------   alpha is the variable used in CDPM2U to evaluate contribution of compressive stresses to principal stress tensor  <Input>
c     totstrain(6)             --------------   array containing strain increments {xx,yy,zz,xy,yz,xz} <Input>
c     deltaTime                 --------------    time step  length <Input>
c     oldStrainRate                 --------------    strain rate  <Input/Output>
c     rateFactor                --------------   rateFactor used to incorporate impact effects on the constitutive law <Input>
c     fc                        --------------   concrete compressive strength <Input>
c     tol                       --------------   tolerance used to identify whether it is a tensile or compressive strain state <Input>
c     princDir(3,3)             --------------   matrix containing the principal directions of the strain rate tensor. Eigenvectors stored columnwise
c     princStrainRate(6)        --------------   array containing the 3 eigenvalues of the strain rate tensor. {rate1,rate2,rate3}
c     max                       --------------   max(tensile) principal strain rate
c     min                       --------------   min (compressive) principal strain rate
c     ratioT                    --------------   max strainrate/sig0Tension (see MC90)
c     ratioC                    --------------   max strainrate/sig0Compression (see MC90)
c     rateFactorTension         --------------   contribution of tensile strain rates to the rate Factor
c     rateFactorCompression     --------------   contribution of compressive strain rates to the rate Factor
c     strainRate                 --------------    strain rate  
      rateFactorTension=1.
      rateFactorCompression=1.
      
      call cdpm2u_computePrincValues(strainRate,princStrainRate,1,
     $     princDir)
      
      max= -1.e-20
      min = 1.e20

      do k=1,3
         if (max .lt. princStrainRate(k)) then
            max=princStrainRate(k)           
         end if
         if (min .gt. princStrainRate(k)) then
            min=princStrainRate(k)
         end if
      enddo
      
      if ( 1. - alpha .gt.tol ) then 
         rate = max
      else
         rate =  min
      endif

      ratioT= rate/1.e-6

      if ( rate .lt. 1.e-6 ) then
         rateFactorTension = 1.
      else if ( 1.e-6 .lt. rate .and. rate .lt. 10.) then
         rateFactorTension = ratioT**0.018
      else 
         rateFactorTension = 0.0062 * (ratioT **(1./3.))
      endif

      
      ratioC= rate/(-30.e-6)

      if ( rate .gt. -30.e-6 ) then
         rateFactorCompression = 1.
      else if (-30.e-6 .gt. rate .and. rate .gt. -30) then
         rateFactorCompression = ratioC**(0.014)        
      else 
         rateFactorCompression =  0.012*(ratioC**(1./3.))
      endif

      
      rateFactor = ( 1. - alpha ) * rateFactorTension + 
     $     alpha * rateFactorCompression
      cdpm2u_computeRateFactor=rateFactor

      return
      end

      subroutine cdpm2u_computeDamage(omegaC,omegaT,strainRate,
     $     rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,kappaDT2,
     $     kappaDC,kappaDC1,kappaDC2,stress,deltaPlasticStrainNormT,
     $     tempKappaP,len,stressOld,oldAlpha,epsilon)
c     Subroutine to perform the damage return. Both compressive and tensile damage variables are calculated.
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)

      real omegaC,omegaT,tempRateFactor,rateFactor,epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     stress(6),deltaPlasticStrainNormT,deltaPlasticStrainNormC,
     $     tempKappaP,len,omegaOldC,omegaOldT,alpha,sigElastic,
     $     rhoElastic, thetaElastic,pHelp,help,tempEquivStrain,qhelp,
     $     tempEquivStrainT,tempEquivStrainC,fsT,fsC,ducMeas,Rs,
     $     strainRate(6),yieldTolDamage,cdpm2u_computeRateFactor,
     $     deltaPlasticStrainNormTNew,
     $     cdpm2u_computeDeltaPlasticStrainNormT,
     $     deltaPlasticStrainNormCNew,
     $     cdpm2u_computeDeltaPlasticStrainNormC,
     $     alphaZero,e0,rFunction,
     $     cdpm2u_computeDamageT,cdpm2u_computeDamageC,
     $     cdpm2u_computeEquivalentStrain,stressOld(6),
     $     minEquivStrain,oldAlpha,epsilon
      integer unloadingFlag,step1flag
c     omegaC          --------------- compressive damage variable                        <Input/Output>
c     omegaT          --------------- tensile damage variable                            <Input/Output>
c     totstrain(6)    --------------- array containing the rate of the strain tensor     <Input>
c     strainRate      ---------------  strain rate     <Input/Output>
c     deltaTime       ---------------  time step length     <Input>
c     rateFactor      --------------- variable to incorporate impact effects on the constitutive law  <Input/Output>
c     alpha           --------------- variable showing the contribution of damage on the strains      <Input/Output>
c     epsilonT        --------------- tensile equivalent strain                <Input/Output>
c     epsilonC        --------------- compressive equivalent strain            <Input/Output>
c     epsilon         ---------------  old equivalent strain (no rate influnce)            <Input/Output>
c     kappaDT         --------------- history parameter kappaDT                <Input/Output>
c     kappaDT1        --------------- history parameter kappaDT1               <Input/Output>
c     kappaDT2        --------------- history parameter kappaDT2               <Input/Output>
c     kappaDC         --------------- history parameter kappaDC                <Input/Output>
c     kappaDC1        --------------- history parameter kappaDC1               <Input/Output>
c     kappaDC2        --------------- history parameter kappaDC2               <Input/Output>
c     stress(6)       --------------- stress tensor                            <Input>
c     deltaPlasticStrainNormT ------- norm of the increment of plastic strains <Input>
c     tempKappaP      --------------- cummulative plastic strain (kappaP)      <Input>
c     deltaElStrain(6)  ---------------  increment of the elastic strain vector  <Input>
c     oldAlpha        ---------------  old alpha  <Input>
c     len             --------------- characteristic length used to combine damage law with the crack-band approach <Input>
c     deltaPlasticStrainNormTNew ---- norm of the increment of the cummulative plastic strain used in the damage algorithm for tension
c     cdpm2u_computeDeltaPlasticStrainNormT- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the tensile damage variable
c     deltaPlasticStrainNormCNew ---- norm of the increment of the cummulative plastic strains used in the damage algorithm for compression
c     cdpm2u_computeDeltaPlasticStrainNormC- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the compressive damage variable
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     alphaZero       --------------- parameter used for the calculation of the damage ductility measure based on eqs.(56-57) of the IJSS paper by Grassl et al.
c     e0              --------------- parameter equal to ft/E
c     cdpm2u_computeDamageT,cdpm2u_computeDamageC - functions used to calculate tensile and compressive damage variables respectively
c     yieldTolDamage  --------------- tolerance used in the damage algorithm for cased where Hp=0
c     Rs              --------------- variable used in the calculation of the ductility measure for damage given in eq.(57) of the IJSS paper by Grassl et al.
c     xs              --------------- ductility measure for damage given in eq.(56) of the IJSS paper by Grassl et al.
c     fsT,fsC         --------------- loading function for the tensile and compressive damage variables
c     cdpm2u_computeRateFactor ------------- function used for the calculation of the rate factor to include impact effects on the constitutive law
c     omegaOldT,omegaOldC ----------- old (previous step) tensile and compressive damage variables
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression)
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     minEquivStrain    -------------- when unloading is occuring corresponds to the minEquivStrain before reloading occurs
c     step1flag    -------------- flag indicating whether we are at the first step
c                                        =1 it is analysis 1st step
c                                        =0  it not analysis 1st step

      e0=ft/ym
      yieldTolDamage=gTol*10.0
      omegaOldC=omegaC
      omegaOldT=omegaT
      deltaPlasticStrainNormC=deltaPlasticStrainNormT
      if (rateFactor .eq. 0.) then
         step1flag=1
         rateFactor=1.
      else
         step1flag=0
      end if
      call cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,tempEquivStrain,epsilon,ym,pr,
     $     gtol)
      
c-------------------Compute tensile and compressive equivalent strains-------------------------------------
      if (strrateflg .gt. 0.0 .and. omegaC .eq. 0. .and. 
     $     omegaT.eq. 0. ) then
         tempRateFactor=cdpm2u_computeRateFactor(alpha,strainrate,
     $        fc,gTol,strrateflg)
      else 
         tempRateFactor=rateFactor
      end if
      tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/
     $     tempRateFactor
      
      if (unloadingFlag .eq. 0) then
         tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $        tempRateFactor
      else
          tempEquivStrainC=epsilonC+ oldAlpha*(minEquivStrain-epsilon)/
     $        tempRateFactor + alpha*(tempEquivStrain-minEquivStrain)/
     $         tempRateFactor
      end if
   
c     Note rate factor is calculated only once at the onset of damage
      if ( ( tempEquivStrainT .gt. e0 .or. tempEquivStrainC .gt. e0
     $     ) .and. ( ( omegaT .eq. 0. ) .and. (omegaC .eq. 0. ) ).and.
     $     strrateflg .gt. 0.0 .and. step1flag .ne. 1) then
         tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/rateFactor
         if (unloadingFlag .eq. 0) then
            tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $           rateFactor
         else
            tempEquivStrainC=epsilonC+oldAlpha*(minEquivStrain-
     $           epsilon)/rateFactor+(tempEquivStrain-minEquivStrain)*
     $           alpha/rateFactor 
         end if
      else
         rateFactor = tempRateFactor 
      endif
      
      fsT = (tempEquivStrainT - kappaDT)/e0
      fsC = (tempEquivStrainC - kappaDC)/e0

      epsilon=tempEquivStrain
      epsilonT=tempEquivStrainT
      epsilonC=tempEquivStrainC
c -------------------- Compute Ductility Measure Damage --------------------------------------------------
      call cdpm2u_computeTrialCoordinates(stress, sigElastic,
     $rhoElastic,thetaElastic)
      Rs = 0.
      alphaZero= 1./sqrt(6.) 
      if ( sigElastic .lt. 0. ) then
         if ( rhoElastic .gt. 1.e-16 ) then
            Rs = -sigElastic /(alphaZero*rhoElastic)
         else 
            Rs = -sigElastic * 1.e16 / alphaZero
         endif
      else 
         Rs = 0.
      endif
      ducMeas = 1. + ( as - 1. ) * Rs ** bs
      
c----- Check which damage surfaces (tensile/compressive) are active --------------------------------------
      if (fsT .lt. -yieldTolDamage .and. 
     $     fsC .lt. -yieldTolDamage) then
c     no increase of the damage variables required
      else if (  fsT .ge. -yieldTolDamage .and. 
     $        fsC .lt. -yieldTolDamage) then
c     only tensile damage surface active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT, deltaPlasticStrainNormT,kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / ducMeas / rateFactor
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / 
     $        ducMeas
         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT,rateFactor)
      else if (  fsT .lt.  -yieldTolDamage .and. 
     $        fsC .ge.  -yieldTolDamage) then
c     only compressive damage surface active
         deltaPlasticStrainNormCNew =
     $        cdpm2u_computeDeltaPlasticStrainNormC(alpha, 
     $        tempEquivStrainC,
     $        deltaPlasticStrainNormC, kappaDC,rhoElastic,
     $        tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas*rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC

         omegaC = 
     $        cdpm2u_computeDamageC(kappaDC, 
     $        kappaDC1, kappaDC2, omegaOldC,rateFactor)
      else if (  fsT .ge.-yieldTolDamage  .and. 
     $        fsC .ge.-yieldTolDamage ) then
c Both compressive and tensile damage surfaces are active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT,deltaPlasticStrainNormT, kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / (ducMeas * rateFactor)
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / ducMeas

         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT,rateFactor)
c     only compressive damage surface active
         deltaPlasticStrainNormCNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormC(
     $        alpha,tempEquivStrainC,deltaPlasticStrainNormC, 
     $        kappaDC,rhoElastic,tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas * rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC
         omegaC = cdpm2u_computeDamageC(kappaDC, kappaDC1, 
     $        kappaDC2, omegaOldC,rateFactor)
      endif
    
      return
      end

      subroutine cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,equivStrainNew,equivStrainOld,
     $     ym,pr,gtol)
c     Function to check if there is unloading and reloading occuring within one step (usually happens during cyclic loading). If such a process is happening the algorithm returns a flag and an approximate value of the minimum equivalent strain. Moreover it returns the current equivalent strain.

      real stress(6),minEquivStrain,stress1(6),ym,pr,
     $     deltaStress(6),stressPlus(6),stressMinus(6),
     $     equivStrainOld,equivStrainNew,equivStrain1,stressOld(6),
     $     equivStrainPlus,equivStrainMinus,sigV,rho,theta,
     $     cdpm2u_computeEquivalentStrain,gtol
      integer i,j,unloadingFlag
c     stress(6)   ------------- effective stress tensor(no damage included) <Input>
c     stressOld(6) ---------- effective stress vector in  previous step <Input>
c     equivStrainOld ---------- equivalent strain of the previous loading step <input>
c     ym             ---------- Young's modulus <Input>
c     pr             ---------- Poisson's ratio <Input>
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression) <Output>
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     equivStrainNew ---------- equivalent strain of the current loading step  <Output>
c     minEquivStrain   -------- minimum equivalent strain <Output>
c     stressPlus    ----------   stress vector representing previous step's elastic strain vector plus 0.99*deltaStrain
c     stressMinus   ----------   stress vector representing current step's elastic strain vector minus 0.01*deltaStrain
c     sigV          ---------- volumentric stress
c     rho           ---------- deviatoric stress
c     theta         ---------- Lode angle
c     equivStrain1  ---------- equivalent strain of various strain vectors
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al

      call cdpm2u_computeTrialCoordinates(stress, sigV, rho,theta)
      equivStrainNew= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
    
      do i=1,6
        deltaStress(i)=stress(i)-stressOld(i)
      enddo
      do i=1,6
         stressPlus(i)=stressOld(i)+0.01*deltaStress(i)
         stressMinus(i)=stressOld(i)+0.99*deltaStress(i)
      enddo

      call cdpm2u_computeTrialCoordinates(stressPlus, sigV, rho,theta)
      equivStrainPlus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      call cdpm2u_computeTrialCoordinates(stressMinus, sigV, rho,theta)
      equivStrainMinus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      unloadingFlag=0
      minEquivStrain=equivStrainOld
      if ( (equivStrainPlus .lt. equivStrainOld .and.
     $     equivStrainMinus .lt. equivStrainNew) .and. 
     $     (abs(equivStrainPlus - equivStrainOld) .gt. gtol/10. .and. 
     $     abs(equivStrainMinus - equivStrainNew) .gt. gtol/10.)) then
         unloadingFlag=1
c        Unloading and reloading occurs within a single step. Subincrementation performed' 
         do i=1,100
            do j=1,6
               stress1(j)=stressOld(j)+deltaStress(j)*FLOAT(i)/100.0 
            enddo
            call cdpm2u_computeTrialCoordinates(stress1,sigV,rho,theta)
            equivStrain1= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
            if( equivStrain1 .le. minEquivStrain) then
               minEquivStrain=equivStrain1
            else
               goto 625
            end if
         enddo
      end if
      
 625  continue
      return
      end
      
      real function  cdpm2u_computeEquivalentStrain(sigElastic,
     $     rhoElastic, thetaElastic)
c     Function to calculate equivalent strain used in the damage part according to eq. 37 of the IJSS paper by Grassl et al.

      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real  answer, rFunction,thetaElastic,rhoElastic,
     $     sigElastic, pHelp,qHelp,help,e0
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     pHelp,help,qhelp -------------- variables used in the calculation of the equivalent strain
c     rFunction       --------------- parameter given by eq.(19) of the IJSS paper by Grassl et al.      
c     e0              --------------- parameter equal to ft/E
      e0=ft/ym
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(thetaElastic)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(thetaElastic) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     $     cos(thetaElastic)**2.+ 5. * ecc**2. - 4. * ecc) )
      pHelp = -m0 * ( rhoElastic * rFunction / ( sqrt(6.) * fc ) + 
     $     sigElastic / fc )
      qHelp = -3. / 2. * (rhoElastic** 2.) / (fc**2.)
      help = -0.5 * pHelp + sqrt((pHelp** 2.) / 4. - qHelp)
c     negative help values are not of interest and create problems since we compute the square root 
      answer = 0.
      if ( help .gt. 0. ) then
         answer= help * e0
      endif
      cdpm2u_computeEquivalentStrain=answer
      return 
      end
      
      real function  cdpm2u_computeDamageT(kappa,kappaOne,kappaTwo,le, 
     $      omegaOld,rateFactor)
c     Function to calculate damage due to tension.

      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real kappa, kappaOne, kappaTwo, omegaOld,residual,le,
     $     residualDerivative,omega,tol,help,e0,yieldTolDamage,
     $     rateFactor,wf1Mod,wfMod,residualStrength
      integer iter,newtonIter
c     kappa           ----------------------- damage history parameter kappaDT <input>
c     kappaOne        ----------------------- damage history parameter kappaDT1 <input>
c     kappaTwo        ----------------------- damage history parameter kappaDT2 <input
c     omegaOld        ----------------------- old damage variable (previous step) omegaT <input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     residualDerivative--------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaT
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaT for exponential softening law
c     help            -----------------------  help variable used in the billinear softening law
c     iter            ----------------------- counter of performed N-R iterations in the exponential softening law
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaT in the exponential softening law
c     yieldTolDamage  ----------------------- tolerance used in the damage algorithm if Hp=0

      newtonIter=100
      tol=gTol/100.0
      e0=ft/ym
      yieldTolDamage=gTol*10.0
      residualStrength=1.e-3*ft
      wfMod = wf
      wf1Mod = wf1
      ftTemp = ft*(1.-yieldTolDamage)
      if(strrateflg .eq. 2) then
         wfMod = wf/rateFactor
         wf1Mod = wf1/rateFactor
      else if(strrateflg .eq. 3) then
         wfMod = wf/(rateFactor*rateFactor)
         wf1Mod = wf1/(rateFactor*rateFactor)
      end if


      if ( kappa  .gt. e0*(1-yieldTolDamage) ) then
        if ( type .eq. 0. ) then
c     Linear damage law
           omega = ( ym * kappa * wfMod - ftTemp * wfMod +
     $      ftTemp * kappaOne * le ) /
     $          ( ym * kappa * wfMod - ftTemp * le * kappaTwo )
           help = le * kappaOne + le * omega * kappaTwo
           if ( help .ge. 0. .and. help .lt. wfMod .and. 
     $          (1-omega)*ym*kappa .gt. residualStrength) then
              goto 185
           endif
c           write(*,*) 'Residual stress computed\n'
           omega = 1.-1.e-3*ft/(ym*kappa)            
        else if ( type .eq. 1. ) then
c     Bilinear damage law
           omega = ( ym * kappa * wf1Mod - ftTemp * wf1Mod - ( ft1 -
     $      ftTemp ) * 
     $          kappaOne * le ) /( ym * kappa * wf1Mod + 
     $          ( ft1 - ftTemp ) * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo
            if ( help .ge. 0. .and. help .lt. wf1Mod .and. 
     $           (1-omega)*ym*kappa .gt. residualStrength) then
               goto 185
            endif
            
            omega = ( ym * kappa * ( wfMod - wf1Mod ) -
     $       ft1 * ( wfMod - wf1Mod ) +
     $           ft1 * kappaOne * le  - ft1 * wf1Mod ) / ( ym * kappa * 
     $           ( wfMod - wf1Mod )  - ft1 * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo

            if ( help .gt. wf1Mod .and. help .lt. wfMod  .and. 
     $           (1-omega)*ym*kappa .gt. residualStrength) then
               goto 185
            endif
c            write(*,*) 'Residual stress computed\n'
            omega = 1.-1.e-3*ftTemp/(ym*kappa)
            
         else if ( type .eq. 2. ) then
c     Exponential: Iterative solution with N-R procedure
            omega = 1.
            residual=0.
            residualDerivative = 0.
            iter = 0

 135        continue
            iter=iter+1
            residual = ( 1 - omega ) * ym * kappa - ftTemp *
     $           exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod)
            residualDerivative = -ym * kappa + ftTemp * le * 
     $           kappaTwo / wfMod * exp(-le * ( omega * kappaTwo + 
     $           kappaOne ) / wfMod)
            omega =omega- residual / residualDerivative;
            if ( iter .gt. newtonIter ) then
               write(*,*) '*** Algorithm for tensile damage-
     $No convergence reached after 100 iterations ***'
               stop
            end if
            if (abs(residual/ftTemp) .ge. 1.e-8) then
               goto 135
            end if
         end if
      else 
         omega = 0.;
      endif
      
       
       if ( omega .gt. 1. ) then
          omega = 1.
       endif
       
       if ( omega .lt. 0. .or. omega .lt. omegaOld) then
          omega=omegaOld
       endif
 185   cdpm2u_computeDamageT= omega
       return 
       end  

      real function cdpm2u_computeDamageC(kappa, kappaOne, 
     $     kappaTwo, omegaOld, rateFactor)
c     Function to calculate damage due to compression
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real kappa, kappaOne, kappaTwo, omegaOld,residual,dResidualDOmega,
     $     exponent,omega,tol,kappaDC,e0,yieldTolDamage,rateFactor,
     $     efcMod
      integer nite,newtonIter
c     kappa           ----------------------- damage history parameter kappaDC <Input>
c     kappaOne        ----------------------- damage history parameter kappaDC1 <Input>
c     kappaTwo        ----------------------- damage history parameter kappaDC2 <Input
c     omegaOld        ----------------------- old damage variable (previous step) omegaC <Input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     dResidualDOmega ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     exponent        ----------------------- exponent used in the formulation of the damage law (in the current version assumed =1.)
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaC
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaC
c     nite            ----------------------- counter of performed N-R iterations
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaC
      omega=1.
      nite = 0
      residual = 0.
      dResidualDOmega = 0.
      exponent = 1.
      e0=ft/ym
      
      newtonIter=100
      tol=gTol/2.0
      yieldTolDamage=gTol*10.
      
      ftTemp = ft*(1-yieldTolDamage)
      efcMod = efc
      if(strrateflg .eq. 2) then
         efcMod = efc/rateFactor
      else if(strrateflg .eq. 3) then
         efcMod = efc/(rateFactor*rateFactor)
      end if

      
      if (damageflag .eq. 1.0) then
         omega=0.0
         kappaOne=0.0
         kappaTwo=0.0
         kappa=0.0
         goto 188
      end if
      

      if ( kappa .gt. e0*(1-yieldTolDamage)) then
 187     continue 
         nite=nite+1;
         residual =  ( 1. - omega ) * ym * kappa - ftTemp * exp( - ( 
     $        kappaOne + omega * kappaTwo ) / efcMod )
         dResidualDOmega =-ym * kappa + ftTemp *
     $   kappaTwo / efcMod * exp( -( 
     $        kappaOne + omega * kappaTwo ) / efcMod )
         omega = omega- residual / dResidualDOmega
         if ( nite .gt. newtonIter ) then
            write(*,*) '*** Algorithm for compressive damage-
     $No convergence reached.'
c Set omega=omegaOld =', omegaOld,' ***'
c            omega = omegaOld
            stop
         endif
         if( abs(residual/ftTemp) .ge. tol )   goto 187
         
       else 
          omega = 0.
       endif

       if ( omega .gt. 1 ) then
          omega = 1.
       endif
       
       if ( omega .lt. 0. .or. omega .lt. omegaOld ) then
          omega = omegaOld
       endif
       
 188   cdpm2u_computeDamageC=omega
       return
       end

      
      real function  cdpm2u_computeDeltaPlasticStrainNormT(tempKappaD,
     $     plastStrNorm, kappaD,e0,yieldTolDamage)      
c     Function returning the norm of the increment of the plastic strain tensor.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch

      real e0,yieldTolDamage,answer,tempKappaD,plastStrNorm, kappaD,
     $     factor
c     tempKappaD      -------------------    temporary (current) KappaDt                   <Input>
c     kappaD          -------------------    old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------    plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------    variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------    tolerance used when Hp=0                      <Input>
c     factor          -------------------    factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     answer          -------------------    calculated norm             
      factor = 0.
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0  * ( 1. - yieldTolDamage )) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer = plastStrNorm*factor
      else 
         answer=plastStrNorm 
      end if
      cdpm2u_computeDeltaPlasticStrainNormT=answer
      return 
      end

      real function  cdpm2u_computeDeltaPlasticStrainNormC( 
     $     alpha,tempKappaD,
     $     plastStrNorm, kappaD,rho,tempKappa,yieldTolDamage)
c     Function returning the norm of the increment of the plastic strain tensor multiplied by alphaC and betaC according to eq. (48) of the IJSS paper by Grassl et al.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch      
      real ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,eratetype,
     1     type,bs,wf,wf1,efc,ft1,sratetype,failflg,m0,gTol,
     1     damageflag,printflag
c$omp threadprivate (/cdpmc/)
      real   tempKappaD,plastStrNorm, kappaD,factor,qh2,alpha,rho,e0,
     $     cdpm2u_qh2fun,extraFactor,tempKappa,kappa,
     $     yieldTolDamage,answer

c     tempKappaD      -------------------  temporary (current) KappaDt                   <Input>
c     kappaD          -------------------  old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------  plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------  variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------  tolerance used when Hp=0                      <Input>
c     rho             -------------------  deviatoric stress                             <Input>
c     alpha           -------------------  variable calculated based on eq. 46 of the IJSS paper by P. Grassl et al.    <Input>
c     tempKappa       -------------------  temporary (current step) cummulative plastic strain (kappaP)   <Input>
c     factor          -------------------  factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     qh2             -------------------  variables containing the results of the hardening function 
c     cdpm2u_qh2fun          -------------------  functions to calculate the hardening functions  given in eq. (31) in IJSS paper by Grassl et al.
c     extraFactor     -------------------  variable betaC given in eq. (50) of the IJSS paper by Grassl et al.
c     answer          -------------------  calculated answer         
      e0=ft/ym
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer=plastStrNorm*factor
      else 
         answer=plastStrNorm
      end if
      
      qh2=cdpm2u_qh2fun(tempKappa,hp)
      
      if (rho<1.e-16) then
         extraFactor =ft * qh2 * sqrt(2. / 3.) / 1.e-16 / sqrt( 1. + 
     $        2.*  (df** 2.) )
      else 
         extraFactor =ft * qh2 * sqrt(2. / 3.) / rho / sqrt( 1. + 
     $        2.* (df** 2.) )
      endif
      answer=answer*extraFactor*alpha
      cdpm2u_computeDeltaPlasticStrainNormC=answer
      return
      end
      

c ---------------------------------------------General functions--------------------------------

      subroutine cdpm2u_computePrincValues(ax,w,flag,z)
c     Subroutine to solve the eigenvalues and eigenvectors of real symmetric matrix by jacobi method.
c     Eigenvalues are stored in descending order (from max to min) in array w and eigenvectors
c     are stored columnwise in matrix z.

      real ax,w,z,matrix
      integer flag
      dimension ax(6),w(3), z(3,3),matrix(3,3)
      real ssum, aa, co, si, tt, tol, sum, aij, aji;
      integer ite, i, j, k, ih
      
      real swap
      integer ii,jj,kk
c     ax(6)    ---------------------            Input tensor. It is a 6x1 array and components are given 
c                                               in order xx,yy,zz,xy,yz,xz  <Input>
c     flag     ---------------------            Flag denoting whether it is a stress or a strain tensor in 
c                                                order to convert voigt shear strains(gammas) to tensorial
c                                                shear strains(epsilons).
c                                                = 0: Stress tensor
c                                                = 1: Strain tensor
c     w(3)     --------------------             3x1 Array containing eigenvectors of ax stored in 
c                                               descending order(max to min)  <Output>
c     z(3,3)   --------------------             Matrix containing all eigenvectors stored columnwise <Output>
c     matrix(3,3) -----------------             Stress/Strain tensor that is being processed
c     tol      --------------------             Tolerance of the algorithm (default 10 figures,1e-10)
c     ite,i,j,k,ih ----------------             Counters used in various loops
c     sum,ssum --------------------             Summation variables used in various addition procedures
c     aa,co,si,tt,aij,aji ---------             Variables used in the algorithm 

      tol=1.e-10
c     Reconstruct tensor based on given array
      matrix(1,1)=ax(1)
      matrix(2,2)=ax(2)
      matrix(3,3)=ax(3)
      if (flag .eq. 0) then
         matrix(2,1)=ax(4)
         matrix(1,2)=ax(4)
         matrix(3,1)=ax(6)
         matrix(1,3)=ax(6)
         matrix(3,2)=ax(5)
         matrix(2,3)=ax(5)
      else
         matrix(2,1)=ax(4)/2.
         matrix(1,2)=ax(4)/2.
         matrix(3,1)=ax(6)/2.
         matrix(1,3)=ax(6)/2.
         matrix(3,2)=ax(5)/2.
         matrix(2,3)=ax(5)/2.
      endif
      
c     Initialise w,z and check if zero stress state
      do i=1,3
         w(i) = matrix(i, i)
      enddo
      sum=0.
      do i=1,3
         do j=1,3
            sum =sum+ abs( matrix(i, j) )
            z(i, j) = 0.0
         enddo
         z(i, i) = 1.0
      enddo
      if ( sum .le. 0.0 ) then
         goto 900
      endif
            
c     Reduce to matrix diagonal
      ite=0

 272  continue
      ssum = 0.0
      do j=2,3
         ih = j - 1
         do i=1,ih
            if ( abs( matrix(i, j) ) / sum  .gt. tol ) then
               ssum =ssum+ abs( matrix(i, j) )
c     CALCULATE ROTATION ANGLE
               aa = atan2( matrix(i, j) * 2.0, w(i) - w(j) ) /  2.0
               si = sin(aa)
               co = cos(aa)
               
c     MODIFY "I" AND "J" COLUMNS OF "matrix" and "z"
               do k=1,i-1
                  tt = matrix(k, i)
                  matrix(k, i) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (i,i)
               tt = w(i)
               w(i) = co * tt + si *matrix(i, j)
               aij = -si * tt + co *matrix(i, j)
               tt = z(i, i)
               z(i, i) = co * tt + si *z(i, j)
               z(i, j) = -si * tt + co *z(i, j)
               
               do k=i+1,j-1
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (j,j)
               tt = matrix(i, j)
               aji = co * tt + si *w(j)
               w(j) = -si * tt + co *w(j)
               
               tt = z(j, i)
               z(j, i) = co * tt + si *z(j, j)
               z(j, j) = -si * tt + co *z(j, j)
               
               do k=j+1,3
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(j, k)
                  matrix(j, k) = -si * tt + co *matrix(j, k)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     MODIFY DIAGONAL TERMS
               w(i) = co * w(i) + si * aji
               w(j) = -si * aij + co *w(j)
               matrix(i, j) = 0.0
            else 
c     matrix(I,J) MADE ZERO BY ROTATION
            endif
         enddo
      enddo
      
      ite=ite+1
c     CHECK FOR CONVERGENCE
      if ( ite .gt. 50 ) then
         write(*,*)  '*** Compute principal values.
     $Too many iterations! ***'
         stop
      endif
      
      if ( abs(ssum) / sum .gt. tol  ) goto 272
      
      do ii=1,2
         do jj=1,2
            if ( w(jj + 1) > w(jj) ) then
c     swap eigenvalues and eigenvectors
               swap = w(jj + 1);
               w(jj + 1) = w(jj);
               w(jj) = swap;
               do kk=1,3
                  swap = z(kk, jj + 1);
                  z(kk, jj + 1) = z(kk, jj);
                  z(kk, jj) = swap;
               enddo
            endif
         enddo
      enddo
 900  continue
      return
      end
      
      subroutine cdpm2u_computeInverseJac(jacinv,jac,error)
c     Subroutine to calculate the inverse of a 4x4 matrix
      real jacinv(4,4),jac(4,4),tmp(4,4)
      real  piv, linkomb,dtol
      integer i,j,k,error
c     jac(4,4)    ---------  matrix whose inverse will be calculated <Input>
c     jacinv(4,4) ---------  inverse matrix  <Output>
c     error       --------- integer showing whether solution has converged <Output>
c                            = 0 solution has converged
c                            =-1 solution has not converged
c     tmp(4,4)    --------- temporary matrix
c     piv,lincomb --------- variables used in gaussian elimination
c     i,j,k       --------- counters
c     dtol        --------- tolerance of the algorithm

      dtol=1.e-20           
      do i=1,4
         do j=1,4
            tmp(i,j)=jac(i,j)
            jacinv(i,j)=0.
         enddo
         jacinv(i,i)=1.
      enddo
      
      do i=1,3
         piv = tmp(i, i)
         if (abs(piv) .lt. dtol) then
            error=-1
            goto 452
         endif

         do j=i+1,4
            linkomb = tmp(j, i) / tmp(i, i)
            do k=i,4
               tmp(j, k) = tmp(j, k) - tmp(i, k) * linkomb
            enddo
            do k=1,4
               jacinv(j, k) =jacinv(j, k)-jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=4,2,-1
         piv = tmp(i, i)
         do j=i-1,1,-1
            linkomb = tmp(j, i) / piv
            do k=i,1,-1
               tmp(j, k) =  tmp(j, k) -tmp(i, k) * linkomb
            enddo
            do k=4,1,-1
               jacinv(j, k) =jacinv(j, k) - jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=1,4
         do j=1,4
            jacinv(i, j) = jacinv(i, j) / tmp(i, i)
         enddo
      enddo 
      error=0
 452  continue
      return
      end
      
      
      
      subroutine cdpm2u_transformStressVectorTo(stress,princDir,
     $     princStress)
c     Rotates principal stress tensor back to the original Coordinate system. It calculates transformation 
c     matrix  and then multiplies with principal stress tensor.

      real stress(6),princDir(3,3),princStress(3),transpose(3,3),
     $     transformationMtrx(6,6),sum,princStressTensor(6)
      integer i,j
c     stress (6)             ----------------- stress vector <Output>
c     princDir(3,3)          ----------------- matrix containing eigenvectors stored columnwise <Input>
c     princStress(3)         ----------------- stress vector containing principal stresses =
c                                              {sigma1,sigma2,sigma3}  <Input>
c     princStressTensor(6)   ----------------- stress vector at principal axis coordinate system 
c                                              ={sigma1,sigma2,sigma3,0,0,0}
c     transpose(3,3)         ----------------- matrix containing eigenvectors stored rowise
c     sum                    ----------------- variable used in summation during matrix multiplication
c     transformationMtrx(6,6)----------------- transformation matrix used to transform principal stress
c                                              vector to stress vector in original CS
c     i and j                -----------------  integers used as counters


      do i=1,3
         princStressTensor(i)= princStress(i)
         princStressTensor(i+3)= 0.
         do j=1,3
            transpose(i,j)=princDir(j,i)
         enddo
      enddo
      
      
      transformationMtrx(1,1)=transpose(1,1)*transpose(1,1)
      transformationMtrx(1,2)=transpose(2,1)*transpose(2,1)
      transformationMtrx(1,3)=transpose(3,1)*transpose(3,1)
      transformationMtrx(1,4)=2.*transpose(1,1)*transpose(2,1)
      transformationMtrx(1,5)=2.*transpose(2,1)*transpose(3,1)
      transformationMtrx(1,6)=2.*transpose(1,1)*transpose(3,1)

      transformationMtrx(2,1)=transpose(1,2)*transpose(1,2)
      transformationMtrx(2,2)=transpose(2,2)*transpose(2,2)
      transformationMtrx(2,3)=transpose(3,2)*transpose(3,2)
      transformationMtrx(2,4)=2.*transpose(1,2)*transpose(2,2)
      transformationMtrx(2,5)=2.*transpose(2,2)*transpose(3,2)
      transformationMtrx(2,6)=2.*transpose(1,2)*transpose(3,2)

      transformationMtrx(3,1)=transpose(1,3)*transpose(1,3)
      transformationMtrx(3,2)=transpose(2,3)*transpose(2,3)
      transformationMtrx(3,3)=transpose(3,3)*transpose(3,3)
      transformationMtrx(3,4)=2.*transpose(1,3)*transpose(2,3)
      transformationMtrx(3,5)=2.*transpose(2,3)*transpose(3,3)
      transformationMtrx(3,6)=2.*transpose(1,3)*transpose(3,3)

      transformationMtrx(4,1)=transpose(1,1)*transpose(1,2)
      transformationMtrx(4,2)=transpose(2,1)*transpose(2,2)
      transformationMtrx(4,3)=transpose(3,1)*transpose(3,2)
      transformationMtrx(4,4)=transpose(1,1)*transpose(2,2)+
     $     transpose(2,1)*transpose(1,2)
      transformationMtrx(4,5)=transpose(2,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(2,2)
      transformationMtrx(4,6)=transpose(1,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(1,2)

      transformationMtrx(5,1)=transpose(1,2)*transpose(1,3)
      transformationMtrx(5,2)=transpose(2,2)*transpose(2,3)
      transformationMtrx(5,3)=transpose(3,2)*transpose(3,3)
      transformationMtrx(5,4)=transpose(1,2)*transpose(2,3)+
     $     transpose(2,2)*transpose(1,3)
      transformationMtrx(5,5)=transpose(2,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(2,3)
      transformationMtrx(5,6)=transpose(1,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(1,3)
    

      transformationMtrx(6,1)=transpose(1,1)*transpose(1,3)
      transformationMtrx(6,2)=transpose(2,1)*transpose(2,3)
      transformationMtrx(6,3)=transpose(3,1)*transpose(3,3)
      transformationMtrx(6,4)=transpose(1,1)*transpose(2,3)+
     $     transpose(2,1)*transpose(1,3)
      transformationMtrx(6,5)=transpose(2,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(2,3)
      transformationMtrx(6,6)=transpose(1,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(1,3)

 
      do i=1,6
         sum=0.
         do j=1,6
            sum=sum+  transformationMtrx(i,j)*princStressTensor(j)
         enddo
         stress(i)=sum
      enddo
c
      return
      end


      subroutine utan50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,failels,cma,qmat,nlq,lq1)
c
c*****************************************************************
c Implementation of the stiffness for CDPM2 in UMAT50V in LS-DYNA
c Support page: http://petergrassl.com/Research/DamagePlasticity/CDPMLSDYNA/index.html
c
C Only elastic stiffness so far.
c
c      Authors of this subroutine: Peter Grassl
c
c Last updated: 20 April 2016
c
c*****************************************************************

      implicit double precision(a-h,o-z)
c      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      logical failels(*)
      character*5 etype
      
      real es(6,6)
      real ym,pr
c
c Define elastic stiffness to start with. Needs to be improved later.
c
c
      do k=1,6
         do j=1,6
            es(j,k) =0.
         enddo
      enddo

      ym=cm(1)
      pr=cm(2)      

      const = ym/((1.+pr)*(1.-2.*pr))
      es(1,1)= const*(1.-pr)
      es(2,2)= es(1,1)
      es(3,3) = es(1,1)
      es(4,4) = const*(1.-2.*pr)/2.
      es(5,5) = es(4,4)
      es(6,6) = es(4,4)
      es(1,2) = const*pr
      es(1,3) = es(1,2)
      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(2,3) = es(1,2)
      es(3,2) = es(2,3)

      do i=lft,llt
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo

      return

      end 
