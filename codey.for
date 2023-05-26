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
