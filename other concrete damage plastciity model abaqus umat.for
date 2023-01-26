c     see Lee,  J.,  and G. L. Fenves, “Plastic-Damage Model 
c     for Cyclic Loading of Concrete Structures,” Journal of 
c     Engineering Mechanics, vol. 124, no. 8, pp. 892–900, 1998. 
      SUBROUTINE UMAT(sig,statev,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 stran,dstran,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION sig(NTENS),statev(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 stran(NTENS),dstran(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     5 JSTEP(4),cc(ntens,ntens),ps(3),an(3,3),dd1(3,6),dd2(3,6),
     $ sig_pr(3),sde_pr(3),dj2_ds2(3,3),cin(3,3),uni(3),dx(6),
     $ x(6),ress(6),pjac(6,6),dm(3),dm_ds(3,3),dh(2),dh_ds(2,3),
     $drr_ds(3),pjac2(6,6),dx6(6),dx66(6),dx5(5),px55(5,5),
     $px33(3,3),px23(2,3),sig_tr(6),sig_pr1(3),e_tr(3),e_tr1(3),
     $df_ds(3),df_dk(2),ress1(6),dx55(5),yx33(3,3),dir63(6,3)
      PARAMETER(zr=0.0D0,on=1.0D0,tw=2.0D0,thr=3.0D0,fr=4.0D0,
     %tol=1.0d-6,tol2=10.0d0)
      pi=atan(on)*fr;
      Em=PROPS(1);emu=PROPS(2);fco=props(3);fto=props(4);
      fb_fc=props(5);ecc=props(6);omega=props(7);pkc=props(8)
      rec_c=props(11);rec_t=props(12);omega=omega*pi/180.0;
      tan_o=tan(omega);alpha=(fb_fc-on)/(fb_fc*tw-on);
      gamma=thr*(on-pkc)/(tw*pkc-on);
      ept=statev(1);epc=statev(2);dd=statev(3);
      p1=Em/(on+emu)/(on-tw*emu);cc=zr;
      cc(1:3,1:3)=p1*emu;
      cc(1,1)=(on-emu)*p1;cc(2,2)=(on-emu)*p1;
      cc(3,3)=(on-emu)*p1;Gm=Em/(on+emu)/tw;
      cc(4,4)=Gm;cc(5,5)=Gm;cc(6,6)=Gm;ddsdde=cc*(on-dd);
      sig_tr=sig/(on-dd)+matmul(cc,dstran);
      
      
      sig_tr=zr;sig_tr(1)=-12.5;ept=zr;epc=zr!0.002-fco/Em+0.0001!!!!
      write(*,*) 'epc1',epc
      call sprind(sig_tr,ps,an,1,ndi,nshr)
      call korder(ps,an,dd1,dd2,ndi)
      pi1=ps(1)+ps(2)+ps(3);sig_pr=ps;sde_pr=sig_pr-pi1/thr
      pj2=(sde_pr(1)**2+sde_pr(2)**2+sde_pr(3)**2)/tw
      call kaxialy(fc_eff,ft_eff,fc,ft,hc,ht,hc_eff,ht_eff,
     $ ddt_dept,ddc_depc,epc,ept,dc,dt,fto,fco,Em,tol,tol2)
      beta=fc_eff/ft_eff*(on-alpha)-(on+alpha)
      ff=(alpha*pi1+sqrt(thr*pj2)+beta*(sig_pr(1)+abs(sig_pr(1)))
     %/tw-gamma*(abs(sig_pr(1))-sig_pr(1))/tw)/(on-alpha)-fc_eff
      if ((ff/Em)<(tol*tw)) then
            sig=(on-dd)*sig_tr
      else
            dj2_ds2=-on;dj2_ds2(1,1)=tw;dj2_ds2(2,2)=tw;
            dj2_ds2(3,3)=tw;dj2_ds2=dj2_ds2/thr;
            cin=-emu/Em;cin(1,1)=on/Em;cin(2,2)=on/Em;
            cin(3,3)=on/Em;dx=zr;x=zr;x(1:3)=sig_pr;
            x(4)=ept;x(5)=epc;sig_pr1=sig_pr;ept1=ept;epc1=epc
            e_tr=matmul(cin,sig_pr);e_tr1=e_tr;ress=zr;
            ress(6)=ff;res_n=abs(ress(6)/Em/tol2);res_n1=on;
            iter=0;uni=zr;uni(1)=on;
            do while(res_n1>=tol/tol2/tw)
                  par=tol/tol2
                  write(*,*) 'ff121',res_n,par,ress
                  !write(*,*) 'ff121',x
                  !write(*,*) 'ff1221',sig_pr
                  iter=iter+1;pjac=zr;
                  sqrt1=sqrt((ecc*fto*tan_o)**2+thr*pj2)
                  dm=tan_o/thr+thr/tw*sde_pr/sqrt1
                  do iii=1,3
                        dm_ds(:,iii)=thr/tw*dj2_ds2(:,iii)/sqrt1
     %-thr*thr/fr*sde_pr*sde_pr(iii)/sqrt1**3
                  end do
                  pi1_n=abs(sig_pr(1))+abs(sig_pr(2))+
     %abs(sig_pr(3))
                  rr=(pi1+pi1_n)/pi1_n/tw
                  dh(1)=rr*dm(1);dh(2)=(rr-on)*dm(3)
                  par=tw*pi1_n**2;
                  do iii=1,3
                        if (sig_pr(iii)<zr) then
                              drr_ds(iii)=(pi1_n+pi1)/par
                        else
                              drr_ds(iii)=(pi1_n-pi1)/par
                        end if
                  end do
                  dh_ds(1,:)=drr_ds*dm(1)+rr*dm_ds(1,:)
                  dh_ds(2,:)=drr_ds*dm(3)+(rr-on)*dm_ds(3,:)
                  pjac(1:3,1:3)=cin+x(6)*dm_ds;pjac(1:3,6)=dm
                  pjac(4:5,1:3)=-x(6)*dh_ds;pjac(4:5,6)=-dh;
                  pjac(4,4)=on;pjac(5,5)=on;
                  if (sig_pr(1)>zr) then
                  df_ds=(alpha+thr/tw*sde_pr/sqrt(thr*pj2)
     %+beta*uni)/(on-alpha)
                  df_dk(1)=-ht_eff*fc_eff/ft_eff**2*sig_pr(1)
                  df_dk(2)=hc_eff/ft_eff*sig_pr(1)-hc_eff
                  else
                  df_ds=(alpha+thr/tw*sde_pr/sqrt(thr*pj2)
     %+gamma*uni)/(on-alpha)
                  df_dk(1)=zr;df_dk(2)=-hc_eff
                  end if      
                  pjac(6,1:3)=df_ds;pjac(6,4:5)=df_dk
                  if (res_n<tol/tol2) exit
                  pjac2=pjac;ress1=ress;
                  call gauss_2(pjac2,-on*ress1,dx,6) 
                  pjac2=pjac
                  !call kinverse(pjac2,pjac2,6)
                  !dx=-matmul(pjac2,ress)
                  !write(*,*) 'pjac2',pjac2
                  x=x+dx;sig_pr=x(1:3);
                  ress(1:3)=matmul(cin,x(1:3))-e_tr1+x(6)*dm
                  if(x(4)<ept1+tol/tol2**5) x(4)=ept1
                  if(x(5)<epc1+tol/tol2**5) x(5)=epc1
                  ress(4)=x(4)-ept1-x(6)*dh(1);ept=x(4)
                  ress(5)=x(5)-epc1-x(6)*dh(2);epc=x(5)
                  !write(*,*) 'ress',ress
                  !write(*,*) 'dx',dx
                  pi1=sig_pr(1)+sig_pr(2)+sig_pr(3)
                  sde_pr=sig_pr-pi1/thr
                  pj2=(sde_pr(1)**2+sde_pr(2)**2
     %+sde_pr(3)**2)/tw
                  call kaxialy(fc_eff,ft_eff,fc,ft,hc,ht,hc_eff,
     $ ht_eff,ddt_dept,ddc_depc,epc,ept,dc,dt,fto,fco,Em,tol,tol2)
                  beta=fc_eff/ft_eff*(on-alpha)-(on+alpha);
                  ff=(alpha*pi1+sqrt(thr*pj2)+beta*(sig_pr(1)
     %+abs(sig_pr(1)))/tw-gamma*(abs(sig_pr(1))-sig_pr(1))/tw)
     $/(on-alpha)-fc_eff
                  ress(6)=ff;res_n=sqrt(ress(1)**2+ress(2)**2+
     $ress(3)**2+ress(4)**2+ress(5)**2+(ress(6)/Em/tol2)**2)
                  res_n1=sqrt(ress(1)**2+ress(2)**2+ress(3)**2)+on
                  write(*,*) 'pjc',pjac
                  write(*,*) 'ress',ress
                  write(*,*) 'dx',dx
                  write(*,*) 'x',x
                  !write(*,*) 'dm',dm
                  !write(*,*) 'dm_ds',dm_ds
                  !write(*,*) 'dh',dh
                  !write(*,*) 'dh_ds',dh_ds
                  if (iter==200) then
                        write(*,*) 'trials alot'
                        !exit
                        call xit
                  end if
            end do
            write(*,*) 'ff11',res_n,sig_pr
            !pjac=zr;pjac(1:3,1:3)=cin;
            !pjac(1:3,6)=dm;pjac(4:5,6)=-dh;
            !pjac(6,1:3)=df_ds;pjac(6,4:5)=df_dk
            !pjac(4,4)=on;pjac(5,5)=on;
            !pjac(1:3,6)=df_ds;pjac(4:5,6)=df_dk;
            pjac2=pjac;
            write(*,*) 'pjac',pjac
            call kinverse(pjac2(1:5,1:5),px55,5)
            pjac2=pjac;
            call kinverse(pjac2,pjac2,6)
            write(*,*) 'px55',px55
            dx5=matmul(px55,pjac(1:5,6))
            dx55=matmul(transpose(px55),pjac(6,1:5))
            par=dx5(1)*pjac(6,1)+dx5(2)*pjac(6,2)+
     $dx5(3)*pjac(6,3)+dx5(4)*pjac(6,4)+dx5(5)*pjac(6,5)
            do ii=1,3
              px33(ii,:)=px55(ii,1:3)-dx5(ii)*dx55(1:3)/par
            end do
            do ii=1,2
              px23(ii,:)=px55(ii+3,1:3)-dx5(ii+3)*dx55(1:3)/par
            end do
            !px23=px23*-on
            write(*,*) 'dc',dc,dt,ddc_depc,ddt_dept
            !write(*,*) 'epc2',epc,fc_eff,fc
            write(*,*) 'px223',px23
            dd=on-(on-dc)*(on-dt)
            ddc=ddc_depc*(on-dt);ddt=ddt_dept*(on-dc)
            do i=1,3 
            do j=1,3
             yx33(i,j)=sig_pr(i)*(px23(1,j)*ddt+px23(2,j)*ddc)
            end do
            end do
            
            sig_tr=matmul(transpose(dd2),sig_pr)
            sig=(on-dd)*sig_tr;
            write(*,*) 'px33a',px33
            write(*,*) 'px33b',yx33
            px33=(on-dd)*px33-yx33
            write(*,*) 'px33',px33
            write(*,*) 'pj',pjac2(1:3,1:3)
            write(*,*) 'pj2',pjac2(4:5,1:3)
            write(*,*) 'ptotj',pjac2
            ddsdde=matmul(matmul(transpose(dd2),px33),dd2)
            iii=0;e_tr=matmul(cin,sig_pr)
            !!!!!!!!!!!!!!!!!!!!!!!!
            call kinverse(px33,px33,3)
            e_tr=matmul(px33,sig_pr)
            !!!!!!!!!!!!!!!!!!!
            do i=1,2
                  do j=i+1,3
                   iii=iii+1
                   dir63(1,iii)=2.0*an(i,1)*an(j,1)
                   dir63(2,iii)=2.0*an(i,2)*an(j,2)
                   dir63(3,iii)=2.0*an(i,3)*an(j,3)
                   dir63(4,iii)=an(i,1)*an(j,2)+an(i,2)*an(j,1)
                   dir63(5,iii)=an(i,1)*an(j,3)+an(i,3)*an(j,1)
                   dir63(6,iii)=an(i,3)*an(j,2)+an(i,2)*an(j,3)
                   if (abs(e_tr(i)-e_tr(j))<tol) then
                     rati=cc(i,i)-cc(i,j)
                   else
                     rati=(sig_pr(i)-sig_pr(j))/(e_tr(i)-e_tr(j))
                   end if 
                   do ii=1,6
                    do jj=1,6
                       ddsdde(ii,jj)=ddsdde(ii,jj)
     $+0.5*(on-dd)*rati*dir63(ii,iii)*dir63(jj,iii)
                    end do
                   end do
                  end do
            end do
      end if
      write(*,*) 'ddsdde',ddsdde
      statev(1)=ept;statev(2)=epc;statev(3)=dd;
      call xit
      return
      end
      subroutine korder(ps,an,dd1,dd2,ndi)
            real*8 ps(3),an(3,3),dd1(3,6),dd2(3,6),temp11(3),temp1
            integer ndi
            do ii=1,ndi
                  if (ps(ii)>ps(1)) then
                        temp1=ps(1);ps(1)=ps(ii);ps(ii)=temp1
                        temp11=an(1,:);an(1,:)=an(ii,:);
                        an(ii,:)=temp11
                  end if
                  if (ps(ii)<ps(ndi)) then
                        temp1=ps(ndi);ps(ndi)=ps(ii);ps(ii)=temp1
                        temp11=an(ndi,:);an(ndi,:)=an(ii,:);
                        an(ii,:)=temp11
                  end if
            end do
            do ii=1,ndi
                  dd1(ii,1)=an(ii,1)*an(ii,1);
                  dd1(ii,2)=an(ii,2)*an(ii,2);
                  dd1(ii,3)=an(ii,3)*an(ii,3);
                  dd1(ii,4)=2.*an(ii,1)*an(ii,2);
                  dd1(ii,5)=2.*an(ii,1)*an(ii,3);
                  dd1(ii,6)=2.*an(ii,2)*an(ii,3);
            end do
            dd2=dd1;dd2(:,4:6)=0.5*dd1(:,4:6)
            !ds_pr=dd1*ds     ds=transpose(dd2)*ds_pr
       return
       end
      subroutine kaxialy(fc_eff,ft_eff,fc,ft,hc,ht,hc_eff,ht_eff,
     $ ddt_dept,ddc_depc,epc,ept,dc,dt,fto,fco,Em,tol,tol2)
            real*8 fc_eff,ft_eff,fc,ft,hc,ht,hc_eff,ht_eff,
     $ ddt_dept,ddc_depc,epc,ept,dc,dt,fto,fco,Em,faci,tol,tol2
            pbb=0.002-fco/Em;pcc=0.01;faci=100.0
            if (abs(epc)<tol/tol2**5) then
                  fc_eff=0.4*fco;fc=fc_eff
                  hc_eff=Em*faci;hc=Em*faci;
                  ddc_depc=0.0;dc=0.0;
            elseif (epc<pbb) then
                  paa=epc/pbb;fc_eff=0.4+(0.6-0.01)*paa**3+
     $(0.01-0.6)*3.0*paa**2+(0.6*3.0-0.01*2.0)*paa;
                  fc_eff=fc_eff*fco;fc=fc_eff;
                  hc_eff=(0.6-0.01)*3.*paa**2+(0.01-0.6)*6.0*paa
     $+(0.6*3.0-0.01*2.0);hc_eff=hc_eff*fco/pbb;hc=hc_eff;
                  ddc_depc=0.0;dc=0.0
            elseif (epc>=pbb.and.epc<pbb*5.0) then
                  paa=epc/pbb;fc_eff=(1.0+0.01*(paa-1.0))*fco;
                  fc=(1.0-0.95/4.0*(paa-1))*fco;
                  hc_eff=0.01*fco/pbb;hc=-0.95/4.0*fco/pbb;
                  ddc_depc=-hc/fc_eff;
                  dc=1.0-fc/fc_eff;
                  !write(*,*) 'hc',hc,dc,ddc_depc,pbb,paa,fc
            else
                  paa=epc/pbb;fc_eff=(1.0+0.01*(paa-1.0))*fco;
                  fc=0.05*fco;
                  hc_eff=0.01*fco/pbb;hc=0.0;
                  ddc_depc=0.0;dc=1.0-fc/fc_eff;
            end if
            if (abs(ept)<tol/tol2**5) then
                  ft_eff=fto;ft=ft_eff;
                  ht_eff=Em*faci;ht=Em*faci;
                  ddt_dept=0.0;dt=0.0;
            elseif (ept<pcc*0.9) then
                  ft_eff=fto+0.01*fco/pbb*ept;
                  ft=fto*(pcc-ept)/pcc;
                  ht_eff=0.01*fco/pbb;ht=-1.0*fto/pcc;
                  ddt_dept=-ht/fto;dt=1.0-ft/ft_eff;
            else
                  ft_eff=fto+0.01*fco/pbb*ept;
                  ft=fto*0.1;
                  ht_eff=0.01*fco/pbb;ht=0.0;
                  ddt_dept=0.0;dt=1.0-ft/ft_eff;
            end if
      return
      end
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
                  ! b(n)   - array of the right hand coefficients b
                  ! n      - number of equations (size of matrix A)
                  ! output ...
                  ! x(n)   - solutions
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
            do i=k,n                       ! loop over rows
            s(i) = 0.0
            do j=k,n                    ! loop over elements of row i
                  s(i) = max(s(i),abs(a(i,j)))
            end do
            end do
            
            ! step 3: "pivoting 1" 
            ! find a row with the largest pivoting element
            pivot = abs(a(k,k)/s(k))
            l = k
            do j=k+1,n
            if(abs(a(j,k)/s(j)) > pivot) then
                  pivot = abs(a(j,k)/s(j))
                  l = j
            end if
            end do
            
            ! Check if the system has a sigular matrix
            if(pivot == 0.0) then
            write(*,*) ' The matrix is sigular '
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
      
      end subroutine gauss_2




            
