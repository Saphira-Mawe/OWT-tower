 %==== ce programme calcule les  deplacement, slope moment and shear force 
%===et frequences propres dune wind turbine subjected=============
% to a  wind load by tranfer matrix method=========================%
%==============donnees de la foundatioN  DS ================================
clear all
format long
tic
%syms om
global E  rok m g K
g=9.8; 
freqrad=[];
vpo=[];

%======================================================================
%=======parameters to identify the wave force===================== ng====
% H=6; dw=20; lam=58; T=10; cd=1; zet=0; 
%======================================================================
%=============Inertia FI and Drag (FD) of a wave force=================

fondation='AF';
%=========FONDATION========================================
switch fondation
% 
case 'AF'
  a=1; b=2; epsilon=1e-6;
 %  1-2.5= 1ere freq, 8-12-2eme freq, 18-21-3e freq
%===example1=======
% rok=8500; E=2.1*1e11; m=350*1e3; J=2.54*1e7;
% Df=6.2132;  Lf=17.5; L=125.1;thf=0.05986; Mf=20; 
% Dp=6;Lp=20+Lf;  thp=0.06;            Mp=20; 
% Db=6; Dt=3.87; tht=0.035; Mt=40;
% a=1; b=2; eps_abs=1e-4;eps_step=1e-2; epsilon=1e-6;
% ================exemple 2-yibo==========
m=218.28E3; E=2.1*1e11; rok=7850; J=2.5*1e7;
 Df=5.5;  Lf=27;  thf=75e-5; Lt=68.67;
 Dp=5.5;Lp=28; thp=75E-3;     Mp=20;  Mf=20;
 Db=5.5; Dt=3.3160; tht=75E-3;  L=Lt+Lp+Lf;   Mt=40;  
om=a;
[Hf]=transferfoundationAF2(Df, Lf, L, thf, Mf, om);
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,om)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,om);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,om);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om.^2  1 0 ; m*om.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpA=HT(3,3)*HT(4,4)-HT(3,4)*HT(4,3);
%============================
om=b;
%=========================
[Hf]=transferfoundationAF2(Df, Lf, L, thf, Mf, om);
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,om)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,om);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,om);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om.^2  1 0 ; m*om.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpB=HT(3,3)*HT(4,4)-HT(3,4)*HT(4,3);
%=============================================
if vpA==0
    freq=a;
elseif (vpB==0)
        freq=b;
elseif ( vpA*vpB > 0 )
        error( 'f(a) and f(b) do not have opposite signs' );
  return
else
err=abs(b-a);
% fprintf('\n x0 \t\t x1 \t\t f(x0) \t\t f(x1)');
% fprintf('\n %0.4f \t %0.4f \t %0.4f \t %0.4f',x0,x1,f(x0),f(x1));
while err > epsilon
    c=(a*vpB-b*vpA)/(vpB-vpA);
    a=b;        
    b=c;
 [Hf]=transferfoundationAF2(Df, Lf, L, thf, Mf, a);
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,a)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,a);
% %=============================TOWER======================jm===
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,a);
%=============================================================
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*a.^2  1 0 ; m*a.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpA=HT(3,3)*HT(4,4)-HT(3,4)*HT(4,3);
%===================
[Hf]=transferfoundationAF2(Df, Lf, L, thf, Mf, b);
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,b)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,b);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,b);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*b.^2  1 0 ; m*b.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpB=HT(3,3)*HT(4,4)-HT(3,4)*HT(4,3);
 err=abs(b-a);
    freq=c; 
    fprintf('\n %0.4f \t %0.4f \t %0.4f \t %0.4f',a,b,vpA,vpB);
  end
end
 freqhertzAF=freq/6.28
 %=====================================
case  'CS'
    a=1; b=2;  epsilon=1e-6;
% % % %  %  1-2.5= 1ere freq, 8-12-2eme freq, 18-21-3e freq
 KL=2.58*1e9; KLR=-2.26*1e10;  KR=2.64*1e11;
  rok=8500; E=2.1*1e11; m=350*1e3; J=2.54*1e7;
  Lf=0; L=107.6;
Dp=6;Lp=20+Lf;  thp=0.06;            Mp=20; 
Db=6; Dt=3.87; tht=0.035;            Mt=40;
% %======================yibo===================
 %kum=-3.08E9; kof=kum; kuf=2.06e8; kom=8.33e10;
 % KL=2.057*1e8; KR=8.3296*1e10; KLR=-3.0769*1e9;  %CS% HOMOGENOUS SOI
 % m=218.28E3; E=2.1*1e11; rok=7850; J=2.5*1e7;
 % Lf=1e-4; L=106.67; 
 % Dp=5.5;Lp=28; thp=75E-3;     Mp=20; % valeur max 70 pour convergence       
 % Db=5.5; Dt=3.3160; tht=75E-3;     Mt=40;  % maleur max a 55 pour convergence
% %====================================================
om=a;
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,om)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,om);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,om);

%=====MATRICE GLOBALE= CS doesnt have foundation ==============================
     HG=Ht*Hp;
% %=========================== boundary at tower top ========================================
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om^2  1 0 ; m*om^2 0 0 1 ];%; 0 0 0 0 1];
  %=====BOTTOM BOUNDARY==============================
  % Up1=[1 0 0 0 ; 0 1 0 0 ;kof kom 1 0 ; -kuf -kum 0 1 ]; %;; 0 0 0 0 1];
    Up1=[1 0 0 0 ; 0 1 0 0; KLR KR 1 0 ; -KL -KLR 0 1 ];
  %===total transfer matrix================
 HT=Upn*(HG*Up1);

 vpA=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
 %==========================================
om=b;
%======================================
%===================MONOPILE==============================================
%function[Hp]=transferpileAF(Db,Lf,Lp, L,thp, M,om)
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,om);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,om);

%=====MATRICE GLOBALE= CS doesnt have foundation ==============================
     HG=Ht*Hp;
% %=========================== boundary at tower top ========================================
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om^2  1 0 ; m*om^2 0 0 1 ];%; 0 0 0 0 1];
  %=====BOTTOM BOUNDARY==============================
  % Up1=[1 0 0 0 ; 0 1 0 0 ;kof kom 1 0 ; -kuf -kum 0 1 ]; %;; 0 0 0 0 1];
    Up1=[1 0 0 0 ; 0 1 0 0; KLR KR 1 0 ; -KL -KLR 0 1 ];
   % KL=2.057*1e8; KR=8.3296*1e10; KLR=-3.0769*1e9;  %CS% HOMOGENOUS SOIL
    % kum=-3.08E9; kof=kum; kuf=2.06e8; kom=8.33e10;
  %===total transfer matrix================
 HT=Upn*(HG*Up1);
 vpB=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
if (vpA==0)
    freq=a;
elseif (vpB==0)
        freq=b;
elseif ( vpA* vpB > 0 )
        error( 'f(a) and f(b) do not have opposite signs' );
  return
else
err=abs(b-a);
% fprintf('\n x0 \t\t x1 \t\t f(x0) \t\t f(x1)');
% fprintf('\n %0.4f \t %0.4f \t %0.4f \t %0.4f',x0,x1,f(x0),f(x1));
while err > epsilon
    c=(a*vpB-b*vpA)/(vpB-vpA);
    a=b;        
    b=c;
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,a);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,a);

%=====MATRICE GLOBALE= CS doesnt have foundation ==============================
     HG=Ht*Hp;
% %=========================== boundary at tower top ========================================
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*a^2  1 0 ; m*a^2 0 0 1 ];%; 0 0 0 0 1];
  %=====BOTTOM BOUNDARY==============================
  % Up1=[1 0 0 0 ; 0 1 0 0 ;kof kom 1 0 ; -kuf -kum 0 1 ]; %;; 0 0 0 0 1];
    Up1=[1 0 0 0 ; 0 1 0 0; KLR KR 1 0 ; -KL -KLR 0 1 ];
  %===total transfer matrix================
 HT=Upn*(HG*Up1);
vpA=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
%================================================
[Hp]=transferpileAF2( Dp, Lf, Lp, L,thp, Mp,b);
% %=============================TOWER=========================================
%function[Ht]=transfertowerAF(Db,Dt,Lp,L, tht,M,om)
[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,Mt,b);
%=====MATRICE GLOBALE= CS doesnt have foundation ==============================
     HG=Ht*Hp;
% %=========================== boundary at tower top ========================================
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*b^2  1 0 ; m*b^2 0 0 1 ];%; 0 0 0 0 1];
  %=====BOTTOM BOUNDARY==============================
  % Up1=[1 0 0 0 ; 0 1 0 0 ;kof kom 1 0 ; -kuf -kum 0 1 ]; %;; 0 0 0 0 1];
    Up1=[1 0 0 0 ; 0 1 0 0; KLR KR 1 0 ; -KL -KLR 0 1 ];
  %===total transfer matrix================
 HT=Upn*(HG*Up1);
vpB=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
err=abs(b-a);
    freq=c; 
     fprintf('\n %0.4f \t %0.4f \t %0.4f \t %0.4f',a,b,vpA,vpB);
  end
end
    freqhertzCS=freq/6.28

    case  'DS'
%=====================================================================
%=======================================================================
 a=8; b=12; epsilon=1e-6;
 rok=8500; E=2.1*1e11; m=350*1e3; J=2.54*1e7; g=9.81;
Df=6;  Lf=36;   Lt=87.6; thf=0.06; Mf=36;
Dp=6;Lp=20+Lf;  thp=0.06;            Mp=200; 
Db=6; Dt=3.87; tht=(0.035+0.025)/2;  L=Lt+Lp; Mt=200;
% %=============================================
 %====example-yibo==================    
% m=218.28*1e3; E=210*1e9; rok=7850;J=2.5E7;
% Df=5.5; Lf=55; thf=75*1e-3;   Mf=55;
% Dp=5.5;Lp=28+Lf; thp=75*1e-3;  Mp=20;
% Db=5.5; Dt=3.3160; Lt=78.67; L=Lp+Lt; tht=75*1e-3; g=9.81;  Mt=40;
%==========================================================
om
=a;
[Hf]=transferfoundationDS2(Df, Lf,L, thf,Mf,om);
%===================MONOPILE==============================================
%function[Hp]=psmpile2( Dp, Lf,Lp, L, thp, omga, Mp,N)
[Hp]=transferpileAF2( Dp, Lf, Lp,L, thp,Mp,om);
% %=============================TOWER=========================================
%function[Ht]=psmtower2(Db,Dt,Lp, L, tht, omga, M,N)
[Ht]=transfertowerAF2(Db,Dt,Lp, L, tht, Mt, om);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om^2  1 0 ; m*om^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpA=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
%=======================
om=b;
%==========
[Hf]=transferfoundationDS2(Df, Lf,L, thf,Mf,om);
%===================MONOPILE==============================================
%function[Hp]=psmpile2( Dp, Lf,Lp, L, thp, omga, Mp,N)
[Hp]=transferpileAF2( Dp, Lf, Lp,L, thp, Mp,om);
% %=============================TOWER=========================================
%function[Ht]=psmtower2(Db,Dt,Lp, L, tht, omga, M,N)
[Ht]=transfertowerAF2(Db,Dt,Lp, L, tht,  Mt,om);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*om^2  1 0 ; m^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpB=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
if (vpA==0)
    freq=a;
elseif (vpB==0)
        freq=b;
elseif ( vpA* vpB > 0 )
   disp('Enter valid interval !!!');
    return
else
err=abs(b-a);
while err > epsilon
    c=(a*vpB-b*vpA)/(vpB-vpA);
    a=b;        
    b=c;
   [Hf]=transferfoundationDS2(Df, Lf,L, thf,Mf,a);
%===================MONOPILE==============================================
%function[Hp]=psmpile2( Dp, Lf,Lp, L, thp, omga, Mp,N) 
[Hp]=transferpileAF2( Dp, Lf, Lp,L, thp, Mp,a);
% %=============================TOWER=========================================
%function[Ht]=psmtower2(Db,Dt,Lp, L, tht, omga, M,N)
[Ht]=transfertowerAF2(Db,Dt,Lp, L, tht, Mt,a);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*a.^2  1 0 ; m*a.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpA=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
%==================================
[Hf]=transferfoundationDS2(Df, Lf,L, thf,Mf,b);
%===================MONOPILE==============================================
%function[Hp]=psmpile2( Dp, Lf,Lp, L, thp, omga, Mp,N)
[Hp]=transferpileAF2( Dp, Lf, Lp,L, thp, Mp,b);
% %=============================TOWER=========================================
%function[Ht]=psmtower2(Db,Dt,Lp, L, tht, omga, M,N)
[Ht]=transfertowerAF2(Db,Dt,Lp, L, tht, Mt,b);
%===============
%=====MATRICE GLOBALE===================================================
     HG=Ht*(Hp*Hf);
 %=========================== boundary at tower top =====================
 % Upn=[1 0 0 0 0; 0 1 0 0 0; 0 -J*om.^2  1 0 0; m*om.^2 0 0 1 0; 0 0 0 0 1];
  Upn=[1 0 0 0 ; 0 1 0 0 ; 0 -J*b.^2  1 0 ; m*b.^2 0 0 1]; 
%===total transfer matrix===============================================m
  HT=Upn*HG;
%=======boundary====condition== AF==U(Z=0)=Teta(z=0)=0 ===and M(L)=Q(L)=0
vpB=HT(3,1)*HT(4,2)-HT(3,2)*HT(4,1);
    err=abs(b-a);
    freq=c; 
   fprintf('\n %0.4f \t %0.4f \t %0.4f \t %0.4f',a,b,vpA,vpB);
end
    % fprintf('\n The root is %4.3f ',root);
    % fprintf('\n The error is %4.3f ',err);
end
     
freqhertzDS=freq/6.28

end
 
toc


