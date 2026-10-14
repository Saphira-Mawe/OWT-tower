function[Hp]=transferpileAF2(Db,Lf,Lp, L,thp, M,om)
global  E rok m g  
%=================================
 Hp=eye(4);
 cm=2; row=1000;
Ik=pi*(Db^4-(Db-2*thp)^4)/64;
Ak=pi*(Db^2-(Db-2*thp)^2)/4;
 Ma=0*cm*row*Ak;
pk=@(z)m*g+g*Ak*rok*(L-z);
l=linspace(Lf,Lp,M);
for i=1:M-1
        lu=(l(i+1)+l(i))/2;
    p=pk(lu);
%===================================================
s11a=-p/(2*E*Ik)+(p.^2./(2*E*Ik)^2+(rok*Ak+Ma)*om.^2./(E*Ik)).^0.5;
s11=(s11a).^0.5;
s12=s11;
s21a=-p/(2*E*Ik)-(p.^2./(2*E*Ik)^2+(rok*Ak+Ma)*om.^2./(E*Ik)).^0.5;
s21=(abs(s21a)).^0.5;
s22=s21;
%======
z=l(i);
%================================
t11=cosh(s11.*z);    t21=s11.*sinh(s11.*z);  t31=(p+E*Ik*s11.^2).*cosh(s11.*z);  t41=(p+E*Ik*s11.^2).*s11.*sinh(s11*z);
t12=sinh(s12.*z);    t22=s12.*cosh(s12.*z);  t32=(p+E*Ik*s12.^2).*sinh(s12.*z); t42=(p+E*Ik*s12.^2).*s12.*cosh(s12.*z);
t13=cos(s21.*z);     t23=-s21.*sin(s21.*z);  t33=(p-E*Ik*s21.^2).*cos(s21.*z);    t43=-(p-E*Ik*s21.^2).*s21.*sin(s21.*z);
t14=sin(s22.*z);    t24=s22.*cos(s22.*z);    t34=(p-E*Ik*s22.^2).*sin(s22.*z);    t44=(p-E*Ik*s22.^2).*s22.*cos(s22.*z);
 t15=0;              t25=0;              t35=0;              t45=0;
 

Tpi=[t11 t12 t13 t14 ;
        t21 t22 t23 t24 ;
        t31 t32 t33 t34; 
        t41 t42 t43 t44];
%===============================
z=l(i+1);
%==================
t11=cosh(s11.*z);    t21=s11.*sinh(s11.*z);  t31=(p+E*Ik*s11.^2).*cosh(s11.*z);  t41=(p+E*Ik*s11.^2).*s11.*sinh(s11*z);
t12=sinh(s12.*z);    t22=s12.*cosh(s12.*z);  t32=(p+E*Ik*s12.^2).*sinh(s12.*z); t42=(p+E*Ik*s12.^2).*s12.*cosh(s12.*z);
t13=cos(s21.*z);     t23=-s21.*sin(s21.*z);  t33=(p-E*Ik*s21.^2).*cos(s21.*z);    t43=-(p-E*Ik*s21.^2).*s21.*sin(s21.*z);
t14=sin(s22.*z);    t24=s22.*cos(s22.*z);    t34=(p-E*Ik*s22.^2).*sin(s22.*z);    t44=(p-E*Ik*s22.^2).*s22.*cos(s22.*z);
 t15=0;              t25=0;              t35=0;              t45=0;
 %====
Tpi1=[t11 t12 t13 t14 ;
        t21 t22 t23 t24 ;
        t31 t32 t33 t34; 
        t41 t42 t43 t44];
  Th=Tpi1*inv(Tpi) ;
  Hp=Th*Hp;
  % multiplication elemement per element 
end
end


