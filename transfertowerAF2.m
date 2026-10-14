\function[Ht]=transfertowerAF2(Db,Dt,Lp,L,tht,M,om)
global   E  rok   m  g
Ht=eye(4);
Ma=0; 
%==============================================
 q=Dt/Db;
 afa1=(q*Lp-L)/(Lp-L);
 bet1=(1-q)/(Lp-L);
 Ik1=pi*(Db^4-(Db-2*tht)^4)/64;
 Ak1=pi*(Db^2-(Db-2*tht)^2)/4;
 Ikk=@(z)Ik1*(afa1^4+4*(afa1^3*bet1)*z+6*(afa1^2*bet1^2)*z.^2+4*afa1*bet1^3*z.^3+bet1^4*z.^4); 
 Akk=@(z)Ak1*(afa1^2+2*afa1*bet1*z+bet1^2*z^2);
 pk=@(z)m*g+rok*g*Ak1*(afa1^2*(L-z)+(afa1*bet1)*(L^2-z^2)+(1/3)*(bet1^2)*(L^3-z^3));
 l=linspace(Lp,L,M);
 for i=1:M-1
    lu=(l(i)+l(i+1))/2;
   p=pk(lu);
  Ak=Akk(lu);
  Ik=Ikk(lu);
 %===================================================
s11a=-p/(2*E*Ik)+(p.^2./(2*E*Ik)^2+(rok*Ak+Ma)*om.^2./(E*Ik)).^0.5;
s11=(s11a).^0.5;
s12=s11;
s21a=-p/(2*E*Ik)-(p.^2./(2*E*Ik)^2+(rok*Ak+Ma)*om.^2./(E*Ik)).^0.5;
%s21=(s21a).^0.5;
s21=(abs(s21a)).^0.5;
s22=s21;
%======================================
z=l(i);
t11=cosh(s11.*z);    t21=s11.*sinh(s11.*z);  t31=(p+E*Ik*s11.^2).*cosh(s11.*z);  t41=(p+E*Ik*s11.^2).*s11.*sinh(s11*z);
t12=sinh(s12.*z);    t22=s12.*cosh(s12.*z);  t32=(p+E*Ik*s12.^2).*sinh(s12.*z); t42=(p+E*Ik*s12.^2).*s12.*cosh(s12.*z);
t13=cos(s21.*z);     t23=-s21.*sin(s21.*z);  t33=(p-E*Ik*s21.^2).*cos(s21.*z);    t43=-(p-E*Ik*s21.^2).*s21.*sin(s21.*z);
t14=sin(s22.*z);    t24=s22.*cos(s22.*z);    t34=(p-E*Ik*s22.^2).*sin(s22.*z);    t44=(p-E*Ik*s22.^2).*s22.*cos(s22.*z);
 t15=0;              t25=0;              t35=0;              t45=0;

%  t11
% t12
Tpi=[t11 t12 t13 t14 ;
        t21 t22 t23 t24 ;
           t31 t32 t33 t34; 
        t41 t42 t43 t44];
z=l(i+1);
t11=cosh(s11.*z);    t21=s11.*sinh(s11.*z);  t31=(p+E*Ik*s11.^2).*cosh(s11.*z);  t41=(p+E*Ik*s11.^2).*s11.*sinh(s11*z);
t12=sinh(s12.*z);    t22=s12.*cosh(s12.*z);  t32=(p+E*Ik*s12.^2).*sinh(s12.*z); t42=(p+E*Ik*s12.^2).*s12.*cosh(s12.*z);
t13=cos(s21.*z);     t23=-s21.*sin(s21.*z);  t33=(p-E*Ik*s21.^2).*cos(s21.*z);    t43=-(p-E*Ik*s21.^2).*s21.*sin(s21.*z);
t14=sin(s22.*z);    t24=s22.*cos(s22.*z);    t34=(p-E*Ik*s22.^2).*sin(s22.*z);    t44=(p-E*Ik*s22.^2).*s22.*cos(s22.*z);
 t15=0;              t25=0;              t35=0;              t45=0;

Tpi1=[t11 t12 t13 t14 ;
        t21 t22 t23 t24 ;
        t31 t32 t33 t34; 
        t41 t42 t43 t44];
 Th=Tpi1*inv(Tpi);
  Ht=Th*Ht;  % multiplication elemement per element 
end
end



