function[Hf]=transferfoundationDS2(Db,Lf,L,thf,M,om)
global E  rok   m  g  

%sympref('floatingPointOutput',false)
Hf=eye(4);   Ma=0;

%  K1=[1.168*1e9 1.165*1e9 1.130*1e9  1.100*1e9 1.060*1e9 1.020*1e9...
% 9.890*1e8  9.540*1e8 9.190*1e8]; 
%  K2=[8.830*1e8 8.48*1e8 8.130*1e8 7.770*1e8 7.42*1e8 7.070*1e8 6.71*1e8 6.36*1e8 6.01*1e8 5.66*1e8];
% K3=[5.300*1e8 4.95*1e8 4.600*1e8 3.85*1e8 3.05*1e8 2.80*1e8 2.54*1e8...
% 2.280*1e8 2.080*1e8 1.72*1e8]; 
% K4=[1.430*1e8 1.140*1e8 8.018*1e7 5.223*1e7 3.556*1e7 2.091*1e7 9*1e6  1.156*1e6]; 
% K2i=[K1,K2,K3,K4];

 %  K2m=cumsum(K2i);
 % K=K2m(37);
%==YIBO========================
K1=[1.16E9 9.72E8 9.54E8 9.36E8 9.18E8 1.15E6 6.37E8...
    6.24E8 6.11E8 5.98E8 5.85E8  5.72E8 5.59E8 5.46E8 5.33E8 5.20E8];
K2=[5.07E8 4.94E8 4.81E8 4.65E8 4.55E8 4.42E8 4.29E8 4.16E8 4.03E8 3.90E8...
    9.14E5 9.14E5 9.14E5 9.14E5 9.14E5 1.43E6 1.43E6 1.98E8 1.89E8];
K3=[1.8E8 1.71E8 1.62E8 1.53E8 1.44E8 1.35E8 1.26E8 1.17E8 4.06E5 7.72E5 7.31E5...
    6.91E5 6.5E5 6.09E5 5.69E5 5.28E5 4.87E5 4.47E5 4.06E5 2E6 1E6  ];
K2i=[K1,K2,K3];
% K2m=cumsum(K2i);
% K=K2m(56);

 Ik=pi*(Db^4-(Db-2*thf)^4)/64;
Ak=pi*(Db^2-(Db-2*thf)^2)/4;
pk=@(z)m*g+g*Ak*rok*(L-z);
l=linspace(0,Lf,M);
for i=1:M-1
    K=K2i(i);
  lu=(l(i)+l(i+1))/2;
    p=pk(lu);

   %=========================
z=l(i);

 c=p/sqrt(4*E*Ik*(K-rok*Ak*om^2));
 d=((K-rok*Ak*om^2)/(4*E*Ik))^0.25;
a=sqrt(1-c)*d;
b=sqrt(1+c)*d;
r=p+(a^2-b^2)*E*Ik;
%=====================coefficients=====================================
t11=cosh(a*z)*cos(b*z);
t12=cosh(a*z)*sin(b*z);
t13=sinh(a*z)*cos(b*z);
t14=sinh(a*z)*sin(b*z);
t21=a*t14-b*t12;
t22=a*t13+b*t11;
t23=a*t11-b*t14;
t24=a*t12+b*t13 ;
t31=r*t11-2*b*a*E*Ik*t14;  
 t32=r*t12+2*b*a*E*Ik*t13;  
  t33=r*t13-2*b*a*E*Ik*t12; 
   t34=r*t14+2*b*a*E*Ik*t11;        
t15=0;                      t25=0;                  t35=0;
   
t41=(r*a-2*a*b^2*E*Ik)*t13-(r*b+2*a^2*b*E*Ik)*t12;
t42=(r*a-2*a*b^2*E*Ik)*t14+(r*b+2*a^2*b*E*Ik)*t11;
t43=(r*a-2*a*b^2*E*Ik)*t11-(r*b+2*a^2*b*E*Ik)*t14;
t44=(r*a-2*a*b^2*E*Ik)*t12+(r*b+2*a^2*b*E*Ik)*t13;
t15=0;
%=========================================================================
Tpi=[t11,t12 t13 t14 ; t21 t22 t23 t24 ; t31 t32 t33 t34 ; t41 t42 t43 t44];
cc=pinv(Tpi);
% %=========================================================================
K=K2i(i+1);
 z=l(i+1);
 c=p/sqrt(4*E*Ik*(K-rok*Ak*om^2));
 d=((K-rok*Ak*om^2)/(4*E*Ik))^0.25;
a=sqrt(1-c)*d;
b=sqrt(1+c)*d;
r=p+(a^2-b^2)*E*Ik ;
%=====================coefficients=====================================
t11=cosh(a*z)*cos(b*z); t21=a*t14-b*t12; t31=r*t11-2*b*a*E*Ik*t14 ; 
t12=cosh(a*z)*sin(b*z); t22=a*t13+b*t11; t32=r*t12+2*b*a*E*Ik*t13 ;
t13=sinh(a*z)*cos(b*z); t23=a*t11-b*t14;  t33=r*t13-2*b*a*E*Ik*t12 ;
t14=sinh(a*z)*sin(b*z); t24=a*t12+b*t13;   t34=r*t14+2*b*a*E*Ik*t11; 
% exp(a*z)
% cos(b*z)
% 
% t12
% t13
% t14


t15=0;                      t25=0;                  t35=0;
   
t41=(r*a-2*a*b^2*E*Ik)*t13-(r*b+2*a^2*b*E*Ik)*t12;
t42=(r*a-2*a*b^2*E*Ik)*t14+(r*b+2*a^2*b*E*Ik)*t11;
t43=(r*a-2*a*b^2*E*Ik)*t11-(r*b+2*a^2*b*E*Ik)*t14;
t44=(r*a-2*a*b^2*E*Ik)*t12+(r*b+2*a^2*b*E*Ik)*t13; 
t15=0;
%=========================================================================
Tpi1=[t11,t12 t13 t14 ; t21 t22 t23 t24 ; t31 t32 t33 t34 ; t41 t42 t43 t44];
  %Th=Tpi1*inv(Tpi)
  Th=Tpi1*cc;
     Hf=Th*Hf; % multiplication elemement per element 

end
end

