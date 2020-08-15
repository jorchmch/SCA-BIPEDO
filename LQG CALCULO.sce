//Cargamos el archivo guardado de "linearizing"
load("PLANTABIPEDO.sod","X","U","sys")

////Ponemos las matrices obtenidas de Scilab///
A=sys.A
B=sys.B
C=sys.C
D=sys.D

C=C(1:2,:);
D=zeros(2,2);


/////Controllability and Observability/////
///Matrices Scilab///
//controlabilidad
Cc = cont_mat(A,B)
rankCc=rank(Cc)
//observabilidad
O = obsv_mat(A, C)
rankO=rank(O)

///Matrices Analiticas///
//controlabilidad
//Cc1 = cont_mat(A1,B1)
//rankCc=rank(Cc1)
//observabilidad
//O1 = obsv_mat(A1,C1)
//rankO=rank(O1)

/////Plotear valores singulares//////
///Valores matrices Scilab///
G = syslin('c', A, B, C, D);
tr = trzeros(G)
w = logspace(-3,3);
sv = svplot(G,w);

//ploteo valores singulares ambos modelos//
//primer ploteo para las matrices de Scilab
scf(1);
plot2d("ln", w, [20*log(sv')/log(10)])
xgrid(12)
xtitle("Valores Singulares de la Planta","Frequency (rad/s)", "Amplitude (dB)");


////Obtencion de las funciones de transferencia////
//MatricesScilab//
[h]=ss2tf(sys)
//MatricesAnaliticas//
//[h1]=ss2tf(G1)

////Obtencion de polos y ceros/////
//Matrices Scilab//
scf(2)
plzr(h);
xtitle("Polos y zeros de la Planta")
//Matrices analiticas
//scf(5)
//plzr(h1);
//xtitle("Polos y zeros matrices analiticas")

// Escalonamiento a la planta //
su = diag( [0.9614, 0.2753] )
sx = diag( [3.157, 11.47, 3.157, 11.47] )
sy = diag( [3.157 3.157] )
 
ap_ = sx*A*inv(sx)
bp_ = sx*B*inv(su)
cp_ = sy*C*inv(sx)
dp_ = sy*D*inv(su)

Gs_= syslin("c",ap_, bp_, cp_, dp_)

// Valores singulares de la planta escalonada //
sv1 = svplot(Gs_,w);
scf(3)
plot2d("ln", w, [20*log(sv1')/log(10)])
xgrid(12)
xtitle("Valores Singulares de Planta Escalonada","Frequency (rad/s)", "Amplitude (dB)");

// Planta aumentada  con integradores antes del proyecto de controlador

[ns,nc] = size(bp_);                     //ns = número de entradas;  
                                        //nc = número de controles;   
a_1 = [ap_            bp_ ; 
       0*ones(nc,ns)  0*ones(nc,nc) ];
b_1 = [0*ones(ns,nc); eye(nc,nc)];
c_1 = [cp_      0*ones(nc,nc)];
d_1 = 0*ones(nc,nc)

Gs_1= syslin("c",a_1, b_1, c_1, d_1)

// Valores singulares de la planta escalonada con  el  integrador //
sv2 = svplot(Gs_1,w);
scf(5)
plot2d("ln", w, [20*log(sv2')/log(10)])
xgrid(12)
xtitle("Valores Singulares de Planta Escalonada con Integradores","Frequency (rad/s)", "Amplitude (dB)");

// LQR controller calculation
// Recuperar Target Loop resolviendo un problema de LQR barato
q = c_1'*c_1;          //Matriz de ponderación del estado
rho = 1e-9;        //Parámetro de recuperación de control barato

r = rho*eye(nc,nc)                 //Matriz de ponderación de control
         //how we calculate B
B=b_1*inv(r)*b_1';
A=a_1;
        //Solv the ricatti equation
X=riccati(A,B,q,'c','eigen');
        //matriz de ganacia
G_1=inv(r)*b_1'*X;

////// PREGUNTA 7 //////
//calculate observer Kalman Filter
ll =  inv(cp_*inv(-ap_)*bp_ + dp_);     
lh = -inv(ap_)*bp_*ll;
l = [lh                        //ll, lh - Para la conformación de 
     ll];                      //bucles de baja y alta frecuencia.

Gs_2= syslin("c",a_1, l, c_1, d_1)

// Valores singulares del filtro de bucle Abierto // 
sv3 = svplot(Gs_2,w);
scf(6)
plot2d("ln", w, [20*log(sv3')/log(10)])
xgrid(12)
xtitle("Valores Singulares de Filtro Abierto","Frequency (rad/s)", "Amplitude (dB)");

// Filtro de Kalman
pnint=eye(nc,nc)               //Proceso  de matriz de intensidad de ruido
mu=0.01;                       //Medicion de la intensidad de ruido

mnint=mu*eye(nc,nc)            //Matriz de intensidad de ruido  de medicion

Ch=l*l';                       //Forma de Ch para "riccati" segun  Scilab
Ah=a_1';                       //Forma de Ah para "riccati" segun  Scilab

Bh=c_1'*inv(mnint)*c_1;

Xh=riccati(Ah,Bh,Ch,'c','eigen');

                           //ganacia H
H_1=(inv(mnint)*c_1*Xh)';

Gs_3= syslin("c",a_1, H_1, c_1, d_1)

// Valores singulares del observador Filtro  Kalman
sv4 = svplot(Gs_3,w);
scf(7)
plot2d("ln", w, [20*log(sv4')/log(10)])
xgrid(12)
xtitle("Valores Singulares de Filtro de Kalman","Frequency (rad/s)", "Amplitude (dB)");


//ACTIVAR ESTA PARTE PARA LOS POLOS Y CEROS DEL FILTRO DE KALMAN GANACIA H 
[h2]=ss2tf(Gs_3)
scf(8)
plzr(h2);
xtitle("Polos y Zeros del Filtro Kalman, Ganancia H")

/////// PREGUNTA 8 ////////
// COMPENSADOR K(S) DE LA FORMA DEL PORF. RODRIGUEZ //
ak = [ a_1-b_1*G_1-H_1*c_1  0*ones(ns+nc,nc)
       G_1                  0*ones(nc,nc) ];
bk = [ H_1
       0*ones(nc,nc) ];
ck = [0*ones(nc, ns+nc) eye(nc,nc) ];
dk = [0*eye(nc,nc)];

Gs_4= syslin("c",ak, bk, ck, dk)
// Valores singulares del compensador "K(s)" //
sv4 = svplot(Gs_4,w);
scf(9)
plot2d("ln", w, [20*log(sv4')/log(10)])
xgrid(12)
xtitle("Valores Singulares del compensador Ks","Frequency (rad/s)", "Amplitude (dB)");

/////// PREGUNTA 9 ///////
//SENSIBILIDAD "S" Y SENSIBILIDAD COMPLEMENTARIA "T"
//Analisis en bucle abierto
al = [ ap_                     bp_*ck
       0*ones(ns+nc+nc,ns)     ak    ];
bl = [ 0*ones(ns,nc)
       bk ];
cl = [ cp_  0*ones(nc,ns+nc+nc) ];
dl = [0*eye(nc,nc)];

Gs_5= syslin("c",al, bl, cl, dl)

// Valores Singulares de bucle abierto//
sv5 = svplot(Gs_5,w);
scf(10)
plot2d("ln", w, [20*log(sv5')/log(10)])
xgrid(12)
xtitle("Valores Singulares del bucle abierto","Frequency (rad/s)", "Amplitude (dB)");

// Valores singulares de sensibilidad S //
Gs_6= syslin("c",al-bl*cl, bl, -cl, eye(nc,nc))
sv6 = svplot(Gs_6,w);
scf(11)
plot2d("ln", w, [20*log(sv6')/log(10)])
xgrid(12)
xtitle("Ploteo de la Sensibilidad","Frequency (rad/s)", "Amplitude (dB)");

// Valores singulares de sensibilidad complementaria T //
Gs_7= syslin('c',al-bl*cl, bl, cl, dl)
sv7 = svplot(Gs_7,w);
scf(12)
plot2d("ln", w, [20*log(sv7')/log(10)])
xgrid(12)
xtitle("Ploteo de la Sensibilidad Complementaria","Frequency (rad/s)", "Amplitude (dB)");

// Valorers Singulares de S y T juntos // 
scf(13)
plot2d("ln", w, [20*log(sv6')/log(10)])
plot2d("ln", w, [20*log(sv7')/log(10)])
xgrid(12)
xtitle("Sensibilidad y Sensibilidad Complementaria","Frequency (rad/s)", "Amplitude (dB)");

