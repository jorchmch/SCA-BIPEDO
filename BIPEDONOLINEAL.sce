clc;clear;
function xdot=BIPEDO(u1,u2,u3,u4,u5,u6)
    
    //Abrimos los parametros del Robot "ParametrosBipedo.sce"
    exec('ParametrosBipedo.sce', -1);
    
    //colocamos mientros estados
q1=u1;  //angulo
q1d=u3; //velocidad angular
q2=u2;   //posicion
q2d=u4;  //velocidad
tau1=u5;
tau2=u6;   //control-fuerza

//%matriz M
M11=m*a^2+m*l^2+mh*l^2;
M21=-m*l*b*cos(u1-u2);
M12=M21;
M22=m*b^2;
M=[M11 M12; M21 M22];

//matriz C
C11=0;
C21=m*l*b*u3*sin(u1-u2);
C12=-m*l*b*u4*sin(u1-u2);
C22=0;
C=[C11 C12; C21 C22];

//matriz g 
G11=-m*g*a*sin(u1)-m*g*l*sin(u1)-mh*g*l*sin(u1);
G21=m*g*b*sin(u2);
G=[G11;G21];

tau=[u5;u6];
qd=[u3;u4];

//espacio de estados obtenidos del item 1.
xdot=[qd;
    inv(M)*(tau-C*qd-G)];
      
endfunction
