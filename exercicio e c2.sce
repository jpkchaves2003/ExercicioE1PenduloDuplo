//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


//Cenário C2

clear; // Limpeza de variáveis
close(winsid()) // Fecha as janelas gráficas abertas anteriormente


//variáveis utilizadas:

lambda = 1.7 //L 1/L 2 = 1+(1+5+8+0)/20
wp = 1 //rad/s  -  wp = (g/L 1)^(1/2)

//vetor de estados do estado inicial

theta_1_ini = 0   //rad
theta_2_ini = 0    //rad
theta_1_p_ini = 0.1      //rad/s
theta_2_p_ini = 0       //rad/s


X0 = [theta_1_ini; theta_2_ini; theta_1_p_ini; theta_2_p_ini] //vetor de espaços do estado inicial


// intervalos de tempo

T = 60 //s - tempo total
Np = 10000 // numero de passos
t = linspace(0,T,Np); //vetor de tempo

//vetores para armazenar os valores de thetas e thetas ponto

theta_1= zeros(Np) //theta 1
theta_2= zeros(Np) //theta 2
theta_1_p= zeros(Np) //theta ponto 1
theta_2_p= zeros(Np) //theta ponto 2


//função não linear
function dy = penduloduplo(t,y)  //integração das equações diferenciais
    a11 = (lambda/3 + 1)*lambda
    a12 = cos(y(1)-y(2))/2
    a21 = lambda*cos(y(1)-y(2))/2
    a22 = 1/3
    b1 = -y(3)*y(4)/2 * sin(y(1)-y(2)) + (lambda/2 + 1)*wp^2*sin(y(1))
    b2 = ((-lambda*y(3)^2)/2)*(sin(y(1)-y(2))) + ((wp^2)/2)*lambda*sin(y(2))
    
    
    dy(1)=y(3);
    dy(2)= y(4);
    dy(3)=(-a12*b2-a22*b1)/(-a12*a21+a11*a22);
    dy(4)= (a11*b2-a21*b1)/(-a12*a21+a11*a22);
endfunction


//resolvendo a EDO
y = ode(X0,0,t,penduloduplo);

theta_1 = y(1,:); //Theta 1 (rad)
theta_2 = y(2,:); //Theta 2 (rad)
theta_1_p = y(3,:); //Theta 1 ponto (rad/s)
theta_2_p = y(4,:); //Theta 2 ponto (rad/s)


scf(1);
plot2d(t,theta_1,2);
title("Theta 1 não linear", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);

scf(2);
plot2d(t,theta_2,3);
title("Theta 2 não linear", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);

scf(3);
plot2d(t,theta_1_p,4);
title("Theta 1 ponto não linear", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);

scf(4);
plot2d(t,theta_2_p,5);
title("Theta 2 ponto não linear", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);

// modelo linear

//vetores para armazenar os valores de thetas e thetas ponto

theta_1z = zeros(Np); //Theta 1 (rad)
theta_2z = zeros(Np); //Theta 2 (rad)
theta_1_pz = zeros(Np); //Theta 1 ponto (rad/s)
theta_2_pz = zeros(Np); //Theta 2 ponto (rad/s)


function dz = penduloduplolinear(t,z)  //integração das equações diferenciais
    A11 = (lambda/3 + 1)*lambda
    A12 = 1/2
    A21 = lambda/2
    A22 = 1/3
    B1 = (lambda/2 + 1)*wp^2*lambda*z(1)
    B2 = ((wp^2)/2)*lambda*z(2)
    
    
    dz(1)=z(3);
    dz(2)= z(4);
    dz(3)=(-A12*B2+A22*B1)/(-A12*A21+A11*A22);
    dz(4)= (A11*B2-A21*B1)/(-A12*A21+A11*A22);
endfunction

//Resolvendo a EDO Linear

z = ode(X0,0,t,penduloduplolinear);

theta_1z = z(1,:); //Theta 1 (rad)
theta_2z = z(2,:); //Theta 2 (rad)
theta_1_pz = z(3,:); //Theta 1 ponto (rad/s)
theta_2_pz = z(4,:); //Theta 2 ponto (rad/s)


scf(5);
plot2d(t,theta_1z,2);
title("Theta 1 linear" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);

scf(6);
plot2d(t,theta_2z,3);
title("Theta 2 linear" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);

scf(7);
plot2d(t,theta_1_pz,4);
title("Theta 1 ponto linear" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);

scf(8);
plot2d(t,theta_2_pz,5);
title("Theta 2 ponto linear" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);


//Comparando linear e não linear

dif_theta_1 = theta_1 - theta_1z
dif_theta_2 = theta_2 - theta_2z
dif_theta_1_p = theta_1_p - theta_1_pz
dif_theta_2_p = theta_2_p - theta_2_pz

//dif theta 1
scf(9);
plot2d(t,dif_theta_1,5);
title("diferença de theta 1 (nlin - lin)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta1(rad)", 'fontsize', 3);

//dif theta 2
scf(10);
plot2d(t,dif_theta_2,6);
title("diferença de theta 2 (nlin - lin)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta2(rad)", 'fontsize', 3);

//dif theta 1 ponto
scf(11);
plot2d(t,dif_theta_1_p,7);
title("diferença de theta 1 ponto (nlin - lin)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta1 p(rad/s)", 'fontsize', 3);

//dif theta 2 ponto
scf(12);
plot2d(t,dif_theta_2_p,6);
title("diferença de theta 2 ponto (nlin - lin)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta2 p(rad/s)", 'fontsize', 3);


