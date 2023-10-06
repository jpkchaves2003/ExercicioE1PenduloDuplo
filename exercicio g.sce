//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


//Cenário C1

//Método 1: Adams
//Método 2: Runge-Kutta

clear; // Limpeza de variáveis
close(winsid()) // Fecha as janelas gráficas abertas anteriormente


//variáveis utilizadas:

lambda = 1.7 //L 1/L 2 = 1+(1+5+8+0)/20
wp = 1 //rad/s  -  wp = (g/L 1)^(1/2)

//vetor de estados do estado inicial

theta_1_ini = 1.0e-15   //rad
theta_2_ini = 1.0e-15    //rad
theta_1_p_ini = 0      //rad/s
theta_2_p_ini = 0       //rad/s


X0 = [theta_1_ini; theta_2_ini; theta_1_p_ini; theta_2_p_ini] //vetor de espaços do estado inicial


// intervalos de tempo

T = 40 //s - tempo total
Np = 10000 // numero de passos
t = linspace(0,T,Np); //vetor de tempo

//vetores para armazenar os valores de thetas e thetas ponto

//Não linear método 1
theta_1= zeros(Np) //Theta 1 (rad)
theta_2= zeros(Np) //Theta 2 ponto (rad/s)
theta_1_p= zeros(Np) //Theta 2 ponto (rad/s)
theta_2_p= zeros(Np) //Theta 2 ponto (rad/s)

//Não linear método 2
theta_1z = zeros(Np); //Theta 1 (rad)
theta_2z = zeros(Np); //Theta 2 ponto (rad/s)
theta_1_pz = zeros(Np); //Theta 2 ponto (rad/s)
theta_2_pz = zeros(Np); //Theta 2 ponto (rad/s)

//Não linear método 1
theta_1k= zeros(Np) //Theta 1 (rad)
theta_2k= zeros(Np) //Theta 2 ponto (rad/s)
theta_1_pk= zeros(Np) //Theta 2 ponto (rad/s)
theta_2_pk= zeros(Np) //Theta 2 ponto (rad/s)

//Não linear método 2
theta_1w = zeros(Np); //Theta 1 (rad)
theta_2w = zeros(Np); //Theta 2 ponto (rad/s)
theta_1_pw = zeros(Np); //Theta 2 ponto (rad/s)
theta_2_pw = zeros(Np); //Theta 2 ponto (rad/s)



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

//função linear
function dL = penduloduplolinear(t,L)  //integração das equações diferenciais
    A11 = (lambda/3 + 1)*lambda
    A12 = 1/2
    A21 = lambda/2
    A22 = 1/3
    B1 = (lambda/2 + 1)*wp^2*lambda*z(1)
    B2 = ((wp^2)/2)*lambda*z(2)
    
    
    dL(1)=L(3);
    dL(2)= L(4);
    dL(3)=(-A12*B2+A22*B1)/(-A12*A21+A11*A22);
    dL(4)= (A11*B2-A21*B1)/(-A12*A21+A11*A22);
endfunction


//resolvendo a EDO

//Não linear e método 1
y = ode("adams",X0,0,t,penduloduplo);

theta_1 = y(1,:); //Theta 1 (rad)
theta_2 = y(2,:); //Theta 2 (rad)
theta_1_p = y(3,:); //Theta 1 ponto (rad/s)
theta_2_p = y(4,:); //Theta 2 ponto (rad/s)

//Não Linear e método 2
z = ode("rk",X0,0,t,penduloduplo);

theta_1z = z(1,:); //Theta 1 (rad)
theta_2z = z(2,:); //Theta 2 (rad)
theta_1_pz = z(3,:); //Theta 1 ponto (rad/s)
theta_2_pz = z(4,:); //Theta 2 ponto (rad/s)

//Linear e método 1
k = ode("adams",X0,0,t,penduloduplolinear);

theta_1k = k(1,:); //Theta 1 (rad)
theta_2k = k(2,:); //Theta 2 (rad)
theta_1_pk = k(3,:); //Theta 1 ponto (rad/s)
theta_2_pk = k(4,:); //Theta 2 ponto (rad/s)

//Linear e método 2
w = ode("rk",X0,0,t,penduloduplolinear);

theta_1w = w(1,:); //Theta 1 (rad)
theta_2w = w(2,:); //Theta 2 (rad)
theta_1_pw = w(3,:); //Theta 1 ponto (rad/s)
theta_2_pw = w(4,:); //Theta 2 ponto (rad/s)


//"nl" ->(normal, log) pros eixos x,y ; abs() é a função de módulo do valor
scf(1);
plot2d("nl",t,[abs(theta_1'),abs(theta_1z'), abs(theta_1k'), abs(theta_1w')],[2,3,4,5]);
title("Theta 1 em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("theta 1 (rad)", 'fontsize', 3);
legend(['M1 N-Lin','M2 N-Lin', 'M1 Lin','M2 Lin'],4);

scf(2);
plot2d("nl",t,[abs(theta_2'),abs(theta_2z'), abs(theta_2k'), abs(theta_2w')],[2,3,4,5]);
title("Theta 2 em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("theta 2 (rad)", 'fontsize', 3);
legend(['M1 N-Lin','M2 N-Lin', 'M1 Lin','M2 Lin'],4);

scf(3);
plot2d("nl",t,[abs(theta_1_p'),abs(theta_1_pz'), abs(theta_1_pk'), abs(theta_1_pw')],[2,3,4,5]);
title("Theta 1p em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("theta 1p (rad/s)", 'fontsize', 3);
legend(['M1 N-Lin','M2 N-Lin', 'M1 Lin','M2 Lin'],4);


scf(4);
plot2d("nl",t,[abs(theta_2_p'),abs(theta_2_pz'), abs(theta_2_pk'), abs(theta_2_pw')],[2,3,4,5]);
title("Theta 2p em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("theta 2p (rad/s)", 'fontsize', 3);
legend(['M1 N-Lin','M2 N-Lin', 'M1 Lin','M2 Lin'],4);






