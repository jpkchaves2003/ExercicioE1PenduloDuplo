//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


//Cenário C2 - Métodos M1 e M2 - Modelo não linear

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


theta_1z = zeros(Np); //Theta 1 (rad)
theta_2z = zeros(Np); //Theta 2 (rad)
theta_1_pz = zeros(Np); //Theta 1 ponto (rad/s)
theta_2_pz = zeros(Np); //Theta 2 ponto (rad/s)


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

//Método 1: método de Adams

y = ode("adams",X0,0,t,penduloduplo);

theta_1 = y(1,:); //Theta 1 (rad)
theta_2 = y(2,:); //Theta 2 (rad)
theta_1_p = y(3,:); //Theta 1 ponto (rad/s)
theta_2_p = y(4,:); //Theta 2 ponto (rad/s)




// Metodo 2: método de Runge-Kutta

//vetores para armazenar os valores de thetas e thetas ponto


z = ode("rk",X0,0,t,penduloduplo);

theta_1z = z(1,:); //Theta 1 (rad)
theta_2z = z(2,:); //Theta 2 (rad)
theta_1_pz = z(3,:); //Theta 1 ponto (rad/s)
theta_2_pz = z(4,:); //Theta 2 ponto (rad/s)


//Energia mecânica

//Em = T + V
// Ec (energia cinética)= T = mi * L2 ^3 *(  1/6*lambda^3*theta_1_p^2 + 1/6*theta_2_p^2 + 1/2*lambda^2*theta_1_p^2 + lambda*theta_1_p*theta_2_p*cos(theta_1 - theta_2)  )
// Ep (energia potencial) = V = -mi*L2^3 *(  (1/2*lambda^2 + lambda)*cos(theta_1) + 1/2*lambda*cos(theta_2)  )

Em1 = zeros(Np)
Em2 = zeros(Np)
Ec1 = zeros(Np)
Ec2 = zeros(Np)
Ep1 = zeros(Np)
Ep2 = zeros(Np)

//definição de valores de mi e L2:
mi = 1 //kg/m, massa específica linear das barras
L2 = 1 //m, comprimento da barra 2



for n = 1 : Np
    Ec1(n) = mi * L2 ^3 *(  1/6*lambda^3*theta_1_p(n)^2 + 1/6*theta_2_p(n)^2 + 1/2*lambda^2*theta_1_p(n)^2 + lambda*theta_1_p(n)*theta_2_p(n)*cos(theta_1(n) - theta_2(n))  )
    Ep1(n) = -mi*L2^3 *(  (1/2*lambda^2 + lambda)*cos(theta_1(n)) + 1/2*lambda*cos(theta_2(n))  )
    Em1(n) = Ec1(n) + Ep1(n)
    
    Ec2(n) = mi * L2 ^3 *(  1/6*lambda^3*theta_1_pz(n)^2 + 1/6*theta_2_pz(n)^2 + 1/2*lambda^2*theta_1_pz(n)^2 + lambda*theta_1_pz(n)*theta_2_pz(n)*cos(theta_1z(n) - theta_2z(n))  )
    Ep2(n) = -mi*L2^3 *(  (1/2*lambda^2 + lambda)*cos(theta_1z(n)) + 1/2*lambda*cos(theta_2z(n))  )
    Em2(n) = Ec2(n) + Ep2(n)
end



//Comparando método 1 e 2

dif_Em = Em1 - Em2

//dif Em
scf(1);
plot2d(t,dif_Em,5);
title("diferença de Energia mecânica (met 1 - met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d Em(J)", 'fontsize', 3);


//E1
scf(2);
plot2d(t,Em1,5);
title("Energia mecânica (met 1)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Em(J)", 'fontsize', 3);


//dif Em
scf(3);
plot2d(t,Em2,5);
title("Energia mecânica (met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Em(J)", 'fontsize', 3);


