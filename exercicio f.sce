//Daniel Marins Silva D´Oliveira Martins : 12554571
//David Lopes F. de Carvalho :             12610305
//João Pedro Kawachi Chaves :              12745188
//Leonardo Fortunato de Carvalho :         9837910


//Cenário aleatório - Métodos M1 e M2

clear; // Limpeza de variáveis
close(winsid()) // Fecha as janelas gráficas abertas anteriormente


//variáveis utilizadas:

lambda = 1.7 //L 1/L 2 = 1+(1+5+8+0)/20
wp = 1 //rad/s  -  wp = (g/L 1)^(1/2)

//vetor de estados do estado inicial

theta_1_ini = 0.01   //rad
theta_2_ini = 0.01    //rad
theta_1_p_ini = -0.01      //rad/s
theta_2_p_ini = 0       //rad/s


X0 = [theta_1_ini; theta_2_ini; theta_1_p_ini; theta_2_p_ini] //vetor de espaços do estado inicial


// intervalos de tempo

T = 100 //s - tempo total
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
tic();
y = ode("adams",X0,0,t,penduloduplo);

tempo_adams = toc();

theta_1 = y(1,:); //Theta 1 (rad)
theta_2 = y(2,:); //Theta 2 (rad)
theta_1_p = y(3,:); //Theta 1 ponto (rad/s)
theta_2_p = y(4,:); //Theta 2 ponto (rad/s)



// Metodo 2: método de Runge-Kutta



tic();
z = ode("rk",X0,0,t,penduloduplo);

tempo_rk = toc();

theta_1z = z(1,:); //Theta 1 (rad)
theta_2z = z(2,:); //Theta 2 (rad)
theta_1_pz = z(3,:); //Theta 1 ponto (rad/s)
theta_2_pz = z(4,:); //Theta 2 ponto (rad/s)


scf(5);
plot2d(t,[theta_1',theta_1z'],[2,3]);
title("Theta 1 em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);
legend(['Método 1';'Método 2'],4);

scf(6);
plot2d(t,[theta_2',theta_2z'],[2,3]);
title("Theta 2 em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Posição angular(rad)", 'fontsize', 3);
legend(['Método 1';'Método 2'],4);


scf(7);
plot2d(t,[theta_1_p',theta_1_pz'],[2,3]);
title("Theta 1 ponto em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);
legend(['Método 1';'Método 2'],4);

scf(8);
plot2d(t,[theta_2_p',theta_2_pz'],[2,3]);
title("Theta 2 ponto em função do tempo" ,'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("Velocidade angular(rad/s)", 'fontsize', 3);
legend(['Método 1';'Método 2'],4);


//Comparando método 1 e 2

dif_theta_1 = theta_1 - theta_1z
dif_theta_2 = theta_2 - theta_2z
dif_theta_1_p = theta_1_p - theta_1_pz
dif_theta_2_p = theta_2_p - theta_2_pz

//dif theta 1
scf(9);
plot2d(t,dif_theta_1,5);
title("diferença de theta 1 (met 1 - met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta1(rad)", 'fontsize', 3);

//dif theta 2
scf(10);
plot2d(t,dif_theta_2,6);
title("diferença de theta 2 (met 1 - met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta2(rad)", 'fontsize', 3);

//dif theta 1 ponto
scf(11);
plot2d(t,dif_theta_1_p,7);
title("diferença de theta 1 ponto (met 1 - met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta1 p(rad/s)", 'fontsize', 3);

//dif theta 2 ponto
scf(12);
plot2d(t,dif_theta_2_p,6);
title("diferença de theta 2 ponto (met 1 - met 2)", 'fontsize', 5);
xlabel("t(s)", 'fontsize', 3);
ylabel("d theta2 p(rad/s)", 'fontsize', 3);

disp("Tempo adams:",tempo_adams, "; Tempo Runge-Jutta:",tempo_rk)


