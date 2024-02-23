format("v")
//________________________________ déclaration drs données du problème  ____________________________________

AU = 149597870690; // 1 unité astronomique en m
r0 = AU; // Rayon de l'orbite initiale en m
rf = 1.5 * AU; // Rayon de l'orbite cible finale en m
m0 = 1000; // Masse satellite à t0 en kg
T = 0.1:0.1:0.6; // Poussée du moteur en N (crée un vecteur de 0.1 à 0.6 avec un pas de 0.01)
g0 = 9.80665; // Gravité terrestre en ms^-2
Isp = 3000; // Impulsion spécifique du moteur en s
mu_body = 1.32712440018e+20; // constante gravitationnelle du soleil en m^3s^-2

r0=r0;
u0=0;
v0=(mu_body/r0)^0.5

rf=rf;
uf=0;
vf=(mu_body/rf)^0.5


//on doit chnager la précision du solveur car r(tf)-rf) à un ordre de grandeur de 12 et on doit avoir 10 chiffres significatif à la fin donc en tout 22 chiffres significatif. Cela est impossible

//________________________________ Pourquoi normaliser ?  ____________________________________

//toujours se poser la question de l'ordre de grandeur
// on vas donc normaliser pour éviter le pb de chiffre significatif
// pour r0 on choisis DU=r0; cad r0=r0/DU => r0=1
//VU = v0 
//MU =m0

DU=r0; //Distance unitaire
VU=v0; // Vitesse unitaire
MU=m0; // Masse unitaire
TU=DU/VU; // Temps unitaire
FU=(MU*DU)/TU^2 // Force unitaire




//  ____________________________________ definition de la strucutre param  ____________________________________

param.DU = DU;
param.VU = VU;
param.MU = MU;
param.TU = TU;
param.FU = FU;

// pour tester dynpol 

/* param.mu_boddy=1
param.T = 0.01
param.g0_Isp = 1
param.m0 = 1
param.t0 = 0
param.t=0  */

// pour gnultmin

param.mu_body = mu_body/param.DU^3*param.TU^2;
param.T = T(1)/ param.FU;
param.g0_Isp = (g0*Isp)/param.VU;
param.m0 = m0/param.MU;
param.t0 = 0;
param.x_t0 = [r0/param.DU 
                u0/param.VU 
                v0/param.VU];

param.x_tf = [rf/param.DU 
                uf/param.VU
                vf/param.VU]; 

//________________________________ Décrire la dynamique du système : état + état adjoint  ____________________________________

//x= [r; u; v; l_r; l_u; l_v;]

function dx = dynpol(t, x)
    mu = param.mu_body;
    T = param.T;
    g0_Isp = param.g0_Isp;
    m0 = param.m0;
    t0 = param.t0;
    dm = -T / g0_Isp;
    phi = atan(x(5), x(6));
    dx=x;
    dx(1) = x(2);
    dx(2) = (x(3)^2) / x(1) - mu / (x(1)^2) + (T / (m0 + dm * t)) * sin(phi);
    dx(3) = -x(2) * x(3) / x(1) + (T / (m0 + dm * t)) * cos(phi);
    dx(4) = -x(5) * ((-x(3)^2) / (x(1)^2) + 2 * mu / (x(1)^3)) - x(6) * (x(2) * x(3) / (x(1)^2));
    dx(5) = -x(4) + x(6) * (x(3) / x(1));
    dx(6) = -x(5) * (2 * x(3) / x(1)) + x(6) * (x(2) / x(1));
endfunction


//x_test=[1;2;3;4;-5;-6]



//dxTest=dynpol(param.t0,x_test)

/*disp(dxTest(1))
disp(dxTest(2))
disp(dxTest(3))
disp(dxTest(4))
disp(dxTest(5))
disp(dxTest(6)) */





function gnmt = gnultmin(p)
    t0 = param.t0;
    r0 = param.x_t0(1);
    u0 = param.x_t0(2);
    v0 = param.x_t0(3);
    rf = param.x_tf(1);
    uf = param.x_tf(2);
    vf = param.x_tf(3);
    y0=[r0 u0 v0 p(2) p(3) p(4)];
    tf=p(1)
    atol=1e-10
    rtol=1e-10
    y=ode("rk",y0,t0,tf,dynpol, atol, rtol)
    gnmt=p;
    gnmt(1)=y(1)- rf ;
    gnmt(2)=y(2)-uf;
    gnmt(3)=y(3)-sqrt(param.mu_body/rf);
    gnmt(4)=p(2)^2+p(3)^2+p(4)^2-1;
endfunction



p0=[3;0.6;0.3;0.7]
test = gnultmin(p0);

disp(test(1))
disp(test(2))
disp(test(3))
disp(test(4))
disp(param);

// on résoud gnultmin(p)=0



for i = 1:6
    param.T = T(i)/ param.FU;
    [x, v, info] = fsolve(p0, gnultmin);
    disp(info);
    disp(v);
    InitialCond = [param.x_t0(1); param.x_t0(2); param.x_t0(3); x(2); x(3); x(4)];
    t = 0:0.01:x(1);
    sol = ode(InitialCond, param.t0, t, dynpol);
    timeInDays = t * param.TU / 86400;
    AngleInDegrees = atan(sol(5,:), sol(6,:))*360/2/%pi;
    figure;
    plot(timeInDays, AngleInDegrees);
    xlabel('Time (days)');
    ylabel('Control in deg');
    powerInNewton=param.T * param.FU
    xgrid
    h = gca();
    h.background=color('white');
   title('Control history for '+ string(powerInNewton) + 'N');   
end







