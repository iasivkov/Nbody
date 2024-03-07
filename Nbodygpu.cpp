#include "Nbodygpu.h"
#include "Global_constants.h"

//массы компонент
const float PI = 3.1415926;

float Md = 1.0*1.0;
float Mh = 1.0*5.8;
float Mb = 1.0*0.3333;
float Ms = 0.1*Md;


//характерные размеры диска
float h = 1.0;

float z0 = 0.2;
float r0 = 1.0/3.5;

float ksi = PI/6.0;
float xs = 6.0*h;
float ys = 6.0*h;
float zs = sqrt(xs*xs + ys*ys) * sin(ksi)/cos(ksi); 
float4 rs = {xs, ys, zs, 0.0 };
float rabs = sqrtf(xs*xs + ys*ys)/cos(ksi); 

float orb;

//-----------------------------------------------------------------------------


float Qt =  1.5;
float rho0= Md/(4.0*PI*h*h*z0);
//хар-ки гало
float gam = 1.0;
float rc = 10.0;//радиус обрезания
const float q = gam/rc;
float alpha =1.0/(1 - sqrt(PI)*q*exp(q*q)*(1 - erf(q)));
//ха-рки балджа
float a = 0.2;
float c = 0.1;

float4 rref = {sqrt(3.0)*h, sqrt(3.0)*h, 0.0, 0.0};

// радиус системы
float  R = 25.0;

float Rs = 5.0;

// параметр точности
float  theta = 0.25;

// шаг по времени
float  TimeStep = 0.08;

// гравитационная постоянная
float	G = 1.0;



float4 *X;
float4 *A;
float4 *V;
float Bin[10000];
float Bin2[10000];
float Bin3[10000];
float Bin4[1000];
float Bin5[100000];
float Bin6[30000];
float BinD[128];
float BinDz[128];

alglib::hqrndstate state;
//-----------------------------------------------------------------------------

float  random()
{	
	return alglib::hqrnduniformr(state);
	//return  (float)rand() / (float)RAND_MAX ;
}

double r_for_potential2;
void potential(double k, double xminusa, double bminusx, double &y, void *ptr) 
{

	double r = r_for_potential2;
  y = (Mb*a/(2.0*PI*r)/((r + a)*(r + a)*(r + a)) + Mh*alpha*sqrt(PI)/2.0/PI/PI/rc*exp(-r*r/rc/rc)/(r*r + gam*gam))*r;
}

void potential2(double k, double xminusa, double bminusx, double &y, void *ptr) 
{

	double r = r_for_potential2;
  y = (Mb*a/(2.0*PI*r)/((r + a)*(r + a)*(r + a)) + Mh*alpha*sqrt(PI)/2.0/PI/PI/rc*exp(-r*r/rc/rc)/(r*r + gam*gam))*r*r;
}


float Phi(float r)
{
		alglib::autogkstate s;
 		double v;
		alglib::autogkreport rep;
		r_for_potential2=r;
		
		alglib:: autogksmooth(r, 100*h, s);
		autogkintegrate(s,potential);
		alglib:: autogkresults(s, v, rep);

		alglib::autogkstate s2;
		double g;
		alglib::autogkreport rep2;
		alglib:: autogksmooth(0, r, s2);
		autogkintegrate(s2,potential2);
		alglib:: autogkresults(s2, g, rep2);

		return -4.0*PI*G*(g/r+v);
}
//-----------------------------------------------------------------------------
float Q(float r)
{			
	
			//гало, балдж, диск все компоненты
			//if(r<0.01)return -r;
			int L=r/R*10000;
			if(r==R)L=9999;
			r=(L+1.0)*R/10000.0;
			float sum = 0;
			/*for(int i = 0; i<=L; i++)
				{
					sum += Bin[i];
				}*/
			return -G*Bin[L]/r;
		
	

	
}

float Q2(float r)
{			
	//гало, балдж
		float step = R*0.0001;
			//if(r<0.01)return -r;
			int L=r/step;
			//if(r == R)L=9999;
			//r=(L+1.0)*R/10000.0;
			//float sum = 0;
			/*for(int i = 0; i<=L; i++)
				{
					sum += Bin3[i];
				}*/
			return -G*Bin2[L]/((L+1.0)*step);

	
}

float Q3(float r)
{			
	
			float step = R*0.0001;
			//if(r<0.01)return -r;
			int L=r/step;
			//if(r == R)L=9999;
			//r=(L+1.0)*R/10000.0;
			//float sum = 0;
			/*for(int i = 0; i<=L; i++)
				{
					sum += Bin3[i];
				}*/
			return -G*Bin3[L]/((L+1.0)*step);
		

	
}

float Qs(float r)
{
			//if(sqrtr<0.01)return -r;
			int L=r /Rs*1000;
			if(r == Rs)L=999;
			r=(L+1.0)*Rs/1000.0;
			float sum = 0;
			for(int i = 0; i<=L; i++)
				{
					sum += Bin4[i];
				}
			return -G*sum/r;

}

float Q5(float r)
{
			//if(sqrtr<0.01)return -r;
			int L=r /R*100000;
			//if(r == R)L=5999;
			r=(L+1.0)*R/100000.0;

			
			/*for(int i = 0; i<=L; i++)
				{
					sum += Bin5[i];
				}*/
			return -G*Bin5[L]/r;

}

float Q6(float r)
{
			//if(sqrtr<0.01)return -r;
			int L=r /R*30000;
			if(r == R)L=29999;
			r=(L+1.0)*R/30000;

			float sum = Bin6[L]*Mh/(double)NforHalo;

			/*for(int i = 0; i<=L; i++)
				{
					sum += Bin6[i];
				}
				*/
			return -G*sum/r;

}



//экспоненциальный интеграл
float Ei(float x)
{
	return	alglib::exponentialintegralei(x);
}
//-----------------------------------------------------------------------------
//куммулятивная масса диска
float Mu(float r,float z)
{
	if(z<0)z=z*(-1.0);
	return (1 - 1 * exp(-r/h) - exp(-r/h)*r/h)*sinh(z/z0)/cosh(z/z0);
	
	
}

float Surfd(float r)
	{
	return Md/(2.0*PI*h*h)*exp(-r/h);
}

//Hernquist
float	Rhod(float r )
{
	return Surfd(r)*2.0*PI*r;
}

//Sellwood
/*float	Rhod_1(float r )
{
	return 2*Md/(3.0*PI*h*h)*(exp(-r/h)-exp(-2*r/h))*2.0*PI*r;
}
*/

float Sech(float z)
{
	return 4.0/(exp(z/z0)+ exp(-z/z0))/(exp(z/z0)+ exp(-z/z0))/2.0/z0;
}


//вспомогательные константы для интегрирования
double r_for_potential;
double z_for_potential;

//потенциал диска
void disk_potential(double k, double xminusa, double bminusx, double &y, void *ptr) 
{
  y = -2.0*PI*G*rho0*h*h*2.0*z0 * alglib::besselj0(k*r_for_potential) * exp(-k*fabs(z_for_potential)) * sqrt(1 + k*k * h*h)/((1 + k*k * h*h)*(1 + k*k * h*h));
}

void disk_potential_derivation(double k, double xminusa, double bminusx, double &y, void *ptr) 
{
   y = 2.0*PI*G*rho0*h*h*2.0*z0 * alglib::besselj1(k*r_for_potential) * k * exp(-k*fabs(z_for_potential)) * sqrt(1 + k*k * h*h)/((1 + k*k * h*h)*(1 + k*k * h*h));
}
void disk_potential_second_derivation(double k, double xminusa, double bminusx, double &y, void *ptr) 
{
   y = 2.0*PI*G*rho0*h*h*2.0*z0 * (alglib::besselj0(k*r_for_potential) - alglib::besselj1(k*r_for_potential)/(k*r_for_potential)) * k*k * exp(-k*fabs(z_for_potential)) * sqrt(1 + k*k * h*h)/((1 + k*k * h*h)*(1 + k*k * h*h));
}


double Qd(double r,double z)
{
  
   
		alglib::autogkstate s;
 		double v;
		alglib::autogkreport rep;
		r_for_potential=r;
		z_for_potential=z;

		alglib:: autogksmooth(0, 25*h, s);
		autogkintegrate(s,disk_potential);
		alglib:: autogkresults(s, v, rep);

    return v;
}

float dQ(float r, float z)
{
	
    alglib::autogkstate s;
 	double v;
    alglib::autogkreport rep;
	r_for_potential=r;
	z_for_potential=z;

    alglib:: autogksmooth(0, 25*h, s);
	autogkintegrate(s,disk_potential_derivation);
    alglib:: autogkresults(s, v, rep);
	if(!is_valid(v))printf("\n integration %f: %i,\t%i,\t%i",r, rep.terminationtype, rep.nfev, rep.nintervals);
	
    return v;
}

float ddQ(float r, float z)
{
	
    alglib::autogkstate s;
 	double v;
    alglib::autogkreport rep;
	r_for_potential=r;
	z_for_potential=z;

    alglib:: autogksmooth(0, 25*h, s);
	autogkintegrate(s,disk_potential_second_derivation);
    alglib:: autogkresults(s, v, rep);

	if(!is_valid(v))printf("\n integration %f: %i,\t%i,\t%i",r, rep.terminationtype, rep.nfev, rep.nintervals);

    return v;
}








void int_function_1_func(double r, double xminusa, double bminusx, double &y, void *ptr) 
{
    
  y = Mh*alpha*sqrt(PI)/2.0/(PI*PI)/rc*exp(-r*r/rc/rc)/(r*r + gam*gam)*PI*4.0*r*r; //y = x*x*exp(-x*x)/(x*x + q*q);
}


float Muh2(float r)
{

	
	float a1 = 0.0;
	float b1 = r;
	alglib::autogkstate s;
    double v;
	alglib::autogkreport rep;

	alglib::autogksmooth(a1, b1, s);
	alglib::autogkintegrate(s, int_function_1_func);
	alglib::autogkresults(s, v, rep);

	return v;
	
	/*int l6;	
	l6=(int)(r/R*30000.0);
	
	if(NforHalo==0)return 0.0;
	else return Bin6[l6]*Mh/NforHalo;
	*/
	
}

float Muh(float r)
{

	
	/*float a1 = 0.0;
	float b1 = r;
	alglib::autogkstate s;
    double v;
	alglib::autogkreport rep;

	alglib::autogksmooth(a1, b1, s);
	alglib::autogkintegrate(s, int_function_1_func);
	alglib::autogkresults(s, v, rep);

	return v;
	*/
	int l6;	
	l6=(int)(r/R*30000.0);
	
	if(NforHalo==0)return 0.0;
	else return Bin6[l6]*Mh/NforHalo;
	
	
}

float Mub(float r)
{
	/*int l5;	
	l5=(int)(r/R*6000.0);
	
	if(NforBulge==0)return 0.0;
	else return Bin5[l5]*Mb/NforBulge;
	*/
	return Mb*r*r /(r + a)/(r + a); 
}

float Mub2(float r)
{
	int l5;	
	l5=(int)(r/R*6000.0);
	
	if(NforBulge==0)return 0.0;
	else return Bin5[l5]*Mb/NforBulge;
	
	//return Mb*r*r /(r + a)/(r + a); 
}


float Qh(float r)
{

	return -G*Muh(r)/r + G*Mh*alpha/rc/sqrt(PI) * Ei(-r*r/(rc*rc) - q*q);
}


float dQh(float r)
{
	return 2.0*Mh*alpha*G/sqrt(PI) * ( r / (rc*rc*rc) /(r*r/(rc*rc) + q*q) * (exp(-r*r/(rc*rc) - q*q) - exp(-r*r/(rc*rc)))) + G*Muh(r)/(r*r);
}

float ddQh(float r)
{
	return 2.0*Mh*alpha*G*exp(-(r*r)/(rc*rc))/(rc*rc*rc)/(r*r/(rc*rc) + q*q)/sqrt(PI)*(2.0*r*r/(rc*rc) + 2.0*r*r/(rc*rc)/(r*r/(rc*rc) + q*q) + exp(-q*q) - 2.0*r*r/(rc*rc*exp(q*q)) - 2.0*r*r/(rc*rc*exp(q*q)*(r*r/(rc*rc) + q*q))) - 2.0*Muh(r)/(r*r*r);
}


float Qb(float r)
{
	//return -G*Mub(r)*(r + a)/(r*r);
	
	return -G*Mb/(r + a);
}

float dQb(float r)
{
	//return G*Mub(r)/(r*r);
	return G*Mb/(r + a)/(r + a);
}

float ddQb(float r)
{
	//return -2.0*Mub(r)/(r*r*(r + a));
	return -2.0*G*Mb/(r + a)/(r + a)/(r + a);
}

float Qtot(float r)
{
	return Qh(r) + Qb(r);
}

float dQtot(float r)
{
	return dQh(r) + dQb(r);
}

float ddQtot(float r)
{
	return ddQh(r) + ddQb(r);
}

//куммулятивная масса балджа
/*float Mub(float m)
{
	return Mb * m * m /(1 + m)/(1 + m);
}
float Rhob(float x,float y,float z)
{
	float	m = pow((x*x+y*y)/a/a + z*z/c/c,0.5f);
	return Mb/2.0/PI/a/c/m/(1 + m)/(1 + m)/(1 + m);
}
float Rhob_1(float x,float y,float z)
{
	float m = pow((x*x+y*y)/a/a + z*z/c/c,0.5f);
	return Mb/2.0/PI/a/c/m/(1 + m)/(1 + m)/(1 + m)*PI*4.0*pow((pow(x*y,1.6075f) + pow(x*z,1.6075f) + pow(y*z,1.6075f))/3.0,1.0/1.6075);
}*/


float Fun(float mu)
{
	return (mu+sqrt(mu*Mb))*a/(Mb-mu);
}


float Rhob(float r)
{
	return Mb*a/(2.0*PI*r)/((r + a)*(r + a)*(r + a))*PI*4.0*r*r;
}
float Rhob_1(float r)
{
	return Mb*a/(2.0*PI*r)/((r + a)*(r + a)*(r + a));
}

//распределеие гаусса


float F(float v,float sig )
{
	return  4*PI*pow(2.0*PI*sig*sig,-3.0/2.0)*v*v*exp(-v*v/2.0/sig/sig);
	
	//return 4.0*PI*pow(2.0*PI*sig*sig,-3.0/2.0)*(-1.0*v*exp(-v*v/2.0/sig/sig)*sig*sig + 1.25331413*sig*sig*sig*erf(v/sqrt(2.0)/sig)  ) ;
	//return (0.7978845606 * (-1.0 * v * exp(-0.5 * v * v /sig/sig) * sig*sig + 1.253314137 * sig * sig * sig * erf(0.7071067812 * v / sig))) / pow(sig,-3.0);
	//return  0.5*(1 + erf(v/sig));
}
float F2(float v,float sig )
{
	return 1.0/(sig*sqrt(2*PI))*exp(-v*v/(2*sig*sig));	
//return  0.5*(1 + erf(v/sqrt(2.0)/sig));
}


float veldistr( float Vesc, float sig)
{
 float v=random()*Vesc;
 float maxf=0.5870506255/sig;
 float fr=random()*maxf;
 while(fr>F(v,sig))
 {
	v=random()*Vesc;
	fr=random()*maxf;
 }
  return v;
}
 
float gauss(float sig)
{
	
	float v=(2.0*random() - 1.0)*3.0*sig;
	//if(v < 0.0)t=-1.0;
	float maxf=1.0/(sig*sqrt(2*PI));
	float fr=random()*maxf;
	while(fr>F2(v,sig))
	{
		v= random() *3.0*sig;
		fr=random() * maxf;
		
	}
	return v;
    
}




//функция плотности гало
float Rhoh(float r)
{
	return Mh*alpha*sqrt(PI)/2.0/PI/PI/rc*exp(-r*r/rc/rc)/(r*r + gam*gam);
}
float Rhoh_1(float r)
{
	return Mh*alpha*sqrt(PI)/2.0/(PI*PI)/rc*exp(-r*r/rc/rc)/(r*r + gam*gam)*PI*4.0*r*r;
}






//распределение плотности спутника (Яффе)
float Rhos(float r)
{
	return Ms / (4 * PI * r0 * r0 * r0) * (r0 / r) * (r0 / r) / (1 + r/r0 ) / (1 + r/r0 );
}

//куммулятивная масса спутника
float Mus(float r)
{
	return Ms*r / (r0 + r); 
}


void int_sig_h(double r, double xminusa, double bminusx, double &y, void *ptr)
{
	if(r<=R)y =  Rhoh(r)  * (-Q3(r))/r;
	else y =  Rhoh(r)  * (-Q3(R))/r;
}
void int_sig_b(double r, double xminusa, double bminusx, double &y, void *ptr)
{
	if(r<=R)y = Rhob(r)  * (-Q5(r))/r;
	else y = Rhob(r)  * (-Q5(R))/r;
	
}

float sigma2_h(float r,int k)
{	
	
	alglib::autogkstate s;
 	double v;
    alglib::autogkreport rep;
	
	

    alglib:: autogksmooth(r, 50.0*h, s);
	autogkintegrate(s,int_sig_h);
    alglib:: autogkresults(s, v, rep);

  
	//printf("\n\n%f\n",1.0/Rhoh(r) * v);

	return 1.0/Rhoh(r) * v ;
	
	
}
float sigma2_b(float r,int k)
{	
	
	alglib::autogkstate s;
 	double v;
    alglib::autogkreport rep;
	
	
    alglib:: autogksmooth(r, 50*h, s);
	autogkintegrate(s,int_sig_b);
    alglib:: autogkresults(s, v, rep);

  
//printf("\n\n%f\n",1.0/Rhoh(r) * v);

	return 1.0/Rhob(r) * v ;
	
	
	
}

//sigma^2
/*
float sigma2(float r,int p)
{	
	float sum=0;
	if(p>R)p=R;

	//if(r < 4.0*h/3.0 && p==1) sum =( Rhoh(r) * G * fabsf(Q3(r))/r + Rhoh(R) * G * fabsf(Q3(R))/R)/2.0; 
	if(p==1) sum =( Rhoh(r)  * (-Q3(r))/r + Rhoh(R) *  (-Q3(R))/R)/2.0;
	//if(p==1) sum =( Rhoh(r)  * (dQ(sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z))/r + Rhoh(R) *  (dQ(R,X[k].z))/2.0;
    if(p==2) sum =( Rhob(r)  * (-Q(r))/r + Rhob(R) *  (-Q(R))/R)/2.0;
	if(p==3) sum = ( Rhos(r)  * (-Qs(r))/r + Rhos(Rs) *  (-Qs(Rs))/Rs)/2.0;
	
	int n = 100;
	
	for(int i=1; i<n;i++)
	{
		//if(r < 4.0*h/3.0 && p==1) sum += Rhoh(i*(R-r)/n + r) * G * fabsf(Q3(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==1) sum += Rhoh(i*(R-r)/n + r)  * (-Q3(i*(R-r)/n + r))/(i*(R-r)/n + r);
		//if(p==1) sum += Rhoh(i*(R-r)/n + r)  * dQ(i*(R-sqrt(X[k].x*X[k].x+X[k].y*X[k].y))/n + sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z);
		if(p==2) sum += Rhob(i*(R-r)/n + r) * (-Q(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==3) sum += Rhos(i*(Rs-r)/n + r) * (-Qs(i*(Rs-r)/n + r))/(i*(Rs-r)/n + r);
	}
		//printf("\n\n%f\n",1.0/Rhoh(r) * sum * (R-r)/n);
	if(p == 1) return 1.0/Rhoh(r) * sum * (R-r)/n;
	if(p == 2)  return 1.0/Rhob(r) * sum * (R-r)/n;
	if(p == 3) return 1.0/Rhos(r) * sum * (Rs-r)/n;
	
}
*/
float sigma2(float r,int p)
{	
	float sum=0;
	

	//if(r < 4.0*h/3.0 && p==1) sum =( Rhoh(r) * G * fabsf(Q3(r))/r + Rhoh(R) * G * fabsf(Q3(R))/R)/2.0; 
	if(p==1) sum =( Rhoh(r)  * (-Q3(r))/r + Rhoh(10.0*R) *  (-Q3(R))/10.0/R)/2.0;
	
	if(p==2) sum =( Rhob_1(r)  * (-Q3(r))/r + Rhob_1(10.0*R) *  (-Q3(R))/10.0/R)/2.0;
	if(p==3) sum = ( Rhos(r)  * (-Qs(r))/r + Rhos(Rs) *  (-Qs(Rs))/Rs)/2.0;
	
	int n;
	if(p==1)n = 10000;
	if(p==2)n = 10000;
	
	for(int i=1; i<n;i++)
	{
		//if(r < 4.0*h/3.0 && p==1) sum += Rhoh(i*(R-r)/n + r) * G * fabsf(Q3(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==1)
		{	
		
			if(i*(10*R-r)/n + r > R) sum += Rhoh(i*(10.0*R-r)/n + r)  * (-Q3(R)/(i*(10.0*R-r)/n + r));
			else { sum += Rhoh(i*(10.0*R-r)/n + r)  * (-Q3(i*(10.0*R-r)/n + r))/(i*(10.0*R-r)/n + r);}
		}
		if(p==2)
		{	
			
			if(i*(10*R-r)/n + r > R){sum += Rhob_1(i*(10.0*R-r)/n + r)  * (-Q3(R)/(i*(10.0*R-r)/n + r));}
			else { sum += Rhob_1(i*(10.0*R-r)/n + r)  * (-Q3(i*(10.0*R-r)/n + r))/(i*(10.0*R-r)/n + r);}
			
			 //sum += Rhob(i*(10.0*R-r)/n + r)  * (dQb(i*(10.0*R-r)/n + r));
		}
		if(p==3) sum += Rhos(i*(Rs-r)/n + r) * (-Qs(i*(Rs-r)/n + r))/(i*(Rs-r)/n + r);
	}
		//printf("\n\n%f\n",1.0/Rhoh(r) * sum * (R-r)/n);
	if(p == 1) {
					//printf("sigma= %f\n ",1.0/Rhoh(r) * sum * (10*R-r)/n); 
					return 1.0/Rhoh(r) * sum * (10.0*R-r)/n;
				}
	else{
		if(p == 2)  return 1.0/Rhob_1(r) * sum * (10.0*R-r)/n;
		else{
			if(p == 3) return 1.0/Rhos(r) * sum * (Rs-r)/n;
			else {
					std::invalid_argument("Wrong variant");
					}
				}
			}
	
	
}

float sigma2(float r,int p, int k)
{	
	float sum=0;
	

	//if(r < 4.0*h/3.0 && p==1) sum =( Rhoh(r) * G * fabsf(Q3(r))/r + Rhoh(R) * G * fabsf(Q3(R))/R)/2.0; 
	//if(p==1) sum =( Rhoh(r)  * (-Q3(r))/r + Rhoh(R) *  (-Q3(R))/R)/2.0;
	if(p==1) sum =( Rhoh(r)  * (dQ(sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z))/r + Rhoh(R) *  (dQ(R,X[k].z)))/2.0 + ( Rhoh(r)  * dQtot(r) + Rhoh(R) *  dQtot(R))/2.0;
    if(p==2) sum =( Rhob(r)  * (dQ(sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z))/r + Rhob(R) *  (dQ(R,X[k].z)))/2.0 + ( Rhoh(r)  * dQtot(r) + Rhoh(R) *  dQtot(R))/2.0;
	if(p==3) sum = ( Rhos(r)  * (-Qs(r))/r + Rhos(Rs) *  (-Qs(Rs))/Rs)/2.0;
	
	int n = 60;
	for(int i=1; i<n;i++)
	{
		//if(r < 4.0*h/3.0 && p==1) sum += Rhoh(i*(R-r)/n + r) * G * fabsf(Q3(i*(R-r)/n + r))/(i*(R-r)/n + r);
		//if(p==1) sum += Rhoh(i*(R-r)/n + r)  * (-Q3(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==1) sum += Rhoh(i*(R-r)/n + r)  * dQ(i*(R-sqrt(X[k].x*X[k].x+X[k].y*X[k].y))/n + sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z) + Rhoh(i*(R-r)/n + r)  * (dQtot(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==2)  sum += Rhob(i*(R-r)/n + r)  * dQ(i*(R-sqrt(X[k].x*X[k].x+X[k].y*X[k].y))/n + sqrt(X[k].x*X[k].x+X[k].y*X[k].y),X[k].z) + Rhob(i*(R-r)/n + r)  * (dQtot(i*(R-r)/n + r))/(i*(R-r)/n + r);
		if(p==3) sum += Rhos(i*(Rs-r)/n + r) * (-Qs(i*(Rs-r)/n + r))/(i*(Rs-r)/n + r);
	}

			
	if(p == 1) return 1.0/Rhoh(r) * sum * (R-r)/n;
	else{ 
		if(p == 2)  return 1.0/Rhob(r) * sum * (R-r)/n;
			else{ 
				if(p == 3) return 1.0/Rhos(r) * sum * (Rs-r)/n;
				else std::invalid_argument("Wrong variant");
			}
		}
	
}






// начальное распределение частиц диска, гало, балджа
void    InitParticles()
{	
	char FileName4[32];
	sprintf(FileName4,"Vr.txt");
	FILE * out4 = fopen(FileName4, "w+");
	srand(time(NULL));
	X = new float4[NParticles];
	A = new float4[NParticles];
	V = new float4[NParticles];
	
	for(int i=0; i<10000 ;i++)
	{
		Bin[i] = 0;
	}
	for(int i=0; i<1000 ;i++)
	{
		Bin2[i] = 0;
	}
	for(int i=0; i<10000 ;i++)
	{
		Bin3[i] = 0;
	}
	for(int i=0; i<1000 ;i++)
	{
		Bin4[i] = 0;
	}
	for(int i=0; i<100000 ;i++)
	{
		Bin5[i] = 0;
	}
	for(int i=0; i<30000 ;i++)
	{
		Bin6[i] = 0;
	}
	for(int i=0; i<128 ;i++)
	{
		BinD[i] = 0;
	}
	for(int i=0; i<128 ;i++)
	{
		BinDz[i] = 0;
	}


	float Mumh =Mh/5.8;//0.6401888889/5.8*Mh;
	float Mumb = 0;
	float Mums = Mus(Rs);
	float r = 0;
	float phi = 0;
	float teta = 0;
	float x = 0;
	float y = 0;
	float z = 0;
	float m = 0;

	printf("alpha %f\n", alpha);
    
	float Mur = 0;
	float Murz = 0;
    //задаем координаты частиц методом отбора-отказа(фон Неймана)
	char FileName7[32];
	sprintf(FileName7,"VrH.txt");
	FILE * out7 = fopen(FileName7, "w+");
	//printf("%f\n",Rhob(0.099999998529215414));

    for (int k=0; k<NParticles; k++)
    {	
		//printf("%i\t%i\t%i",k, NforDisc, NParticles);
		if(k<NforDisc)
		{	
			float Mum2 = 0.367*Md;
			//float Mum2 = 3.0/8.0*Md;
			phi = random()*2.0*PI;
			r = random()*R*h;
			z = (random()*2.0 - 1.0) * R;

			
		/*	Mur = random()* Mum2;
		
			while( Mur > Rhod(r)*Sech(z) )
			{
				//phi = random()*2*PI;		
				r = random()*R;
				z = (random()*2.0 - 1.0) * R;
				Mur = random() * Mum2;
			
			}
		*/
			
			Mur = random()* Mum2;
		
			while( Mur > Rhod(r) )
			{
					
				
				r = random()*R*h;
				
				Mur = random() * Mum2;
			//	printf("rho_%i %f\t%f\t%f\n",k,r, Mur,Rhod(r));
			
			}
			
			float Mumd=Sech(0.0);
			Mur = random()* Mumd;
			while( Mur > Sech(z) )
			{
				z = (random()*2.0 - 1.0) * R;
				Mur = random() * Mumd;
				//printf("rho_%i %f\t%f\t%f\n",k,r, Mur,Rhod(r));
			}
			
			
			X[k].x = r * cos(phi); 
			A[k].x = 0.0;
	
			X[k].y = r * sin(phi);
			A[k].y = 0.0;

			X[k].z = z;
			A[k].z = 0.0;

			X[k].w = Md * 1.0/NforDisc;
			A[k].w = 0.0;
			
			
			//printf("%f\t%f\t%f\t%f", Mumd, Mur, r, z);
			
			int l = r/R*10000;



			int l2=r*0.35/R*10000;
			
			//Bin3[l2] += X[k].w;
			//if(r==R)l=999;
			if(10000-l<=3)
			{
				for(int i=-3;i<10000-l;i++)
				{	
				int l=(r/R*10000  + i);
				
				 Bin[l] +=X[k].w/(3.0 + 10000.0 - l);

				//if(r>0.35*h*4.0)Bin3[l] += X[k].w/(3.0 + 10000.0 - l);
				Bin3[l2] += X[k].w/(3.0 + 10000.0 - l);
				
				}
			}
			else if(l < 3)
			{
			
				for(int i=-l;i<=3;i++)
				{	
					int l=(r/R*10000  + i);
					Bin[l] += X[k].w/(l + 1.0 + 3.0 );
					//if(r>0.35*h*4.0) Bin3[l] += X[k].w/(l + 1.0 + 3.0 );
					Bin3[l2] += X[k].w/(l + 1.0 + 3.0 );
					
				}
			}
			else
			{
				for(int i=-3;i<=3;i++)
				{	
					int l=(r/R*10000  + i);
					Bin[l] +=  X[k].w/7.0;

					//if(r>0.35*h*4.0)Bin3[l] += X[k].w/7.0;
					Bin3[l2] += X[k].w/7.0;
				}
			}
			
			//int l3=fabsf(r-rabs)/Rs*1000;
			//if(l3==1000)l3=999;
			//Bin5[l]+=X[k].w;
			//Bin[l] += X[k].w;
			//if(l3<1000)Bin4[l3] += X[k].w;
			
		}
		if(k >= NforDisc && k < (NforHalo + NforDisc)&& NforHalo != 0)
		{
			r = random()*R;
			Mur = random()* Mumh;
			while(Mur > Rhoh_1(r))
			{
			  r = random()*R;
			  Mur = random()*Mumh;

			}
			
			 //r=r/2.0;

			phi = phi = random()*2.0*PI;;
			teta = asin(random()*2.0 - 1.0);
			X[k].x = r * cos(phi)*cos(teta); 
			A[k].x = 0.0;

			X[k].y = r * sin(phi)*cos(teta);
			A[k].y = 0.0;

			X[k].z = r * sin(teta);
			A[k].z = 0.0;

			X[k].w = Mh * 1.0/NforHalo;
			A[k].w = 0.0;
			
			int l=r/R*10000;
			int lh=r/25.0*128;

			int l6;
			l6=(int)(r/R*30000.0);
			//printf("halo: %f\t%i",r/R*30000.0,l6);
			Bin6[l6] +=1.0;


			if(r==R)l=999;
			//printf("Mumb %d\t%f\t%f\n",l, r , r/R);
			//int l3=fabsf(r-rabs)/Rs*1000;
			//if(l3==1000)l3=999;
			
			if(10000-l<=8)
			{
				for(int i=-8;i<10000-l;i++)
				{	
				int l=(r/R*10000  + i);
			

				Bin[l] += X[k].w/(8.0 + 10000.0 - l);
				Bin2[l] += X[k].w/(8.0 + 10000.0 - l);
				Bin3[l] += X[k].w/(8.0 + 10000.0 - l);
			

				}
			}
			else if(l<8)
			{
			
				for(int i=-l;i<=8;i++)
				{	
					int l=(r/R*10000  + i);
				

					Bin[l] += X[k].w/(l + 1.0 + 8.0 );
					Bin2[l] += X[k].w/(l + 1.0 + 8.0 );
					Bin3[l] += X[k].w/(l + 1.0 + 8.0 );
				
				}
			}
			else
			{
				for(int i=-8;i<=8;i++)
				{	
					int l=(r/R*10000  + i);
			

					Bin[l] += X[k].w/17.0;
					Bin2[l] += X[k].w/17.0;
					Bin3[l] += X[k].w/17.0;
			
				}
			}


			/*
				if(6000-l6<=8)
			{
				for(int i=-8;i<6000-l6;i++)
				{	
				int l5 = (r/R*6000  + i);

				Bin5[l5] += X[k].w/(8.0 + 6000.0 - l6);
				}
			}
			else if(l6<8)
			{
			
				for(int i=-l6;i<=8;i++)
				{	
					
					int l6 = (r/R*6000  + i);

					Bin6[l6] += X[k].w/(l6 + 1.0 + 8.0 );
				}
			}
			else
			{
				for(int i=-8;i<=8;i++)
				{	
					
					int l6 = (r/R*6000  + i);

				
					Bin6[l6] +=  X[k].w/17.0;
				}
			}*/
			
			//Bin[l] += X[k].w;
			//Bin2[l] += X[k].w;
			//Bin3[l] += X[k].w; 

		//	if(l!=0)Bin6[lh] += X[k].w/(4.0/3.0*PI*pow((lh+1)*25.0/128.0,3.0)-4.0/3.0*PI*pow((lh)*25.0/128.0,3.0)); 
		//	else Bin6[lh] += X[k].w/(4.0/3.0*PI*pow(25.0/128,3.0)); 
			//Bin6[l]+= X[k].w;
			
		}
		if(k >= (NforHalo + NforDisc) && k < (NforDisc + NforHalo + NforBulge) && NforBulge!=0)
		{
			/*Mumb = Mub(R/0.1);
			r = random()*R;
			phi =  random()*2.0*PI;
			teta = asin(random()*2.0 -1.0);
			float x = r * cos(phi)*cos(teta);
			float y = r * sin(phi)*cos(teta);
			float z = r * sin(teta);
			m = pow((x*x+y*y)/a/a + z*z/c/c,0.5f);
			Mur = random()* Mumb;
			
			while(Mur < Mub(m))
			{
				r = random()*R;
				phi =  random()*2.0*PI;
				teta = asin(random()*2.0 -1.0);
				x = r * cos(phi)*cos(teta);
				y = r * sin(phi)*cos(teta);
				z = r * sin(teta);
				m = pow((x*x+y*y)/a/a + z*z/c/c,0.5f);
				Mur = random()* Mumb;
			}*/
			Mumb = 0.3/a*Mb;
			
			r = random()*R;
			Mur = random()* Mumb;
			while(Mur > Rhob(r))
			{
			  r = random()*R;
			  Mur = random()*Mumb;

			}
			
		
			

		//	r=r/2.0;
			phi = phi = random()*2.0*PI;;
			teta = asin(random()*2.0 - 1.0);

		
			X[k].x = r * cos(phi)*cos(teta); 
			A[k].x = 0.0;
			
			
        
			X[k].y = r * sin(phi)*cos(teta);
			A[k].y = 0.0;
           
			
			
			X[k].z = r * sin(teta);
			A[k].z = 0.0;
			
			
			
			
			X[k].w = Mb * 1.0/NforBulge;
			A[k].w = 0.0;
			
			int l=r/R*10000;
			int l5 = r/R*100000;

			Bin5[l5] += X[k].w;

			
			if(10000-l<=2)
			{
				for(int i=-2;i<10000-l;i++)
				{	
				int l=(r/R*10000  + i);
			

				Bin[l] += X[k].w/(2.0 + 10000.0 - l);
				Bin2[l] += X[k].w/(2.0 + 10000.0 - l);
				Bin3[l] += X[k].w/(2.0 + 10000.0 - l);
			
				}
			}
			else if(l<2)
			{
			
				for(int i=-l;i<=2;i++)
				{	
					int l=(r/R*10000  + i);
			

					Bin[l] += X[k].w/(l + 1.0 + 2.0 );
					Bin2[l] += X[k].w/(l + 1.0 + 2.0 );
					Bin3[l] += X[k].w/(l + 1.0 + 2.0 );
				
				}
			}
			else
			{
				for(int i=-2;i<=2;i++)
				{	
					int l=(r/R*10000 + i);
					
					Bin[l] += X[k].w/5.0;
					Bin2[l] += X[k].w/5.0;
					Bin3[l] += X[k].w/5.0;
			
				}
			}

			/*
			if(6000-l5<=2)
			{
				for(int i=-2;i<6000-l5;i++)
				{	
				int l5 = (r/R*6000  + i);

				Bin5[l5] += X[k].w/(2.0 + 6000.0 - l5);
				}
			}
			else if(l5<2)
			{
			
				for(int i=-l5;i<=2;i++)
				{	
					
					int l5 = (r/R*6000  + i);

					Bin5[l5] += X[k].w/(l5 + 1.0 + 2.0 );
				}
			}
			else
			{
				for(int i=-2;i<=2;i++)
				{	
					
					int l5 = (r/R*6000  + i);

				
					Bin5[l5] +=  X[k].w/5.0;
				}
			}*/
			
			//Bin[l] += X[k].w;
			//Bin2[l] += X[k].w;
			//if(l3<1000)Bin4[l3] += X[k].w;
			
			//Bin3[l] +=  X[k].w;
			//Bin5[l] +=  X[k].w;
		}
		/*else if(0)
		{
			
			r = random()*Rs;
			phi = random()*2.0*PI;
			teta = asin(random()*2.0 -1.0);
			x = xs + r * cos(phi)*cos(teta);
			y = ys + r * sin(phi)*cos(teta);
			z = zs + r * sin(teta);
			
			Mur = random()* Mums;
			while(Mur <= Mus(r))
			{
				r = random()*Rs;
				
				x = xs + r * cos(phi)*cos(teta);
				y = ys + r * sin(phi)*cos(teta);
				z = zs + r * sin(teta);
				Mur = random()* Mums;
			}
			//printf("Mumb %f\t%f\t%f\t%f\n", r , Mums,Mus(r),Mur);
			
			

		
			X[k].x = x;
			A[k].x = 0.0;
			
			
        
			X[k].y = y;
			A[k].y = 0.0;
           
			
			
			X[k].z = z;
			A[k].z = 0.0;
			
			
			
			
			X[k].w = Ms * 1.0/NforSatellite;
			A[k].w = 0.0;
			
			int l=r/R*1000;
			

			Bin4[l] +=  X[k].w;
		}*/
	//	fprintf(out7,"%f\t%f\t%f\t\n", V[k].x, V[k].y, V[k].z  );
		
            
	}
	fclose(out7);
	
   

   	system("pause");
	float ref = 0.0;
	int numref = 0;

	/*for(int k=0; k<NforDisc; k++)
	{
		float r = sqrt(X[k].x*X[k].x + X[k].y*X[k].y);
		if(sqrt((sqrt(3.0) - X[k].x)*(sqrt(3.0) - X[k].x) + (sqrt(3.0) - X[k].y)*(sqrt(3.0) - X[k].y)) < h/5.0  )
			{	
				float t = sqrt(X[k].x*X[k].x + X[k].y*X[k].y + X[k].z*X[k].z);
				//float O = sqrt(dQ(r,0)/r);//sqrt((dQ(r,0) + G*(Mub(m)+Muh(t))/r/r)/r);
				float K =sqrt(3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r)));
				//float K =sqrt(3.0*(dQ(r,0) - Q2(r)/r)/r + (ddQ(r,0) + 2.0*Q2(r)/r/r));
				//float K =sqrt((3.0*(dQ(r,0) + G*(Mub(m)+Muh(t))/r/r)/r) + (ddQ(r,0)-2.0*G*(Mub(m)+Muh(t))/r/r/r));
				ref += 3.36 * G * Md/(2.0*PI*h*h)* exp(-r/h)/K;
					
					numref++;
			}
	}*/

	for(int k=0; k<NforDisc; k++)
	{
		float r = sqrt(X[k].x*X[k].x + X[k].y*X[k].y);
		if(2.43 - r < h/15.0 && 2.43 - r> -h/15.0 )
			{	
				//float t = sqrt(X[k].x*X[k].x + X[k].y*X[k].y + X[k].z*X[k].z);
				//float O = sqrt(dQ(r,0)/r);//sqrt((dQ(r,0) + G*(Mub(m)+Muh(t))/r/r)/r);
				//float K =sqrt(3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r)));
				float K =sqrt(3.0*(dQ(r,0) - Q2(r)/r)/r + (ddQ(r,0) + 2.0*Q2(r)/r/r));
				//float K =sqrt((3.0*(dQ(r,0) + G*(Mub(m)+Muh(t))/r/r)/r) + (ddQ(r,0)-2.0*G*(Mub(m)+Muh(t))/r/r/r));
				ref += 3.36 * G * Md/(2.0*PI*h*h)* exp(-r/h)/K;
					
					numref++;
			}
	}


	//ref *= 3.36 * G * (float)numref*Md/NforDisc/(PI*(2.43 + h/30.0)*(2.43 + h/30.0) - PI*(2.43 - h/30.0)*(2.43 - h/30.0));
	float Scrit = ref/(float)numref;


	float c0 =Qt*Scrit/exp(-2.43/(2.0*h));
	

    //printf("halo %f\n%i\n%f\n%f\n",c0,numref, Md/(2.0*PI*h*h) *exp(-2.45/h),(float)numref*Md/NforDisc/(PI*(2.45 + h/30.0)*(2.45 + h/30.0) - PI*(2.45 - h/30.0)*(2.45 - h/30.0)) );
	//system("pause");

	//задаем начальые скорости частиц через БУБ
	float vteta = 0;
	float vphi = 0;
	float vr = 0;
	

	for(int i = 1; i<10000; i++)
			{
				
				Bin[i]+=Bin[i-1];
			}
	for(int i = 1; i<10000; i++)
			{
				
				Bin3[i]+=Bin3[i-1];
			}
	for(int i = 1; i<100000; i++)
			{
				
				Bin5[i]+=Bin5[i-1];
				
			}
		for(int i = 1; i<30000; i++)
			{
				
				
				Bin6[i]+=Bin6[i-1];
			}

		printf("sigma %f\n",  sigma2_b(0.005784,2)- sigma2(0.005784,2));
	for (int k=0; k<NParticles; k++)
	{
		//printf("%i\n",k);
		//диск
		if(k<NforDisc)
		{	
			
			r = sqrt(X[k].x * X[k].x + X[k].y * X[k].y);
			float t = sqrt(X[k].x * X[k].x + X[k].y * X[k].y + X[k].z*X[k].z);
				if(r == 0.0)
		{

			V[k].x = 0.0;
			V[k].y = 0.0;
			V[k].z = 0.0;
			V[k].w = 0.0;
			continue;
		}
			
			float pot=0.0;
			float surf=0.0;
		/*	for (int j=0; j<NParticles; j++)
		{
			float rad = sqrt((X[j].x - X[k].x)*(X[j].x - X[k].x) + (X[j].y - X[k].y)*(X[j].y - X[k].y) + (X[j].z - X[k].z)*(X[j].z - X[k].z));
			if(rad==0.0)continue;
			pot+=-X[j].w/rad;
			if(rad<h/50)surf+=1.0/NforDisc;
			
		}
			surf=surf/(PI*h*h/50.0/50.0);*/
			/*float K = sqrt(3.0/r*(dQ(r,0) - Q2(r)/r) + (ddQ(r,0)+2.0*Q2(r)/r/r));//sqrt(fabs(3.0*(dQ(r,0) - Q2(r)/r)/r + (ddQ(r,0)+2.0*Q2(r)/r/r)));//sqrt(fabs(3.0*( - Q(r)/r)/r + 2.0*Q(r)/r/r));
			float O =sqrt(1/r*(dQ(r,0) - Q2(r)/r));//sqrt(1/r*(dQ(r,0) - Q2(r)/r));//sqrt(Vc2(r))/r;//sqrt(dQ(r,0)/r - Q2(r)/r/r);//sqrt(-Q(r))/r;//sqrt(dQ(r,0)/r - Q2(r)/r/r);
			float VR2 = c0*c0*exp(-sqrt(r*r +h*h/8.0)/h);//c0*c0*exp(-sqrt(r*r + h*h/8.0)/h); //pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/K * Qt, 2.0);//c0*c0*exp(-sqrt(r*r + h*h/8.0)/h);//c0*c0*exp(-sqrt(r*r + h*h/8.0)/h);//3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/K*Qt;
			float Vphi2=0.0f;
			//if(NforHalo + NforBulge==0) Vphi2 = pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/2.0/O , 2.0);//VR2*K*K/O/O/4;//pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/2.0/O , 2.0);
			  Vphi2 =VR2*K*K/O/O/4.0;
			
			float Vz2 = PI*G*rho0*2.0*z0*exp(-sqrt(r*r + h*h/8.0 )/h)*z0;//PI*G*rho0*2.0*z0*exp(-sqrt(r*r)/h)*z0;//PI*G*rho0*2.0*z0*exp(-sqrt(r*r + h*h/8.0)/h)*z0;
			float Vphis=0.0;
			//if(dQ(r,z)*r - Q2(r) - VR2 * (1 - K * K /(4.0 * O * O) - 2*sqrt(X[k].x * X[k].x + X[k].y * X[k].y)/h)<0.0)Vphis= 0.0f;
			//else  
			Vphis = sqrt(dQ(r,0)*r - Q2(r) - VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h) );//sqrt(fabs(Vc2(r) - VR2 * (1 - K * K /(4.0 * O * O) - 2*sqrt(X[k].x * X[k].x + X[k].y * X[k].y)/h)) );//sqrt(fabs((dQ(r,0) + G*(Mub(m)+Muh(t))/r/r)*r - VR2 * (1 - K * K /(4.0 * O * O) - 2*sqrt(X[k].x * X[k].x + X[k].y * X[k].y)/h)) );
			float Vz = gauss(sqrt(Vz2));
			float VR = gauss(sqrt(VR2));
			float Vphi = gauss(sqrt(Vphi2));
			*/
			
			//float K =sqrt(fabs(3.0*(dQ2(r,0) + dQtot(r))/r + ddQ2(r,0)+ ddQtot(r)));//sqrt(fabs(3.0*( - Q(r)/r)/r + 2.0*Q(r)/r/r));//sqrt(3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r)));//sqrt(fabs(3.0*(dQ(r,0) - Q2(r)/r)/r + (ddQ(r,0)+2.0*Q2(r)/r/r)));//sqrt(fabs(3.0*( - Q(r)/r)/r + 2.0*Q(r)/r/r));
			//if(3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r))<0.0)K=0.0;
			
			float O = sqrt(dQ(r,0)*r - Q2(r))/r;//sqrt(dQ2(r,0)*r)/r; //sqrt(dQ2(r,0)*r + dQtot(r))/r;//PI*sqrt(G*3.0*Md/(32.0*PI*h*h)/(8.0));//sqrt(1/r*(dQ(r,0) + dQtot(r)));//sqrt(1/r*(dQ(r,0) - Q2(r)/r));//sqrt(Vc2(r))/r;//sqrt(dQ(r,0)/r - Q2(r)/r/r);//sqrt(-Q(r))/r;//sqrt(dQ(r,0)/r - Q2(r)/r/r);
			float K = sqrt(3.0*(dQ(r,0) - Q2(r)/r)/r + (ddQ(r,0) + 2.0*Q2(r)/r/r));
			//float K =sqrt(2.0)*O;
									
			float VR2 = c0*c0*exp(-sqrt(r*r  + h*h/8)/h);//pow(1.5f*3.36 * G * Md/(2.0*PI*h*h) * exp(-r/h)/K,2.0);//pow(3.36 * G * surf/K,2.0);//(0.1*O*4.0)*(0.1*O*4.0);////c0*c0*exp(-sqrt(r*r + h*h/8.0)/h); //pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/K * Qt, 2.0);//c0*c0*exp(-sqrt(r*r + h*h/8.0)/h);//c0*c0*exp(-sqrt(r*r + h*h/8.0)/h);//3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/K*Qt;
			
			float Vphi2 = 0.0f;
			//if(NforHalo + NforBulge==0) Vphi2 = pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/2.0/O , 2.0);//VR2*K*K/O/O/4;//pow(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h)/2.0/O , 2.0);
			  //Vphi2 =VR2*K*K/O/O/4.0;
			Vphi2 =VR2*K*K/O/O/4;//pow(3.36 * G * Md/(2.0*PI*h*h)* exp(-r/h)/2.0/O,2.0);//pow(3.36 * G *surf/2.0/O,2.0);//(0.1*O*4.0)*(0.1*O*4.0);

			float Vz2 = PI*G*rho0*2.0*z0*exp(-sqrt(r*r + h*h/8)/h)*z0;//0.25*Vphi2;// PI*G*rho0*2.0*z0*exp(-sqrt(r*r + h/8.0  )/h)*z0;//PI*G*rho0*2.0*z0*exp(-sqrt(r*r)/h)*z0;//PI*G*rho0*2.0*z0*exp(-sqrt(r*r + h*h/8.0)/h)*z0;
			float Vphis=0.0;
			if(O*r*O*r + VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h)<0.0)Vphis= 0.0f;
			else  Vphis =sqrt(O*r*O*r + VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h));
			//sqrt((dQ(r,0) + dQtot(r))*r + VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h) );
			//Vphis = sqrt(O*r*O*r);//sqrt(fabs(dQ2(r,0)*r));
			//if(O*r*O*r + dQtot(r)*r+ VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h)<0.0)Vphis=0.0f;
			float Vz = gauss(sqrt(Vz2));
			float VR = gauss(sqrt(VR2));
			float Vphi = gauss(sqrt(Vphi2));

			//printf("while %f %f %f\n",Vphis,VR,Vphi);



		//	if(O*r*O*r + VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h)<0.0)printf("test %f\t%f\n",r,O*r*O*r + VR2 * (1 - K * K /(4.0 * O * O) - 2*r/h));
			//if(3.0*(dQ2(r,0) + dQtot(r))/r + (ddQ2(r,0) + ddQtot(r))<0.0)printf("K %f\t%f\n",r, 3.0*(dQ2(r,0) + dQtot(r))/r + (ddQ2(r,0) + ddQtot(r)));
			//if(_isnan(Vphis))printf("O %f\t%f\n",r,dQ2(r,0)*r);
			
			

			x = -(Vphis + Vphi)*X[k].y/r - VR*X[k].x/r ;
			y = (Vphis + Vphi)*X[k].x/r - VR*X[k].y/r ;
			z = Vz;
			//float vm = sqrtf((Vphis + Vphi) *(Vphis + Vphi)  + Vz * Vz + VR * VR);
			float v = sqrt(x * x + y * y + z * z); 
			//if(_isnan(v))printf("while%i %f \n",k,v);
		
		
			//printf("test %f\t%f\t%f\n",r ,v,Q2(r));
			//printf("test %f\t%f\n",r ,dQ(14.5,0));
			//sqrt(VR2)*K/(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h))
			if (!is_valid(x) && !is_valid(y))
			{
				printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", r,O, K,dQ(r,0),dQtot(r),ddQ(r,0),ddQtot(r),3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r)));}
			fprintf(out4,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", r,sqrt(VR2), sqrt(Vz2),sqrt(Vphi2),O,K,sqrt(VR2)*K/(3.36 * G * Md/(2.0*PI*h*h)* exp(-r/h)),Qd(r,0));
			V[k].x = x ;
			V[k].y = y ;
			V[k].z = z;
			V[k].w = 0.0;
			
				if(k== NforDisc-1 )
		{ 
			fclose(out4); 
			//system("pause");
		}
		}
	
		//гало
	 if(k >= NforDisc && k < (NforHalo + NforDisc)&& NforHalo != 0 )
		{
			
			r = sqrtf(X[k].x * X[k].x + X[k].y * X[k].y + X[k].z * X[k].z);
		//	printf("halo %f\t%f\t%f\t%f\n", X[k].x, X[k].y, X[k].z, r );
					if(fabs(r) == 0.0)
		{

			V[k].x = 0.0;
			V[k].y = 0.0;
			V[k].z = 0.0;
			V[k].w = 0.0;
			continue;
		}
			
			float VR2 = sigma2(r,1);
			
			float Vesc = sqrtf(fabs(2.0*Phi(r))); 
			
			vr = veldistr(4.0*sqrtf(VR2),sqrtf(VR2));
		
		
			while(vr > 0.95 * Vesc)
			{
								
				vr = veldistr(4.0*sqrtf(VR2),sqrtf(VR2));								
			
			}
	
			phi  = random()*2.0*PI;
			teta = asin(random()*2.0 - 1.0);
			
			//vr=0.0;
			printf("halo %f\t%f\t%f\t%f\n",r, vr, Phi(r), Q3(r));
			x = vr * cos(phi)*cos(teta); 
			y = vr * sin(phi)*cos(teta);
			z = vr * sin(teta);
		
			V[k].x = x;
			V[k].y = y;
			V[k].z = z;
			V[k].w = 0;
			
		}
		//балдж
		if(k >= NforDisc + NforHalo && k <( NforDisc + NforHalo + NforBulge) && NforBulge != 0)
		{r = sqrtf(X[k].x * X[k].x + X[k].y * X[k].y + X[k].z * X[k].z);
		 
		
		
			if(r == 0.0)
				{
					V[k].x = 0.0;
					V[k].y = 0.0;
					V[k].z = 0.0;
					V[k].w = 0.0;
					continue;
				}
									
			float VR2 = sigma2(r,2);			
					
		
				
			float Vesc = sqrtf(fabs(2.0*Q5(r)));
			vr = veldistr(4.0*sqrtf(VR2),sqrtf(VR2));
			
			while(vr > 0.95 * Vesc)
			{
				
				vr = veldistr(4.0*sqrtf(VR2),sqrtf(VR2));
				
			}
			printf("Vr = %f\t%f\n",r,VR2);
			phi  = random()*2.0*PI;
			teta = asin(random()*2.0 - 1.0);
			//vr=0.0;
			x = vr * cos(phi)*cos(teta); 
			y = vr * sin(phi)*cos(teta);
			z = vr * sin(teta);
		
			V[k].x = x;
			V[k].y = y;
			V[k].z = z;
			V[k].w = 0;
			
		}
		/*else if(NforSatellite!=0)
		{
			r = sqrtf((X[k].x-xs) * (X[k].x-xs) + (X[k].y-ys) * (X[k].y-ys) + (X[k].z - zs) * (X[k].z - zs));
			//float rab = sqrt(xs*xs + ys*ys + zs*zs); 
			if(r == 0.0)
		{
			V[k].x = 0.0 + sqrtf(-Q(rabs)) * cosf(ksi) * (-ys/rabs);
			V[k].y = 0.0 + sqrtf(-Q(rabs)) * cosf(ksi) * (-xs/rabs);
			V[k].z = 0.0 + sqrtf(-Q(rabs)) * sinf(ksi) * (-ys/rabs);
			V[k].w = 0.0;
			continue;
		}
				if(r < 0.001)
				{
					V[k].x = 1.0/sqrtf(3.0)* sqrtf(2*fabsf(Q(r))) + 1.0/sqrtf(3.0)*gauss(0.5*sqrtf(2*fabsf(Q(r)))) + sqrtf(-Q(rabs)) * cosf(ksi) * ys/rabs;
					V[k].y = 1.0/sqrtf(3.0)* sqrtf(2*fabsf(Q(r))) + 1.0/sqrtf(3.0)*gauss(0.5*sqrtf(2*fabsf(Q(r)))) + sqrtf(-Q(rabs)) * cosf(ksi) * (-xs/rabs);
					V[k].z = 1.0/sqrtf(3.0)* sqrtf(2*fabsf(Q(r))) + 1.0/sqrtf(3.0)*gauss(0.5*sqrtf(2*fabsf(Q(r)))) + sqrtf(-Q(rabs)) * sinf(ksi) * (-ys/rabs);
					V[k].w = 0.0;
					continue;
				}

			//vteta = gauss(sqrtf(sigma2(r,3)));
			//vphi = gauss(sqrtf(sigma2(r,3)));
			float VR2 =0.0;// sigma2(r,3); 


			//phi = acos(X[k].x/r);
			//teta = asin(X[k].z/r);
			//float Xs = X[k].x - xs;
			//float Ys = X[k].y - ys;
			//float Zs = X[k].z - zs;
			//x = vr * (-Xs/sqrtf(r*r - Zs*Zs)) * sqrtf(r*r - Zs*Zs)/r + vphi * (Ys/sqrtf(r*r - Zs*Zs)) + vteta * fabsf(Zs/r) * (-Xs/sqrtf(r*r - Zs*Zs)); 
			//y = vr * (-Ys/sqrtf(r*r - Zs*Zs)) * sqrtf(r*r - Zs*Zs)/r + vphi * -Xs/sqrtf(r*r - Zs*Zs) + vteta * fabsf(Zs/r) * (-Ys/sqrtf(r*r - Zs*Zs)); 
			//z = vr * (-Zs/r)  + vteta * sqrtf(r*r - Zs*Zs)/r;  
			 
				
		//	float v = sqrtf(x * x + y * y + z * z); 
			//float Vesc = sqrtf(2*fabsf(Qs(r)));
			float Vesc = sqrtf(2*fabsf(Qs(r)));
			float vr = veldistr(10.0,sqrtf(VR2));
		
			while(vr > 0.95*Vesc)
			{
				
				//vteta = gauss(sqrtf(sigma2(r,3)));
				//vphi = gauss(sqrtf(sigma2(r,3)));
				vr = veldistr(10.0,sqrtf(VR2));
				//vteta =vteta;//*random()*0.8;
				//vphi =vphi;//*random()*0.8;;
				//vr = vr;//*random()*0.8;;


				//x = vr * (-Xs/sqrtf(r*r - Zs*Zs)) * sqrtf(r*r - Zs*Zs)/r + vphi * (Ys/sqrtf(r*r - Zs*Zs)) + vteta * (Zs/r) * (-Xs/sqrtf(r*r - Zs*Zs)); 
				//	y = vr * (-Ys/sqrtf(r*r - Zs*Zs)) * sqrtf(r*r - Zs*Zs)/r + vphi * (-Xs/sqrtf(r*r - Zs*Zs)) + vteta * (Zs/r) * (-Ys/sqrtf(r*r - Zs*Zs)); 
				//z = vr * (-Zs/r)  + vteta * sqrtf(r*r - Zs*Zs)/r;  
				//v = sqrtf(x * x + y * y + z * z); 
				//printf("bulge %f\t%f\t%f\t%f\t%f\t%f\n",X[k].x,X[k].y,X[k].z, r,vr ,Vesc);
				//system("pause");
			}
			phi  = random()*2.0*PI;
			teta = asin(random()*2.0 - 1.0);
			x = vr * cos(phi)*cos(teta); 
			y = vr * sin(phi)*cos(teta);
			z = vr * sin(teta);

			
			V[k].x = x + sqrtf(-Q(rabs))  * (-ys/rabs) * orb;
			V[k].y = y + sqrtf(-Q(rabs))  * (xs/rabs) * orb;
			V[k].z = z ;
			V[k].w = 0;
			
			//printf("satellite %f\t%f\t%f\n",(-ys/rabs)  ,orb, sqrtf(-Q(rabs)) );
		}*/
	}
	
  
}
//-----------------------------------------------------------------------------

// запись координат частиц в файл c именем out/iteration_(i).txt

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

int  main(int argc, char *argv[])
{
	/*FILE *file;

	file = fopen("initialparam.txt","r");

	int t=0;
	float input[8];
	char result_sting[20];
	while(fgets(result_sting,sizeof(result_sting),file))
        {

            input[t]=atof( result_sting  );                 
            t++;   
        }
 
        fclose(file);

		for (int i=0; i<8; i++)
    {
        
		
		printf(" %f\n", input[i]);
        
    }
	*/
	//параметры спутника
	/*psi = input[0]*PI/180.0;
	rabs = input[1];
	xs = 0.0;
	zs = rabs*sin(psi);
	ys = rabs*cos(psi);
	orb = input[6];
	r0 = input[2];
	//параметры галактики
	Qt = input[3];
	//if(input[7]== 0.0)NforBulge=0;
	Mb = input[4];
	a =0.05;// input[5];
*/
	alglib::hqrndseed(time(NULL),time(NULL),state);

    InitParticles();
	char FileName[32];
	
	
    sprintf(FileName, "test.txt");
	
    FILE * out = fopen(FileName, "w+");

	for (int i=0; i<1000; i++)
    {
        
		fprintf(out," %f\n", Bin[i]);
        
    }
	
    fclose(out);

	NbodyCu::Run(TimeStep, X, V, A);
	system("pause");
    delete [] X;
	delete [] A;
	delete [] V;

	return 0;
}

