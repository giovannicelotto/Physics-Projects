#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

double v0, L,N,h;		//valore del potenziale dentro la buca, larghezza della buca, N grid spacing, h spaziatura
double result;			//double usato nella funzione successiva
double prop(int N, double v2,double v1,double v0, double epsilon, double fn1, double fn0)	//propagatore di Numerov
{
result=(2*(1-5*(epsilon-v1)/(6*pow(N,2)))*fn1-(1+(epsilon-v0)/(6*N*N))*fn0)/(1+(epsilon-v2)/(6*N*N));
return result;
}


int main()
{

cout<<"Inserire valore di L, larghezza buca."<<endl;
cin>>L;
cout<<"Inserire N numero di punti nell'intervallo [0,L]"<<endl;
cin>>N;

h=L/N;	
														//spaziatura grid spacing
vector <double> fna;												//tre vettori utilizzati per la covnergenza nello shooting
vector <double> fnb;
vector <double> fnc;
vector <double> v;

fna.push_back(0);												//primo elemento imposto dalle condizioni al contorno
fnb.push_back(0);
fnc.push_back(0);

fna.push_back(2/N);												//secondo elemento determinato dalla scelta che la derivata prima nel primo punto valga 1
fnb.push_back(2/N);
fnc.push_back(2/N);

for(int i=2;i<N+1;i++)												//riempio momentaneamente con zero i vettori
{fna.push_back(0);
fnb.push_back(0);
fnc.push_back(0);}


int numero;													//variabile nello switch
cout<<"Che tipo di potenziale deve esserci all'interno della buca infinita? \n 1=costante nullo \n 2=lineare \n 3=cosinusoidale \n 4=gaussiano centrato a metà buca"<<endl;
cin>>numero;
switch (numero){
case 1:
	{
	//cout<<"Inserire valore di v, potenziale all'interno della buca nella regione [0, L]"<<endl;
	//cin>>v0;
	for(int i=0;i<N+1;i++)
		{
		v.push_back(0);
		}
		break;
	}

case 2:
	{
	cout<<"Inserire valore di k tale che il potenziale all'interno della buca nella regione [0, L] sia V=k*x/L"<<endl;
	double k;
	cin>>k;
	for(int i=0;i<N+1;i++)
		{
		v.push_back(k*i/N);
		}
		break;
	}

case 3:
	{
	cout<<"Inserire valore di k tale che il potenziale all'interno della buca nella regione [0, L] sia V=k*cos((pi*x)/L)"<<endl;
	double k;
	cin>>k;
	for(int i=0;i<N+1;i++)
		{
		v.push_back(k*cos(M_PI*i/N));
		}
		break;
	}
	
case 4:
	{
	cout<<"Inserire valore di k tale che il potenziale all'interno della buca nella regione [0, L] sia V=k*exp(-(x/L-1/2)^2/(sigma)^2)"<<endl;
	double k;
	cin>>k;
	cout<<"Inserire valore di sigma tale che il potenziale all'interno della buca nella regione [0, L] sia V=k*exp(-(x/L-1/2)^2/(2*(sigma)^2))"<<endl;
	double sigma;
	cin>>sigma;
	for(int i=0;i<N+1;i++)
		{
		v.push_back(k*exp(-pow((i/N-0.5),2)/(2*pow(sigma,2))));
		}
		break;
	
	
	}

}

ofstream scrittura("potenziale.txt");
for(int i=0;i<N+1;i++)
{
scrittura<<i*L/N<<"	"<<v.at(i)<<endl;
}

cout<<"Scegliere due guess iniziali per energie epsilon_a e epsilon_b"<<endl;		//epsilon è uguale all'autovalore E_n diviso (h_plankridotta/(m*L^2))
double epsa;
double epsb;
double epsc;
double temp;
cin>>epsa;
cout<<endl;
cin>>epsb;
cout<<endl;

if(epsa>epsb)
{
temp=epsb;
epsb=epsa;
epsa=temp;

}


cout<<"Definisci accuratezza richiesta"<<endl;
double eta;
cin>>eta;

double diff;
diff=pow(epsa-epsb,2);
diff=sqrt(diff);


while (diff>eta)
{
diff=pow(epsa-epsb,2);
diff=sqrt(diff);
cout.precision((9));
cout<<epsa<<"	"<<epsb<<endl;
if(diff<eta) break;

for(int i=2;i<N+1;i++)
	{
	fna.at(i)=prop(N,v.at(i),v.at(i-1),v.at(i-2),epsa,fna.at(i-1),fna.at(i-2));
	fnb.at(i)=prop(N,v.at(i),v.at(i-1),v.at(i-2),epsb,fnb.at(i-1),fnb.at(i-2));
	}

if(fna.at(N)*fnb.at(N)>0)
	{
	epsb+=1;
	continue;
	}
if(fna.at(N)*fnb.at(N)<0)
	{
	epsc=(epsa+epsb)/2;
	for(int i=2;i<N+1;i++)
	{
	fnc.at(i)=prop(N,v.at(i),v.at(i-1),v.at(i-2),epsc,fnc.at(i-1),fnc.at(i-2));
	}

if(fna.at(N)*fnc.at(N)<0) 
	{
	epsb=epsc;
	continue;
	}
if(fnb.at(N)*fnc.at(N)<0) 
	{
	epsa=epsc;
	continue;}
	}

}

cout<<"Il corrispondente autovalore in unità atomiche (h_planck_ridotta=1 e m=1) risulta "<<epsa/pow(L,2)<<endl;
ofstream output("autofunzione.txt");


for(int i=2;i<N+1;i++)
{
fna.at(i)=prop(N,v.at(i),v.at(i-1),v.at(i-2),epsa,fna.at(i-1),fna.at(i-2));
}


double somma=0;				//normalizzo il modulo quadro
for(int i=0;i<N+1;i++)
{
somma+=pow(fna.at(i),2)*h;
}


for(int i=0;i<N+1;i++)
{
output<<i*h<<"	"<<fna.at(i)/sqrt(somma)<<endl;}
cout<<"è stata stampata in \"autofunzione.txt\" l'autofunzione corrispondente all'autovalore: "<<epsa/pow(L,2)<<endl;


return 0;
}
