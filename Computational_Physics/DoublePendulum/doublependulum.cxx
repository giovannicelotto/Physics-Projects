#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
//number of points, spacing, inital conditions of hamiltonian coordinates
int n;
double h, l, fi1_0 ,fi2_0,p1_0,p2_0;
double g=9.81;
vector <double> t;

//funzioni per R-G

double fi1point(vector <double> y)
{
double result=(y.at(2)-y.at(3)*cos(y.at(0)-y.at(1)))/(pow(l,2)*(1+pow(sin(y.at(0)-y.at(1)),2)));
return result;
}

double fi2point(vector <double> y)
{
double result=(2*y.at(3)-y.at(2)*cos(y.at(0)-y.at(1)))/(pow(l,2)*(1+pow(sin(y.at(0)-y.at(1)),2)));
return result;
}

double p1point(vector <double> y)
{
double result=1/(l*l*(1+pow(sin(y.at(0)-y.at(1)),2)))*(-y.at(2)*y.at(3)*sin(y.at(0)-y.at(1))+(y.at(2)*y.at(2)+2*y.at(3)*y.at(3)-2*y.at(2)*y.at(3)*cos(y.at(0)-y.at(1)))/(1+pow(sin(y.at(0)-y.at(1)),2))*sin(y.at(0)-y.at(1))*cos(y.at(0)-y.at(1)))-2*g*l*sin(y.at(0));
return result;
}


double p2point(vector <double> y)
{
double result=1/(l*l*(1+pow(sin(y.at(0)-y.at(1)),2)))*(y.at(2)*y.at(3)*sin(y.at(0)-y.at(1))-(y.at(2)*y.at(2)+2*y.at(3)*y.at(3)-2*y.at(2)*y.at(3)*cos(y.at(0)-y.at(1)))/(1+pow(sin(y.at(0)-y.at(1)),2))*sin(y.at(0)-y.at(1))*cos(y.at(0)-y.at(1)))-g*l*sin(y.at(1));
return result;
}


int main()
{
cout<<"insert spacing h"<<endl;
cin>>h;
cout<<"insert number of points"<<endl;
cin>>n;

for (int i=0; i<n;i++)
{t.push_back(i*h);}

//finora ho discretizzato l'asse temporale in N punti equispaziati

cout<<"insert lengths of the two pendulum";
cin>>l;
cout<<"Insert fi1_0 (angle in rad of first pendulum with respect to the vertical line)";
cin>>fi1_0;
cout<<"Insert fi2_0 (angle in rad of second pendulum with respect to the vertical line)";
cin>>fi2_0;
cout<<"Insert p1_0 generalized momentum of first pendulum";
cin>>p1_0;
cout<<"IInsert p2_0 generalized momentum of second pendulum";
cin>>p2_0;

vector <double> f1;
vector <double> f2;
vector <double> p1;
vector <double> p2;

vector <double> Y1;
vector <double> Y2;
vector <double> Y3;
vector <double> Y4;

for(int i=0;i<4;i++)
{Y1.push_back(0);
Y2.push_back(0);
Y3.push_back(0);
Y4.push_back(0);
}
//initial conditions
f1.push_back(fi1_0);
f2.push_back(fi2_0);
p1.push_back(p1_0);
p2.push_back(p2_0);

for(int i=0;i<n;i++)
{
//first point Y1
Y1.at(0)=f1.at(i);
Y1.at(1)=f2.at(i);
Y1.at(2)=p1.at(i);
Y1.at(3)=p2.at(i);
//second point Y2
Y2.at(0)=f1.at(i)+h/2*fi1point(Y1);
Y2.at(1)=f2.at(i)+h/2*fi2point(Y1);
Y2.at(2)=p1.at(i)+h/2*p1point(Y1);
Y2.at(3)=p2.at(i)+h/2*p2point(Y1);

//third point Y3
Y3.at(0)+=h/2*fi1point(Y2);
Y3.at(1)+=h/2*fi2point(Y2);
Y3.at(2)+=h/2*p1point(Y2);
Y3.at(3)+=h/2*p2point(Y2);
//fourth point Y4
Y4.at(0)=f1.at(i)+h*fi1point(Y3);
Y4.at(1)=f2.at(i)+h*fi2point(Y3);
Y4.at(2)=p1.at(i)+h*p1point(Y3);
Y4.at(3)=p2.at(i)+h*p2point(Y3);
//compute of y_(n+1)
f1.push_back(f1.at(i)+h/6*(fi1point(Y1)+2*fi1point(Y2)+2*fi1point(Y3)+fi1point(Y4)));
f2.push_back(f2.at(i)+h/6*(fi2point(Y1)+2*fi2point(Y2)+2*fi2point(Y3)+fi2point(Y4)));
p1.push_back(p1.at(i)+h/6*(p1point(Y1)+2*p1point(Y2)+2*p1point(Y3)+p1point(Y4)));
p2.push_back(p2.at(i)+h/6*(p2point(Y1)+2*p2point(Y2)+2*p2point(Y3)+p2point(Y4)));
}

ofstream output("doublependulum.dat");
for(int i=0;i<n;i++)
{output<<l*sin(f1.at(i))+l*sin(f2.at(i))<<"	"<<-l*cos(f1.at(i))-cos(f2.at(i))<<endl;}

return 0;
}
