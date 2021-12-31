/*
//firstly load functions.C

this macro provides a Monte Carlo Simulation based on the Glauber model
that estimates the number of collisions between nucleons and number of nucleons
participants between in a scattering between two nuclei of Pb-208
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
void mcg(int entries=10000,int binpart=450, int binb=250){

//CONSTANTS
int A=208;					//number of nucleons per nucleus
double R=6.624;				//radius of fermi distribution
double a=0.549;				//diffusivity of fermi distribution
double sigma=6.76;			//cross section in fm**2
double d=0.4;				//distance required betweeen two nucleons of the same nucles	
double maxrad=100;			//maximum of radial coordinate

//functions used to generate random numbers
TF1* f=new TF1("f","pow(x,2)/(1+exp((x-[0])/[1]))",0,maxrad);	//radial distribution
TF1* f2=new TF1("f2","sin(x)",0,TMath::Pi());					//polar angle distribution
TF1* fb= new TF1("fb","x",0,20);									//impact parameter distribution
TF1* fphi=new TF1("fphi","1",0,1);

f->SetParameters(R,a);
f->SetNpx(5000);
f2->SetNpx(5000);
fb->SetNpx(5000);
fphi->SetNpx(5000);


//coordinate radiali e cartesiane, vettore con coordinate dei nucleoniA e vettore con coordinate nuncleoniB
double r,theta,phi,x,y,z;
vector <double> xi,yi,zi,xi2,yi2,zi2;

vector <bool> Npart;		//boolean vector for nucleons false=not partecipating, true =partecipanting
TNtuple* ntuple = new TNtuple("ntuple","Demo ntuple","Npart:Ncoll:b");

//parametri e counter
double b;			//parametero di impatto
int p=0;			//counter partecipanti
int count=0;			//counter collisioni
bool flag=false;		//nucleone tocca o non tocca i precedenti

//OUTPUT

TGraph* gPart=new TGraph();												//TGraph: partecipanti vs b
TGraph* gColl=new TGraph();												//TGraph: collisioni vs b
TH1F* hcoll=new TH1F("h","N_{events} vs N_{coll};N_{coll};N_{events}",1250,0,2500);				//eventi vs collisioni
TH1F* hpart=new TH1F("hpart","N_{events} vs N_{part};N_{part};N_{events}",binpart,0,450);				//eventi vs partecipanti
TH1F* hb=new TH1F("hb","N_{events}with at least one collision;b [fm]; Counts (normalized)",binb,0,20);				//counts vs b

//histograms to plot distributions of nucleons
TH1F* hr=new TH1F("hr","radial distribution; Radius [fm];Counts (normalized)",binb,0,20);
TH1F* hphi=new TH1F("hphi","azimutal angle #phi distribution;[rad];Counts (normalized)",62800,0,2*TMath::Pi());
TH1F* htheta=new TH1F("htheta","polar angle #theta distribution;[rad];Counts (normalized)",1570,0,TMath::Pi());

beauty_histo(hr,1.4);
beauty_histo(hphi,1.4);
beauty_histo(htheta,1.4);


//Ricentramento
double xcm1,ycm1,zcm1,xcm2,ycm2,zcm2;
//variabili di stato per mostrare percentuale del processo
int stato;
int st_output=0.;
cout<<st_output<<"\% \r";
cout.flush();
int t=0;		//counter of the event t-th
int nn=0;		//mumber of events without collisions
	

	while(nn!=entries)
	{
		
		xcm1=0;ycm1=0;zcm1=0;xcm2=0;ycm2=0;zcm2=0;
	//verifica in output lo stato	
	{
	stato=(nn*100)/(entries);
	if(stato!=st_output)
	{st_output=stato;
	cout<<st_output<<"\%\r";
	cout.flush();}
	}
	b=fb->GetRandom(0,20);
	for(int i=0;i<2*A;i++)Npart.push_back(false);
	
	//NUCLEO A:	genero 208 nucleoni
	for(int i=0;i<A;i++)
			{
			r=f->GetRandom(0,maxrad);
			theta = f2->GetRandom(0,TMath::Pi());
			phi = 2*TMath::Pi()*fphi->GetRandom(0,1);

			x=r*cos(phi)*sin(theta);
			y=r*sin(phi)*sin(theta);
			z=r*cos(theta);

//check if nucleon touches the other nucleons (and changes value of the boolean and decreases of 1 unit numbe of nucleons generated)
			for (int j=0;j<i;j++)
				{
				if(sqrt(pow(x-xi.at(j),2) + pow(y-yi.at(j),2) + pow(z-zi.at(j),2)  )<d)	
					{
					i=i-1;
					flag=true;
					break;
					}
				}
	//if the nucleons touches one nucleon already generated breaks the cycle
				if(flag==true)
				{
					flag=false;
					continue;
				}
	//otherwise save it
				xi.push_back(x);
				yi.push_back(y);
				zi.push_back(z);
				hr->Fill(r);
				htheta->Fill(theta);
				hphi->Fill(phi);
			}

	//NUCLEO B
	for(int i=0;i<A;i++)
			{
		r=f->GetRandom(0,maxrad);
		theta = f2->GetRandom(0,TMath::Pi());
		phi = 2*TMath::Pi()*fphi->GetRandom(0,1);

		x=b+r*cos(phi)*sin(theta);
		y=r*sin(phi)*sin(theta);
		z=20+r*cos(theta);


	//check if nucleon touches the other nucleons (and changes value of the boolean and decreases of 1 unit numbe of nucleons generated)
		for (int j=0;j<i;j++)
			{
			if(sqrt(pow(x-xi2.at(j),2) + pow(y-yi2.at(j),2) + pow(z-zi2.at(j),2)     )<d)
				{
				i=i-1;
				flag=true;
				break;
				}
			}
	//if the nucleons touches one nucleon already generated breaks the cycle		
	if(flag==true)
		{
		flag=false;
		continue;
		}
	//otherwise save it
		xi2.push_back(x);
		yi2.push_back(y);
		zi2.push_back(z);
		}



//center the nuclei with respect to its center of mass
for (int i=0;i<A;i++)
{
	xcm1=xcm1+xi.at(i);
	ycm1=ycm1+yi.at(i);

	xcm2=xcm2+xi2.at(i);
	ycm2=ycm2+yi2.at(i);

}

//move all the nucleons to have a center of mass centered in the wanted position
for (int i=0;i<A;i++)
{
		xi.at(i)=xi.at(i)-xcm1/208;
		yi.at(i)=yi.at(i)-ycm1/208;

		xi2.at(i)=xi2.at(i)-xcm2/208+b;
		yi2.at(i)=yi2.at(i)-ycm2/208;

	
}



//counting collisions
for(int i=0;i<A;i++)	//nucleon of nucleus A fixed
{
	for(int j=0;j<A;j++)	//reading all nucleons of nucleus B
	{
	
	if	(sqrt(pow(xi.at(i)-xi2.at(j),2) + pow(yi.at(i)-yi2.at(j),2))<sqrt(sigma/TMath::Pi()))
	{
		count++;
		Npart.at(i)=true;
		Npart.at(j+A)=true;
		}
	}
}

//counting participants
for(int i=0;i<Npart.size();i++){if (Npart.at(i)==true)p++;}
		
/*So far:
number of collisions stored in "count"
number of participants is stored in "p"
*/		
if(p!=0)
{
	gColl->SetPoint(nn,b,count);
	gPart->SetPoint(nn,b,p);
	nn++;		//increasing number of events with at least one collision
	ntuple->Fill(p,count,b);

}


hcoll->Fill(count);
hpart->Fill(p);
p=0;
count=0;
xi.clear();
yi.clear();
zi.clear();
xi2.clear();
yi2.clear();
zi2.clear();	
Npart.clear();
t++;			//number of total cycles
}


		//GRAPHIC
	
TCanvas* c1 = new TCanvas("c1","mycanvas",1500,800);
TCanvas* c2= new TCanvas("c2","canvas consistenza",1200,800);
TCanvas* c3= new TCanvas("c3","canvas b!=0",1080,800);


c1->Divide(2,2);
c2->Divide(2,2);
c2->cd(1);
gPad->SetRightMargin(0.001);
gPad->SetLeftMargin(0.2);
double scale=hr->Integral();
hr->DrawNormalized();
c2->cd(2);
gPad->SetLeftMargin(0.2);
htheta->DrawNormalized();
c2->cd(3);
gPad->SetLeftMargin(0.2);
hphi->DrawNormalized();
c1->cd(2);
gColl->SetTitle("N_{coll} vs b;b[fm];N_{coll}");
gColl->SetMarkerStyle(1);
beauty_graph(gColl,0.65,0.07);
//gColl->GetYaxis()->SetTicks("9,5,0");
gColl->DrawClone("ap");
c1->cd(1);
gPad->SetRightMargin(0.05);
gPart->SetTitle("N_{part} vs b;b[fm];N_{part}");
gPart->SetMarkerStyle(1);
beauty_graph(gPart,0.65,0.07);
gPart->DrawClone("ap");


c1->cd(4);
beauty_histo(hcoll,0.65,0.07);
hcoll->DrawNormalized();		//#events_vs_coll
gPad->SetLogy();
c1->cd(3);
gPad->SetRightMargin(0.05);
//#events_vs_part
beauty_histo(hpart,0.65,0.07);
hpart->DrawNormalized();
gPad->SetLogy();		

vector <double> fasciacoll;
vector <double> fasciapart;
vector <double> collnonnulle;
vector <double> partnonnulli;
vector <double> b_ord;
vector <double> b_ord2;


//elimino gli eventi con zero collisioni
for(int i=0;i<nn;i++)
{

	collnonnulle.push_back(gColl->GetY()[i]);
	b_ord.push_back(gColl->GetX()[i]);

	partnonnulli.push_back(gPart->GetY()[i]);
	b_ord2.push_back(gPart->GetX()[i]);


}


cout<<"*************************"<<endl<<endl;


//INIZIO ORDINAMENTO del vettore dei b e in modo parallelo i vettori di Ncoll e Npart

double inf=0;
double infcoll=0;
double infpart=0;
int sign=0;


for(int j=0;j<b_ord.size();j++)
{
inf=b_ord.at(j);
for(int i=j;i<collnonnulle.size();i++)
{
	if (b_ord.at(i)<inf)
	{
		inf=b_ord.at(i);
		infcoll=collnonnulle.at(i);
		infpart=partnonnulli.at(i);
		sign=i;
	}
}
b_ord.at(sign)=b_ord.at(j);
b_ord.at(j)=inf;

collnonnulle.at(sign)=collnonnulle.at(j);
collnonnulle.at(j)=infcoll;

partnonnulli.at(sign)=partnonnulli.at(j);
partnonnulli.at(j)=infpart;
}






for(int i=0;i<b_ord.size();i++)
{
	hb->Fill(b_ord.at(i));
}
c3->cd();
auto *p1 = new TPad("p3_1","p3_1",0.,0.5,0.49,1); p1->Draw();
	p1->SetTopMargin(0.1);
	p1->SetRightMargin(0.1);
	p1->SetBottomMargin(0.1);
	
auto *p2 = new TPad("p3_2","p3_2",0.51,0.5,1.,1.);p2->Draw();
	p2->SetTopMargin(0.1);
	p2->SetLeftMargin(0.08);
	p2->SetBottomMargin(0.1);
	
	auto *p3 = new TPad("p3_3","p3_3",0.,0.,1.,0.5);p3->Draw();
	p3->SetBottomMargin(0.1);
	p3->SetTopMargin(0.005);
	
p3->cd();
beauty_histo(hb);
hb->DrawNormalized();


double max=collnonnulle.size()/20;
double min=0;
/*
TH1F* hb2=new TH1F("hb2","p ",binb,0,20);


for(int j=0;j<2;j++)
{
if((j%2)!=0) hb2->SetFillColor(2);
if((j%2)==0) hb2->SetFillColor(17);

for(int i=0;i<b_ord.size();i++)
{
	if (b_ord.at(i)<b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hb2->Fill(b_ord.at(i));}
}
//for(int i=0;i<199;i++) hpart2->SetBinContent(i, 1.*hpart2->GetBinContent(i)*(bincontent) );
	hb2->Scale(1./hb->Integral());
	hb2->DrawClone("same hist");
	hb2->Reset("ICES");
	min=min+collnonnulle.size()/20;
	max=max+collnonnulle.size()/20;
}

min=collnonnulle.size()/10;	
max=collnonnulle.size()/5;


for(int j=1;j<10;j++)
{
if((j%2)!=0) hb2->SetFillColor(17);
if((j%2)==0) hb2->SetFillColor(2);
for(int i=0;i<b_ord.size();i++)
{
	if (b_ord.at(i)<=b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hb2->Fill(b_ord.at(i));}
}
//for(int i=0;i<199;i++) hpart2->SetBinContent(i, 1.*hpart2->GetBinContent(i)*(bincontent) );
	hb2->Scale(1./hb->Integral());
	hb2->DrawClone("same hist");
	hb2->Reset("ICES");
	min=min+collnonnulle.size()/10;
	max=max+collnonnulle.size()/10;
}
*/

//histo  of participants with class of centrality


p1->cd();
gPad->SetLogy();
TH1F* hpart2= new TH1F("N_{part} divided in class of centrality","N_{part} divided in class of centrality;N_{part};Counts",binpart,0,450);
hpart2->SetFillColor(10);
beauty_histo(hpart2);
for(int i=0;i<b_ord.size();i++)
{
hpart2->Fill(partnonnulli.at(i));
}
hpart2->SetFillStyle(3003);
hpart2->DrawClone("hist");
hpart2->Reset("ices");
max=collnonnulle.size()/20;
min=0;

for(int j=0;j<2;j++)
{
if((j%2)!=0) hpart2->SetFillColor(2);
if((j%2)==0) hpart2->SetFillColor(17);
for(int i=0;i<b_ord.size();i++)
{
	if (b_ord.at(i)<=b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hpart2->Fill(partnonnulli.at(i));}
}
//for(int i=0;i<199;i++) hpart2->SetBinContent(i, 1.*hpart2->GetBinContent(i)*(bincontent) );
	hpart2->DrawClone("same hist");
	hpart2->Reset("ICES");
	min=min+collnonnulle.size()/20;
	max=max+collnonnulle.size()/20;
}
min=collnonnulle.size()/10;
max=2*collnonnulle.size()/10;
for(int j=1;j<10;j++)
{
if((j%2)!=0) hpart2->SetFillColor(17);
if((j%2)==0) hpart2->SetFillColor(2);
for(int i=0;i<b_ord.size();i++)
{
	
	if (b_ord.at(i)<=b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hpart2->Fill(partnonnulli.at(i));}
}
//for(int i=0;i<199;i++) hpart2->SetBinContent(i,1.*hpart2->GetBinContent(i)*(bincontent) );
	hpart2->DrawClone("same hist");
	hpart2->Reset("ICES");
	min=min+collnonnulle.size()/10;
	max=max+collnonnulle.size()/10;
}




//histo of collisions with class of centrality

p2->cd();
gPad->SetLogy();
TH1F* hcoll2= new TH1F("N_{coll} divided in class of centrality","N_{coll} divided in class of centrality;N_{coll};Counts",1250,0,2500);
hcoll2->SetFillColor(10);
beauty_histo(hcoll2);
for(int i=0;i<b_ord.size();i++)
{
hcoll2->Fill(collnonnulle.at(i));
}
hcoll2->SetFillStyle(3003);
hcoll2->DrawClone("hist");
hcoll2->Reset("ices");
max=collnonnulle.size()/20;
min=0;
for(int j=0;j<2;j++)
{
if((j%2)!=0) hcoll2->SetFillColor(2);
if((j%2)==0) hcoll2->SetFillColor(17);
for(int i=0;i<b_ord.size();i++)
{
	if (b_ord.at(i)<=b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hcoll2->Fill(collnonnulle.at(i));}
}
	hcoll2->DrawClone("same hist");
	hcoll2->Reset("ICES");
	min=min+collnonnulle.size()/20;
	max=max+collnonnulle.size()/20;
}
min=collnonnulle.size()/10;
max=2*collnonnulle.size()/10;
for(int j=1;j<10;j++)
{
if((j%2)!=0) hcoll2->SetFillColor(17);
if((j%2)==0) hcoll2->SetFillColor(2);
for(int i=0;i<b_ord.size();i++)
{	
	if (b_ord.at(i)<=b_ord.at(max-1) and b_ord.at(i)>=b_ord.at(min)) 
	{hcoll2->Fill(collnonnulle.at(i));}
}
//for(int i=0;i<199;i++) hcoll2->SetBinContent(i,1.*hcoll2->GetBinContent(i)*(bincontent) );
	hcoll2->DrawClone("same hist");
	hcoll2->Reset("ICES");
	min=min+collnonnulle.size()/10;
	max=max+collnonnulle.size()/10;
}



//OUTPUT



ofstream output("output.txt");
output.precision(6);
max=collnonnulle.size()/20;
min=0;
output<<"Number of events with at least one collision	"<<collnonnulle.size()<<" out of total: "<<t<<" "<<endl;

//20 classes with 5% of total events

output<<"cent	"<<"b_min	b_max	N_part	RMS	N_coll	RMS"<<endl;
for(int j=0;j<20;j++)
{
	for(int i=min;i<max;i++)
	{	
		fasciacoll.push_back(collnonnulle.at(i));
		fasciapart.push_back(partnonnulli.at(i));
	}
	output<<j*5<<"\%-"<<(j+1)*5<<"\%:  	"<<b_ord.at(min)<<"	"<<b_ord.at(max-1)<<"	"<<med(fasciapart)<<"	"<<dev(fasciapart)<<"	"<<med(fasciacoll)<<"	"<<dev(fasciacoll)<<endl;
	min=min+collnonnulle.size()/20;
	max=max+collnonnulle.size()/20;
	fasciacoll.clear();
	fasciapart.clear();
}


output<<endl<<"\n\n************************************************\n\n"<<endl;

//first 5 classes with 1% each

max=collnonnulle.size()/100;
min=0;
for(int j=0;j<5;j++)
{
	for(int i=min;i<max;i++)
	{	
		fasciacoll.push_back(collnonnulle.at(i));
		fasciapart.push_back(partnonnulli.at(i));
	}
	output<<j<<"-"<<(j+1)<<"\%	"<<b_ord.at(min)<<"	"<<b_ord.at(max-1)<<"	"<<med(fasciapart)<<"	"<<dev(fasciapart)<<"	"<<med(fasciacoll)<<"	"<<dev(fasciacoll)<<endl;
	min=min+collnonnulle.size()/100;
	max=max+collnonnulle.size()/100;
	fasciacoll.clear();
	fasciapart.clear();

}
output<<endl<<"\n\n************************************************\n\n"<<endl;
//classes with 10%
max=collnonnulle.size()/10;
min=0;
for(int j=0;j<10;j++)
{
	for(int i=min;i<max;i++)
	{	
		fasciacoll.push_back(collnonnulle.at(i));
		fasciapart.push_back(partnonnulli.at(i));
	}
	output<<j*10<<"-"<<(j+1)*10<<"\%	"<<b_ord.at(min)<<"	"<<b_ord.at(max-1)<<"	"<<med(fasciapart)<<"	"<<dev(fasciapart)<<"	"<<med(fasciacoll)<<"	"<<dev(fasciacoll)<<endl;
	min=min+collnonnulle.size()/10;
	max=max+collnonnulle.size()/10;
	fasciacoll.clear();
	fasciapart.clear();

}

output<<endl<<"\n\n************************************************\n\n"<<endl;
//classes with 20%
max=collnonnulle.size()/5;
min=0;
for(int j=0;j<5;j++)
{
	for(int i=min;i<max;i++)
	{	
		fasciacoll.push_back(collnonnulle.at(i));
		fasciapart.push_back(partnonnulli.at(i));
	}
	output<<j*20<<"-"<<(j+1)*20<<"\%	"<<b_ord.at(min)<<"	"<<b_ord.at(max-1)<<"	"<<med(fasciapart)<<"	"<<dev(fasciapart)<<"	"<<med(fasciacoll)<<"	"<<dev(fasciacoll)<<endl;
	min=min+collnonnulle.size()/5;
	max=max+collnonnulle.size()/5;
	fasciacoll.clear();
	fasciapart.clear();

}
max=0.5*collnonnulle.size();
min=0.3*collnonnulle.size();
for(int i=min;i<max;i++)
	{	
		fasciacoll.push_back(collnonnulle.at(i));
		fasciapart.push_back(partnonnulli.at(i));
	}
	output<<30<<"-"<<50<<"\%	"<<b_ord.at(min)<<"	"<<b_ord.at(max-1)<<"	"<<med(fasciapart)<<"	"<<dev(fasciapart)<<"	"<<med(fasciacoll)<<"	"<<dev(fasciacoll)<<endl;
cout<<"\n\n Saved data in \"output.txt\" with classes of centrlaties\n\n Saved a TNtuple file as a root file \"CollPartB.root\""<<endl;

TFile* hfile = new TFile("CollPartB.root","RECREATE","Demo ROOT file");
ntuple->Write();
hfile->Close();
gStyle->SetOptStat(1110);
}

