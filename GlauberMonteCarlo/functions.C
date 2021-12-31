//functions used in the code

//boolean used to sort a vector
bool wayToSort(int i, int j) { return i > j; }

void beauty_histo(TH1F* h, double Yoffset=0.65,double Ytitlesize=0.05)
{
	h->GetYaxis()->SetLabelSize(0.05);
h->GetYaxis()->SetTitleSize(Ytitlesize);
h->GetYaxis()->SetTitleOffset(Yoffset);
h->GetXaxis()->SetLabelSize(0.05);
h->GetXaxis()->SetTitleSize(0.05);
h->GetXaxis()->SetTitleOffset(0.85);
	
}

void beauty_graph(TGraph* g, double Yoffset=0.65,double Ytitlesize=0.05)
{
	g->GetYaxis()->SetLabelSize(0.05);
g->GetYaxis()->SetTitleSize(Ytitlesize);
g->GetYaxis()->SetTitleOffset(Yoffset);
g->GetXaxis()->SetLabelSize(0.05);
g->GetXaxis()->SetTitleSize(0.05);
g->GetXaxis()->SetTitleOffset(0.85);
	
}
//sum
double sum(vector <double> a)
{
double value=0;
for (int i=0;i<a.size();i++)
{value=value+a.at(i);}
return value;
}

//somma quadratica di un vettore
double sumq(vector <double> a)
{
double value=0;
for (int i=0;i<a.size();i++)
{value=value+pow(a.at(i),2);}
return value;
}

//MEDIA
double med(vector <double> a)
{
double value=0;
for (int i=0;i<a.size();i++)
{
value=value+a.at(i);
}
double media=value/a.size();
return media;
}

//DEVIAZIONE STANDARD
double dev(vector <double> a)
{
double value=0;
for (int i=0;i<a.size();i++)
{
	value=value+(pow((a.at(i)-med(a)),2))/(a.size()-1);

}
return sqrt(value);
}
//ERRORE DELLA MEDIA
double err(vector <double> a)
{
double value=dev(a)/sqrt(a.size());
return value;
}
