//---------------------------------------------------------------------------
#задачка "перепад давления"
#include <vcl.h>
#pragma hdrstop
// иванов2лаба.cpp:
#include <iostream>
#include <math.h>
#include "conio.h"
#include <fstream>

#include "ololo.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "iComponent"
#pragma link "iCustomComponent"
#pragma link "iPlotComponent"
#pragma link "iVCLComponent"
#pragma link "iXYPlot"
#pragma link "iAnalogOutput"
#pragma link "iEditCustom"
#pragma resource "*.dfm"
using namespace std;
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
{

        // Объявление переменных
	const int N=1000000;
	double x0=0.1, t1, gamma=3.4, k=3.5, q;

	double ron1=11.777,  ron2=1.1777, pn1=2000000.0, pn2=100000.0;
	double h=0.001, tau=0.00000001, alfa, delta/*, D=320.0*/;

	double *ro1, *p1, *E1, *m1, *x1, *u1;
	double *ro2, *p2, *E2,  *x2, *u2, *X, *ron, *q1;
	double t, x;
	int i, j, ikr, iter;
	ro1=(double *)malloc((N+3)*sizeof(double));
	ro2=(double *)malloc((N+3)*sizeof(double));
	ron=(double *)malloc((N+3)*sizeof(double));
	p1=(double *)malloc((N+3)*sizeof(double));
	p2=(double *)malloc((N+3)*sizeof(double));
	E1=(double *)malloc((N+3)*sizeof(double));
	E2=(double *)malloc((N+3)*sizeof(double));
	m1=(double *)malloc((N+3)*sizeof(double));
	X=(double *)malloc((N+3)*sizeof(double));
	u1=(double *)malloc((N+3)*sizeof(double));
	u2=(double *)malloc((N+3)*sizeof(double));
	x1=(double *)malloc((N+3)*sizeof(double));
	x2=(double *)malloc((N+3)*sizeof(double));
        q1=(double *)malloc((N+3)*sizeof(double));
        t1= iAnalogOutput1->Value;
        j=0;
	for (i=1, x=0.; x<=1.; x+=h, i++)
	{
		u1[i]=0.;
		if (x>=0. && x<=x0)
                {
		   p1[i]=pn1;
		   ro1[i]=ron1;
		}
		else
		{
		p1[i]=pn2;
		ro1[i]=ron2;
		}
		ron[i]=ro1[i];
		E1[i]=(1./(gamma-1.)*(p1[i]/ro1[i]));
		x1[i]=x;
		m1[i]=ro1[i]*h;
		ikr=i;

	}


	for (j=1, t=tau; t<=t1; t+=tau, j++)
	{ u1[0]=0.;
   	  u2[0]=0.;



		for (i=1, x=0; x<=1.; x+=h, i+=1)
			{
				p1[0]=p1[1];
				u2[i]=u1[i]-(tau/m1[i])*(p1[i]-p1[i-1]); //_1. полусумма, тк аш меняется
				x2[i]=x1[i]+u2[i]*tau;                                 //_2_
				ikr=i;
			}
		u2[ikr+1]=u2[ikr];
		for (i=1, x=0; x<=1.; x+=h, i++)
			{
				ro2[i]=1./((1./ro1[i])+(tau/m1[i])*(u2[i+1]-u2[i]));//_3_
                                if (/*ro2[i]>ro2[i+1]&&*/u2[i]>u2[i+1])
                                {
                                       alfa = (gamma+1)*pow(k,2)*pow(h,2)/(2*pow(M_PI,2));
                                       q = alfa*ro2[i]*pow((u2[i+1]-u2[i])/h,2);
                                       q1[i] = q;
                                }
                                else q = 0.;
                                E2[i]=(E1[i]-(tau*q/m1[i])*(u2[i+1]-u2[i]))/(1.+(tau/m1[i])*(gamma-1.)*ro2[i]*(u2[i+1]-u2[i]));    //_4_  (берем на j слое давление)
			}
			for (i=1, x=0; x<=1.; x+=h, i++)
			{
                                if (/*ro2[i]>ro2[i+1]&&*/u2[i]>u2[i+1])
                                {
                                        alfa = (gamma+1)*pow(k,2)*pow(h,2)/(2*pow(M_PI,2));
                                        q = alfa*ro2[i]*pow((u2[i+1]-u2[i])/h,2);
                                }
                                else q = 0.;
                                p2[i]=(gamma-1.)*ro2[i]*E2[i]+q;
			}

		for (i=1, x=0; x<=1.; x+=h, i++)  //обновление массива
			{
				p1[i]=p2[i];
				ro1[i]=ro2[i];
				E1[i]=E2[i];
				u1[i]=u2[i];
				x1[i]=x2[i];
				X[i]=x;
				ikr=i;
		        }

	}


					for (i=1; i<=ikr; i+=1)
					{
                                                iXYPlot1->Channel[0]->AddXY(X[i],p1[i+1]);
                                                iXYPlot2->Channel[0]->AddXY(X[i],ro1[i+1]);
                                                iXYPlot3->Channel[0]->AddXY(X[i],E1[i+1]);
                                                iXYPlot4->Channel[0]->AddXY(X[i],q1[i+1]);

        				}
        delete [] ro1;
        delete [] ro2;
        delete [] ron;
        delete [] p1;
        delete [] p2;
        delete [] E1;
        delete [] E2;
        delete [] m1;
        delete [] X;
        delete [] u1;
        delete [] u2;
        delete [] x1;
        delete [] x2;
        delete [] q1;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Gohome1Click(TObject *Sender)
{
exit(0);        
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender)
{
        iXYPlot1->ClearAllData();
        iXYPlot2->ClearAllData();
        iXYPlot3->ClearAllData();
        iXYPlot4->ClearAllData();
}
//---------------------------------------------------------------------------
