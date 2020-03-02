//00ベクトルのみの減衰率分布を作成するマクロ

#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<math.h>
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

// 00データを読み込んでカウンター配列に格納する関数
void read00( string data00, Int_t counter[16][16] )
{
  Int_t x,y=0;
  Int_t i,j;
  ifstream data(data00.c_str());

  for(i=0; i<16; i++)
  {
    for(j=0; j<16; j++)
    {
      counter[i][j] = 0;
    }
  }

	while(!data.eof())
	{
		data >> x >> y;
		counter[x][y-16]++;
	}
  x = y = 0;

  return;
}

Int_t i, j;
Int_t counter[16][16];
Double_t BGtime=724;                      //要変更
Double_t OBtime=337;                      //要変更
Double_t mtime_ratio = OBtime/BGtime;

void ana00(){

  // 00データからoututfileを作るか、アウトプットからヒストを作るかを選択する
  string data00;
  Int_t choice;
  Int_t select;

  cout << "Select 1 or 2" << endl << "1: Make a output file from 00data." <<
    endl << "2: Make 2D histgrams from output data." << endl << "INPUT>>" ;
  cin >> choice;

// choice=1 00_data.txt から 00_output.txtを作成する。
  if(choice == 1)
  {
  cout << "Select 1:BG or 2:OB" << endl << "INPUT>>";
  cin >> select;

  if (select == 1)
  {
    data00 = "BG_00_data.txt";
    ofstream BG00("BG_00_output.txt");
    read00( data00, counter);
    cout << ">> Loding BG_00_data.txt" << endl;
    for(i=0; i <16; i++)
    {
      for(j=0; j<16; j++)
      {
        BG00 << i << " " << j << " " << counter[i][j] << endl;
      }
    }
    cout << ">> Implimantation is done." << endl;
  }
  else if (select == 2)
  {
    data00 = "OB_00_data.txt";
    ofstream OB00("OB_00_output.txt");
    read00( data00, counter);
    cout << ">> Loading OB_00_data.txt" << endl;
    for(i=0; i <16; i++)
    {
      for(j=0; j<16; j++)
      {
        OB00 << i << " " << j << " " << counter[i][j] << endl;
      }
    }
    cout << ">> Implimantation is done." << endl;
  }
  else if (select  > 2){cout << "Please select 1 or 2";}
  } //choice == 1


// choice=2 00_output.txtを読み込んで減衰率分布を作る
  else if ( choice == 2)
  {

   ifstream BGoutput("BG_00_output.txt");
   ifstream OBoutput("OB_00_output.txt");
   Double_t BGcounter[16][16];
   Double_t OBcounter[16][16];
   Int_t i,j;
   Int_t x,y;

  for(i=0; i <16; i++)
   {
     for(j=0; j<16; j++)
     {
        BGoutput >> x >> y >> BGcounter[i][j];
        cout << x << " " << y << " " << BGcounter[i][j] << endl;
     }
   }
   for(i=0; i <16; i++)
   {
     for(j=0; j<16; j++)
     {
        OBoutput >> x >> y >> OBcounter[i][j];
        cout << x << " " << y << " " << OBcounter[i][j] << endl;
     }
   }

	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPalette(1);

	Int_t x_bin_width=16;
	Int_t y_bin_width=16;
	Double_t xmin=-0.5;
	Double_t xmax=15.5;
	Double_t ymin=-0.5;
	Double_t ymax=15.5;

	TCanvas *cvs1 = new TCanvas("cvs1", "cvs1", 900, 900 );
	TH2D *hist1 =new TH2D("hist1","hist1",x_bin_width,xmin,xmax,y_bin_width,ymin,ymax);
	hist1->SetStats(0);
	hist1->SetTitle("BG 00_event");
	hist1->GetXaxis()->SetTitle("sheet1 ch");
	hist1->GetXaxis()->CenterTitle();
	hist1->GetYaxis()->SetTitle("sheet2 ch");
	hist1->GetYaxis()->CenterTitle();

	for(i=0; i<16; i++)
	{
		for(j=0; j<16; j++)
		{
			hist1->Fill(i,j,BGcounter[i][j]);
		}
	}

	hist1->Draw("colz");

	TCanvas *cvs2 = new TCanvas("cvs2", "cvs2", 900, 900 );
	TH2D *hist2 =new TH2D("hist2","hist2",x_bin_width,xmin,xmax,y_bin_width,ymin,ymax);
	hist2->SetStats(0);
	hist2->SetTitle("OB 00_event");
	hist2->GetXaxis()->SetTitle("sheet1 ch");
	hist2->GetXaxis()->CenterTitle();
	hist2->GetYaxis()->SetTitle("sheet2 ch");
	hist2->GetYaxis()->CenterTitle();

	for(i=0; i<16; i++)
	{
		for(j=0; j<16; j++)
		{
			hist2->Fill(i,j,OBcounter[i][j]);
		}
	}

	hist2->Draw("colz");

  TCanvas *cvs3=new TCanvas("cvs3","cvs3",900,900);
  TH2D *hist3 =new TH2D("hist3","hist3",x_bin_width,xmin,xmax,y_bin_width,ymin,ymax);
  hist3-> SetTitle("substraction");
  hist3->SetStats(0);
  hist3->Add(hist1,hist2,1,-mtime_ratio );
  hist3->GetXaxis()->SetTitle("sheet1 ch");
  hist3->GetXaxis()->CenterTitle();
  hist3->GetYaxis()->SetTitle("sheet2 ch");
  hist3->GetYaxis()->CenterTitle();
  hist3->Draw("colz");
  cvs3->Print("0_Substraction.gif");


  TCanvas *cvs4=new TCanvas("cvs4","cvs4",900,900);
  TH2D *hist4 =new TH2D("hist4","hist4",x_bin_width,xmin,xmax,y_bin_width,ymin,ymax);
  hist4->Divide(hist3,hist1);
  hist4->SetTitle("Attenuation");
  hist4->GetXaxis()->SetTitle("sheet1 ch");
  hist4->GetXaxis()->CenterTitle();
  hist4->GetYaxis()->SetTitle("sheet2 ch");
  hist4->GetYaxis()->CenterTitle();

  hist4->Draw("colz");
  //  attenuation->SetMaximum(1);
  // attenuation->SetMinimum(0);
  hist4->SetStats(0);


  TCanvas *cvs5=new TCanvas("cvs5","cvs5",900,900);
  TH2D *hist5 =new TH2D("hist5","hist5",x_bin_width,xmin,xmax,y_bin_width,ymin,ymax);

  Double_t minius_reg;
  for(Int_t i=0;i<32;i++)
   {
     for(Int_t j=0;j<32;j++)
       {

	 minius_reg = hist4->GetBinContent(i,j);
	 if(minius_reg<=0)
	   {
	     minius_reg=0.001;
	     hist5->Fill(i-16,j-16,minius_reg);
	   }
	 else if(minius_reg<1)
	   {
	     hist5->Fill(i-16,j-16,minius_reg);
     }
	 else if(minius_reg<100)
	   {
	     minius_reg=0.001;
	     hist5->Fill(i-16,j-16,minius_reg);
	   }
       }
   }

  hist5->SetTitle("Attenuation_custom");
  hist5->GetXaxis()->SetTitle("sheet1 ch");
  hist5->GetXaxis()->CenterTitle();
  hist5->GetYaxis()->SetTitle("sheet2 ch");
  hist5->GetYaxis()->CenterTitle();
  hist5->SetMaximum(0.2);
  hist5->SetMinimum(0);
  hist5->Draw("colz");
  hist5->SetStats(0);
  cvs5->Print("0_atten_custom.gif");
  } //choice 2
  return;
}
