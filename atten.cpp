#include <sstream>
#include <fstream>
#include <math.h>
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"

void atten()
{

  TFile *f = new TFile("out.root", "RECREATE");

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPalette(1); //ヒストグラムの色を指定する。

  string input_file_name_1 = "1_Background_output.txt"; //バックグラウンドのデータの名前
  string input_file_name_2 = "2_Object_output.txt";     //試験片計測のデータの名前
  Double_t background_time = 0.;                        //バックグラウンドの計測時間
  Double_t object_time = 0.;                            //試験片の計測時間

  Int_t x_bin_width = 31;
  Int_t y_bin_width = 31;
  Double_t xmin = -15.5;
  Double_t xmax = 15.5;
  Double_t ymin = -15.5;
  Double_t ymax = 15.5; //ヒストグラムの諸条件

  Double_t BG[31][31];
  Double_t OB[31][31];
  Double_t OB_error[31][31];
  Double_t BG_error[31][31];
  Double_t attenError[31][31];

  std::ifstream file1(input_file_name_1); //ファイルを読み込む
  std::ifstream file2(input_file_name_2); //ファイルを読み込む

  Int_t i, j;
  Double_t k;
  file1 >> i >> j >> k;
  background_time = k;
  cout << "BG_time = " << background_time << endl;
  file2 >> i >> j >> k;
  object_time = k;
  cout << "OB_time = " << object_time << endl;

  //BG_MEASUREMENT///////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist1 = new TH2D("hist1", "hist1", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax); //ヒストグラムを定義
  hist1->SetStats(0);                                                                         //統計ボックスの非表示
  hist1->SetTitle("Background_measurement");
  hist1->GetXaxis()->SetTitle("deltaX");
  hist1->GetXaxis()->CenterTitle();
  hist1->GetYaxis()->SetTitle("deltaY");
  hist1->GetYaxis()->CenterTitle();
  while (!file1.eof())
  {
    Double_t j, k, l;
    file1 >> j >> k >> l; //ファイルの3つの要素を変数に格納
    hist1->Fill(j, k, l); //ヒストグラムに値を入れ込む
    BG[(Int_t)j + 15][(Int_t)k + 15] = l;
    BG_error[(Int_t)j + 15][(Int_t)k + 15] = sqrt(l) / l;
  }

  //BG_ERROR////////////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist11 = new TH2D("hist11", "hist11", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax); //ヒストグラムを定義
  hist11->SetStats(0);
  hist11->SetTitle("Background_measurement_error");
  hist11->GetXaxis()->SetTitle("deltaX");
  hist11->GetXaxis()->CenterTitle();
  hist11->GetYaxis()->SetTitle("deltaY");
  hist11->GetYaxis()->CenterTitle();
  for (i = -15; i < 16; i++)
  {
    for (j = -15; j < 16; j++)
    {
      hist11->Fill(i, j, BG_error[i + 15][j + 15]);
    }
  }

  //OB_MEASUREMENT//////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist2 = new TH2D("hist2", "hist2", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax); //ヒストグラムを定義
  hist2->SetTitle("Object_measurement");
  hist2->GetXaxis()->SetTitle("deltaX");
  hist2->GetXaxis()->CenterTitle();
  hist2->GetYaxis()->SetTitle("deltaY");
  hist2->GetYaxis()->CenterTitle();
  hist2->SetStats(0);

  while (!file2.eof())
  {
    Double_t j, k, l;
    file2 >> j >> k >> l;
    hist2->Fill(j, k, l);
    OB[(Int_t)j + 15][(Int_t)k + 15] = l;
    OB_error[(Int_t)j + 15][(Int_t)k + 15] = sqrt(l) / l;
  }

  //OB_ERROR//////////////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist21 = new TH2D("hist21", "hist21", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax); //ヒストグラムを定義
  hist21->SetStats(0);                                                                           //統計ボックスの非表示
  hist21->SetTitle("Object_measurement_error");
  hist21->GetXaxis()->SetTitle("deltaX");
  hist21->GetXaxis()->CenterTitle();
  hist21->GetYaxis()->SetTitle("deltaY");
  hist21->GetYaxis()->CenterTitle();
  for (i = -15; i < 16; i++)
  {
    for (j = -15; j < 16; j++)
    {
      hist21->Fill(i, j, OB_error[i + 15][j + 15]);
    }
  }

  //SUBSTRACT//////////////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist3 = new TH2D("hist3", "hist3", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  hist3->SetTitle("substraction");
  hist3->SetStats(0);                                              //統計ボックスの非表示
  hist3->Add(hist1, hist2, 1 / background_time, -1 / object_time); //両データを計測時間で割ってデータの引き算をする
  hist3->GetXaxis()->SetTitle("deltaX");
  hist3->GetXaxis()->CenterTitle();
  hist3->GetYaxis()->SetTitle("deltaY");
  hist3->GetYaxis()->CenterTitle();

  //ATTENUATION/////////////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist4 = new TH2D("hist4", "hist4", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  hist4->Divide(hist3, hist1, 1, 1 / background_time); //データの割り算
  hist4->SetTitle("Attenuation");
  hist4->GetXaxis()->SetTitle("deltaX");
  hist4->GetXaxis()->CenterTitle();
  hist4->GetYaxis()->SetTitle("deltaY");
  hist4->GetYaxis()->CenterTitle();
  // hist4->SetMaximum(0.3);
  // hist4->SetMinimum(0);
  hist4->SetStats(0);

  //ATTENUATION ERROR///////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist41 = new TH2D("hist41", "hist41", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  Double_t obj, bg;

  for (Int_t i = -15; i < 16; i++)
  {
    for (Int_t j = -15; j < 16; j++)
    {
      obj = BG[i + 15][j + 15];
      bg = OB[i + 15][j + 15];
      attenError[i + 15][j + 15] = sqrt(obj / pow(bg, 2.0) + pow(obj, 2.0) / pow(bg, 3.0));
      hist41->Fill(i, j, attenError[i + 15][j + 15]);
    }
  }
  hist41->SetTitle("Attenuation Error");
  hist41->GetXaxis()->SetTitle("deltaX");
  hist41->GetXaxis()->CenterTitle();
  hist41->GetYaxis()->SetTitle("deltaY");
  hist41->GetYaxis()->CenterTitle();
  hist41->SetMaximum(0.1);
  // hist4->SetMinimum(0);
  hist41->SetStats(0);

  //ATRENUATION CUSTOM//////////////////////////////////////////////////////////////////////////////////////////
  TH2D *hist5 = new TH2D("hist5", "hist5", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  Double_t minius_reg;
  for (Int_t i = 0; i < 32; i++)
  {
    for (Int_t j = 0; j < 32; j++)
    {
      minius_reg = hist4->GetBinContent(i, j);
      if (minius_reg <= 0)
      {
        minius_reg = 0.001;
        hist5->Fill(i - 16, j - 16, minius_reg);
      }
      else if (minius_reg < 1)
      {
        hist5->Fill(i - 16, j - 16, minius_reg);
      }
      else if (minius_reg < 100)
      {
        minius_reg = 0.001;
        hist5->Fill(i - 16, j - 16, minius_reg);
      }
    }
  }
  hist5->SetTitle("Attenuation_custom");
  hist5->GetXaxis()->SetTitle("deltaX");
  hist5->GetXaxis()->CenterTitle();
  hist5->GetYaxis()->SetTitle("deltaY");
  hist5->GetYaxis()->CenterTitle();
  hist5->SetMaximum(0.35);
  hist5->SetMinimum(0);
  hist5->SetStats(0);

  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// BG  cph ////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  TH2D *hist6 = new TH2D("hist6", "hist6", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  hist6->Add(hist1, 1 / background_time); //データの割り算
  hist6->SetTitle("BG_cph");
  hist6->GetXaxis()->SetTitle("deltaX");
  hist6->GetXaxis()->CenterTitle();
  hist6->GetYaxis()->SetTitle("deltaY");
  hist6->GetYaxis()->CenterTitle();
  hist6->SetMaximum(130);
  // attenuation->SetMinimum(0);
  hist6->SetStats(0);

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// OBJ cph//////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  TH2D *hist7 = new TH2D("hist7", "hist7", x_bin_width, xmin, xmax, y_bin_width, ymin, ymax);
  hist7->Add(hist2, 1 / object_time); //データの割り算
  hist7->SetTitle("OBJ_cph");
  hist7->GetXaxis()->SetTitle("deltaX");
  hist7->GetXaxis()->CenterTitle();
  hist7->GetYaxis()->SetTitle("deltaY");
  hist7->GetYaxis()->CenterTitle();
  hist7->SetMaximum(130);
  // attenuation->SetMinimum(0);
  hist7->SetStats(0);

  //DRAW HIST//////////////////////////////////////////////////////////////////////
  TCanvas *c = new TCanvas("c", "c", 900, 900);
  hist1->Draw("Colz");
  c->Print("1_Background_measurement.pdf");
  hist11->Draw("Colz");
  c->Print("11_Background_measurement_error.pdf");
  hist2->Draw("Colz");
  c->Print("2_Object_measurement.pdf");
  hist21->Draw("Colz");
  c->Print("21_Object_measurement_error.pdf");
  hist4->Draw("Colz");
  c->Print("4_Attenuation.pdf");
  hist41->Draw("Colz");
  c->Print("41_Attenuation_error.pdf");
  hist5->Draw("Colz");
  c->Print("5_Attenuation_custom.pdf");
  hist6->Draw("Colz");
  c->Print("6_BG_cph.pdf");
  hist7->Draw("Colz");
  c->Print("7_OB_cph.pdf");

  c->Close();

  f->Write();
  return;
}
