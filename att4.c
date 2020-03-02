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
#include "TGaxis.h"
#include "TMath.h"
#include "TGraph.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <algorithm>
#include <dirent.h>
#include <stdlib.h>

bool StrString(const char*s1,const char *s2);
time_t  get_caltime(char *tstr, char *tstr2);

typedef struct
{
  double n_pressure;
  time_t n_time;
}pres_t;


void att4(){
const int timespan=3600.; //[sec]分割ファイルの時間間隔
const double b_pressure=1013.;//基準とする気圧[hpa]
const double b_temp=10.;//基準となる気温[℃]
//const double factor_pre=0.124; //補正の値[%/hPa]
const double factor_pre=0.1234; //補正の値[cph/hPa]

const int pix=16;
const int dive=1;

int pixs=pix/dive;


float psd_b[16][16][31][31];
float psd_o[16][16][31][31];
float psd_att[16][16][31][31];
float psd_err[16][16][31][31];
float psd_thick[16][16][31][31];
float back_sum[31][31];



float psd_b_dive[16][16];
float psd_o_dive[16][16];
float psd_att_dive[16][16];
float psd_thick_dive[16][16];
float psd_err_dive[16][16];

float vec_o[31][31];
float vec_b[31][31];
float vec_att[31][31];
float vec_err[31][31];
float vec_thick[31][31];


double count_n[400];
double count_np[400];
double pres_n[400];
int time[400];

for(int i=0; i<400; i++)
{
count_n[i]=0.;
count_np[i]=0.;
pres_n[i]=0.;
time[i]=0;
}

ofstream backd("back_data_vec.txt");

for(int i =0; i<pix; i++)
	{
	for(int j =0; j<pix; j++)
		{
			psd_b_dive[i][j]=0;
			psd_o_dive[i][j]=0;
			vec_thick[i][j]=0.;
			psd_thick_dive[i][j]=0;

			psd_att_dive[i][j]=0;
		}
	}



int x_vec=0;
int y_vec=0;



for(int l=0; l<31; l++)
{
	for(int k=0; k<31; k++)
	{
			vec_o[l][k]=0;
			vec_b[l][k]=0;
			vec_att[l][k]=0;
			vec_err[l][k]=0;
			back_sum[l][k]=0.0;

		for (int i=0; i<16; i++ )
			{
				for (int j=0; j<16; j++)
				{
						psd_b[i][j][k][l]=0.0;
						psd_o[i][j][k][l]=0.0;
						psd_att[i][j][k][l]=0.0;
				}
			}
	}
}
  gStyle->SetPalette(1);


	ifstream pre_file("pressure.txt");//行数を調べたい(配列のため)
	if(!pre_file.is_open())
{
	cout<<"no pressure file"<<endl;
	return;
}
	int n_line=0;
		 while(!pre_file.eof())
	   {
	     string dum;
	     getline(pre_file, dum);
	     n_line++;
	   }
	pre_file.close();
	////////////////大気圧データ保存
	pres_t *pres;
	size_t elm_size=sizeof(pres_t);
	pres=(pres_t*) malloc(n_line * elm_size);//配列の動的確保
	////////////////

	char buff[100];
	char buff2[100];
		time_t  cal_time;
	ifstream pre_file2("pressure.txt");
	int pre_n=0;
	while(!pre_file2.eof())
	{
//	  pre_file2>>buff>>buff2>>pres[pre_n].n_pressure;
pre_file2>>buff>>buff2>>pres[pre_n].n_pressure;//ここでの太陽磁場はまだ中性子[cpm]
//cout<<buff<<"   "<<buff2<<"  "<<endl;
	  cal_time = get_caltime(buff, buff2);
//cout<<ctime(&cal_time);
	  if(cal_time == -1)
	           printf("変換失敗\n");
	   else  pres[pre_n].n_time=cal_time;//時間に変換
	  pre_n++;
	}
	pre_file2.close();

/////読み込み準備
	stringstream ss1;
	const char* path_o="./object/";//文字列型
//	const char* path_o="./background/";//文字列型
	DIR *dp_o;//ファイル一覧を取得
	dirent* entry_o;//direntとは実行が成功したら0,失敗したら-1を返す
	string filename_o[4096];
	int filenum_o=0;
	dp_o = opendir(path_o);//path以下のやつを全部開く

	while((entry_o=readdir(dp_o))!=NULL)
		{
string tmp;
tmp=entry_o -> d_name;//読み込んだものの名前を記憶(おそらく拡張子以外のところ)
cout<<tmp<<endl;
if(StrString(tmp.c_str(),".dat"))
	{
		filename_o[filenum_o]=tmp;

		cout<<tmp<<" "<<filenum_o<<" "<<filename_o[filenum_o]<<endl;
		filenum_o++;
	}
		}

float o_time=((float)filenum_o+1.)*timespan;


////////////////background
		stringstream ss2;
		const char* path_b="./background/";//文字列型
//		const char* path_b="./background_aaa/";//文字列型
		DIR *dp_b;//ファイル一覧を取得
		dirent* entry_b;//direntとは実行が成功したら0,失敗したら-1を返す
		string filename_b[4096];
		int filenum_b=0;
		dp_b = opendir(path_b);//path以下のやつを全部開く

		while((entry_b=readdir(dp_b))!=NULL)
			{
	string tmp_b;
	tmp_b=entry_b -> d_name;//読み込んだものの名前を記憶(おそらく拡張子以外のところ)
	cout<<tmp_b<<endl;
	if(StrString(tmp_b.c_str(),".dat"))
		{
			filename_b[filenum_b]=tmp_b;

			cout<<tmp_b<<" "<<filenum_b<<" "<<filename_b[filenum_b]<<endl;
			filenum_b++;
		}
			}
float b_time=((float)filenum_b+1.)*timespan;


//////////////////////////ファイル読み込み開始
///object
			for(Int_t iiii=0;iiii<filenum_o;iiii++)
				{
		ss1.str("");

		ss1<<path_o<<filename_o[iiii];
		std::ifstream hoge(ss1.str().c_str());
		if(!hoge.is_open())
			{
				cout<<"no file"<<endl;
				return;
			}
		else
			{
					cout<<iiii<<": "<<ss1.str()<<" loading "<<endl;
			}

string dummy;
			time_t re_time;//ファイルの更新時間
			getline(hoge, dummy);//1行飛ばし

			hoge>>re_time;

			////////////ここで気圧の補正
						double correct_pre=0;

			      bool flagt=false;
			      time_t diftime=re_time-1800;//更新時間が分割ファイルの最後のため1800秒引く
			      for(int j=0; j<n_line-1; j++)
			        {
			          if(difftime(diftime,pres[j].n_time)<0 && flagt==false)
			            {
			              flagt=true;
			              if(abs(difftime(diftime,pres[j].n_time))>abs(difftime(diftime,pres[j-1].n_time)))//時間が近いほうの圧力を使用
			                {
												correct_pre=pres[j-1].n_pressure;
			                }
			              else if(abs(difftime(diftime,pres[j].n_time))<abs(difftime(diftime,pres[j-1].n_time)))//上に同じ
			                {
												correct_pre=pres[j].n_pressure;
			                }
			                else
			                cout<<"       correct error      "<<endl;
							cout<<correct_pre<<endl;
			            }
			            
//			            else {cout<<j<<"error"<<diftime<<"    "<<pres[j].n_time<<"   "<<flagt<<endl;};
			        }
			///////////
			while(!hoge.eof())
				{
					int mu1x_o,mu1y_o,mu2x_o,mu2y_o;
				//	back>>mu1x_b>>mu1y_b>>mu2y_b>>mu2x_b;
					hoge>>mu1x_o>>mu1y_o>>mu2x_o>>mu2y_o;
//								y_vec = (63-mu2y_o) - (31-mu1y_o) + 15;
//								x_vec = (mu2x_o-32) - mu1x_o + 15;
//								y_vec = (31-mu1y_o) - (63-mu2y_o) + 15;
	//							x_vec = mu1x_o - (mu2x_o-32) + 15;

//								y_vec = (31-mu1y_o) - (63-mu2y_o) + 15;
	//							x_vec = mu1x_o - (mu2x_o-32) + 15;

							  x_vec = (47-mu2x_o) - (15-mu1x_o) + 15;
					         y_vec = (63-mu2y_o) - (31-mu1y_o) + 15;


//					         x_vec = (mu2x_o-32) - mu1x_o + 15;
	//				         y_vec = (63-mu2y_o) - (31-mu1y_o) + 15;



//					psd_o[mu1x_o][mu1y_o-16][x_vec][y_vec]+=1.*(1.+((correct_pre-b_pressure)*factor_pre/100.))*(1.+((correct_temp-b_temp)*factor_temp/100.));
					psd_o[mu1x_o][mu1y_o-16][x_vec][y_vec]+=1.*(1.+((correct_pre-b_pressure)*factor_pre/100.));
//					psd_o[mu1x_o][mu1y_o-16][x_vec][y_vec]+=1.;
//					vec_o[x_vec][y_vec]+=1.*(1.+((correct_pre-b_pressure)*factor_pre/100.))*(1.+((correct_temp-b_temp)*factor_temp/100.));
					vec_o[x_vec][y_vec]+=1.*(1.+((correct_pre-b_pressure)*factor_pre/100.));
//					vec_o[x_vec][y_vec]+=1.;

					
				}

			}



/////////////////////backgrounfd

for(Int_t jjjj=0;jjjj<filenum_b;jjjj++)
	{
ss2.str("");

ss2<<path_b<<filename_b[jjjj];
std::ifstream hoge2(ss2.str().c_str());
if(!hoge2.is_open())
{
	cout<<"no file"<<endl;
	return;
}
else
{
		cout<<jjjj<<": "<<ss2.str()<<" loading "<<endl;
}

string dummy2;
time_t re_time2;//ファイルの更新時間
getline(hoge2, dummy2);//1行飛ばし

hoge2>>re_time2;

////////////ここで気圧の補正
			double correct_pre2=0;
			bool flagt2=false;
			time_t diftime2=re_time2-1800;//更新時間が分割ファイルの最後のため1800秒引く
			for(int j=0; j<n_line-1; j++)
				{
//						cout<<j<<"   "<<ctime(&pres[j].n_time);
					if(difftime(diftime2,pres[j].n_time)<0 && flagt2==false)
						{
							flagt2=true;
							if(abs(difftime(diftime2,pres[j].n_time))>abs(difftime(diftime2,pres[j-1].n_time)))//時間が近いほうの圧力を使用
								{
									correct_pre2=pres[j-1].n_pressure;
								}
							else if(abs(difftime(diftime2,pres[j].n_time))<abs(difftime(diftime2,pres[j-1].n_time)))//上に同じ
								{
									correct_pre2=pres[j].n_pressure;
								}
								else
								cout<<"       correct error  2    "<<endl;
											cout<<correct_pre2<<endl;
											time[jjjj]=pres[j].n_time;

						}						
				}
				
				
				
///////////
while(!hoge2.eof())
	{
		int mu1x_b,mu1y_b,mu2x_b,mu2y_b;
	//	back>>mu1x_b>>mu1y_b>>mu2y_b>>mu2x_b;
		hoge2>>mu1x_b>>mu1y_b>>mu2x_b>>mu2y_b;
//					y_vec = (63-mu2y_b) - (31-mu1y_b) + 15;
//					x_vec = (mu2x_b-32) - mu1x_b + 15;

//								y_vec = (31-mu1y_b) - (63-mu2y_b) + 15;
	//							x_vec = mu1x_b - (mu2x_b-32) + 15;

							  x_vec = (47-mu2x_b) - (15-mu1x_b) + 15;
					         y_vec = (63-mu2y_b) - (31-mu1y_b) + 15;


//					         x_vec = (mu2x_b-32) - mu1x_b + 15;
	//				         y_vec = (63-mu2y_b) - (31-mu1y_b) + 15;



//		psd_b[mu1x_b][mu1y_b-16][x_vec][y_vec]+=1.*(1.+((correct_pre2-b_pressure)*factor_pre/100.))*(1.+((correct_temp2-b_temp)*factor_temp/100.));
		psd_b[mu1x_b][mu1y_b-16][x_vec][y_vec]+=1.*(1.+((correct_pre2-b_pressure)*factor_pre/100.));
//		psd_b[mu1x_b][mu1y_b-16][x_vec][y_vec]+=1.;
//		vec_b[x_vec][y_vec]+=1.*(1.+((correct_pre2-b_pressure)*factor_pre/100.))*(1.+((correct_temp2-b_temp)*factor_temp/100.));
		vec_b[x_vec][y_vec]+=1.*(1.+((correct_pre2-b_pressure)*factor_pre/100.));
//		vec_b[x_vec][y_vec]+=1.;

count_n[jjjj]++;
count_np[jjjj]+=1.*(1.+((correct_pre2-b_pressure)*factor_pre/100.));
	}
pres_n[jjjj]=correct_pre2;
}


for(int i=0; i<31; i++)
	{
	for(int j=0; j<31; j++)
		{
		backd<<i<<"  "<<j<<"  "<<vec_b[i][j]<<endl;
		}
	}

double sumn=0.;
double sumnp=0.;

for(int i=0; i<filenum_b ; i++)
{
sumn+=count_n[i];
sumnp+=count_np[i];
}

double average_n=sumn/(double)filenum_b;
double average_np=sumnp/(double)filenum_b;

sumn=0.;
sumnp=0.;




for(int l=0; l<pix; l++)
{
	for(int k=0 ; k<pix; k++)
	{
		for(int i=0 ; i<31; i++)
		{
			for(int j=0 ; j<31; j++)
			{
			back_sum[i][j]+=psd_b[k][l][i][j];
			}
		}
	}
}

for(int i=0 ; i<31; i++)
	{
		for(int j=0 ; j<31; j++)
		{
		if(i>=15)
		{
			if(j>=15)
			{
				for(int k=0+abs(i-15) ; k<pix; k++)
				{
					for(int l=0+abs(j-15) ; l<pix; l++)
					{
					psd_b[k][l][i][j]=back_sum[i][j]/(((float)pix-(float)abs(i-15))*((float)pix-(float)abs(j-15)));
					}
				}
				}
			else
			{
				for(int k=0+abs(i-15) ; k<pix; k++)
				{
					for(int l=0 ; l<pix-abs(j-15); l++)
					{
					psd_b[k][l][i][j]=back_sum[i][j]/(((float)pix-(float)abs(i-15))*((float)pix-(float)abs(j-15)));
					}
				}			
			}
		}
		else
		{
			if(j>=15)
			{
				for(int k=0 ; k<pix-abs(i-15); k++)
				{
					for(int l=0+abs(j-15) ; l<pix; l++)
					{
					psd_b[k][l][i][j]=back_sum[i][j]/(((float)pix-(float)abs(i-15))*((float)pix-(float)abs(j-15)));
					}
				}
				}
			else
			{
				for(int k=0 ; k<pix-abs(i-15); k++)
				{
					for(int l=0 ; l<pix-abs(j-15); l++)
					{
					psd_b[k][l][i][j]=back_sum[i][j]/(((float)pix-(float)abs(i-15))*((float)pix-(float)abs(j-15)));
					}
				}
			
			}
		
		
		
		}


		}
		
	}


ofstream data("pres.txt");

for(int i=0; i<filenum_b ; i++)
{
sumn+=(count_n[i]-average_n)*(count_n[i]-average_n);
sumnp+=(count_np[i]-average_np)*(count_np[i]-average_np);

data<<i<<"   "<<time[i]<<"   "<<pres_n[i]<<"   "<<count_n[i]<<"   "<<count_np[i]<<endl;

}


cout<<" average_n "<<average_n <<"  sigme_n  "<<sqrt(sumn/(double)filenum_b)<<" average_np "<<average_np <<"  sigme_np  "<<sqrt(sumnp/(double)filenum_b)<<endl;


				free(pres);

//00だけ
TH2D *bg_map =new TH2D("bg_map","",16,-0.5,15.5,16,15.5,31.5);
TH2D *ob_map =new TH2D("ob_map","",16,-0.5,15.5,16,15.5,31.5);
TH2D *att_map =new TH2D("att_map","",16,-0.5,15.5,16,15.5,31.5);
TH2D *att_map_err =new TH2D("att_map_err","",16,-0.5,15.5,16,15.5,31.5);
TH2D *att_map_thick =new TH2D("att_map_thick","",16,-0.5,15.5,16,15.5,31.5);

//


///

TH2D *att_vec =new TH2D("att_vec","",31,-0.5,30.5,31,-0.5,30.5);//ベクトルベクトル
TH2D *att_bg =new TH2D("att_bg","",31,-0.5,30.5,31,-0.5,30.5);
TH2D *att_ob =new TH2D("att_ob","",31,-0.5,30.5,31,-0.5,30.5);
TH2D *att_err =new TH2D("att_err","",31,-0.5,30.5,31,-0.5,30.5);//誤差
TH2D *att_thick =new TH2D("att_thick","",31,-0.5,30.5,31,-0.5,30.5);//ベクトルベクトル







 TCanvas *cvs3=new TCanvas("attenuation","attenuation",900,900);

for(int i=0; i<16; i++)
	{
	for(int j=0; j<16; j++)
		{
				psd_att[i][j][15][15]=1-((psd_o[i][j][15][15]/o_time)/(psd_b[i][j][15][15]/b_time));
//				psd_err[i][j][15][15]=sqrt(pow(sqrt(psd_o[i][j][15][15])/psd_b[i][j][15][15],2)+pow((psd_o[i][j][15][15]*sqrt(psd_b[i][j][15][15]))/pow(psd_b[i][j][15][15],2),2));

				psd_err[i][j][15][15]=sqrt(pow(sqrt(psd_o[i][j][15][15])/psd_b[i][j][15][15],2)+pow((psd_o[i][j][15][15]*sqrt(psd_b[i][j][15][15]))/pow(psd_b[i][j][15][15],2),2))/((psd_o[i][j][15][15]/o_time)/(psd_b[i][j][15][15]/b_time));

				if(psd_att[i][j][15][15]<=0) psd_att[i][j][15][15]=0.0000001;
				psd_thick[i][j][15][15]=(1.14e-3)*pow(psd_att[i][j][15][15]*100.,3)-(4.56e-2)*pow(psd_att[i][j][15][15]*100.,2)+2.79*psd_att[i][j][15][15]*100.;

//cout<<psd_o[i][j][0][0]<<"   "<<psd_b[i][j][0][0]<<"  "<<psd_att[i][j][0][0]<<endl;
				bg_map->Fill(i,j+16,psd_b[i][j][15][15]);
				ob_map->Fill(i,j+16,psd_o[i][j][15][15]);
				att_map->Fill(i,j+16,psd_att[i][j][15][15]);
				att_map_err->Fill(i,j+16,psd_err[i][j][15][15]);
				att_map_thick->Fill(i,j+16,psd_thick[i][j][15][15]);

		}
	}
	att_map->Draw("colz");
 TCanvas *cvs10=new TCanvas("attenuation_err","attenuation_err",900,900);
	att_map_err->Draw("colz");

	TCanvas *cvs1=new TCanvas("bg","bg",900,900);
	bg_map->Draw("colz");

	TCanvas *cvs2=new TCanvas("ob","ob",900,900);
	ob_map->Draw("colz");

	TCanvas *cvs130=new TCanvas("thick","thick",900,900);
	att_map_thick->Draw("colz");

//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////

 TCanvas *cvs4=new TCanvas("attenuation_s","attenuation_s",900,900);

//cout<<psd_b[1][1][15][15]<<"  "<<psd_b[1][1][14][15]<<"  "<<psd_b[1][1][16][15]<<"  "<<psd_b[1][1][15][16]<<"  "<<psd_b[1][1][15][14]<<endl;

TH2D *bg_map_s =new TH2D("bg_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *ob_map_s =new TH2D("ob_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s =new TH2D("att_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_err =new TH2D("att_map_s_err","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_thick =new TH2D("att_map_s_thick","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);

/*
if(dive==1)
{
for(int i=0; i<pixs; i++)
{
	for(int j=0; j<pixs; j++)
	{
		for(int m=-1; m<2; m++)
			{
			for(int n=-1; n<2; n++)
				{
				if((abs(m)+abs(n))!=2)//(1,1????)
				{
				psd_b_dive[i][j]+=psd_b[i][j][15+n][15+m];
				psd_o_dive[i][j]+=psd_o[i][j][15+n][15+m];
				}
				{
				psd_b_dive[i][j]+=psd_b[i][j][15+n][15+m];
				psd_o_dive[i][j]+=psd_o[i][j][15+n][15+m];
				}

				}
			}
	}
}
}
*/

if(dive==1)
{

for(int i=0; i<pixs; i++)
{
	for(int j=0; j<pixs; j++)
	{
//				psd_b_dive[i][j]+=psd_b[i][j][15][12];
	//			psd_o_dive[i][j]+=psd_o[i][j][15][12];
				psd_b_dive[i][j]+=psd_b[i][j][15][15];
				psd_o_dive[i][j]+=psd_o[i][j][15][15];
				psd_b_dive[i][j]+=psd_b[i][j][15][14];
				psd_o_dive[i][j]+=psd_o[i][j][15][14];
		//		psd_b_dive[i][j]+=psd_b[i][j][15][13];
			//	psd_o_dive[i][j]+=psd_o[i][j][15][13];

//				psd_b_dive[i][j]+=psd_b[i][j][15][16];
	//			psd_o_dive[i][j]+=psd_o[i][j][15][16];

/*
				psd_b_dive[i][j]+=psd_b[i][j][15][15];
				psd_o_dive[i][j]+=psd_o[i][j][15][15];
				psd_b_dive[i][j]+=psd_b[i][j][14][15];
				psd_o_dive[i][j]+=psd_o[i][j][14][15];
				psd_b_dive[i][j]+=psd_b[i][j][16][15];
				psd_o_dive[i][j]+=psd_o[i][j][16][15];
				psd_b_dive[i][j]+=psd_b[i][j][13][15];
				psd_o_dive[i][j]+=psd_o[i][j][13][15];
*/
/*				psd_b_dive[i][j]+=psd_b[i][j][14][15];
				psd_o_dive[i][j]+=psd_o[i][j][14][15];
				psd_b_dive[i][j]+=psd_b[i][j][15][15];
				psd_o_dive[i][j]+=psd_o[i][j][15][15];
				psd_b_dive[i][j]+=psd_b[i][j][16][15];
				psd_o_dive[i][j]+=psd_o[i][j][16][15];
//				psd_b_dive[i][j]+=psd_b[i][j][15][17];
//				psd_o_dive[i][j]+=psd_o[i][j][15][17];*/
	}
}
}

else
{
for(int i=0; i<pixs; i++)
{
	for(int j=0; j<pixs; j++)
	{
		for(int k=0; k<dive; k++)
		{
			for(int l=0; l<dive; l++)
			{
				for(int m=-1; m<2; m++)
					{
					for(int n=-1; n<2; n++)
						{
						if((abs(m)+abs(n))!=2)
						{
						psd_b_dive[i][j]+=psd_b[i*2+k][j*2+l][15+n][15+m];
						psd_o_dive[i][j]+=psd_o[i*2+k][j*2+l][15+n][15+m];
						}
						}
					}
			}
		}
	}
}
}



for(int i=0; i<pixs; i++)

	{

	for(int j=0; j<pixs; j++)
		{
				psd_att_dive[i][j]=1-((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
//				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2));

				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2))/((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));

				if(psd_att_dive[i][j] <= 0) psd_att_dive[i][j]=0.0000001;

				psd_thick_dive[i][j]=(1.14e-3)*pow(100.*psd_att_dive[i][j],3)-(4.56e-2)*pow(100.*psd_att_dive[i][j],2)+2.79*psd_att_dive[i][j]*100.;


				att_map_s_err->Fill(i,j+pixs,psd_err_dive[i][j]);
				att_map_s->Fill(i,j+pixs,psd_att_dive[i][j]);
				ob_map_s->Fill(i,j+pixs,psd_o_dive[i][j]);
				bg_map_s->Fill(i,j+pixs,psd_b_dive[i][j]);
				att_map_s_thick->Fill(i,j+pixs,psd_thick_dive[i][j]);

		}

}

///まとめたりしたやつ


//////////////////////////xy
/*
TH2D *bg_map_s =new TH2D("bg_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *ob_map_s =new TH2D("ob_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s =new TH2D("att_map_s","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_err =new TH2D("att_map_s_err","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_thick =new TH2D("att_map_s_thick","",pixs,-0.5,pixs-0.5,pixs,pixs-0.5,2*pixs-0.5);


for(int i=0; i<pixs; i++)
{
	for(int j=0; j<pixs; j++)
	{
		for(int k=0; k<dive; k++)
		{
			for(int l=0; l<dive; l++)
			{
			psd_b_dive[i][j]+=psd_b[i*2+k][j*2+l][15][15];
			psd_o_dive[i][j]+=psd_o[i*2+k][j*2+l][15][15];
			}
		}
	}
}



for(int i=0; i<pixs; i++)

	{

	for(int j=0; j<pixs; j++)
		{
				psd_att_dive[i][j]=1-((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
//				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2));

				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2))/((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));

				if(psd_att_dive[i][j] <= 0) psd_att_dive[i][j]=0.0000001;

				psd_thick_dive[i][j]=(1.14e-3)*pow(100.*psd_att_dive[i][j],3)-(4.56e-2)*pow(100.*psd_att_dive[i][j],2)+2.79*psd_att_dive[i][j]*100.;


				att_map_s_err->Fill(i,j+pixs,psd_err_dive[i][j]);
				att_map_s->Fill(i,j+pixs,psd_att_dive[i][j]);
				ob_map_s->Fill(i,j+pixs,psd_o_dive[i][j]);
				bg_map_s->Fill(i,j+pixs,psd_b_dive[i][j]);
				att_map_s_thick->Fill(i,j+pixs,psd_thick_dive[i][j]);

		}
	}

/////////////////////x
/*
TH2D *bg_map_s =new TH2D("bg_map_s","",pixs,-0.5,pixs-0.5,16,15.5,31.5);
TH2D *ob_map_s =new TH2D("ob_map_s","",pixs,-0.5,pixs-0.5,16,15.5,31.5);
TH2D *att_map_s =new TH2D("att_map_s","",pixs,-0.5,pixs-0.5,16,15.5,31.5);
TH2D *att_map_s_err =new TH2D("att_map_s_err","",pixs,-0.5,pixs-0.5,16,15.5,31.5);
TH2D *att_map_s_thick =new TH2D("att_map_s_thick","",pixs,-0.5,pixs-0.5,16,15.5,31.5);


for(int i=0; i<pixs; i++)
{
	for(int j=0; j<pix; j++)
	{
		for(int k=0; k<dive; k++)
		{
			psd_b_dive[i][j]+=psd_b[i*2+k][j][15][15];
			psd_o_dive[i][j]+=psd_o[i*2+k][j][15][15];
		}
	cout<<psd_b_dive[i][j]<<endl;
	}
}



for(int i=0; i<pixs; i++)
	{
	for(int j=0; j<pix; j++)
		{
				psd_att_dive[i][j]=1-((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
//				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2));
				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2))/((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
				if(psd_att_dive[i][j] <= 0) psd_att_dive[i][j]=0.0000001;
				psd_thick_dive[i][j]=(1.14e-3)*pow(100.*psd_att_dive[i][j],3)-(4.56e-2)*pow(100.*psd_att_dive[i][j],2)+2.79*psd_att_dive[i][j]*100.;
				att_map_s_err->Fill(i,j+pix,psd_err_dive[i][j]);
				att_map_s->Fill(i,j+pix,psd_att_dive[i][j]);
				ob_map_s->Fill(i,j+pix,psd_o_dive[i][j]);
				bg_map_s->Fill(i,j+pix,psd_b_dive[i][j]);
				att_map_s_thick->Fill(i,j+pix,psd_thick_dive[i][j]);
		}
	}

*/
/////////////////////////y
/*
TH2D *bg_map_s =new TH2D("bg_map_s","",16,-0.5,15.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *ob_map_s =new TH2D("ob_map_s","",16,-0.5,15.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s =new TH2D("att_map_s","",16,-0.5,15.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_err =new TH2D("att_map_s_err","",16,-0.5,15.5,pixs,pixs-0.5,2*pixs-0.5);
TH2D *att_map_s_thick =new TH2D("att_map_s_thick","",16,-0.5,15.5,pixs,pixs-0.5,2*pixs-0.5);



for(int i=0; i<16; i++)
{
	for(int j=0; j<pixs; j++)
	{
			for(int l=0; l<dive; l++)
			{
			psd_b_dive[i][j]+=psd_b[i][j*2+l][15][15];
			psd_o_dive[i][j]+=psd_o[i][j*2+l][15][15];
			}
	}
}


//cout<<"a"<<endl;
for(int i=0; i<pix; i++)
	{
	for(int j=0; j<pixs; j++)
		{
				psd_att_dive[i][j]=1-((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
//				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2));
				psd_err_dive[i][j]=sqrt(pow(sqrt(psd_o_dive[i][j])/psd_b_dive[i][j],2)+pow((psd_o_dive[i][j]*sqrt(psd_b_dive[i][j]))/pow(psd_b_dive[i][j],2),2))/((psd_o_dive[i][j]/o_time)/(psd_b_dive[i][j]/b_time));
				if(psd_att_dive[i][j] <= 0) psd_att_dive[i][j]=0.0000001;
				psd_thick_dive[i][j]=(1.14e-3)*pow(100.*psd_att_dive[i][j],3)-(4.56e-2)*pow(100.*psd_att_dive[i][j],2)+2.79*psd_att_dive[i][j]*100.;
				att_map_s_err->Fill(i,j+pixs,psd_err_dive[i][j]);
				att_map_s->Fill(i,j+pixs,psd_att_dive[i][j]);
				ob_map_s->Fill(i,j+pixs,psd_o_dive[i][j]);
				bg_map_s->Fill(i,j+pixs,psd_b_dive[i][j]);
				att_map_s_thick->Fill(i,j+pixs,psd_thick_dive[i][j]);
		}
	}
*/
//////////////////////////////////


	att_map_s->Draw("colz");

 TCanvas *cvs5=new TCanvas("ob_s","ob_s",900,900);

	ob_map_s->Draw("colz");
 TCanvas *cvs6=new TCanvas("bg_s","bg_s",900,900);

	bg_map_s->Draw("colz");

 TCanvas *cvs11=new TCanvas("att_s_err","att_s_err",900,900);

	att_map_s_err->Draw("colz");

 TCanvas *cvs110=new TCanvas("att_s_thick","att_s_thick",900,900);

	att_map_s_thick->Draw("colz");


///////////////////////?x?N?g???x?N?g??
for(int i=0; i<31; i++)
	{
	for(int j=0; j<31; j++)

		{
				vec_att[i][j]=1-((vec_o[i][j]/o_time)/(vec_b[i][j]/b_time));
				vec_thick[i][j]=(1.14e-3)*pow(100.*vec_att[i][j],3)-(4.56e-2)*pow(100.*vec_att[i][j],2)+2.79*vec_att[i][j]*100.;
				vec_err[i][j]=sqrt(pow(sqrt(vec_o[i][j])/vec_b[i][j],2)+pow((vec_o[i][j]*sqrt(vec_b[i][j]))/pow(vec_b[i][j],2),2));
//				vec_err[i][j]=sqrt(pow(sqrt(vec_o[i][j])/vec_b[i][j],2)+pow((vec_o[i][j]*sqrt(vec_b[i][j]))/pow(vec_b[i][j],2),2))/((vec_o[i][j]/o_time)/(vec_b[i][j]/b_time));

				if(vec_att[i][j] <= 0) vec_att[i][j]=0.0000001;
				att_vec->Fill(i,j,vec_att[i][j]);
				att_err->Fill(i,j,vec_err[i][j]);
				att_bg->Fill(i,j,vec_b[i][j]);
				att_ob->Fill(i,j,vec_o[i][j]);
				att_thick->Fill(i,j,vec_thick[i][j]);

		}
	}
 TCanvas *cvs7=new TCanvas("att_vec","att_vec",900,900);
	att_vec->Draw("colz");
 TCanvas *cvs8=new TCanvas("att_err","att_err",900,900);
	att_err->Draw("colz");

 TCanvas *cvs9=new TCanvas("att_bg","att_bg",900,900);
	att_bg->Draw("colz");
 TCanvas *cvs20=new TCanvas("att_ob","att_ob",900,900);
	att_ob->Draw("colz");

 TCanvas *cvs120=new TCanvas("att_thick","att_thick",900,900);
	att_thick->Draw("colz");

//////////////////////
cout<<" bg  "<<filenum_b<<endl;
cout<<" ob  "<<filenum_o<<endl;


   return;

}



time_t  get_caltime(char *tstr, char *tstr2)//時間変換
{
        char    buf[100];
        char    buf2[100];
        time_t  cal_time;

        struct tm
                work_tm;

/*        if(strlen(tstr) != 14)
                return(-1);*/
        strncpy(buf, tstr, 4);
        buf[4] = '\0';
        work_tm.tm_year = atoi(buf) - 1900;//年
        strncpy(buf, tstr+5, 2);
        buf[2] = '\0';
        work_tm.tm_mon = atoi(buf) - 1;//月
        strncpy(buf, tstr+5+3, 2);
        work_tm.tm_mday = atoi(buf);//日

        if(strlen(tstr2) == 4)
          {
            strncpy(buf2, tstr2, 1);
            work_tm.tm_hour = atoi(buf2);
            strncpy(buf2, tstr2+2, 2);
            work_tm.tm_min = atoi(buf2);
          }
          else if(strlen(tstr2) == 5)
            {
              strncpy(buf2, tstr2, 2);
              work_tm.tm_hour = atoi(buf2);
              strncpy(buf2, tstr2+3, 2);
              work_tm.tm_min = atoi(buf2);
            }

            work_tm.tm_sec = 00;
            work_tm.tm_isdst = -1;
        if((cal_time = mktime(&work_tm)) == -1){
                return(-1);
        }
        return(cal_time);
}


bool StrString(const char*s1,const char *s2)//サブルーチンみたいなもの
{
  int n;
  n=strlen(s2);//strlenは文字列の長さを返す

  while(1){
    s1=strchr(s1,s2[0]);//s1の最初からs2[0]の文字を探す。その文字がなければNULLを返す

    if(s1==NULL)

      return false;

    if(strncmp(s1,s2,n)==0)//s1とs2を前からn文字分比べる、s1=s2の場合0を返す
      return true;
    s1++;
  }
}
