#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include <iosfwd>
#include <ios>
#include <ostream>

bool StrString(const char *s1,const char *s2);//文字列を比較する関数

int main()
{
  std::stringstream ss1;//文字列操作のためのクラス stringstream
  const char* path="./";//現在のフォルダのパスを格納
  DIR *dp;//ディレクトリ内のファイルを受け取るための変数型　DIR
  dirent *entry;//ファイル情報を受け取るための変数型
  std::string BGfilename[4096];//先頭に"BG"が書かれているファイル名を受け取る
  std::string OBfilename[4096];//先頭に"OB"が書かれているファイル名を受け取る
  int BGfilenum=0;//ファイルの個数
  int OBfilenum=0;//ファイルの個数

  dp = opendir(path);//dpに現在のディレクトリのファイル情報を渡す
  while((entry=readdir(dp))!=NULL)//一つずつファイル読んでいく
  {
    std::string tmp;
    tmp=entry -> d_name;//ファイル名をtmpに格納

    if(StrString(tmp.c_str(),"BG"))//BGが含まれているならtrueを返す
    {
      BGfilename[BGfilenum]=tmp;
      BGfilenum++;
    }
    else if(StrString(tmp.c_str(),"OB"))//OBが含まれているならtrueを返す
    {
      OBfilename[OBfilenum]=tmp;
      OBfilenum++;
    }
  }
  //ここまでで、filenameにはBG,OBがファイル名に含まれるファイル名が格納されている

  std::cout<<"BGfilenum = "<<BGfilenum<<std::endl;//ファイルの数
  std::cout<<"OBfilenum = "<<OBfilenum<<std::endl;

  std::ofstream BG_of("0_BackGround_output.txt");//書き出すファイル
  std::ofstream OB_of("0_Object_output.txt");

  int i,j,k;
  int x,y = 0;
  int count[31][31];
  int BG[962];
  int OB[962];
  int line=0;
  double BG_time_buf[4096];
  double OB_time_buf[4096];
  double BG_time=0.;//計測時間を足し合わせていく
  double OB_time=0.;

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////配列の初期化
  for(i=0; i<962; i++)//初期化
  {
    BG[i]=0;
    OB[i]=0;
  }

  for(i=0; i<31; i++)//初期化
  {
    for(j=0; j<31; j++)
    {
      count[i][j]=0;
    }
  }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////ファイルの読み込みと足し合わせ
  for(i=0; i<BGfilenum; i++)//BG
  {
    std::string tmpBG;
    tmpBG = BGfilename[i];
    std::cout<<BGfilename[i]<<std::endl;//読み込むファイル名を出力
    std::ifstream BG_if(tmpBG);//ファイルを読み込む
    BG_if >> x >> y >> BG_time_buf[i];//各ファイルの計測時間をbufに格納しておく

    for(j=0; j<31; j++)
    {
      for(k=0; k<31; k++)
      {
        x = 0;
        y = 0;
        BG_if >> x >> y >> count[x+15][y+15];//スペース区切りで入力
        BG[line]+=count[x+15][y+15];//カウント数を足し合わせていく
        line++;//行とlineを対応させてBG[line]に足し合わせていく
      }
    }
    line=0;//初期化
  }
  for(i=0; i<OBfilenum; i++)//OB
  {
    std::string tmpOB;
    tmpOB = OBfilename[i];
    std::cout<<OBfilename[i]<<std::endl;
    std::ifstream OB_if(tmpOB);
    OB_if >> x >> y >> OB_time_buf[i];

    for(j=0; j<31; j++)
    {
      for(k=0; k<31; k++)
      {
        x = 0;
        y = 0;
        OB_if >> x >> y >> count[x+15][y+15];
        OB[line]+=count[x+15][y+15];
        line++;
      }
    }
    line=0;
  }

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////出力ファイルの一行目に時間を書き出す
  for(i=0; i<BGfilenum; i++)
  {
    BG_time+=BG_time_buf[i];
  }
  for(i=0; i<OBfilenum; i++)
  {
    OB_time+=OB_time_buf[i];
  }
  BG_of<<100<<" "<<100<<" "<<BG_time<<std::endl;
  OB_of<<100<<" "<<100<<" "<<OB_time<<std::endl;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////出力ファイルを書き出す
  for(i=-15; i<16; i++)
  {
    for(j=-15; j<16; j++)
    {
      BG_of<<i<<" "<<j<<" "<<BG[line]<<std::endl;//1_BG_output.txtに書き込み
      OB_of<<i<<" "<<j<<" "<<OB[line]<<std::endl;//2_OB_output.txtに書き込み
      line++;
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////文字列の一致を判定する関数
bool StrString(const char*s1,const char *s2)
{
  int n;
  n=strlen(s2);//strlenは文字列の長さを返す
  while(1)
  {
    s1=strchr(s1,s2[0]);//s1の一文字目へのpointa,s1とs2[0]が違うとNULLを返す、拡張子を調べてる。
    if(s1==NULL)
      return false;
    if(strncmp(s1,s2,n)==0)//s1とs2を前からn文字分比べる、s1=s2の場合0を返す
      return true;
    s1++;
  }
}
