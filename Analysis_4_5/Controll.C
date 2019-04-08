#include "Analysis.h"
#include "makeHists.h"
#include <iostream>
#include <fstream>

using namespace std;
int main(int argc, char** argv)
{
//usage output
  if(argc != 3){
    std::cout<<"usage : "<<argv[0]<<" ntuple_file_name type_of_samples  "<<std::endl;
    return 0;
  }

//read in filelist
  ifstream infile;
   infile.open(argv[1],ios::in);
  int count=0;

  cout<<"*************************"<<endl;
  cout<<"*** Analysis Begin ****"<<endl;
  cout<<"*************************"<<endl;
  cout<<endl;
  Analysis Run(argv[1], argv[2]);
  TString InputRoot;
//infile>>InputRoot;
  while(infile>>InputRoot)
   {
     count++;
cout<<" Initial Begin: "<<endl;
     Run.Initial(InputRoot,count);
cout<<" Loop    Begin:"<<endl;
     Run.Loop(argv[2]);
     Run.End(count);
   }

  Run.Output();
  Run.Finish(count);
  cout<<endl;
  cout<<"*************************"<<endl;
  cout<<"**** zeeAna End ****"<<endl;
  cout<<"*************************"<<endl;


  return 1;
}
