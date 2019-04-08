g++ -L/home/yxiang/root_version6/lib -g -std=c++1y  Controll.C makeHists.C Analysis.C -I$ROOTSYS/include `root-config --libs ` -lMinuit -o Analysis.exe
./Analysis.exe Filelist.list MC
#./zeeAna.exe /home/cwang/ui2-clued0/cwang/Charge-Dependent_Scale_Resolution/MC_1.list MC
#mv MC_results.root MCIIa1560.root
#./zeeAna.exe filelist/ui01/MCFileIIa60130.list MC
#mv MC_results.root MCIIa60130.root
#./zeeAna.exe filelist/ui01/MCFileIIa130250.list MC
#mv MC_results.root MCIIa130250.root 
#./zeeAna.exe /home/yxiang/work/filelist/AfterPreskimmed/DataIIb.list data
#mv data_results.root DataIIb.root


#./zeeAna.exe filelist/DataFileIIa.list data
#mv data_results.root data_results_s1.root
#sed -i 's/solenoid_P == -1/solenoid_P == 1/g' zeeAna1.C
#./zeeAna.exe filelist/DataFileIIa.list data
#mv data_results.root data_results_IIa.list
#./zeeAna.exe filelist/DataFile.list data
#mv data_results.root data_results_IIb.list
#gcc -g  Controll.C makeHists.C zeeAna1.C -I$ROOTSYS/include `root-config --libs ` -lMinuit -o zeeAna.exe
#./zeeAna.exe filelist/DataFile.list data
#mv data_results.root data_results_s0.root
