//void MyRunMuon_ZJets(double myflag=20112250, bool isQCD=false, int runflag=0)
void MyRunMuon_ZJets(double myflag=20121005, bool isQCD=false, int runflag=0)
//void MyRunMuon_ZJets(double myflag= -1, bool isQCD=false, int runflag=0)
{
    gSystem->Load("libFWCoreFWLite.so");
    gSystem->Load("libPhysicsToolsUtilities.so");
    gSystem->Load("libPhysicsToolsKinFitter.so");
    gSystem->Load("../../../../../lib/slc5_amd64_gcc462/libMMozerpowhegweight.so");
    gROOT->ProcessLine(".include ../../../../");
    gROOT->ProcessLine(".L Resolution.cc+");
    gROOT->ProcessLine(".L ../../src/METzCalculator.cc+");
    gROOT->ProcessLine(".L ../../src/QGLikelihoodCalculator.C+");
    gROOT->ProcessLine(".L EffTableReader.cc+");
    gROOT->ProcessLine(".L EffTableLoader.cc+");
    gROOT->ProcessLine(".L ClassifierOut/TMVAClassification_withqg_nJ2_mu_BDT.class.C+");
    gROOT->ProcessLine(".L ClassifierOut/TMVAClassification_withqg_nJ3_mu_BDT.class.C+");
    gROOT->ProcessLine(".L ClassifierOut/TMVAClassification_noqg_nJ2_mu_BDT.class.C+");
    gROOT->ProcessLine(".L ClassifierOut/TMVAClassification_noqg_nJ3_mu_BDT.class.C+");
    gROOT->ProcessLine(".L kanamuon_ZJets.cxx+");
    gROOT->ProcessLine("kanamuon runover");
    //Set true/false for isQCD
    char mycmd[500]; sprintf(mycmd,"runover.myana(%.d,%i,%i)",myflag, isQCD,runflag);
    cout << "running :: "<<mycmd << endl;
    gROOT->ProcessLine(mycmd);
}
