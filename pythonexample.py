import lhapdf
import uproot
import numpy as np
import math
import sys

##--NanoAOD file open
#events = uproot.open("/x6/cms/store/mc/RunIISummer19UL18NanoAODv2/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/280000/3FCF21DC-F295-754F-9FFF-557FF2195660.root")
events = uproot.open("/home/taebh/CMSSW_11_0_3/src/test/MiniAnalyzer/For_PDFuncert_M"+sys.argv[1]+"TeV.root")

##--get generator information
id1=events["Events"].array("pdg1")
x1=events["Events"].array("x1")
id2=events["Events"].array("pdg2")
x2=events["Events"].array("x2")
xpdf1=events["Events"].array("xpdf1")
xpdf2=events["Events"].array("xpdf2")
Q=events["Events"].array("scalePDF")

##--set LHAPDF path
lhapdf.pathsPrepend("/home/taebh/2021_Research/share/LHAPDF")

numEvents = len(id1)
nnpdf = lhapdf.mkPDF("NNPDF31_nnlo_as_0118_mc_hessian_pdfas",0)
PDFsets = lhapdf.mkPDFs("CT18NNLO")
numPDF = len(PDFsets)
reweight = np.empty([numPDF])

for j in range(numPDF):
    sumw = 0
    ct18 = PDFsets[j]
    for i in range(numEvents):
        sumw += (ct18.xfxQ(id1[i],x1[i],Q[i])*ct18.xfxQ(id2[i],x2[i],Q[i]))/(nnpdf.xfxQ(id1[i],x1[i],Q[i])*nnpdf.xfxQ(id2[i],x2[i],Q[i]))
    reweight[j] = sumw/numEvents

a = np.square(reweight[1:]-reweight[0])

PDFuncert=math.sqrt(a.sum())
print(PDFuncert)

Alphas1 = lhapdf.mkPDF("CT18NNLO_as_0120",0)
Alphas2 = lhapdf.mkPDF("CT18NNLO_as_0116",0)

sumw = 0
for i in range(numEvents):
    sumw += (Alphas1.xfxQ(id1[i],x1[i],Q[i])*Alphas1.xfxQ(id2[i],x2[i],Q[i]))/(nnpdf.xfxQ(id1[i],x1[i],Q[i])*nnpdf.xfxQ(id2[i],x2[i],Q[i]))
R_Alphas1 = sumw/numEvents

sumw = 0
for i in range(numEvents):
    sumw += (Alphas2.xfxQ(id1[i],x1[i],Q[i])*Alphas2.xfxQ(id2[i],x2[i],Q[i]))/(nnpdf.xfxQ(id1[i],x1[i],Q[i])*nnpdf.xfxQ(id2[i],x2[i],Q[i]))
R_Alphas2 = sumw/numEvents

R_Alphas = (R_Alphas1-R_Alphas2)/2

Alphasuncert = math.sqrt(PDFuncert**2+R_Alphas**2)
print(Alphasuncert)

PDFsets = lhapdf.mkPDFs("CT14nnlo")
numPDF = len(PDFsets)
reweight = np.empty([numPDF])

for j in range(numPDF):
    sumw = 0
    ct14 = PDFsets[j]
    for i in range(numEvents):
        sumw += (ct14.xfxQ(id1[i],x1[i],Q[i])*ct14.xfxQ(id2[i],x2[i],Q[i]))/(nnpdf.xfxQ(id1[i],x1[i],Q[i])*nnpdf.xfxQ(id2[i],x2[i],Q[i]))
    reweight[j] = sumw/numEvents

a = np.square(reweight[1:]-reweight[0])
print(math.sqrt(a.sum()))
