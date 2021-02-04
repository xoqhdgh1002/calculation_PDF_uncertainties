import lhapdf
import uproot
import numpy as np
import math

class cal_pdfuncert:
    #def __init__(self):

    def set_para(self,file_path):

        events = uproot.open(file_path)["MiniAnalyzer"]
        
        ##--get generator information
        self.id1=events["Events"].array("pdg1")
        self.x1=events["Events"].array("x1")
        self.id2=events["Events"].array("pdg2")
        self.x2=events["Events"].array("x2")
        self.xpdf1=events["Events"].array("xpdf1")
        self.xpdf2=events["Events"].array("xpdf2")
        self.Q=events["Events"].array("scalePDF")
        
        self.numEvents = len(self.id1)

    def set_pdf(self,refer,target):

        ##--set LHAPDF path
        lhapdf.pathsPrepend("/home/taebh/2021_Research/share/LHAPDF")
        self.nnpdf = lhapdf.mkPDF(refer,0)
        self.PDFsets = lhapdf.mkPDFs(target)
        self.numPDF = len(self.PDFsets)

    def pdf_uncert(self):

        reweight = np.empty([self.numPDF])
        
        for j in range(self.numPDF):
            sumw = 0
            ct18 = self.PDFsets[j]
            for i in range(self.numEvents):
                sumw += (ct18.xfxQ(self.id1[i],self.x1[i],self.Q[i])*ct18.xfxQ(self.id2[i],self.x2[i],self.Q[i]))/(self.nnpdf.xfxQ(self.id1[i],self.x1[i],self.Q[i])*self.nnpdf.xfxQ(self.id2[i],self.x2[i],self.Q[i]))
            reweight[j] = sumw/self.numEvents
        
        a = np.square(reweight[1:]-reweight[0])
        
        self.PDF_uncert = math.sqrt(a.sum())
        return self.PDF_uncert

    def as_uncert(self,target1,target2):

        Alphas1 = lhapdf.mkPDF(target1,0)
        Alphas2 = lhapdf.mkPDF(target2,0)
        
        sumw = 0
        for i in range(self.numEvents):
            sumw += (Alphas1.xfxQ(self.id1[i],self.x1[i],self.Q[i])*Alphas1.xfxQ(self.id2[i],self.x2[i],self.Q[i]))/(self.nnpdf.xfxQ(self.id1[i],self.x1[i],self.Q[i])*self.nnpdf.xfxQ(self.id2[i],self.x2[i],self.Q[i]))
        R_Alphas1 = sumw/self.numEvents
        
        sumw = 0
        for i in range(self.numEvents):
            sumw += (Alphas2.xfxQ(self.id1[i],self.x1[i],self.Q[i])*Alphas2.xfxQ(self.id2[i],self.x2[i],self.Q[i]))/(self.nnpdf.xfxQ(self.id1[i],self.x1[i],self.Q[i])*self.nnpdf.xfxQ(self.id2[i],self.x2[i],self.Q[i]))
        R_Alphas2 = sumw/self.numEvents
        
        self.As_uncert = (R_Alphas1-R_Alphas2)/2
        return self.As_uncert


    def pdf_as_uncert(self):
        
        self.PDF_As_uncert = math.sqrt(self.PDF_uncert**2+self.As_uncert**2)
        return self.PDF_As_uncert
