import pdf_calculation
import lhapdf

a = pdf_calculation.cal_pdfuncert()

a.set_para("/home/taebh/CMSSW_11_0_3/src/test/MiniAnalyzer/For_PDFuncert_M3TeV.root")
a.set_pdf("NNPDF31_nnlo_as_0118_mc_hessian_pdfas","CT18NNLO")
b = a.pdf_uncert()
print(b)
a.set_pdf("NNPDF31_nnlo_as_0118_mc_hessian_pdfas","CT14nnlo")
b = a.pdf_uncert()
print(b)

