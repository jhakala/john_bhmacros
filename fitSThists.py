from ROOT import *
from sys import argv
STexcComparisons = []
STincComparisons = []
upperExcPads = []
upperIncPads = []
lowerExcPads = []
lowerIncPads = []
PlotsFile = TFile("METvsMHTandSTdists.root")
OutFile = TFile("FormattedFitPlots.root", "RECREATE")
for i in range(2,12):
    stExcMEThist=PlotsFile.Get("stExc%02iHist"%i)
    stExc2METhist=PlotsFile.Get("stExc02Hist")
    stIncMEThist=PlotsFile.Get("stInc%02iHist"%i)

    stExcMHThist=PlotsFile.Get("stExc%02iHistMHT"%i)
    stExc2MHThist=PlotsFile.Get("stExc02HistMHT")
    stIncMHThist=PlotsFile.Get("stInc%02iHistMHT"%i)

    if (i == 2):
        STexcComparisons.append(TCanvas("stExc%02iCanvas"%i, "ST, N=%i"%i, 800, 800))
        STexcComparisons[i-2].cd()
        upperExcPads.append(TPad("Exc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperExcPads[i-2].SetBottomMargin(0)
        upperExcPads[i-2].Draw()
        upperExcPads[i-2].cd()
        upperExcPads[i-2].SetLogy()
#        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
#        stExcMEThist.Draw()

        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMEThist.SetLineColor(kBlue)
        stExcMEThist.Draw()

        stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMHThist.SetLineColor(kRed)
        stExcMHThist.Draw("SAME")

        f1 = TF1("f1", "[0]/([1]*x)**[2]", 1000, 7000)
        f1.SetParameters(2.08e9, 4.28e-3, 6.7)
        for j in range(0, 20):
            stExcMEThist.Fit("f1", "", "", float(argv[1]), float(argv[2]) )
        f2 = TF1("f2", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
        f2.SetParameters(5e3, 2, 0.3, 0.2)
        for j in range(0, 20):
            stExcMEThist.Fit("f2", "", "", float(argv[1]), float(argv[2]) )
        f3 = TF1("f3", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
        f3.SetParameters(5e9, 3, 5e3, 0.6)
        for j in range(0, 20):
            stExcMEThist.Fit("f3", "", "", float(argv[1]), float(argv[2]) )
        f1.SetLineColor(kCyan)
        f1.Draw("SAME")
        f2.SetLineColor(kCyan)
        f2.Draw("SAME")
        f3.SetLineColor(kCyan)
        f3.Draw("SAME")

        f1MHT = TF1("f1MHT", "[0]/([1]*x)**[2]", 1000, 7000)
        f1MHT.SetParameters(2.08e9, 4.28e-3, 6.7)
        for j in range(0, 20):
            stExcMHThist.Fit("f1MHT", "", "", float(argv[1]), float(argv[2]) )
        f2MHT = TF1("f2MHT", "([0]*(1+x)^[1])/(x**([2]+[3]*TMath::Log(x)))", 1000, 7000)
        f2MHT.SetParameters(5e3, 2, 0.3, 0.2)
        for j in range(0, 20):
            stExcMHThist.Fit("f2MHT", "", "", float(argv[1]), float(argv[2]) )
        f3MHT = TF1("f3MHT", "[0]/([1] + [2]*x + x**2)**[3]", 1000, 7000)
        f3MHT.SetParameters(5e9, 3, 5e3, 0.6)
        for j in range(0, 20):
            stExcMHThist.Fit("f3MHT", "", "", float(argv[1]), float(argv[2]) )
        f1MHT.SetLineColor(kMagenta)
        f1MHT.Draw("SAME")
        f2MHT.SetLineColor(kMagenta)
        f2MHT.Draw("SAME")
        f2MHT.SetLineColor(kMagenta)
        f3MHT.Draw("SAME")


        STexcComparisons[i-2].cd()
        lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
        lowerExcPads[i-2].SetTopMargin(0)
        lowerExcPads[i-2].SetBottomMargin(0.2)
        lowerExcPads[i-2].Draw()
        lowerExcPads[i-2].cd()
        lowerExcPads[i-2].Draw()
        STexcComparisons[i-2].Write()

        STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
        STincComparisons[i-2].cd()
        upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperIncPads[i-2].Draw()
        upperIncPads[i-2].cd()
        upperIncPads[i-2].Draw()
        STincComparisons[i-2].cd()
        lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
        lowerIncPads[i-2].SetTopMargin(0)
        lowerIncPads[i-2].SetBottomMargin(0.2)
        lowerIncPads[i-2].Draw()
        STincComparisons[i-2].Write()

    else:
        STexcComparisons.append(TCanvas("stExc%02iCanvas"%i, "ST, N=%i"%i, 800, 800))
        STexcComparisons[i-2].cd()
        upperExcPads.append(TPad("Exc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperExcPads[i-2].SetBottomMargin(0)
        upperExcPads[i-2].Draw()
        upperExcPads[i-2].cd()
        upperExcPads[i-2].SetLogy()
        stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMEThist.Draw()
        stExcMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stExcMHThist.SetLineColor(kRed)
        stExcMHThist.Draw("SAME")

        STexcComparisons[i-2].cd()
        lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
        lowerExcPads[i-2].SetTopMargin(0)
        lowerExcPads[i-2].SetBottomMargin(0.2)
        lowerExcPads[i-2].Draw()
        lowerExcPads[i-2].cd()

        stExcMETRatio = stExcMEThist.Clone("stExcMETRatio")
        stExcMETRatio.Divide(stExc2METhist)
        stExcMETRatio.SetLineColor(kBlue)
        stExcMETRatio.Draw()
        stExcMHTRatio = stExcMHThist.Clone("stExcMHTRatio")
        stExcMHTRatio.Divide(stExc2MHThist)
        stExcMHTRatio.SetLineColor(kRed)
        stExcMHTRatio.Draw("SAME")
        STexcComparisons[i-2].Write()

        STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
        STincComparisons[i-2].cd()
        upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
        upperIncPads[i-2].SetBottomMargin(0)
        upperIncPads[i-2].Draw()
        upperIncPads[i-2].cd()
        upperIncPads[i-2].SetLogy()
        stIncMEThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMEThist.SetLineColor(kBlue)
        stIncMEThist.Draw()
        stIncMHThist.GetXaxis().SetRangeUser(1000, 7000)
        stIncMHThist.SetLineColor(kRed)
        stIncMHThist.Draw("SAME")

        STincComparisons[i-2].cd()
        lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
        lowerIncPads[i-2].SetTopMargin(0)
        lowerIncPads[i-2].SetBottomMargin(0.2)
        lowerIncPads[i-2].Draw()
        lowerIncPads[i-2].cd()

        stIncMETRatio = stIncMEThist.Clone("stIncMETRatio")
        stIncMETRatio.Divide(stExc2METhist)
        stIncMETRatio.SetLineColor(kBlue)
        stIncMETRatio.Draw()
        stIncMHTRatio = stIncMHThist.Clone("stIncMHTRatio")
        stIncMHTRatio.Divide(stExc2MHThist)
        stIncMHTRatio.SetLineColor(kRed)
        stIncMHTRatio.Draw("SAME")
        STincComparisons[i-2].Write()

#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]
