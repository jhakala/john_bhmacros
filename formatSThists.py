from ROOT import *
STexcComparisons = []
STincComparisons = []
upperExcPads = []
upperIncPads = []
lowerExcPads = []
lowerIncPads = []
PlotsFile = TFile("METvsMHTandSTdists.root")
OutFile = TFile("FormattedSTratioPlots.root", "RECREATE")
for i in range(2,12):
  stExcMEThist=PlotsFile.Get("stExc%02iHist"%i)
  stIncMEThist=PlotsFile.Get("stInc%02iHist"%i)
  stExcMHThist=PlotsFile.Get("stExc%02iHistMHT"%i)
  stIncMHThist=PlotsFile.Get("stInc%02iHistMHT"%i)

  STexcComparisons.append(TCanvas("stExc%02iCanvas"%i, "ST, N=%i"%i, 800, 800))
  STexcComparisons[i-2].cd()
  upperExcPads.append(TPad("Exc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
  upperExcPads[i-2].SetBottomMargin(0)
  upperExcPads[i-2].Draw()
  upperExcPads[i-2].cd()
  upperExcPads[i-2].SetLogy()
  stExcMEThist.GetXaxis().SetRangeUser(1000, 7000)
  stExcMEThist.Draw()
  stExcMHThist.SetLineColor(kRed)
  stExcMHThist.Draw("SAMES")

  STexcComparisons[i-2].cd()
  lowerExcPads.append(TPad("Exc%02iratiopad"%i, "ratiopad1", 0, 0.04, 1, 0.3))
  lowerExcPads[i-2].SetTopMargin(0)
  lowerExcPads[i-2].SetBottomMargin(0.2)
  lowerExcPads[i-2].Draw()
  lowerExcPads[i-2].cd()

  stExcRatio = stExcMEThist.Clone("stExcRatio")
  stExcRatio.Divide(stExcMHThist)
  stExcRatio.Draw()
  STexcComparisons[i-2].Write()

  STincComparisons.append(TCanvas("stInc%02iCanvas"%i, "ST, N>=%i"%i, 800, 800))
  upperIncPads.append(TPad("Inc%02ipad"%i, "pad1", 0, 0.3, 1, 1.0))
  upperIncPads[i-2].SetBottomMargin(0)
  upperIncPads[i-2].Draw()
  upperIncPads[i-2].cd()
  upperIncPads[i-2].SetLogy()
  stIncMEThist.GetXaxis().SetRangeUser(1000, 7000)
  stIncMEThist.Draw()
  stIncMHThist.SetLineColor(kRed)
  stIncMHThist.Draw("SAMES")

  STincComparisons[i-2].cd()
  lowerIncPads.append(TPad("Inc%02iratiopad"%i, "incratiopad1", 0, 0.04, 1, 0.3))
  lowerIncPads[i-2].SetTopMargin(0)
  lowerIncPads[i-2].SetBottomMargin(0.2)
  lowerIncPads[i-2].Draw()
  lowerIncPads[i-2].cd()

  stIncRatio = stIncMEThist.Clone("stIncRatio")
  stIncRatio.Divide(stIncMHThist)
  stIncRatio.Draw()

  STincComparisons[i-2].Write()

#gStyle.SetOptTitle("E_{T}^{/mu} / E_{T}^{jet}")
#gStyle.SetOptTitle(0)
#MuonJetC.SetLogy()
#MuonJetC.cd()
#MuonJet1 = isoPlotsFile.Get("MuonJetIso1")
#MuonJet2 = isoPlotsFile.Get("MuonJetIso2")
#MuonJet3 = isoPlotsFile.Get("MuonJetIso3")
#MuonJet4 = isoPlotsFile.Get("MuonJetIso4")
#MuonJet1.Rebin(2)
#MuonJet2.Rebin(2)
#MuonJet3.Rebin(2)
#MuonJet4.Rebin(2s
#MuonJet1.GetXaxis().SetRangeUser(0,1.05)
#MuonJet1.GetXaxis().SetTitle("E_{T}^{#mu}/E_{T}^{jet}")
#MuonJet1.GetYaxis().SetTitle("Events")
#MuonJet1.SetLineColor(kBlue)
#MuonJet2.SetLineColor(kGreen-2)
#MuonJet3.SetLineColor(kRed)
#MuonJet4.SetLineColor(kBlack)
#MuonJet1.Draw()
#gPad.Update()
#st1 = MuonJet1.FindObject("stats")
#st1.SetX1NDC(0.660804)
#st1.SetX2NDC(0.771357)
#st1.SetY1NDC(0.901682)
#st1.SetY2NDC(0.994825)
#MuonJetC.Update()
#MuonJet2.Draw("SAMES")
#gPad.Update()
#st2 = MuonJet2.FindObject("stats")
#st2.SetX1NDC(0.771357)
#st2.SetX2NDC(0.901682)
#st2.SetY1NDC(0.901682)
#st2.SetY2NDC(0.994825)
#MuonJetC.Update()
#MuonJet3.Draw("SAMES")
#gPad.Update()
#st3 = MuonJet3.FindObject("stats")
#st3.SetX1NDC(0.660804)
#st3.SetX2NDC(0.771357)
#st3.SetY1NDC(0.901682)
#st3.SetY2NDC(0.8085391)
#MuonJetC.Update()
#MuonJet4.Draw("SAMES")
#gPad.Update()
#st4 = MuonJet4.FindObject("stats")
#st4.SetX1NDC(0.771357)
#st4.SetX2NDC(0.901682)
#st4.SetY1NDC(0.901682)
#st4.SetY2NDC(0.8085391)
#gPad.Update()
#MuonJetC.Update()
#legend = TLegend(0.5, 0.85, 0.88, 0.65, "Overlapping E_{T}^{#mu} / E_{T}^{jet} for E_{T}^{jet} ranges:", "brNDC")
#legend.SetTextSize(0.02);
#legend.SetLineColor(1);
#legend.SetLineStyle(1);
#legend.SetLineWidth(1);
#legend.SetFillStyle(1001);
#legend.SetFillColor(10);
#legend.AddEntry(MuonJet1,"50 < Jet E_{T} < 150","lp");
#legend.AddEntry(MuonJet2,"150 < Jet E_{T} < 250","lp");
#legend.AddEntry(MuonJet3,"250 < Jet E_{T} < 400","lp");
#legend.AddEntry(MuonJet4,"400 < Jet E_{T}","lp");
#legend.Draw()
#gPad.Update()
#MuonJetC.Update()
#legend.SetX1NDC(0.278894)
#legend.SetX2NDC(0.658291)
#legend.SetY1NDC(0.829237)
#legend.SetY2NDC(0.997413)
#gPad.Update()
#MuonJetC.Update()
#gStyle.SetOptTitle(0)
#gPad.Update()
#MuonJetC.Update()
#
#outfile = TFile("IsoPlotsOutFile.root", "RECREATE")
#MuonJetC.Write()
#MuonJetC.SaveAs("MuonIsoPlot.pdf")

if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
      rep = raw_input( 'enter "q" to quit: ' )
      if 1 < len(rep):
         rep = rep[0]
