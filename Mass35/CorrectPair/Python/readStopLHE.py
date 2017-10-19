from __future__ import print_function
from __future__ import division
import sys
import ROOT as rt
import math
from array import array
from LHEevent import *
from LHEfile import *
import plotTools

if __name__ == '__main__':

    # Name
    fileName = "plots"
    cutName = "cuts"
    cutsOn = False

    # Cuts
    ISO_MAX = 0.1
    LEAD_PT_MIN = 25
    PT_MIN = 5
    LEAD_ETA_MAX = 2.1
    ETA_MAX = 2.4
    M12_MIN = 12
    M4L_MIN = 80
    M4L_MAX = 100

    # Histograms
    MLL_BINS, MLL_LOW, MLL_UP = 48, 0, 120
    M4l_r = rt.TH1F("M4l_r", "M4l_r", 80, 50, 130)
    M12_r = rt.TH1F("M12_r", "M12_r", 100, MLL_LOW, 100)
    M34_r = rt.TH1F("M34_r", "M34_r", 100, MLL_LOW, 100)
    PT_BINS, PT_LOW, PT_UP = 60, 0, 60
    pT1_r = rt.TH1F("pT1_r", "pT1_r", PT_BINS, PT_LOW, PT_UP)
    pT2_r = rt.TH1F("pT2_r", "pT2_r", PT_BINS, PT_LOW, PT_UP)
    pT3_r = rt.TH1F("pT3_r", "pT3_r", PT_BINS, PT_LOW, PT_UP)
    pT4_r = rt.TH1F("pT4_r", "pT4_r", PT_BINS, PT_LOW, PT_UP)
    ETA_BINS, ETA_LOW, ETA_UP = 60, -3, 3
    eta1_r = rt.TH1F("eta1_r", "eta1_r", ETA_BINS, ETA_LOW, ETA_UP)
    eta2_r = rt.TH1F("eta2_r", "eta2_r", ETA_BINS, ETA_LOW, ETA_UP)
    eta3_r = rt.TH1F("eta3_r", "eta3_r", ETA_BINS, ETA_LOW, ETA_UP)
    eta4_r = rt.TH1F("eta4_r", "eta4_r", ETA_BINS, ETA_LOW, ETA_UP)
    PHI_BINS, PHI_LOW, PHI_UP = 64, -3.2, 3.2
    phi1_r = rt.TH1F("phi1_r", "phi1_r", PHI_BINS, PHI_LOW, PHI_UP)
    phi2_r = rt.TH1F("phi2_r", "phi2_r", PHI_BINS, PHI_LOW, PHI_UP)
    phi3_r = rt.TH1F("phi3_r", "phi3_r", PHI_BINS, PHI_LOW, PHI_UP)
    phi4_r = rt.TH1F("phi4_r", "phi4_r", PHI_BINS, PHI_LOW, PHI_UP)
    DELR_BINS, DELR_LOW, DELR_UP = 60, 0, 6
    DELR_BINS2, DELR_LOW2, DELR_UP2 = 60, 3, 9
    delr_r = rt.TH1F("delr_r", "delr_r", DELR_BINS2, DELR_LOW2, DELR_UP2)
    delr12 = rt.TH1F("delr12(Z)", "delr12(Z)", DELR_BINS, DELR_LOW, DELR_UP)
    delr34 = rt.TH1F("delr34(U)", "delr34(U)", DELR_BINS, DELR_LOW, DELR_UP)
    delr13 = rt.TH1F("delr13", "delr13", DELR_BINS, DELR_LOW, DELR_UP)
    delr14 = rt.TH1F("delr14", "delr14", DELR_BINS, DELR_LOW, DELR_UP)
    delr23 = rt.TH1F("delr23", "delr23", DELR_BINS, DELR_LOW, DELR_UP)
    delr24 = rt.TH1F("delr24", "delr24", DELR_BINS, DELR_LOW, DELR_UP)
    DPHI_BINS, DPHI_LOW, DPHI_UP = 64, 0, 3.2
    DPHI_BINS2, DPHI_LOW2, DPHI_UP2 = 100, 3.1, 3.2
    dphi_r = rt.TH1F("dphi_r", "dphi_r", DPHI_BINS2, DPHI_LOW2, DPHI_UP2)
    dphi12 = rt.TH1F("abs dphi12(Z)", "abs dphi12(Z)", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi34 = rt.TH1F("abs dphi34(U)", "abs dphi34(U)", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi13 = rt.TH1F("abs dphi13", "abs dphi13", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi14 = rt.TH1F("abs dphi14", "abs dphi14", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi23 = rt.TH1F("abs dphi23", "abs dphi23", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi24 = rt.TH1F("abs dphi24", "abs dphi24", DPHI_BINS, DPHI_LOW, DPHI_UP)

    p1_mother = rt.TH1I("pair 1 mother", "pair 1 mother", 22, 36, 14)
    p2_mother = rt.TH1I("pair 2 mother", "pair 2 mother", 22, 36, 14)

    # Lists to make TGraphs
    M12_arr_r, M34_arr_r = array('d'), array('d')
    pT1_arr_r, pT2_arr_r, pT3_arr_r, pT4_arr_r = array('d'), array('d'), array('d'), array('d')
    delr12_arr_r, delr34_arr_r = array('d'), array('d')
    dphi12_arr_r, dphi34_arr_r = array('d'), array('d')
    phi12_arr_r, phi34_arr_r = array('d'), array('d')
    gam12_arr_r, gam34_arr_r = array('d'), array('d')

    # Loop over events
    myLHEfile = LHEfile("../unweighted_events.lhe")
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    n_acc = 0
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        accepted = False
        particles, mus = [], []
        Z_mu, U_mu = [], []
        Z_mu_q, U_mu_q = [], []
        mu_mother, mu_q = {}, {}
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            particles.append(p)
        for p in particles:
            if abs(p['ID']) == 13 and p['mIdx'] == p['mIdx2']:
                mu_mother = particles[p['mIdx']]
                if abs(mu_mother['ID']) == 23:
                    Z_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                    Z_mu_q.append(p['ID'])
                elif abs(mu_mother['ID']) == 35:
                    U_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                    U_mu_q.append(p['ID'])
        if Z_mu[0].Pt() < Z_mu[1].Pt():
            Z_mu[0], Z_mu[1] = Z_mu[1], Z_mu[0]
            Z_mu_q[0], Z_mu_q[1] = Z_mu_q[1], Z_mu_q[0]
        if U_mu[0].Pt() < U_mu[1].Pt():
            U_mu[0], U_mu[1] = U_mu[1], U_mu[0]
            U_mu_q[0], U_mu_q[1] = U_mu_q[1], U_mu_q[0]
        mus.append(Z_mu[0])
        mus.append(Z_mu[1])
        mus.append(U_mu[0])
        mus.append(U_mu[1])

        if cutsOn:
            if (mus[0].Pt() > LEAD_PT_MIN and mus[1].Pt() > LEAD_PT_MIN
                and mus[2].Pt() > PT_MIN and mus[3].Pt() > PT_MIN
                and abs(mus[0].Eta()) < LEAD_ETA_MAX and abs(mus[1].Eta()) < LEAD_ETA_MAX
                and abs(mus[2].Eta()) < ETA_MAX and abs(mus[3].Eta()) < ETA_MAX
                and (mus[0] + mus[1]).M() > M12_MIN
                and M4L_MIN < (mus[0] + mus[1] + mus[2] + mus[3]).M() < M4L_MAX
                ):
                accepted = True
        else:
            if (M4L_MIN < (mus[0] + mus[1] + mus[2] + mus[3]).M() < M4L_MAX
                and abs(mus[0].Eta()) < ETA_MAX and abs(mus[1].Eta()) < ETA_MAX
                and abs(mus[2].Eta()) < ETA_MAX and abs(mus[3].Eta()) < ETA_MAX
                ):
                accepted = True

        if accepted:
            n_acc += 1
            M4l_r.Fill((mus[0] + mus[1] + mus[2] + mus[3]).M())
            M12_r.Fill((mus[0] + mus[1]).M())
            M34_r.Fill((mus[2] + mus[3]).M())
            pT1_r.Fill(mus[0].Pt())
            pT2_r.Fill(mus[1].Pt())
            pT3_r.Fill(mus[2].Pt())
            pT4_r.Fill(mus[3].Pt())
            eta1_r.Fill(mus[0].Eta())
            eta2_r.Fill(mus[1].Eta())
            eta3_r.Fill(mus[2].Eta())
            eta4_r.Fill(mus[3].Eta())
            phi1_r.Fill(mus[0].Phi())
            phi2_r.Fill(mus[1].Phi())
            phi3_r.Fill(mus[2].Phi())
            phi4_r.Fill(mus[3].Phi())

            M12_arr_r.append((mus[0] + mus[1]).M())
            M34_arr_r.append((mus[2] + mus[3]).M())
            pT1_arr_r.append(mus[0].Pt())
            pT2_arr_r.append(mus[1].Pt())
            pT3_arr_r.append(mus[2].Pt())
            pT4_arr_r.append(mus[3].Pt())
            delr12_arr_r.append(mus[0].DeltaR(mus[1]))
            delr34_arr_r.append(mus[2].DeltaR(mus[3]))
            dphi12_arr_r.append(abs(mus[0].DeltaPhi(mus[1])))
            dphi34_arr_r.append(abs(mus[2].DeltaPhi(mus[3])))
            phi12_arr_r.append((mus[0] + mus[1]).Phi())
            phi34_arr_r.append((mus[2] + mus[3]).Phi())
            gam12_arr_r.append((mus[0] + mus[1]).Gamma())
            gam34_arr_r.append((mus[2] + mus[3]).Gamma())

            delr_r.Fill((mus[0] + mus[1]).DeltaR(mus[2] + mus[3]))
            dphi_r.Fill((mus[0] + mus[1]).DeltaPhi(mus[2] + mus[3]))
        del oneEvent, myLHEevent

    # Create dummy graph for line
    line_arr = array('d')
    line_arr.append(0)
    line_arr.append(100)
    graph_line = rt.TGraph(2, line_arr, line_arr)
    graph_line.SetLineColor(rt.kBlack)
    graph_line.SetLineWidth(2)

    # Create scatter plots
    mSize = 0.5
    canvas_Mll = rt.TCanvas("Mll", "Mll", 800, 600)
    canvas_Mll.cd()
    graph_Mll_r = rt.TGraph(n_acc, M12_arr_r, M34_arr_r)
    graph_Mll_r.SetName("graph_Mll")
    graph_Mll_r.SetMarkerStyle(rt.kFullCircle)
    graph_Mll_r.SetMarkerSize(mSize)
    graph_Mll_r.Draw("AP")
    graph_line.Draw("L")

    canvas_pT23 = rt.TCanvas("pT23", "pT23", 800, 600)
    canvas_pT23.cd()
    graph_pT23_r = rt.TGraph(n_acc, pT2_arr_r, pT3_arr_r)
    graph_pT23_r.SetName("graph_pT23")
    graph_pT23_r.SetMarkerStyle(rt.kFullCircle)
    graph_pT23_r.SetMarkerSize(mSize)
    graph_pT23_r.Draw("AP")
    graph_line.Draw("L")

    canvas_delr = rt.TCanvas("DeltaR", "DeltaR", 800, 600)
    canvas_delr.cd()
    graph_delr_r = rt.TGraph(n_acc, delr12_arr_r, delr34_arr_r)
    graph_delr_r.SetName("graph_delr")
    graph_delr_r.GetXaxis().SetTitle('#DeltaR_{12}(Z)')
    graph_delr_r.GetYaxis().SetTitle('#DeltaR_{34}(U)')
    graph_delr_r.SetMarkerStyle(rt.kFullCircle)
    graph_delr_r.SetMarkerSize(mSize)
    graph_delr_r.Draw("AP")
    graph_line.Draw("L")

    canvas_dphi = rt.TCanvas("DeltaPhi", "DeltaPhi", 800, 600)
    canvas_dphi.cd()
    graph_dphi_r = rt.TGraph(n_acc, dphi12_arr_r, dphi34_arr_r)
    graph_dphi_r.SetName("graph_dphi")
    graph_dphi_r.GetXaxis().SetTitle('|#Delta#phi_{12}|(Z)')
    graph_dphi_r.GetYaxis().SetTitle('|#Delta#phi_{34}|(U)')
    graph_dphi_r.SetMarkerStyle(rt.kFullCircle)
    graph_dphi_r.SetMarkerSize(mSize)
    graph_dphi_r.Draw("AP")
    graph_line.Draw("L")

    canvas_gamma = rt.TCanvas("Gamma", "Gamma", 800, 600)
    canvas_gamma.cd()
    graph_gamma_r = rt.TGraph(n_acc, gam12_arr_r, gam34_arr_r)
    graph_gamma_r.SetName("graph_gamma")
    graph_gamma_r.GetXaxis().SetTitle('#gamma_{12}(Z)')
    graph_gamma_r.GetYaxis().SetTitle('#gamma_{34}(U)')
    graph_gamma_r.SetMarkerStyle(rt.kFullCircle)
    graph_gamma_r.SetMarkerSize(mSize)
    graph_gamma_r.Draw("AP")
    graph_line.Draw("L")

    # Write hists and graphs
    if cutsOn:
        fileName = fileName + "_" + cutName
    fileName = fileName + ".root"
    outFile = rt.TFile(fileName, "RECREATE")
    outFile.mkdir("Histograms")
    outFile.cd("Histograms")
    p1_mother.Write()
    p2_mother.Write()

    delr_r.Write()
    dphi_r.Write()
    M4l_r.Write()
    M12_r.Write()
    M34_r.Write()
    pT1_r.Write()
    pT2_r.Write()
    pT3_r.Write()
    pT4_r.Write()
    eta1_r.Write()
    eta2_r.Write()
    eta3_r.Write()
    eta4_r.Write()
    phi1_r.Write()
    phi2_r.Write()
    phi3_r.Write()
    phi4_r.Write()

    outFile.cd()
    outFile.mkdir("Scatterplots")
    outFile.cd("Scatterplots")
    canvas_Mll.Write()
    canvas_pT23.Write()
    canvas_delr.Write()
    canvas_dphi.Write()
    canvas_gamma.Write()
#   canvas_phi.Write()

    graph_Mll_r.Write()
    graph_pT23_r.Write()
    graph_delr_r.Write()
    graph_dphi_r.Write()
    graph_gamma_r.Write()

#   canvas_ZpT.Write()
#   canvas_UpT.Write()
#   canvas_pT13.Write()
#   canvas_pT14.Write()
#   canvas_pT23.Write()
#   canvas_pT24.Write()
#   canvas_Zeta.Write()
#   canvas_Ueta.Write()
#   outFile.cd()
#   outFile.mkdir("deltaR")
#   outFile.cd("deltaR")
#   delr12.Write()
#   delr34.Write()
#   delr13.Write()
#   delr14.Write()
#   delr23.Write()
#   delr24.Write()
#   outFile.cd()
#   outFile.mkdir("deltaPhi")
#   outFile.cd("deltaPhi")
#   dphi12.Write()
#   dphi34.Write()
#   dphi13.Write()
#   dphi14.Write()
#   dphi23.Write()
#   dphi24.Write()
    outFile.cd()
    outFile.Close()

    # Print info
    print("Accepted", n_acc, "events (of 100000)")
    print("\nCreated file", fileName)


#   del M4l, M12, M34, pT1, pT2, pT3, pT4, eta1, eta2, eta3, eta4, phi1, phi2, phi3, phi4
#   del delr12, delr34, delr13, delr14, delr23, delr24
#   del graph_Mll, graph_ZpT, graph_UpT, graph_pT13, graph_pT14, graph_pT23, graph_pT24, graph_Zeta, graph_Ueta, graph_delr, graph_dphi
#   del canvas_Mll, canvas_ZpT, canvas_UpT, canvas_pT13, canvas_pT14, canvas_pT23, canvas_pT24, canvas_Zeta, canvas_Ueta, canvas_delr, canvas_dphi
