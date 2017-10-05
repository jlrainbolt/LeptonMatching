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
    cutsOn = False

    # Cuts
    ISO_MAX = 0.1
    LEAD_PT_MIN = 10
    PT_MIN = 5
    LEAD_ETA_MAX = 2.1
    ETA_MAX = 2.4
    M12_MIN = 12
    M4L_MIN = 80
    M4L_MAX = 100

    # Histograms
    MLL_BINS, MLL_LOW, MLL_UP = 48, 0, 120
    M4l = rt.TH1F("M4l", "M4l", 80, 50, 130)
#   MZ = rt.TH1F("MZ", "MZ", 80, 50, 130)
    M12 = rt.TH1F("M12", "M12", 100, MLL_LOW, 100)
    M34 = rt.TH1F("M34", "M34", 80, MLL_LOW, 40)
    PT_BINS, PT_LOW, PT_UP = 60, 0, 60
    pT1 = rt.TH1F("pT1(Z)", "pT1(Z)", PT_BINS, PT_LOW, PT_UP)
    pT2 = rt.TH1F("pT2(Z)", "pT2(Z)", PT_BINS, PT_LOW, PT_UP)
    pT3 = rt.TH1F("pT3(U)", "pT3(U)", PT_BINS, PT_LOW, PT_UP)
    pT4 = rt.TH1F("pT4(U)", "pT4(U)", PT_BINS, PT_LOW, PT_UP)
    ETA_BINS, ETA_LOW, ETA_UP = 60, -3, 3
    eta1 = rt.TH1F("eta1(Z)", "eta1(Z)", ETA_BINS, ETA_LOW, ETA_UP)
    eta2 = rt.TH1F("eta2(Z)", "eta2(Z)", ETA_BINS, ETA_LOW, ETA_UP)
    eta3 = rt.TH1F("eta3(U)", "eta3(U)", ETA_BINS, ETA_LOW, ETA_UP)
    eta4 = rt.TH1F("eta4(U)", "eta4(U)", ETA_BINS, ETA_LOW, ETA_UP)
    PHI_BINS, PHI_LOW, PHI_UP = 64, -3.2, 3.2
    phi1 = rt.TH1F("phi1(Z)", "phi1(Z)", PHI_BINS, PHI_LOW, PHI_UP)
    phi2 = rt.TH1F("phi2(Z)", "phi2(Z)", PHI_BINS, PHI_LOW, PHI_UP)
    phi3 = rt.TH1F("phi3(U)", "phi3(U)", PHI_BINS, PHI_LOW, PHI_UP)
    phi4 = rt.TH1F("phi4(U)", "phi4(U)", PHI_BINS, PHI_LOW, PHI_UP)
    DELR_BINS, DELR_LOW, DELR_UP = 60, 0, 6
    delr12 = rt.TH1F("delr12(Z)", "delr12(Z)", DELR_BINS, DELR_LOW, DELR_UP)
    delr34 = rt.TH1F("delr34(U)", "delr34(U)", DELR_BINS, DELR_LOW, DELR_UP)
    delr13 = rt.TH1F("delr13", "delr13", DELR_BINS, DELR_LOW, DELR_UP)
    delr14 = rt.TH1F("delr14", "delr14", DELR_BINS, DELR_LOW, DELR_UP)
    delr23 = rt.TH1F("delr23", "delr23", DELR_BINS, DELR_LOW, DELR_UP)
    delr24 = rt.TH1F("delr24", "delr24", DELR_BINS, DELR_LOW, DELR_UP)
    DPHI_BINS, DPHI_LOW, DPHI_UP = 64, 0, 3.2
    dphi12 = rt.TH1F("abs dphi12(Z)", "abs dphi12(Z)", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi34 = rt.TH1F("abs dphi34(U)", "abs dphi34(U)", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi13 = rt.TH1F("abs dphi13", "abs dphi13", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi14 = rt.TH1F("abs dphi14", "abs dphi14", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi23 = rt.TH1F("abs dphi23", "abs dphi23", DPHI_BINS, DPHI_LOW, DPHI_UP)
    dphi24 = rt.TH1F("abs dphi24", "abs dphi24", DPHI_BINS, DPHI_LOW, DPHI_UP)

    # Lists to make TGraphs
    M12_arr, M34_arr = array('d'), array('d')
    pT1_arr, pT2_arr, pT3_arr, pT4_arr = array('d'), array('d'), array('d'), array('d')
    eta1_arr, eta2_arr, eta3_arr, eta4_arr = array('d'), array('d'), array('d'), array('d')
    delr12_arr, delr34_arr = array('d'), array('d')
    dphi12_arr, dphi34_arr = array('d'), array('d')

    # Loop over events
    myLHEfile = LHEfile("../unweighted_events.lhe")
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    n_acc, n_guar, n_poss, n_uhoh = 0, 0, 0, 0
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        accepted = False
        particles = []
        Z_mu, U_mu = [], []
        Z_mu_q, U_mu_q = [], []
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            particles.append(p)
        for p in particles:
            if abs(p['ID']) == 13:
                if p['mIdx'] == p['mIdx2']:
                    mu_mother = particles[p['mIdx']]
                    if abs(mu_mother['ID']) == 23:
                        Z_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                        Z_mu_q.append(p['ID'])
                    elif abs(mu_mother['ID']) == 35:
                        U_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                        U_mu_q.append(p['ID'])
#       for p in particles:
#           if abs(p['ID']) == 23:
#               Z = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])

#       Z_mu.sort(key = rt.TLorentzVector.Pt, reverse = True)
#       U_mu.sort(key = rt.TLorentzVector.Pt, reverse = True)

        if Z_mu[0].Pt() < Z_mu[1].Pt():
            Z_mu_q[0], Z_mu_q[1] = Z_mu_q[1], Z_mu_q[0]
            Z_mu[0], Z_mu[1] = Z_mu[1], Z_mu[0]

        if U_mu[0].Pt() < U_mu[1].Pt():
            U_mu_q[0], U_mu_q[1] = U_mu_q[1], U_mu_q[0]
            U_mu[0], U_mu[1] = U_mu[1], U_mu[0]

        if cutsOn:
            if (Z_mu[0].Pt() > LEAD_PT_MIN and Z_mu[1].Pt() > LEAD_PT_MIN
                and U_mu[0].Pt() > PT_MIN and U_mu[1].Pt() > PT_MIN
                and abs(Z_mu[0].Eta()) < LEAD_ETA_MAX and abs(Z_mu[1].Eta()) < LEAD_ETA_MAX
                and abs(U_mu[0].Eta()) < ETA_MAX and abs(U_mu[1].Eta()) < ETA_MAX
                and (Z_mu[0] + Z_mu[1]).M() > M12_MIN
                and M4L_MIN < (Z_mu[0] + Z_mu[1] + U_mu[0] + U_mu[1]).M() < M4L_MAX
                ):
                accepted = True
        else:
            if M4L_MIN < (Z_mu[0] + Z_mu[1] + U_mu[0] + U_mu[1]).M() < M4L_MAX:
                accepted = True

        if accepted:
            n_acc += 1
            M4l.Fill((Z_mu[0] + Z_mu[1] + U_mu[0] + U_mu[1]).M())
#           MZ.Fill(Z.M())
            M12.Fill((Z_mu[0] + Z_mu[1]).M())
            M34.Fill((U_mu[0] + U_mu[1]).M())
            pT1.Fill(Z_mu[0].Pt())
            pT2.Fill(Z_mu[1].Pt())
            pT3.Fill(U_mu[0].Pt())
            pT4.Fill(U_mu[1].Pt())
            eta1.Fill(Z_mu[0].Eta())
            eta2.Fill(Z_mu[1].Eta())
            eta3.Fill(U_mu[0].Eta())
            eta4.Fill(U_mu[1].Eta())
            phi1.Fill(Z_mu[0].Phi())
            phi2.Fill(Z_mu[1].Phi())
            phi3.Fill(U_mu[0].Phi())
            phi4.Fill(U_mu[1].Phi())

            M12_arr.append((Z_mu[0] + Z_mu[1]).M())
            M34_arr.append((U_mu[0] + U_mu[1]).M())
            pT1_arr.append(Z_mu[0].Pt())
            pT2_arr.append(Z_mu[1].Pt())
            pT3_arr.append(U_mu[0].Pt())
            pT4_arr.append(U_mu[1].Pt())
            eta1_arr.append(Z_mu[0].Eta())
            eta2_arr.append(Z_mu[1].Eta())
            eta3_arr.append(U_mu[0].Eta())
            eta4_arr.append(U_mu[1].Eta())

            delr12.Fill(Z_mu[0].DeltaR(Z_mu[1]))
            delr34.Fill(U_mu[0].DeltaR(U_mu[1]))
            delr13.Fill(Z_mu[0].DeltaR(U_mu[0]))
            delr14.Fill(Z_mu[0].DeltaR(U_mu[1]))
            delr23.Fill(Z_mu[1].DeltaR(U_mu[0]))
            delr24.Fill(Z_mu[1].DeltaR(U_mu[1]))
            delr12_arr.append(Z_mu[0].DeltaR(Z_mu[1]))
            delr34_arr.append(U_mu[0].DeltaR(U_mu[1]))

            dphi12.Fill(abs(Z_mu[0].DeltaPhi(Z_mu[1])))
            dphi34.Fill(abs(U_mu[0].DeltaPhi(U_mu[1])))
            dphi13.Fill(abs(Z_mu[0].DeltaPhi(U_mu[0])))
            dphi14.Fill(abs(Z_mu[0].DeltaPhi(U_mu[1])))
            dphi23.Fill(abs(Z_mu[1].DeltaPhi(U_mu[0])))
            dphi24.Fill(abs(Z_mu[1].DeltaPhi(U_mu[1])))
            dphi12_arr.append(abs(Z_mu[0].DeltaPhi(Z_mu[1])))
            dphi34_arr.append(abs(U_mu[0].DeltaPhi(U_mu[1])))

#           if Z_mu[1].Pt() < U_mu[1].Pt and Z_mu_q[1] == U_mu_q[1]:
#               n_mis += 1
#           elif Z_mu[1].Pt() < U_mu[0].Pt and Z_mu_q[1] == U_mu_q[0]:
#               n_mis += 1

            if Z_mu[0].Pt() < U_mu[0].Pt():
                n_uhoh += 1
            if Z_mu[1].Pt() < U_mu[0].Pt():
                n_guar += 1
            if U_mu[1].Pt() < Z_mu[1].Pt() < U_mu[0].Pt():
                n_poss += 1

        del oneEvent, myLHEevent

    # Creat dummy graph for line
    line_arr = array('d')
    line_arr.append(0)
    line_arr.append(100)
    graph_line = rt.TGraph(2, line_arr, line_arr)
    graph_line.SetLineColor(rt.kRed)
    graph_line.SetLineWidth(2)

    # Create scatter plots
    graph_Mll = rt.TGraph(n_acc, M12_arr, M34_arr)
    graph_Mll.GetXaxis().SetTitle('M_{12}(Z)')
    graph_Mll.GetYaxis().SetTitle('M_{34}(U)')
    graph_Mll.SetMarkerStyle(rt.kDot)
#   graph_Mll.SetMarkerSize(0.5)
    canvas_Mll = rt.TCanvas("Mll", "Mll", 800, 600)
    canvas_Mll.cd()
    graph_Mll.Draw("AP")
    graph_line.Draw("L")

    graph_ZpT = rt.TGraph(n_acc, pT1_arr, pT2_arr)
    graph_ZpT.GetXaxis().SetTitle('p_{T1}(Z)')
    graph_ZpT.GetYaxis().SetTitle('p_{T2}(Z)')
    graph_ZpT.SetMarkerStyle(rt.kDot)
#   graph_ZpT.SetMarkerSize(0.5)
    canvas_ZpT = rt.TCanvas("Z muon pT", "Z muon pT", 800, 600)
    canvas_ZpT.cd()
    graph_ZpT.Draw("AP")
    graph_line.Draw("L")

    graph_UpT = rt.TGraph(n_acc, pT3_arr, pT4_arr)
    graph_UpT.GetXaxis().SetTitle('p_{T3}(U)')
    graph_UpT.GetYaxis().SetTitle('p_{T4}(U)')
    graph_UpT.SetMarkerStyle(rt.kDot)
#   graph_UpT.SetMarkerSize(0.5)
    canvas_UpT = rt.TCanvas("U muon pT", "U muon pT", 800, 600)
    canvas_UpT.cd()
    graph_UpT.Draw("AP")
    graph_line.Draw("L")

    graph_pT13 = rt.TGraph(n_acc, pT1_arr, pT3_arr)
    graph_pT13.GetXaxis().SetTitle('p_{T1}(Z)')
    graph_pT13.GetYaxis().SetTitle('p_{T3}(U)')
    graph_pT13.SetMarkerStyle(rt.kDot)
#   graph_pT13.SetMarkerSize(0.5)
    canvas_pT13 = rt.TCanvas("High U vs High Z pT", "High U vs High Z pT", 800, 600)
    canvas_pT13.cd()
    graph_pT13.Draw("AP")
    graph_line.Draw("L")

    graph_pT14 = rt.TGraph(n_acc, pT1_arr, pT4_arr)
    graph_pT14.GetXaxis().SetTitle('p_{T1}(Z)')
    graph_pT14.GetYaxis().SetTitle('p_{T4}(U)')
    graph_pT14.SetMarkerStyle(rt.kDot)
#   graph_pT14.SetMarkerSize(0.5)
    canvas_pT14 = rt.TCanvas("Low U vs High Z pT", "Low U vs High Z pT", 800, 600)
    canvas_pT14.cd()
    graph_pT14.Draw("AP")
    graph_line.Draw("L")

    graph_pT23 = rt.TGraph(n_acc, pT2_arr, pT3_arr)
    graph_pT23.GetXaxis().SetTitle('p_{T2}(Z)')
    graph_pT23.GetYaxis().SetTitle('p_{T3}(U)')
    graph_pT23.SetMarkerStyle(rt.kDot)
#   graph_pT23.SetMarkerSize(0.5)
    canvas_pT23 = rt.TCanvas("High U vs Low Z pT", "High U vs Low Z pT", 800, 600)
    canvas_pT23.cd()
    graph_pT23.Draw("AP")
    graph_line.Draw("L")

    graph_pT24 = rt.TGraph(n_acc, pT2_arr, pT4_arr)
    graph_pT24.GetXaxis().SetTitle('p_{T2}(Z)')
    graph_pT24.GetYaxis().SetTitle('p_{T4}(U)')
    graph_pT24.SetMarkerStyle(rt.kDot)
#   graph_pT24.SetMarkerSize(0.5)
    canvas_pT24 = rt.TCanvas("Low U vs Low Z pT", "Low U vs Low Z pT", 800, 600)
    canvas_pT24.cd()
    graph_pT24.Draw("AP")
    graph_line.Draw("L")

    graph_Zeta = rt.TGraph(n_acc, eta1_arr, eta2_arr)
    graph_Zeta.GetXaxis().SetTitle('#eta_{1}(Z)')
    graph_Zeta.GetYaxis().SetTitle('#eta_{2}(Z)')
    graph_Zeta.SetMarkerStyle(rt.kDot)
#   graph_Zeta.SetMarkerSize(0.5)
    canvas_Zeta = rt.TCanvas("Z muon eta", "Z muon eta", 800, 600)
    canvas_Zeta.cd()
    graph_Zeta.Draw("AP")

    graph_Ueta = rt.TGraph(n_acc, eta3_arr, eta4_arr)
    graph_Ueta.GetXaxis().SetTitle('#eta_{3}(U)')
    graph_Ueta.GetYaxis().SetTitle('#eta_{4}(U)')
    graph_Ueta.SetMarkerStyle(rt.kDot)
#   graph_Ueta.SetMarkerSize(0.5)
    canvas_Ueta = rt.TCanvas("U muon eta", "U muon eta", 800, 600)
    canvas_Ueta.cd()
    graph_Ueta.Draw("AP")

    graph_delr = rt.TGraph(n_acc, delr12_arr, delr34_arr)
    graph_delr.GetXaxis().SetTitle('#DeltaR_{12}(Z)')
    graph_delr.GetYaxis().SetTitle('#DeltaR_{34}(U)')
    graph_delr.SetMarkerStyle(rt.kDot)
#   graph_delr.SetMarkerSize(0.5)
    canvas_delr = rt.TCanvas("DeltaR", "DeltaR", 800, 600)
    canvas_delr.cd()
    graph_delr.Draw("AP")
    graph_line.Draw("L")

    graph_dphi = rt.TGraph(n_acc, dphi12_arr, dphi34_arr)
    graph_dphi.GetXaxis().SetTitle('|#Delta#phi_{12}|(Z)')
    graph_dphi.GetYaxis().SetTitle('|#Delta#phi_{34}|(U)')
    graph_dphi.SetMarkerStyle(rt.kDot)
#   graph_dphi.SetMarkerSize(0.5)
    canvas_dphi = rt.TCanvas("DeltaPhi", "DeltaPhi", 800, 600)
    canvas_dphi.cd()
    graph_dphi.Draw("AP")
    graph_line.Draw("L")

    # Write hists and graphs
    if cutsOn:
        fileName = fileName + "_cuts"
    fileName = fileName + ".root"
    outFile = rt.TFile(fileName, "RECREATE")
    outFile.mkdir("Histograms")
    outFile.cd("Histograms")
    M4l.Write()
#   MZ.Write()
    M12.Write()
    M34.Write()
    pT1.Write()
    pT2.Write()
    pT3.Write()
    pT4.Write()
    eta1.Write()
    eta2.Write()
    eta3.Write()
    eta4.Write()
    phi1.Write()
    phi2.Write()
    phi3.Write()
    phi4.Write()
    outFile.cd()
    outFile.mkdir("Scatterplots")
    outFile.cd("Scatterplots")
    canvas_Mll.Write()
    canvas_ZpT.Write()
    canvas_UpT.Write()
    canvas_pT13.Write()
    canvas_pT14.Write()
    canvas_pT23.Write()
    canvas_pT24.Write()
    canvas_Zeta.Write()
    canvas_Ueta.Write()
    outFile.cd()
    outFile.mkdir("deltaR")
    outFile.cd("deltaR")
    delr12.Write()
    delr34.Write()
    delr13.Write()
    delr14.Write()
    delr23.Write()
    delr24.Write()
    canvas_delr.Write()
    outFile.cd()
    outFile.mkdir("deltaPhi")
    outFile.cd("deltaPhi")
    dphi12.Write()
    dphi34.Write()
    dphi13.Write()
    dphi14.Write()
    dphi23.Write()
    dphi24.Write()
    canvas_dphi.Write()
    outFile.cd()
    outFile.Close()

    # Print info
    print("Accepted", n_acc, "events")
#   print("Found", n_guar, "guaranteed mismatches")
#   print("Found", n_poss, "possible mismatches")
#   print("Found", n_uhoh, "big problems")
    print("Created file", fileName)

    del M4l, M12, M34, pT1, pT2, pT3, pT4, eta1, eta2, eta3, eta4, phi1, phi2, phi3, phi4
    del delr12, delr34, delr13, delr14, delr23, delr24
    del graph_Mll, graph_ZpT, graph_UpT, graph_pT13, graph_pT14, graph_pT23, graph_pT24, graph_Zeta, graph_Ueta, graph_delr, graph_dphi
    del canvas_Mll, canvas_ZpT, canvas_UpT, canvas_pT13, canvas_pT14, canvas_pT23, canvas_pT24, canvas_Zeta, canvas_Ueta, canvas_delr, canvas_dphi
