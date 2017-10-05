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
    fileName = "hists"
    cutsOn = True

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

    # Lists to make TGraphs
    M12_arr, M34_arr = array('d'), array('d')
    pT1_arr, pT2_arr, pT3_arr, pT4_arr = array('d'), array('d'), array('d'), array('d')
    eta1_arr, eta2_arr, eta3_arr, eta4_arr = array('d'), array('d'), array('d'), array('d')

    # Loop over events
    myLHEfile = LHEfile("../unweighted_events.lhe")
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    n_acc =0 
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
            if abs(p['ID']) == 13 and p['mIdx'] == p['mIdx2']:
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

        if Z_mu_q[0] < Z_mu_q[1]:
            Z_mu_q[0], Z_mu_q[1] = Z_mu_q[1], Z_mu_q[0]
            Z_mu[0], Z_mu[1] = Z_mu[1], Z_mu[0]

        if U_mu_q[0] < U_mu_q[1]:
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
            M12_arr.append((Z_mu[0] + Z_mu[1]).M())
            M34_arr.append((U_mu[0] + U_mu[1]).M())
            pT1.Fill(Z_mu[0].Pt())
            pT2.Fill(Z_mu[1].Pt())
            pT3.Fill(U_mu[0].Pt())
            pT4.Fill(U_mu[1].Pt())
            pT1_arr.append(Z_mu[0].Pt())
            pT2_arr.append(Z_mu[1].Pt())
            pT3_arr.append(U_mu[0].Pt())
            pT4_arr.append(U_mu[1].Pt())
            eta1.Fill(Z_mu[0].Eta())
            eta2.Fill(Z_mu[1].Eta())
            eta3.Fill(U_mu[0].Eta())
            eta4.Fill(U_mu[1].Eta())
            eta1_arr.append(Z_mu[0].Eta())
            eta2_arr.append(Z_mu[1].Eta())
            eta3_arr.append(U_mu[0].Eta())
            eta4_arr.append(U_mu[1].Eta())
        del oneEvent, myLHEevent

    # Create scatter plots
    graph_Mll = rt.TGraph(n_acc, M12_arr, M34_arr)
    graph_Mll.GetXaxis().SetTitle('M_{12}(Z)')
    graph_Mll.GetYaxis().SetTitle('M_{34}(U)')
    graph_Mll.SetMarkerStyle(rt.kDot)
#   graph_Mll.SetMarkerSize(0.5)
    canvas_Mll = rt.TCanvas("Mll", "Mll", 800, 600)
    canvas_Mll.cd()
    graph_Mll.Draw("AP")

    graph_ZpT = rt.TGraph(n_acc, pT1_arr, pT2_arr)
    graph_ZpT.GetXaxis().SetTitle('p_{T1}(Z,+)')
    graph_ZpT.GetYaxis().SetTitle('p_{T2}(Z,-)')
    graph_ZpT.SetMarkerStyle(rt.kDot)
#   graph_ZpT.SetMarkerSize(0.5)
    canvas_ZpT = rt.TCanvas("Z muon pT", "Z muon pT", 800, 600)
    canvas_ZpT.cd()
    graph_ZpT.Draw("AP")

    graph_UpT = rt.TGraph(n_acc, pT3_arr, pT4_arr)
    graph_UpT.GetXaxis().SetTitle('p_{T3}(U,+)')
    graph_UpT.GetYaxis().SetTitle('p_{T4}(U,-)')
    graph_UpT.SetMarkerStyle(rt.kDot)
#   graph_UpT.SetMarkerSize(0.5)
    canvas_UpT = rt.TCanvas("U muon pT", "U muon pT", 800, 600)
    canvas_UpT.cd()
    graph_UpT.Draw("AP")

    graph_Zeta = rt.TGraph(n_acc, eta1_arr, eta2_arr)
    graph_Zeta.GetXaxis().SetTitle('#eta_{1}(Z,+)')
    graph_Zeta.GetYaxis().SetTitle('#eta_{2}(Z,-)')
    graph_Zeta.SetMarkerStyle(rt.kDot)
#   graph_Zeta.SetMarkerSize(0.5)
    canvas_Zeta = rt.TCanvas("Z muon eta", "Z muon eta", 800, 600)
    canvas_Zeta.cd()
    graph_Zeta.Draw("AP")

    graph_Ueta = rt.TGraph(n_acc, eta3_arr, eta4_arr)
    graph_Ueta.GetXaxis().SetTitle('#eta_{3}(U,+)')
    graph_Ueta.GetYaxis().SetTitle('#eta_{4}(U,-)')
    graph_Ueta.SetMarkerStyle(rt.kDot)
#   graph_Ueta.SetMarkerSize(0.5)
    canvas_Ueta = rt.TCanvas("U muon eta", "U muon eta", 800, 600)
    canvas_Ueta.cd()
    graph_Ueta.Draw("AP")

    # Write hists and graphs
    if cutsOn:
        fileName = fileName + "_cuts"
    fileName = fileName + ".root"
    histoFILE = rt.TFile(fileName, "RECREATE")
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
    canvas_Mll.Write()
    canvas_ZpT.Write()
    canvas_UpT.Write()
    canvas_Zeta.Write()
    canvas_Ueta.Write()
    histoFILE.Close()

    # Print info
    print("Accepted", n_acc, "events")
    print("Created file", fileName)

    del M4l, M12, M34, pT1, pT2, pT3, pT4, eta1, eta2, eta3, eta4
    del graph_Mll, graph_ZpT, graph_UpT, graph_Zeta, graph_Ueta
    del canvas_Mll, canvas_ZpT, canvas_UpT, canvas_Zeta, canvas_Ueta
