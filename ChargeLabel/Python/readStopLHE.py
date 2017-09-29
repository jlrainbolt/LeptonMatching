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

    # Cuts
    fileName = "hists_CutZM.root"

    # Histograms
    MLL_BINS, MLL_MIN, MLL_MAX = 48, 0, 120
    M4l = rt.TH1F("M4l", "M4l", 80, 50, 130)
    MZ = rt.TH1F("MZ", "MZ", 80, 50, 130)
    M12 = rt.TH1F("M12", "M12", 100, MLL_MIN, 100)
    M34 = rt.TH1F("M34", "M34", 80, MLL_MIN, 40)
    PT_BINS, PT_MIN, PT_MAX = 60, 0, 60
    pT1 = rt.TH1F("pT1", "pT1", PT_BINS, PT_MIN, PT_MAX)
    pT2 = rt.TH1F("pT2", "pT2", PT_BINS, PT_MIN, PT_MAX)
    pT3 = rt.TH1F("pT3", "pT3", PT_BINS, PT_MIN, PT_MAX)
    pT4 = rt.TH1F("pT4", "pT4", PT_BINS, PT_MIN, PT_MAX)
    ETA_BINS, ETA_MIN, ETA_MAX = 60, -3, 3
    eta1 = rt.TH1F("eta1", "eta1", ETA_BINS, ETA_MIN, ETA_MAX)
    eta2 = rt.TH1F("eta2", "eta2", ETA_BINS, ETA_MIN, ETA_MAX)
    eta3 = rt.TH1F("eta3", "eta3", ETA_BINS, ETA_MIN, ETA_MAX)
    eta4 = rt.TH1F("eta4", "eta4", ETA_BINS, ETA_MIN, ETA_MAX)

    # Loop over events
    myLHEfile = LHEfile("unweighted_events.lhe")
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    n_A, n_B, n_C, n_D = 0, 0, 0, 0
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        proc_A = False
        particles = []
        Z_mu, U_mu = [], []
        Z_mu_q, U_mu_q = [], []
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            particles.append(p)

        for p in particles:
            if abs(p['ID']) == 25:
                n_C += 1

        if len(particles) == 7:
            n_D += 1

        # Must have 2 quarks, 1 Z, 1 U, and 4 mu = 8 particles
        if len(particles) == 8:
            # U must be daughter of Z (not quarks)
            for p in particles:
                if abs(p['ID']) == 35:
                    if p['mIdx'] != p['mIdx2']:
                        n_B += 1
                    else:
                        U_mother = particles[p['mIdx']]
                        if abs(U_mother['ID']) == 23:
                            for p in particles:
                                if abs(p['ID']) == 13 and p['mIdx'] == p['mIdx2']:
                                    mu_mother = particles[p['mIdx']]
                                    if abs(mu_mother['ID']) == 23:
                                        Z_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                                        Z_mu_q.append(p['ID'])
                                    elif abs(mu_mother['ID']) == 35:
                                        U_mu.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                                        U_mu_q.append(p['ID'])
                            if len(Z_mu) == 2 and len(U_mu) == 2:
                                if Z_mu_q[0] * Z_mu_q[1] < 0 and U_mu_q[0] * U_mu_q[1] < 0:
                                    proc_A = True
                                    n_A += 1
        if proc_A:
            for p in particles:
                if abs(p['ID']) == 23:
                    Z = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])

            Z_mu.sort(key = rt.TLorentzVector.Pt, reverse = True) 
            U_mu.sort(key = rt.TLorentzVector.Pt, reverse = True)

#           if (Z_mu[0] + Z_mu[1] + U_mu[0] + U_mu[1]).M() < 100:
            if Z.M() < 100:
                M4l.Fill((Z_mu[0] + Z_mu[1] + U_mu[0] + U_mu[1]).M())
                MZ.Fill(Z.M())
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
        del oneEvent, myLHEevent

    # Write hists
    histoFILE = rt.TFile(fileName, "RECREATE")
    M4l.Write()
    MZ.Write()
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
    histoFILE.Close()

    # Print info
    del M4l, MZ, M12, M34, pT1, pT2, pT3, pT4, eta1, eta2, eta3, eta4
    print('Number of proc A events:', n_A)
    print('Number of proc B events:', n_B)
    print('Number of proc C events:', n_C)
    print('Number of proc D events:', n_D)
    print('Total:', n_A + n_B + n_C + n_D)
    print("Created file", fileName)
