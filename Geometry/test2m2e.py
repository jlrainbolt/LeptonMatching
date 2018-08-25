from __future__ import print_function
from __future__ import division
import sys
import ROOT as rt
import math
from array import array
from LHEevent import *
from LHEfile import *
import plotTools


def GetBoosted(p4_, beta):
    p4 = p4_.Clone()
    p4.Boost(beta)
    return p4



if __name__ == '__main__':



    ################
    #  INITIALIZE  #
    ################

    printEvents = False
    printEvery = 100

    lheName = "unweighted_events.lhe"
    n_max = 100000
    nbins = 100


    # Counters
    n_evts, n_2mu2e = 0, 0


    # Momentum of each lepton                                         # looking at 2mu2e events

    # Z REST FRAME
    xlo, xup = 0, 50
    P_mu1 = rt.TH1F("P_mu1", "P_{\mu1}", nbins, xlo, xup);          # leading mu
    P_mu2 = rt.TH1F("P_mu2", "P_{\mu2}", nbins, xlo, xup);          # subleading mu
    P_e1 = rt.TH1F("P_e1", "P_{e1}", nbins, xlo, xup);              # leading e
    P_e2 = rt.TH1F("P_e2", "P_{e2}", nbins, xlo, xup);              # subleading e

    xlo, xup = 0, 100
    P_mm = rt.TH1F("P_mm", "P_{\mu1,\mu2}", nbins, xlo, xup);       # muon pair
    E_mm = rt.TH1F("E_mm", "E_{\mu1,\mu2}", nbins, xlo, xup);
    P_ee = rt.TH1F("P_ee", "P_{e1,e2}", nbins, xlo, xup);           # electron pair
    E_ee = rt.TH1F("E_ee", "E_{e1,e2}", nbins, xlo, xup);
    P_trio = rt.TH1F("P_trio", "P_{\mu2,e1,e2}", nbins, xlo, xup);  # "trailing trio"
    E_trio = rt.TH1F("E_trio", "E_{\mu2,e1,e2}", nbins, xlo, xup);  # "trailing trio"

    xlo, xup = 80, 100
    M_Z_pair = rt.TH1F("M_Z_pair", "E_{\mu1,\mu2} + E_{e1,e2}", nbins, xlo, xup); # reconstruct m_Z
    M_Z_trio = rt.TH1F("M_Z_trio", "E_{\mu1} + E_{\mu2,e1,e2}", nbins, xlo, xup); # reconstruct m_Z




    # Angles between pairs
    xlo, xup = 0, 3.15

    # Z REST FRAME
    theta_Z = rt.TH1F("theta_Z", "theta_{Z}", nbins, xlo, xup);         # muon 1 & trio
    theta_3 = rt.TH1F("theta_3", "theta_{3}", nbins, xlo, xup);         # muon 2 & electron pair
    theta_e = rt.TH1F("theta_e", "theta_{e}", nbins, xlo, xup);         # electron 1 & electron 2

    # "TRAILING TRIO" REST FRAME (prime)
    theta_3_p = rt.TH1F("theta_3_p", "theta'_{3}", nbins, xlo, xup);    # muon 2 & electron pair
    theta_e_p = rt.TH1F("theta_e_p", "theta'_{e}", nbins, xlo, xup);    # electron 1 & electron 2

    # ELECTRON PAIR REST FRAME (double prime)
    theta_e_pp = rt.TH1F("theta_e_pp", "theta''_{e^}", nbins, xlo, xup);# electron 1 & electron 2


    
    # 2D

    ylo, yup = 0, 50

    # muon 2 E vs. theta_3
    P_mu1__theta_3 = rt.TH2F("P_mu1__theta_3", "P_{\mu1} vs. theta_{3}", \
            nbins, xlo, xup, nbins, ylo, yup)
    ylo, yup = -50, 50
    E_diff__theta_3 = rt.TH2F("E_diff__theta_3", "(E_{\mu2} - E_{ee}) vs. theta_{3}", \
            nbins, xlo, xup, nbins, ylo, yup)
    xlo, xup = 0, 50
    P_mu1__cos_3 = rt.TH2F("P_mu1__cos_3", "P_{\mu1} vs. P_{\mu2} cos(theta_{3}/2)", \
            nbins, xlo, xup, nbins, ylo, yup)








    ################
    #  GET EVENTS  #
    ################

    if printEvents:
        print("\nEVENTS (1 IN ", printEvery, ")", sep = "")
        print("=================================================================")
        print("Lookup (l1 px)\t", " First decay", "\t Second decay", sep = "\t")
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    else:
        print("\nRunning...", end = "")
        sys.stdout.flush()

    myLHEfile = LHEfile(lheName)
    myLHEfile.setMax(n_max)
    
    print("Looping over", n_max, "events...", end='')
    sys.stdout.flush()

    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:

        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        accepted = True



        #################
        #  "SELECTION"  #
        #################

        # Get leptons
        particles = myLHEevent.Particles
        elecs, muons = [], []

        for i in range(0, len(particles)):
            p = particles[i]
            if abs(p['ID']) == 11:
                elecs.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
            elif abs(p['ID']) == 13:
                muons.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))


        # Accept only 2mu2e events (for now)
        if len(elecs) != 2 or len(muons) != 2:
            continue;

        muon_pair = muons[0] + muons[1]
        elec_pair = elecs[0] + elecs[1]

        if (muon_pair.M() < elec_pair.M()):
            continue;
        
        n_2mu2e += 1;


        # Check leading/subleading
        muon1 = muons[0]
        muon2 = muons[1]
        if muon1.P() < muon2.P():
            muon1, muon2 = muon2, muon1

        elec1 = elecs[0]
        elec2 = elecs[1]
        if elec1.P() < elec2.P():
            elec1, elec2 = elec2, elec1

        trio = muon2 + elec_pair


        
        ################
        #     BOOST    #
        ################

        boost_p = -trio.BoostVector()

        muon1_p = GetBoosted(muon1, boost_p)
        muon2_p = GetBoosted(muon2, boost_p)
        muon_pair_p = muon1_p + muon2_p

        elec1_p = GetBoosted(elec1, boost_p)
        elec2_p = GetBoosted(elec2, boost_p)
        elec_pair_p = elec1_p + elec2_p


        boost_pp = -elec_pair_p.BoostVector()
        elec1_pp = GetBoosted(elec1_p, boost_pp)
        elec2_pp = GetBoosted(elec2_p, boost_pp)



        ################
        #  FILL HISTS  #
        ################

        P_mu1.Fill(muon1.P())
        P_mu2.Fill(muon2.P())
        P_e1.Fill(elec1.P())
        P_e2.Fill(elec2.P())

        P_mm.Fill(muon_pair.P())
        E_mm.Fill(muon_pair.E())
        P_ee.Fill(elec_pair.P())
        E_ee.Fill(elec_pair.E())
        P_trio.Fill(trio.P())
        E_trio.Fill(trio.E())

        M_Z_pair.Fill(muon_pair.E() + elec_pair.E());
        M_Z_trio.Fill(trio.E() + muon1.E());

        theta_Z.Fill(muon1.Angle(trio.Vect()))
        theta_3.Fill(muon2.Angle(elec_pair.Vect()))
        theta_e.Fill(elec1.Angle(elec2.Vect()))

        theta_3_p.Fill(muon2_p.Angle(elec_pair_p.Vect()))
        theta_e_p.Fill(elec1_p.Angle(elec2_p.Vect()))

        theta_e_pp.Fill(elec1_pp.Angle(elec2_pp.Vect()))

        P_mu1__theta_3.Fill(muon2.Angle(elec_pair.Vect()), muon1.P())
        E_diff__theta_3.Fill(muon2.Angle(elec_pair.Vect()), muon2.E() - elec_pair.E())
        P_mu1__cos_3.Fill((muon2.P() + elec_pair.P()) * math.cos(0.5 * muon2.Angle(elec_pair.Vect())), muon1.P())



        del oneEvent, myLHEevent
    if not printEvents:
        print("done!")

    print(n_2mu2e, "2mu2e decays")


    ################
    # SAVE TO FILE #
    ################

    outFile = rt.TFile("2m2e.root", "RECREATE")

    P_mu1.Write()
    P_mu2.Write()
    P_e1.Write()
    P_e2.Write()

    P_mm.Write()
    E_mm.Write()
    P_ee.Write()
    E_ee.Write()
    P_trio.Write()
    E_trio.Write()

    M_Z_pair.Write()
    M_Z_trio.Write()

    theta_Z.Write()
    theta_3.Write()
    theta_e.Write()
    theta_3_p.Write()
    theta_e_p.Write()
    theta_e_pp.Write()

    P_mu1__theta_3.Write()
    E_diff__theta_3.Write()
    P_mu1__cos_3.Write()

    outFile.Close()
