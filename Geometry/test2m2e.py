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

    lheName = "unweighted_events.lhe"
    n_max = 100000
    nbins = 100


    # Counters
    n_evts, n_2mu2e = 0, 0


    # Momenta and energies                                          # looking at 2mu2e events

    # Z REST FRAME
    xlo, xup = 0, 50
    P_mu1 = rt.TH1D("P_mu1", "P_{\mu1}", nbins, xlo, xup);          # leading mu
    E_mu1 = rt.TH1D("E_mu1", "E_{\mu1}", nbins, xlo, xup);
    P_mu2 = rt.TH1D("P_mu2", "P_{\mu2}", nbins, xlo, xup);          # subleading mu
    E_mu2 = rt.TH1D("E_mu2", "E_{\mu2}", nbins, xlo, xup);
    P_e1 = rt.TH1D("P_e1", "P_{e1}", nbins, xlo, xup);              # leading e
    E_e1 = rt.TH1D("E_e1", "E_{e1}", nbins, xlo, xup);
    P_e2 = rt.TH1D("P_e2", "P_{e2}", nbins, xlo, xup);              # subleading e
    E_e2 = rt.TH1D("E_e2", "E_{e2}", nbins, xlo, xup);

    xlo, xup = 0, 100
    P_mm = rt.TH1D("P_mm", "P_{\mu1,\mu2}", nbins, xlo, xup);       # muon pair
    E_mm = rt.TH1D("E_mm", "E_{\mu1,\mu2}", nbins, xlo, xup);
    M_mm = rt.TH1D("M_mm", "M_{\mu1,\mu2}", nbins, xlo, xup);
    P_ee = rt.TH1D("P_ee", "P_{e1,e2}", nbins, xlo, xup);           # electron pair
    E_ee = rt.TH1D("E_ee", "E_{e1,e2}", nbins, xlo, xup);
    M_ee = rt.TH1D("M_ee", "M_{e1,e2}", nbins, xlo, xup);
    P_trio = rt.TH1D("P_trio", "P_{\mu2,e1,e2}", nbins, xlo, xup);  # "trailing trio"
    E_trio = rt.TH1D("E_trio", "E_{\mu2,e1,e2}", nbins, xlo, xup);
    M_trio = rt.TH1D("M_trio", "M_{\mu2,e1,e2}", nbins, xlo, xup);

    xlo, xup = 80, 100
    M_Z_pair = rt.TH1D("M_Z_pair", "E_{\mu1,\mu2} + E_{e1,e2}", nbins, xlo, xup); # reconstruct m_Z
    M_Z_trio = rt.TH1D("M_Z_trio", "E_{\mu1} + E_{\mu2,e1,e2}", nbins, xlo, xup); # reconstruct m_Z




    # Angles between pairs
    xlo, xup = 0, 3.2

    # Z REST FRAME
    theta_Z = rt.TH1D("theta_Z", "\\theta_{Z}", nbins, xlo, xup);         # muon 1 & trio
    theta_3 = rt.TH1D("theta_3", "\\theta_{3}", nbins, xlo, xup);         # muon 2 & electron pair
    theta_e = rt.TH1D("theta_e", "\\theta_{e}", nbins, xlo, xup);         # electron 1 & electron 2

    theta_mu1_mu2 = rt.TH1D("theta_mu1_mu2", "\\theta_{\mu1,\mu2}", nbins, xlo, xup)
    theta_mu1_e1 = rt.TH1D("theta_mu1_e1", "\\theta_{\mu1,e1}", nbins, xlo, xup)
    theta_mu1_e2 = rt.TH1D("theta_mu1_e2", "\\theta_{\mu1,e2}", nbins, xlo, xup)

    # "TRAILING TRIO" REST FRAME (prime)
    theta_3_p = rt.TH1D("theta_3_p", "\\theta'_{3}", nbins, xlo, xup);    # muon 2 & electron pair
    theta_e_p = rt.TH1D("theta_e_p", "\\theta'_{e}", nbins, xlo, xup);    # electron 1 & electron 2

    # ELECTRON PAIR REST FRAME (double prime)
    theta_e_pp = rt.TH1D("theta_e_pp", "\\theta''_{e^}", nbins, xlo, xup);# electron 1 & electron 2


    
    # Derived quantities
   
    # Fraction of energy to muon 2 or electron pair
    xlo, xup = 0, 2
    P_frac_mu2 = rt.TH1D("P_frac_mu2", "P_{\mu2}/P_{\mu1}", nbins, xlo, xup)
    E_frac_mu2 = rt.TH1D("E_frac_mu2", "E_{\mu2}/E_{\mu1}", nbins, xlo, xup)
    P_frac_ee = rt.TH1D("P_frac_ee", "P_{ee}/P_{\mu1}", nbins, xlo, xup)
    E_frac_ee = rt.TH1D("E_frac_ee", "E_{ee}/E_{\mu1}", nbins, xlo, xup)
    P_frac_e1 = rt.TH1D("P_frac_e1", "P_{e1}/P_{\mu2}", nbins, xlo, xup)
    P_frac_e2 = rt.TH1D("P_frac_e2", "P_{e2}/P_{\mu2}", nbins, xlo, xup)

    # "Difference" between muon 2 electron pair energy
    xlo, xup = -1, 1
    P_diff = rt.TH1D("P_diff", "(P_{\mu2} - P_{ee})/P_{\mu1}", nbins, xlo, xup)
    E_diff = rt.TH1D("E_diff", "(E_{\mu2} - E_{ee})/E_{\mu1}", nbins, xlo, xup)
    

    
    # 2D histograms
    nbins = 50

    # muon 1 E vs. theta_3
    xlo, xup, ylo, yup = 0, 3.15, 20, 50
    E_mu1__theta_3 = rt.TH2D("E_mu1__theta_3", "E_{\mu1} vs. \\theta_{3}", \
            nbins, xlo, xup, nbins, ylo, yup)
    xlo, xup = 0, 1
    E_mu1__cos2_3 = rt.TH2D("E_mu1__cos2_3", "E_{\mu1} vs. cos^{2}(\\theta_{3}/2)", \
            nbins, xlo, xup, nbins, ylo, yup)

    # muon 2 E vs. theta_3
    xlo, xup, ylo, yup = 0, 3.15, 0, 50
    E_mu2__theta_3 = rt.TH2D("E_mu2__theta_3", "E_{\mu2} vs. \\theta_{3}", \
            nbins, xlo, xup, nbins, ylo, yup)
    
    # "Difference" between muon 2 and elec pair energy vs. cos(theta_3)
    ylo, yup = -1, 1
    E_diff__theta_3 = rt.TH2D("E_diff__theta_3", "(E_{\mu2} - E_{ee})/E_{\mu1} vs. \\theta_{3}", \
            nbins, xlo, xup, nbins, ylo, yup)
    xlo, xup = -1, 1
    E_diff__cos_3 = rt.TH2D("E_diff__cos_3", "(E_{\mu2} - E_{ee})/E_{\mu1} vs. cos(\\theta_{3})", \
            nbins, xlo, xup, nbins, ylo, yup)
    E_diff__tan_3 = rt.TH2D("E_diff__tan_3", "(E_{\mu2} - E_{ee})/E_{\mu1} vs. tan(\\theta_{3})", \
            nbins, xlo, xup, nbins, ylo, yup)








    ################
    #  GET EVENTS  #
    ################

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
            continue

        muon_pair = muons[0] + muons[1]
        elec_pair = elecs[0] + elecs[1]

        if muon_pair.M() < elec_pair.M():
            continue
        
        n_2mu2e += 1


        # Check leading/subleading
        muon1 = muons[0]
        muon2 = muons[1]
        if muon1.E() < muon2.E():
            muon1, muon2 = muon2, muon1

        elec1 = elecs[0]
        elec2 = elecs[1]
        if elec1.E() < elec2.E():
            elec1, elec2 = elec2, elec1

        trio = muon2 + elec_pair


        
        ################
        #     BOOST    #
        ################

        # Primed ("trailing trio") frame
        boost_p = -trio.BoostVector()

        muon1_p = GetBoosted(muon1, boost_p)
        muon2_p = GetBoosted(muon2, boost_p)
        muon_pair_p = muon1_p + muon2_p

        elec1_p = GetBoosted(elec1, boost_p)
        elec2_p = GetBoosted(elec2, boost_p)
        elec_pair_p = elec1_p + elec2_p

        
        # Double primed (electron pair) frame
        boost_pp = -elec_pair_p.BoostVector()
        elec1_pp = GetBoosted(elec1_p, boost_pp)
        elec2_pp = GetBoosted(elec2_p, boost_pp)



        ################
        #  FILL HISTS  #
        ################

        # Momentum and energy
        P_mu1.Fill(muon1.P())
        E_mu1.Fill(muon1.E())
        P_mu2.Fill(muon2.P())
        E_mu2.Fill(muon2.E())
        P_e1.Fill(elec1.P())
        E_e1.Fill(elec1.E())
        P_e2.Fill(elec2.P())
        E_e2.Fill(elec2.E())

        P_mm.Fill(muon_pair.P())
        E_mm.Fill(muon_pair.E())
        M_mm.Fill(muon_pair.M())
        P_ee.Fill(elec_pair.P())
        E_ee.Fill(elec_pair.E())
        M_ee.Fill(elec_pair.M())
        P_trio.Fill(trio.P())
        E_trio.Fill(trio.E())
        M_trio.Fill(trio.M())

        # Z mass tests
        M_Z_pair.Fill(muon_pair.E() + elec_pair.E());
        M_Z_trio.Fill(trio.E() + muon1.E());


        # Angles
        angle_Z = muon1.Angle(trio.Vect())
        angle_3 = muon2.Angle(elec_pair.Vect())
        angle_e = elec1.Angle(elec2.Vect())
        angle_mu1_mu2 = muon1.Angle(muon2.Vect())
        angle_mu1_e1 = muon1.Angle(elec1.Vect())
        angle_mu1_e2 = muon1.Angle(elec2.Vect())
        angle_3_p = muon2_p.Angle(elec_pair_p.Vect())
        angle_e_p = elec1_p.Angle(elec2_p.Vect())
        angle_e_pp = elec1_pp.Angle(elec2_pp.Vect())

        theta_Z.Fill(angle_Z)
        theta_3.Fill(angle_3)
        theta_e.Fill(angle_e)
        theta_mu1_mu2.Fill(angle_mu1_mu2)
        theta_mu1_e1.Fill(angle_mu1_e1)
        theta_mu1_e2.Fill(angle_mu1_e2)
        theta_3_p.Fill(angle_3_p)
        theta_e_p.Fill(angle_e_p)
        theta_e_pp.Fill(angle_e_pp)


        # Derived quantities
        E_fraction_mu2 = muon2.E() / muon1.E()
        P_fraction_mu2 = muon2.P() / muon1.P()
        E_fraction_ee = elec_pair.E() / muon1.E()
        P_fraction_ee = elec_pair.P() / muon1.P()
        P_fraction_e1 = elec1.P() / muon2.P()
        P_fraction_e2 = elec2.P() / muon2.P()
        E_difference = (muon2.E() - elec_pair.E()) / muon1.E()
        P_difference = (muon2.P() - elec_pair.P()) / muon1.P()

        E_frac_mu2.Fill(E_fraction_mu2)
        P_frac_mu2.Fill(P_fraction_mu2)
        E_frac_ee.Fill(E_fraction_ee)
        P_frac_ee.Fill(P_fraction_ee)
        P_frac_e1.Fill(P_fraction_e1)
        P_frac_e2.Fill(P_fraction_e2)
        E_diff.Fill(E_difference)
        P_diff.Fill(P_difference)


        # 2D hists
        E_mu1__theta_3.Fill(angle_3, muon1.E())
        E_mu1__cos2_3.Fill(math.cos(angle_3 / 2)**2, muon1.E())
        E_mu2__theta_3.Fill(angle_3, muon2.E())
        E_diff__theta_3.Fill(angle_3, E_difference)
        E_diff__cos_3.Fill(math.cos(angle_3), E_difference)
        E_diff__tan_3.Fill(math.tan(angle_3), E_difference)



        del oneEvent, myLHEevent
    print(n_2mu2e, "2mu2e decays")


    ################
    # SAVE TO FILE #
    ################

    outFile = rt.TFile("2m2e.root", "RECREATE")

    P_mu1.Write()
    E_mu1.Write()
    P_mu2.Write()
    E_mu2.Write()
    P_e1.Write()
    E_e1.Write()
    P_e2.Write()
    E_e2.Write()

    P_mm.Write()
    E_mm.Write()
    M_mm.Write()
    P_ee.Write()
    E_ee.Write()
    M_ee.Write()
    P_trio.Write()
    E_trio.Write()
    M_trio.Write()

    M_Z_pair.Write()
    M_Z_trio.Write()

    theta_Z.Write()
    theta_3.Write()
    theta_e.Write()
    theta_mu1_mu2.Write()
    theta_mu1_e1.Write()
    theta_mu1_e2.Write()
    theta_3_p.Write()
    theta_e_p.Write()
    theta_e_pp.Write()

    E_frac_mu2.Write()
    P_frac_mu2.Write()
    E_frac_ee.Write()
    P_frac_ee.Write()
    P_frac_e1.Write()
    P_frac_e2.Write()
    E_diff.Write()
    P_diff.Write()

    E_mu1__theta_3.Write()
    E_mu1__cos2_3.Write()
    E_mu2__theta_3.Write()
    E_diff__theta_3.Write()
    E_diff__cos_3.Write()
    E_diff__tan_3.Write()

    outFile.Close()
