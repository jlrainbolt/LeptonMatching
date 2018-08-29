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


def GetSum(p4_list):
    p4 = rt.TLorentzVector()
    for p4_ in p4_list:
        p4 += p4_
    return p4


if __name__ == '__main__':



    ################
    #  INITIALIZE  #
    ################

    lheName = "unweighted_events.lhe"
    n_max = 100
    nbins = 100

    pi_tol = 0.0001


    # Counters
    n_4mu = 0


    # Momenta and energies

    # Z REST FRAME

    # Pair "a" leptons (high mass)
    xlo, xup = 0, 50
    P_lep_a1 = rt.TH1D("P_lep_a1", "P_{a1}", nbins, xlo, xup);
    E_lep_a1 = rt.TH1D("E_lep_a1", "E_{a1}", nbins, xlo, xup);
#   P_lep_a2 = rt.TH1D("P_lep_a2", "P_{a2}", nbins, xlo, xup);
#   E_lep_a2 = rt.TH1D("E_lep_a2", "E_{a2}", nbins, xlo, xup);
#   xlo, xup = 0, 100
#   P_pair_a = rt.TH1D("P_pair_a", "P_{a1,a2}", nbins, xlo, xup);
#   E_pair_a = rt.TH1D("E_pair_a", "E_{a1,a2}", nbins, xlo, xup);
#   M_pair_a = rt.TH1D("M_pair_a", "M_{a1,a2}", nbins, xlo, xup);

#   # Pair "b" leptons (low mass)
#   P_lep_b1 = rt.TH1D("P_lep_b1", "P_{b1}", nbins, xlo, xup);
#   E_lep_b1 = rt.TH1D("E_lep_b1", "E_{b1}", nbins, xlo, xup);
#   P_lep_b2 = rt.TH1D("P_lep_b2", "P_{b2}", nbins, xlo, xup);
#   E_lep_b2 = rt.TH1D("E_lep_b2", "E_{b2}", nbins, xlo, xup);
    xlo, xup = 0, 100
#   P_pair_b = rt.TH1D("P_pair_b", "P_{b1,b2}", nbins, xlo, xup);
#   E_pair_b = rt.TH1D("E_pair_b", "E_{b1,b2}", nbins, xlo, xup);
#   M_pair_b = rt.TH1D("M_pair_b", "M_{b1,b2}", nbins, xlo, xup);

    # "Trailing trio"
    P_trio = rt.TH1D("P_trio", "P_{a2,b1,b2}", nbins, xlo, xup);
    E_trio = rt.TH1D("E_trio", "E_{a2,b1,b2}", nbins, xlo, xup);
    M_trio = rt.TH1D("M_trio", "M_{a2,b1,b2}", nbins, xlo, xup);



    # Angles between pairs
    xlo, xup = 0, 3.2

    # Z REST FRAME
    theta_Z = rt.TH1D("theta_Z", "\\theta_{Z}", nbins, xlo, xup);
#   theta_3 = rt.TH1D("theta_3", "\\theta_{3}", nbins, xlo, xup);
#   theta_b = rt.TH1D("theta_b", "\\theta_{b}", nbins, xlo, xup);

#   # "TRAILING TRIO" REST FRAME (prime)
#   theta_3_p = rt.TH1D("theta_3_p", "\\theta'_{3}", nbins, xlo, xup);
#   theta_b_p = rt.TH1D("theta_b_p", "\\theta'_{b}", nbins, xlo, xup);

#   # ELECTRON PAIR REST FRAME (double prime)
#   theta_b_pp = rt.TH1D("theta_b_pp", "\\theta''_{b}", nbins, xlo, xup);


#   
#   # Derived quantities
#  
#   # Fraction of energy to lep b2 or pair b
#   xlo, xup = 0, 1
#   P_frac_lep_a2 = rt.TH1D("P_frac_lep_a2", "P_{a2}/P_{a1}", nbins, xlo, xup)
#   E_frac_lep_a2 = rt.TH1D("E_frac_lep_a2", "E_{a2}/E_{a1}", nbins, xlo, xup)
#   P_frac_pair_b = rt.TH1D("P_frac_pair_b", "P_{b1,b2}/P_{a1}", nbins, xlo, xup)
#   E_frac_pair_b = rt.TH1D("E_frac_pair_b", "E_{b1,b2}/E_{a1}", nbins, xlo, xup)

#   # "Difference" between lep a2 and pair b energy
#   xlo, xup = -1, 1
#   P_diff = rt.TH1D("P_diff", "(P_{a2} - P_{b1,b2})/P_{a1}", nbins, xlo, xup)
#   E_diff = rt.TH1D("E_diff", "(E_{a2} - E_{b1,b2})/E_{a1}", nbins, xlo, xup)
#   

#   
#   # 2D histograms
#   nbins = 50

#   # Lep a1 E vs. theta_3
#   xlo, xup, ylo, yup = 0, 3.15, 20, 50
#   E_lep_a1__theta_3 = rt.TH2D("E_lep_a1__theta_3", "E_{a1} vs. \\theta_{3}", \
#           nbins, xlo, xup, nbins, ylo, yup)
#   xlo, xup = 0, 1
#   E_lep_a1__cos2_3 = rt.TH2D("E_lep_a1__cos2_3", "E_{a1} vs. cos^{2}(\\theta_{3}/2)", \
#           nbins, xlo, xup, nbins, ylo, yup)

#   # Lep a2 E vs. theta_3
#   xlo, xup, ylo, yup = 0, 3.15, 0, 50
#   E_lep_a2__theta_3 = rt.TH2D("E_lep_a2__theta_3", "E_{a2} vs. \\theta_{3}", \
#           nbins, xlo, xup, nbins, ylo, yup)
#   
#   # "Difference" between lep a2 and pair b energy vs. cos(theta_3)
#   ylo, yup = -1, 1
#   E_diff__theta_3 = rt.TH2D("E_diff__theta_3", "(E_{a2} - E_{b1,b2})/E_{a1} vs. \\theta_{3}", \
#           nbins, xlo, xup, nbins, ylo, yup)
#   xlo, xup = -1, 1
#   E_diff__cos_3 = rt.TH2D("E_diff__cos_3", "(E_{a2} - E_{b1,b2})/E_{a1} vs. cos(\\theta_{3})", \
#           nbins, xlo, xup, nbins, ylo, yup)
#   E_diff__tan_3 = rt.TH2D("E_diff__tan_3", "(E_{a2} - E_{b1,b2})/E_{a1} vs. tan(\\theta_{3})", \
#           nbins, xlo, xup, nbins, ylo, yup)








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
        elecs_q, muons_q = {}, {}

        for i in range(0, len(particles)):
            p = particles[i]
            if abs(p['ID']) == 11:
                elecs.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                elecs_q[elecs[-1]] = p['Q']
            elif abs(p['ID']) == 13:
                muons.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                muons_q[muons[-1]] = p['Q']


        # Accept only 4mu events (for now)
        if len(muons) != 4:
            continue
    # elif len(muons) == 4:
        rem_leps = muons
        leps_q = muons_q
        n_4mu += 1


        # Order by energy and assign lep a1
        rem_leps.sort(key=lambda fourvec: fourvec.E(), reverse=True)
        lep_a1 = rem_leps.pop(0)

        
        # Assign "trailing trio", get theta_Z
        trio = GetSum(rem_leps)
        angle_Z = lep_a1.Angle(trio.Vect())

        if (math.fabs(angle_Z - math.pi) > pi_tol):
            print("Problem in S")



        ################
        #     BOOST    #
        ################

        # Primed ("trailing trio") frame
        boost_p = -trio.BoostVector()
        lep_a1_p = GetBoosted(lep_a1, boost_p)

        for i in range(len(rem_leps)):
            # Make sure a2 candidate has opposite charge (remove later as sanity check?)
            if leps_q[rem_leps[i]] != leps_q[lep_a1]:
                rem_leps_ = [GetBoosted(lep, boost_p) for lep in rem_leps]
                lep_a2_ = rem_leps_.pop(i)
                pair_b_ = GetSum(rem_leps_)
#               print(lep_a2_.P(), pair_b_.P(), lep_a2_.Angle(pair_b_.Vect()))
                print(rem_leps[i].E(), rem_leps_[0].Angle(rem_leps_[1].Vect()))
        print("")



#       muon2_p = GetBoosted(muon2, boost_p)
#       muon_pair_p = muon1_p + muon2_p

#       elec1_p = GetBoosted(elec1, boost_p)
#       elec2_p = GetBoosted(elec2, boost_p)
#       elec_pair_p = elec1_p + elec2_p

#       
#       # Double primed (electron pair) frame
#       boost_pp = -elec_pair_p.BoostVector()
#       elec1_pp = GetBoosted(elec1_p, boost_pp)
#       elec2_pp = GetBoosted(elec2_p, boost_pp)



        ################
        #  FILL HISTS  #
        ################

        # Momentum and energy
        P_lep_a1.Fill(lep_a1.P())
        E_lep_a1.Fill(lep_a1.E())
#       P_lep_a2.Fill(muon2.P())
#       E_lep_a2.Fill(muon2.E())
#       P_lep_b1.Fill(elec1.P())
#       E_lep_b1.Fill(elec1.E())
#       P_lep_b2.Fill(elec2.P())
#       E_lep_b2.Fill(elec2.E())

#       P_pair_a.Fill(muon_pair.P())
#       E_pair_a.Fill(muon_pair.E())
#       M_pair_a.Fill(muon_pair.M())
#       P_pair_b.Fill(elec_pair.P())
#       E_pair_b.Fill(elec_pair.E())
#       M_pair_b.Fill(elec_pair.M())
        P_trio.Fill(trio.P())
        E_trio.Fill(trio.E())
        M_trio.Fill(trio.M())


        # Angles
#       angle_3 = muon2.Angle(elec_pair.Vect())
#       angle_b = elec1.Angle(elec2.Vect())
#       angle_3_p = muon2_p.Angle(elec_pair_p.Vect())
#       angle_b_p = elec1_p.Angle(elec2_p.Vect())
#       angle_b_pp = elec1_pp.Angle(elec2_pp.Vect())

        theta_Z.Fill(angle_Z)
#       theta_3.Fill(angle_3)
#       theta_b.Fill(angle_b)
#       theta_3_p.Fill(angle_3_p)
#       theta_b_p.Fill(angle_b_p)
#       theta_b_pp.Fill(angle_b_pp)


#       # Derived quantities
#       E_fraction_lep_a2 = muon2.E() / muon1.E()
#       P_fraction_lep_a2 = muon2.P() / muon1.P()
#       E_fraction_pair_b = elec_pair.E() / muon1.E()
#       P_fraction_pair_b = elec_pair.P() / muon1.P()
#       E_difference = (muon2.E() - elec_pair.E()) / muon1.E()
#       P_difference = (muon2.P() - elec_pair.P()) / muon1.P()

#       E_frac_lep_a2.Fill(E_fraction_lep_a2)
#       P_frac_lep_a2.Fill(P_fraction_lep_a2)
#       E_frac_pair_b.Fill(E_fraction_pair_b)
#       P_frac_pair_b.Fill(P_fraction_pair_b)
#       E_diff.Fill(E_difference)
#       P_diff.Fill(P_difference)


#       # 2D hists
#       E_lep_a1__theta_3.Fill(angle_3, muon1.E())
#       E_lep_a1__cos2_3.Fill(math.cos(angle_3 / 2)**2, muon1.E())
#       E_lep_a2__theta_3.Fill(angle_3, muon2.E())
#       E_diff__theta_3.Fill(angle_3, E_difference)
#       E_diff__cos_3.Fill(math.cos(angle_3), E_difference)
#       E_diff__tan_3.Fill(math.tan(angle_3), E_difference)



        del oneEvent, myLHEevent
    print(n_4mu, "4mu decays")


    ################
    # SAVE TO FILE #
    ################

    outFile = rt.TFile("4m.root", "RECREATE")

    P_lep_a1.Write()
    E_lep_a1.Write()
#   P_lep_a2.Write()
#   E_lep_a2.Write()
#   P_lep_b1.Write()
#   E_lep_b1.Write()
#   P_lep_b2.Write()
#   E_lep_b2.Write()

#   P_pair_a.Write()
#   E_pair_a.Write()
#   M_pair_a.Write()
#   P_pair_b.Write()
#   E_pair_b.Write()
#   M_pair_b.Write()
    P_trio.Write()
    E_trio.Write()
    M_trio.Write()

    theta_Z.Write()
#   theta_3.Write()
#   theta_b.Write()
#   theta_3_p.Write()
#   theta_b_p.Write()
#   theta_b_pp.Write()

#   E_frac_lep_a2.Write()
#   P_frac_lep_a2.Write()
#   E_frac_pair_b.Write()
#   P_frac_pair_b.Write()
#   E_diff.Write()
#   P_diff.Write()

#   E_lep_a1__theta_3.Write()
#   E_lep_a1__cos2_3.Write()
#   E_lep_a2__theta_3.Write()
#   E_diff__theta_3.Write()
#   E_diff__cos_3.Write()
#   E_diff__tan_3.Write()

    outFile.Close()
