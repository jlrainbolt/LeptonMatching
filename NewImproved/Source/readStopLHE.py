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

    ################
    #  INITIALIZE  #
    ################

    printEvents = False
    printEvery = 100

    MZ = 91.2
    gamma = {'ID':22, 'mIdx':-1, 'mIdx2':-1, 'Px':0.0, 'Py':0.0, 'Pz':0.0, 'E':0.0, 'M':0.0, 'Q':0.0}

    n_gen = 0
    n_AA = n_ZA = n_AZ = n_ZZ = n_AU = n_ZU = n_XX = 0
    n_4e = n_2e2mu = n_4mu = 0

#   M_BINS, M_HMIN, M_HMAX = 50, 0, 250
#   MLL_BINS, MLL_HMIN, MLL_HMAX = 50, 0, 200
#   P_BINS, P_HMIN, P_HMAX = 50, 0, 100

#   M4l_gen = rt.TH1F("M4l_gen", "M4l_gen", M_BINS, M_HMIN, M_HMAX)
#   M4l_acc = rt.TH1F("M4l_acc", "M4l_acc", M_BINS, M_HMIN, M_HMAX)
#   M12_gen = rt.TH1F("M12_gen", "M12_gen", MLL_BINS, MLL_HMIN, MLL_HMAX)
#   M12_acc = rt.TH1F("M12_acc", "M12_acc", MLL_BINS, MLL_HMIN, MLL_HMAX)
#   M34_gen = rt.TH1F("M34_gen", "M34_gen", MLL_BINS, MLL_HMIN, MLL_HMAX)
#   M34_acc = rt.TH1F("M34_acc", "M34_acc", MLL_BINS, MLL_HMIN, MLL_HMAX)
#   Mll_gen, Mll_acc = [M12_gen, M34_gen], [M12_acc, M34_acc]
#   pT_gen = [rt.TH1F("pT"+str(ii+1)+"_gen", "pT"+str(ii+1)+"_gen", P_BINS, P_HMIN, P_HMAX) for ii in range(4)]
#   pT_acc = [rt.TH1F("pT"+str(ii+1)+"_acc", "pT"+str(ii+1)+"_acc", P_BINS, P_HMIN, P_HMAX) for ii in range(4)]
#   h_gen, h_acc = [M4l_gen] + Mll_gen + pT_gen, [M4l_acc] + Mll_acc + pT_acc



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

    myLHEfile = LHEfile(sys.argv[1])
    myLHEfile.setMax(100000)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        n_gen += 1
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        accepted = True
        n_Z = n_U = n_e = n_mu = n_X = 0
        isAA = isZA = isAZ = isZZ = isAU = isZU = isXX = False
        is4e = is2e2mu = is4mu = False



        ################
        #  CATEGORIZE  #
        ################

        # Separate propagators and leptons
        particles = myLHEevent.Particles
        props_idx = []      # indices of propagators
        electrons, muons = [], []
        for i in range(0, len(particles)):
            p = particles[i]
            if p['ID'] == 23 or p['ID'] == 99:
                props_idx.append(i)
            elif abs(p['ID']) == 11:
                electrons.append(p)
            elif abs(p['ID']) == 13:
                muons.append(p)
        props = [particles[i] for i in props_idx]
        leptons = electrons + muons
        lookup = "{0:+.10e}".format(leptons[0]['Px'])

        # Determine final state
        if len(electrons) == 4:                         # 4e
            is4e = True
            n_4e += 1
        elif len(electrons) == 2 and len(muons) == 2:   # 2e2mu
            is2e2mu = True
            n_2e2mu += 1
        elif len(muons) == 4:                           # 4mu
            is4mu = True
            n_4mu += 1



        ################
        # SORT MOTHERS #
        ################

        # No boson in event means both are photons
        if len(props) == 0:                             # AA
            isAA = True
            n_AA += 1
            props_idx = (-1, -1)
            props = (gamma, gamma)
            props[0]['mIdx', 'mIdx2'] = 1, 2

        # One boson in event means second is photon
        elif len(props) == 1:
            props_idx.append(-1)
            props.append(gamma)
            # Determine if photon is first or second
            #   by checking for leptons with quarks
            #   (indices 1 & 2) as their "mothers"
            isAX = False
            for l in leptons:
                if l['mIdx'] != l['mIdx2']:    
                    isAX = True
            # If this is the case, swap order of props
            if isAX:
                props_idx[0], props_idx[1] = props_idx[1], props_idx[0]
                props[0], props[1] = props[1], props[0]
                props[0]['mIdx', 'mIdx2'] = 1, 2
                if props[1]['ID'] == 23:                # AZ
                    isAZ = True
                    n_AZ += 1
                elif props[1]['ID'] == 99:              # AU
                    isAU = True
                    n_AU += 1
            # Otherwise, we know we have Z gamma
            elif props[0]['ID'] == 23:                  # ZA
                    isZA = True
                    n_ZA += 1

        # Two bosons in event--easy!
        elif len(props) == 2:
            if props[0]['ID'] == 23:
                if props[1]['ID'] == 99:                # ZU
                    isZU = True
                    n_ZU += 1
                elif props[1]['ID'] == 23:              # ZZ
                    isZZ = True
                    n_ZZ += 1
                    # First Z must be on shell, so its
                    #   mass must be closer to 91.2
                    if abs(props[0]['M'] - MZ) > abs(props[1]['M'] - MZ):   
                        props_idx[0], props_idx[1] = props_idx[1], props_idx[0]
                        props[0], props[1] = props[1], props[0]

        # ~Mystery process~ just in case this Higgs shows up
        if not any([isAA, isZA, isAZ, isZZ, isAU, isZU]):
            isXX = True
            n_XX += 1



        ##################
        # SORT DAUGHTERS #
        ##################

        prop1_kids, prop2_kids = [], []     # daughters of each propagator
        for l in leptons:
            mom1, mom2 = l['mIdx'], l['mIdx2']
            if mom1 == mom2:
                if mom1 == props_idx[0]:
                    prop1_kids.append(l)
                elif mom1 == props_idx[1]:
                    prop2_kids.append(l)
            else:
                prop1_kids.append(l)


        # The above doesn't work for AA and ZA events
        # Pick P1 daughters as pair with largest Mll???
        if len(prop1_kids) > 2:
            prop1_kids, prop2_kids = prop2_kids, prop1_kids
            candidates = []
            for i in range(len(prop2_kids)):
                for j in range(i):
                    l1, l2 = prop2_kids[i], prop2_kids[j]
                    if abs(l1['ID']) == abs(l2['ID']) and l1['ID'] * l2['ID'] < 0:
                        l1_lv = rt.TLorentzVector(l1['Px'], l1['Py'], l1['Pz'], l1['E'])
                        l2_lv = rt.TLorentzVector(l2['Px'], l2['Py'], l2['Pz'], l2['E'])
                        Mll = (l1_lv + l2_lv).M()
                        candidates.append((Mll, i, j))
#                       candidates.append((abs(Mll - MZ), i, j))        # Mll closest to Z mass
                        # This would really be a separate case for ZA events and I would extract the mass of the actual Z
            candidates.sort(reverse = True)
            prop1_kid_indices = [candidates[0][1], candidates[0][2]]
            for i in prop1_kid_indices:
                prop1_kids.append(prop2_kids.pop(i))
        kids1, kids2 = [l['ID'] for l in prop1_kids], [l['ID'] for l in prop2_kids]
        if printEvents and n_gen % printEvery == 0:
            print(lookup, "\t", props[0]['ID'], "->", kids1, "\t", props[1]['ID'], "->", kids2)





#               n_Z += 1
#           elif abs(p['ID']) == 99:
#               n_U += 1
#           elif abs(p['ID']) == 11: 
#               e_lv.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
#               e_q.append(p['Q'])
#               l_lv.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
#               l_q.append(p['Q'])
#               l_id.append(abs(p['ID']))
#           elif abs(p['ID']) == 13:
#               mu_lv.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
#               mu_q.append(p['Q'])
#               l_lv.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
#               l_q.append(p['Q'])
#               l_id.append(abs(p['ID']))
#           else:
#               n_X += 1
#       n_e, n_mu = len(e_lv), len(mu_lv)

#       # Arrange lepton pairs by pT
#       if l_lv[0].Pt() < l_lv[1].Pt():
#           l_lv[0], l_lv[1] = l_lv[1], l_lv[0]
#           l_q[0], l_q[1] = l_q[1], l_q[0]
#           l_id[0], l_id[1] = l_id[1], l_id[0]
#       if l_lv[2].Pt() < l_lv[2].Pt():
#           l_lv[2], l_lv[3] = l_lv[3], l_lv[2]
#           l_q[2], l_q[3] = l_q[3], l_q[2]
#           l_id[2], l_id[3] = l_id[3], l_id[2]

#       if (l_id[0] != l_id[1] or l_id[2] != l_id[3]):
#           print("Flavor mismatch!")


        ################
        #  FILL HISTS  #
        ################

#       # Labels for printing
#       if props[0]['ID'] == 22:V
#           prop1 = "A"
#       elif props[0]['ID'] == 23:
#           prop1 = "Z"
#       else:
#           prop1 = "?"
#       if props[1]['ID'] == 22:
#           prop2 = "A"
#       elif props[1]['ID'] == 23:
#           prop2 = "Z"
#       elif props[1]['ID'] == 99:
#           prop2 = "U"
#       else:
#           prop2 = "?"
#       proc = prop1 + prop2

        del oneEvent, myLHEevent
    if not printEvents:
        print("done!")


    ################
    # PRINT RESULT #
    ################

    print("\n")
    print("DECAY SEQUENCES")
    print("=================================================================")
    print("Total", "AA", "AZ", "ZA", "ZZ", "AU", "ZU", "??", sep = "\t")
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print(n_gen, n_AA, n_AZ, n_ZA, n_ZZ, n_AU, n_ZU, n_XX, sep = "\t")
    print(n_gen/n_gen, n_AA/n_gen, n_AZ/n_gen, n_ZA/n_gen, n_ZZ/n_gen, n_AU/n_gen, n_ZU/n_gen, n_XX/n_gen, sep = "\t")
    print("\n")
    print("SUBPROCESSES")
    print("=================================================================")
    print("Total", "4e", "2e2mu", "4mu", sep = "\t")
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print(n_gen, n_4e, n_2e2mu, n_4mu, sep = "\t")
    print(n_gen/n_gen, n_4e/n_gen, n_2e2mu/n_gen, n_4mu/n_gen, sep = "\t")
    print("\n")
