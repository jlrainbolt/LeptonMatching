from __future__ import print_function
from __future__ import division
import sys
import ROOT as rt
from array import array
from LHEevent import *
from LHEfile import *


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


    # LHE file and number of events
    lheName = "unweighted_events.lhe"
    n_max = 100000


    # Counters
    n_evts = n_4m = n_2m2e = n_2e2m = n_4e = 0


    # Tree and ROOT file
    outName = "Zto4l_decays.root"
    outFile = rt.TFile(outName, "RECREATE")
    tree = rt.TTree("tree", "Z to 4l decays")


    # Create objects for branches
    lep1p4, lep1q, lep1id = rt.TLorentzVector(), array("i", [0]), array("i", [0]) 
    lep2p4, lep2q, lep2id = rt.TLorentzVector(), array("i", [0]), array("i", [0]) 
    lep3p4, lep3q, lep3id = rt.TLorentzVector(), array("i", [0]), array("i", [0]) 
    lep4p4, lep4q, lep4id = rt.TLorentzVector(), array("i", [0]), array("i", [0]) 
    pair1id, pair2id = array("i", [0]), array("i", [0])

    
    # Add branches to tree
    tree.Branch("pair1id", pair1id, "pair1id/I")
    tree.Branch("pair2id", pair2id, "pair2id/I")
    tree.Branch("lep1p4", lep1p4)
    tree.Branch("lep1q", lep1q, "lep1q/I")
    tree.Branch("lep1id", lep1id, "lep1id/I")
    tree.Branch("lep2p4", lep2p4)
    tree.Branch("lep2q", lep2q, "lep2q/I")
    tree.Branch("lep2id", lep2id, "lep2id/I")
    tree.Branch("lep3p4", lep3p4)
    tree.Branch("lep3q", lep3q, "lep3q/I")
    tree.Branch("lep3id", lep3id, "lep3id/I")
    tree.Branch("lep4p4", lep4p4)
    tree.Branch("lep4q", lep4q, "lep4q/I")
    tree.Branch("lep4id", lep4id, "lep4id/I")


    # Apparently this must be done to ensure the TLorentzVectors "live long enough":
    #   https://root-forum.cern.ch/t/tlorentzvector-in-a-ttree/18795
#   tree._lep1p4 = lep1p4
#   tree._lep2p4 = lep2p4
#   tree._lep3p4 = lep3p4
#   tree._lep4p4 = lep4p4



    ################
    #  EVENT LOOP  #
    ################


    myLHEfile = LHEfile(lheName)
    myLHEfile.setMax(n_max)
    
    print("Looping over", n_max, "events...", end='')
    sys.stdout.flush()

    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:

        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_evts += 1




        ################
        #  CATEGORIZE  #
        ################


        # Get lepton info
        particles = myLHEevent.Particles
        elecs, muons = [], []
        leptons_q, leptons_id = {}, {}

        for i in range(0, len(particles)):
            p = particles[i]
            if abs(p['ID']) == 11:
                elecs.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                leptons_q[elecs[-1]] = int(p['Q'])
                leptons_id[elecs[-1]] = 11
            elif abs(p['ID']) == 13:
                muons.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                leptons_q[muons[-1]] = int(p['Q'])
                leptons_id[muons[-1]] = 13


        # Determine event type
        if len(muons) == 4:
            n_4m += 1
            pair1id[0] = 13
            pair2id[0] = 13

        elif len(elecs) == 4:
            n_4e += 1
            pair1id[0] = 11
            pair2id[0] = 11

        elif len(muons) == 2 and len(elecs) == 2:
            elec_pair = GetSum(elecs)
            muon_pair = GetSum(muons)

            if muon_pair.M() > elec_pair.M():
                n_2m2e += 1
                pair1id[0] = 13
                pair2id[0] = 11

            else:
                n_2e2m += 1
                pair1id[0] = 11
                pair2id[0] = 13

        else:
            print("Issue with number of leptons")


        # Combine flavors and sort by energy
        leptons = muons + elecs
        leptons.sort(key=lambda fourvec: fourvec.E(), reverse=True)




        ###############
        #  FILL TREE  #
        ###############


        # Simply writing "lep1p4 = leptons[0]" seems to produce a seg fault :(
        lep1p4.SetPxPyPzE(leptons[0].Px(), leptons[0].Py(), leptons[0].Pz(), leptons[0].E())
        lep2p4.SetPxPyPzE(leptons[1].Px(), leptons[1].Py(), leptons[1].Pz(), leptons[1].E())
        lep3p4.SetPxPyPzE(leptons[2].Px(), leptons[2].Py(), leptons[2].Pz(), leptons[2].E())
        lep4p4.SetPxPyPzE(leptons[3].Px(), leptons[3].Py(), leptons[3].Pz(), leptons[3].E())

        lep1q[0], lep1id[0] = leptons_q[leptons[0]], leptons_id[leptons[0]]
        lep2q[0], lep2id[0] = leptons_q[leptons[1]], leptons_id[leptons[1]]
        lep3q[0], lep3id[0] = leptons_q[leptons[2]], leptons_id[leptons[2]]
        lep4q[0], lep4id[0] = leptons_q[leptons[3]], leptons_id[leptons[3]]

        tree.Fill()


        # For debugging
        if (False):
            print(n_evts)
            print(lep1p4, lep1q[0], lep1id[0])
            print(lep2p4, lep2q[0], lep2id[0])
            print(lep3p4, lep3q[0], lep3id[0])
            print(lep4p4, lep4q[0], lep4id[0])
            print(pair1id[0], pair2id[0])
            print("")
            sys.stdout.flush()


        del oneEvent, myLHEevent


    # Print output message
    print("done!")
    print("\nFound:")
    print("\t", n_4m, "4mu decays")
    print("\t", n_2m2e, "2mu2e decays")
    print("\t", n_2e2m, "2e2mu decays")
    print("\t", n_4e, "4e decays")
    print("\nWrote events to", outName)




    ###############
    #  SAVE FILE  #
    ###############

    outFile.Write()
    outFile.Close()
