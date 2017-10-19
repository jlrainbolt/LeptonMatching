import ROOT as rt

                     
class LHEevent():
    
    def __init__(self):
        self.Particles = []
        self.Model = "NONE"
        
    def fillEvent(self, lheLines):
        # check that this is a good event
        #if lheLines[0].find("<event>") == -1 or lheLines[-1].find("</event>") == -1:
            #print "THIS IS NOT A LHE EVENT"
            #return 0
        for i in range(2,len(lheLines)-1):
            self.Particles.append(self.readParticle(lheLines[i]))
        return 1

    def readParticle(self, lheLine):
        dataIN = lheLine[:-1].split(" ")
        dataINgood = []
        for entry in dataIN:
            if entry != "": dataINgood.append(entry)
        return {'ID': int(dataINgood[0]),
                'mIdx': int(dataINgood[2])-1,
                'mIdx2': int(dataINgood[3])-1,
                'Px' : float(dataINgood[6]),
                'Py' : float(dataINgood[7]),
                'Pz' : float(dataINgood[8]),
                'E' : float(dataINgood[9]),
                'M' : float(dataINgood[10]),
                'Q' : float(dataINgood[12])}
