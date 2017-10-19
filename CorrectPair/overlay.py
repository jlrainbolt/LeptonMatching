from __future__ import print_function
from __future__ import division
import sys
import rootpy.ROOT as rt
import math
from array import array

if __name__ == '__main__':

    # Set cuts, paths, names
    cutsOn = True
    direc = "/tthome/jaj710/MuonMatching/Mass"
    subdirec = "/CorrectPair/"
    masses = ["5", "18", "35", "50"]
    h_names = ["M4l", "pT1", "pT2", "pT3", "pT4", "eta1", "eta2", "eta3", "eta4"]
    g_names = ["Mll", "pT23", "delr", "dphi", "gamma"]

    # Open all files (that exist)
    if cutsOn:
        cuts = "_cuts"
        msize = 0.25
    else:
        cuts = ""
        msize = 0.1
    paths = [direc + m + subdirec + "plots" + cuts + ".root" for m in masses]

    files = []
    for p in paths:
        try:
            files.append(rt.TFile(p))
        except:
            pass
    NFILES = len(files)
    NHISTS = len(h_names)
    NPLOTS = len(g_names)
    lines = range(NFILES)
    lines = [l + 1 for l in lines]
    colors = range(NFILES)
    colors = [c + 1 for c in colors]

    # Get histograms and graphs
    hists, graphs = [], []
    for n in h_names:
        hists.append([f.Get("Histograms/" + n + "_r") for f in files])
    h_legend = rt.TLegend(0.8, 0.8, 1, 1, "", "brNDC")
    h_legend.SetFillColor(0)
    h_legend.SetTextSize(.03)
    for i in range(NFILES):
        for j in range(NHISTS):
            hists[j][i].SetName(h_names[j] + "_" + masses[i])
            hists[j][i].SetTitle(h_names[j])
            hists[j][i].SetStats(0)
            hists[j][i].SetLineWidth(2)
            hists[j][i].SetLineColor(colors[i])
            hists[j][i].SetLineStyle(lines[i])
        h_legend.AddEntry(hists[0][i], "U_{M} = " + masses[i], "L")

    for n in g_names:
        graphs.append([f.Get("Scatterplots/graph_" + n) for f in files])
    g_legend = rt.TLegend(0.8, 0.8, 1, 1, "", "brNDC")
    g_legend.SetFillColor(0)
    g_legend.SetTextSize(.03)
    for i in range(NFILES):
        graphs[0][i].SetFillColor(colors[i])
        g_legend.AddEntry(graphs[0][i], "U_{M} = " + masses[i], "F")
        for j in range(NPLOTS):
            graphs[j][i].SetName(g_names[j] + "_" + masses[i])
            graphs[j][i].SetTitle(g_names[j])
            graphs[j][i].SetMarkerSize(msize)
            graphs[j][i].SetMarkerColor(colors[i])

    # Make canvases
    h_canvases, g_canvases = [], []
    c_width, c_height = 800, 600
    for n in h_names:
        h_canvases.append(rt.TCanvas(n, n, c_width, c_height))
    for j in range(NHISTS):
        h_canvases[j].cd()
        h_ymin, h_ymax = 0, 0
        for h in hists[j]:
            h.Draw()
            if h_canvases[j].GetFrame().GetY2() > h_ymax:
                h_ymax = h_canvases[j].GetFrame().GetY2()
        hists[j][0].GetYaxis().SetRangeUser(h_ymin, h_ymax)
        hists[j][0].Draw()
        for h in hists[j]:
            h.Draw("SAME")
        h_legend.Draw()

    line = array('d')
    line.append(-1000)
    line.append(1000)
    g_line = rt.TGraph(2, line, line)
    g_line.SetLineColor(rt.kBlack)
    g_line.SetLineWidth(2)
    for n in g_names:
        g_canvases.append(rt.TCanvas(n, n, c_width, c_height))
    for j in range(NPLOTS):
        g_canvases[j].cd()
        g_xmin, g_xmax, g_ymin, g_ymax = float('nan'), float('nan'), float('nan'), float('nan')
        for g in graphs[j]:
            g.Draw("AP")
            if not g.GetXaxis().GetXmin() > g_xmin:
                g_xmin = g.GetXaxis().GetXmin()
            if not g.GetXaxis().GetXmax() < g_xmax:
                g_xmax = g.GetXaxis().GetXmax()
            if not g.GetYaxis().GetXmin() > g_ymin:
                g_ymin = g.GetYaxis().GetXmin()
            if not g.GetYaxis().GetXmax() < g_ymax:
                g_ymax = g.GetYaxis().GetXmax()
        graphs[j][0].GetXaxis().SetRangeUser(g_xmin, g_xmax)
        graphs[j][0].GetYaxis().SetRangeUser(g_ymin, g_ymax)
        graphs[j][0].Draw("AP")
        for g in graphs[j]:
            g.Draw("P")
        g_line.Draw()
        g_legend.Draw()

    # Save histograms, graphs, and canvases to root file
    saveFile = rt.TFile("overlay" + cuts + ".root", "RECREATE")
    h_dir = saveFile.mkdir("Histograms")
    h_dir.cd()
    for c in h_canvases:
        c.Write()
    g_dir = saveFile.mkdir("Scatterplots")
    g_dir.cd()
    for c in g_canvases:
        c.Write()
    m_dir = [saveFile.mkdir("U_M = " + m) for m in masses]
    for i in range(NFILES):
        m_dir[i].cd()
        for j in range(NHISTS):
            hists[j][i].Write()
        for j in range(NPLOTS):
            graphs[j][i].Write()
    saveFile.Close()
    print("Created file", saveFile.GetName())

    # Clean house!
    for l in hists:
        for h in l:
            del h
    for l in graphs:
        for g in l:
            del g
    for c in h_canvases:
        del c
    for c in g_canvases:
        del c
