import ROOT
import math

class H1ComparisonPlotter:
    """ Create plots to compare two histograms
    
    Use different functions to create the plots needed, and draw them. All the needed 
    """
    
    DEFAULT_CONFIG_ = { 'log_yaxis_hist': True,
                        'out_fpath': "",
                        'plot_title': "Fit vs Data",
                        'xaxis_title': "Energy [keV]",
                        'yaxis_title': "Counts/keV",
                        'h1_title': "h1",
                        'h2_title': "h2",
                        'fit': None,
                        'hist_stats': 0,
                      }
    
    def __init__(self, h1, h2, cfg={}, no_errorbars_on2=False):
        self.h1_ = h1
        self.h2_ = h2
        self.no_errorbars_on2 = no_errorbars_on2

        self.h1_.SetBinErrorOption(ROOT.TH1.kPoisson)
        self.config_ = self.DEFAULT_CONFIG_.copy()
        self.set_plot_config(cfg)
        
    def set_plot_config(self, cfg):
        for key in cfg:
            if key in self.config_:
                self.config_[key] = cfg[key]
            else:
                # TODO: change not tested
                print("adding new plot config: ", key)
                self.config_[key] = cfg[key]
                
    def init_pullrat_graphs(self):
        self.g_ratio_ = ROOT.TGraph()
        self.ge_sig1_ = ROOT.TGraphErrors()
        self.ge_sig2_ = ROOT.TGraphErrors()
        self.ge_sig3_ = ROOT.TGraphErrors()
        
        for i in range(self.h1_.GetNbinsX()):
            nh1 = self.h1_.GetBinContent(i+1)
            nh2 = self.h2_.GetBinContent(i+1)
            eh1 = self.h1_.GetBinError(i+1)
            eh2 = self.h2_.GetBinError(i+1)

            ratio = nh1/nh2 if nh2 != 0. else 0.
            sigma = ratio * math.sqrt(eh1*eh1/nh1/nh1 + eh2*eh2/nh2/nh2) if nh1 != 0 and nh2 != 0 else 0.
            mid = 1

            self.g_ratio_.SetPoint(i, self.h1_.GetBinCenter(i+1), ratio)
            self.ge_sig1_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig2_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig3_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig1_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 1.*sigma)
            self.ge_sig2_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 2.*sigma)
            self.ge_sig3_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 3.*sigma)

    def init_pull_graphs(self):
        """plot residuals instead of pulls (h1-h2/h2)""" 

        self.ge_ = ROOT.TGraph()
        self.ge_sig1_ = ROOT.TGraphErrors()
        self.ge_sig2_ = ROOT.TGraphErrors()
        self.ge_sig3_ = ROOT.TGraphErrors()
        for i in range(self.h1_.GetNbinsX()):
            nh1 = self.h1_.GetBinContent(i+1)
            nh2 = self.h2_.GetBinContent(i+1)

            if nh1 == 0:
                self.h1_.SetBinError(i+1, 1)

            res = (nh1-nh2)
            eh1 = self.h1_.GetBinErrorUp(i+1) if res < 0 else self.h1_.GetBinErrorLow(i+1)
            eh2 = self.h2_.GetBinError(i+1) #/nh2 if nh2 != 0 else 1

            sigma = math.sqrt(eh1*eh1 + eh2*eh2) 
            mid = 0
            
            self.ge_.SetPoint(i, self.h1_.GetBinCenter(i+1), res/sigma)
            self.ge_sig1_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig2_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig3_.SetPoint(i, self.h1_.GetBinCenter(i+1), mid)
            self.ge_sig1_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 1)
            self.ge_sig2_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 2)
            self.ge_sig3_.SetPointError (i, 0.5*self.h1_.GetBinWidth(i+1), 3)


    def draw_hist_comparison(self, can1, cfg=None):
        """ draw """
        # TODO: maybe returning the canvas might work?
        # Set aliases
        hares = self.h1_
        hg4 = self.h2_

        hpad = can1.GetPad(1)

        stats_xyndc = [
             [ (0.8, 0.9), (0.8, 0.9)],
             [ (0.8, 0.9), (0.7, 0.8)]
        ]
        if 'stats_xyndc' in cfg:
               stats_xyndc = cfg['stats_xyndc']    

        if(cfg['log_yaxis_hist']): ROOT.gPad.SetLogy()

        # TODO: make better, not tested
        if 'yaxis' in cfg and 'xmax' in cfg['yaxis']:
            set_ymax = cfg['yaxis']['xmax']
        else:
            set_ymax = False

        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetBottomMargin(-0.05)

        hares.SetLineColor(ROOT.kBlue+1)
        hares.SetFillColor(ROOT.kBlue+1)
        hares.SetFillStyle(3002)

        hg4.SetLineColor(ROOT.kRed+1)
        hg4.SetFillColor(ROOT.kRed+1)
        hg4.SetFillStyle(3003)

        hist_stats = cfg['hist_stats'] if 'hist_stats' in cfg else 0
        ROOT.gStyle.SetOptStat(hist_stats)

        hares.Draw("he")
        if self.no_errorbars_on2:
            hg4.Draw("h SAMES")
        else:
            hg4.Draw("he SAMES")
        
        # can1.GetPad(1).Update()
        hpad.Update()

        fit = cfg['fit']
        if(fit):
            ROOT.gStyle.SetOptFit(1111)
            fit_str = fit.get('fit_str', 'gaus')
            print(fit_str)
            
            f1 = ROOT.TF1("f1", fit_str, fit['xmin'], fit['xmax'])
            f2 = ROOT.TF1("f2", fit_str, fit['xmin'], fit['xmax'])

            ipar = 0
            for par in fit.get('pars', []):
                if 'val' in par:
                    f1.SetParameter(ipar, par['val'])
                    f2.SetParameter(ipar, par['val'])
                if 'limits' in par:
                    f1.SetParLimits(ipar, par['limits'][0], par['limits'][1])
                    f2.SetParLimits(ipar, par['limits'][0], par['limits'][1])
                if 'name' in par:
                    f1.SetParName(ipar, par['name'])
                    f2.SetParName(ipar, par['name'])
                
                ipar += 1

            hares.Fit("f1", 'R')
            hg4.Fit("f2", 'R')
            ps1 = hares.GetListOfFunctions().FindObject("stats")
            ps2 = hg4.GetListOfFunctions().FindObject("stats")
            ps1.SetX1NDC(0.7)
            ps1.SetX2NDC(0.9)
            ps1.SetY1NDC(0.65)
            ps1.SetY2NDC(0.9)
            ps2.SetX1NDC(0.7)
            ps2.SetX2NDC(0.9)
            ps2.SetY1NDC(0.4)
            ps2.SetY2NDC(0.65)
        elif(hist_stats):
            hares.SetStats(1)
            hg4.SetStats(1)
            ps1 = hares.GetListOfFunctions().FindObject("stats")
            ps2 = hg4.GetListOfFunctions().FindObject("stats")
            ps1.SetX1NDC(stats_xyndc[0][0][0])
            ps1.SetX2NDC(stats_xyndc[0][0][1])
            ps1.SetY1NDC(stats_xyndc[0][1][0])
            ps1.SetY2NDC(stats_xyndc[0][1][1])
            ps2.SetX1NDC(stats_xyndc[1][0][0])
            ps2.SetX2NDC(stats_xyndc[1][0][1])
            ps2.SetY1NDC(stats_xyndc[1][1][0])
            ps2.SetY2NDC(stats_xyndc[1][1][1])
        else:
            hares.SetStats(0)
            hg4.SetStats(0)
            
        hpad.Update()
        
        if set_ymax:
            hares.GetYaxis().SetRangeUser(0, set_ymax)
        # adjust y axis range if it's not log scale
        elif not cfg['log_yaxis_hist']:
            ymax = max(hares.GetMaximum(), hg4.GetMaximum())
            hares.GetYaxis().SetRangeUser(0, ymax*1.1)
        

        hares.GetXaxis().SetTitle(cfg['xaxis_title'])
        hares.GetYaxis().SetTitle(cfg['yaxis_title'])
        hares.GetXaxis().SetLabelSize(0.04)
        hares.GetYaxis().SetLabelSize(0.04)
        hares.GetXaxis().SetTitleSize(0.04)
        hares.GetYaxis().SetTitleSize(0.04)
        hares.GetYaxis().SetTitleOffset(0.5)
        hares.SetTitle(cfg['plot_title'])

        # print(cfg['xaxis_title'])
        # Legend cannot be drawn here...! 
        # Find it where the canvas+hists are being created

    def draw_pulls_plot(self, can, cfg, pad=0):

        can.cd(pad)
        hmeas = self.h1_
        self.init_pull_graphs()

        ge_ = self.ge_
        ge_sig1 = self.ge_sig1_
        ge_sig2 = self.ge_sig2_
        ge_sig3 = self.ge_sig3_


        # ROOT.gStyle.SetFrameFillStyle(0) # May need to draw a line over
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()


        leg_ratio = ROOT.TLegend(0.54,0.66,0.7,0.9)
        ge_.Draw("AP")
        # ROOT.gPad.Range(2000, -10, 2500, 10)
        # mean_ratio_line = ROOT.TLine(2000, 0, 2500, 0)
        # mean_ratio_line.SetLineColor(ROOT.kBlack)
        # mean_ratio_line.SetLineWidth(2)
        # mean_ratio_line.Draw('same')
        # ROOT.gPad.Update()

        ge_sig3.SetFillColor(ROOT.kRed-4)
        ge_sig3.Draw("2same")

        ge_sig2.SetFillColor(ROOT.kOrange)
        ge_sig2.Draw("2same")

        ge_sig1.SetFillColor(ROOT.kCyan)
        ge_sig1.Draw("2same")
        
        # Calculate the Y axis range...
        g_ymin, g_ymax = min(ge_.GetY()), max(ge_.GetY())
        ge_ymax = max(ge_sig3.GetEY())
        
        pad_ymax = max(g_ymax, 0+ge_ymax)
        pad_ymin = min(g_ymin, 0-ge_ymax)
        # extend the range by 10%
        pad_ymin, pad_ymax = 0.5*(-0.1*pad_ymax + 2.1*pad_ymin), 0.5*(2.1*pad_ymax - 0.1*pad_ymin)
        ge_.GetYaxis().SetRangeUser(pad_ymin,pad_ymax)

        xmin, xmax = hmeas.GetBinLowEdge(1), hmeas.GetBinLowEdge(hmeas.GetNbinsX()+1)
        ge_.GetXaxis().SetTitle(cfg['xaxis_title'])
        ge_.GetYaxis().SetTitle(f"({cfg['h1_title']} - {cfg['h2_title']})/#sigma")
        ge_.GetXaxis().SetLabelSize(0.04)
        ge_.GetYaxis().SetLabelSize(0.04)
        ge_.GetXaxis().SetTitleSize(0.04)
        ge_.GetXaxis().SetRangeUser(xmin, xmax)
        ge_.GetYaxis().SetTitleSize(0.04)
        ge_.GetXaxis().SetTitleOffset(0.9)
        ge_.GetYaxis().SetTitleOffset(0.5)

        ge_.SetTitle(f"Normalized Residual (mean={ge_.GetMean(axis=2):.2f}#pm{ge_.GetRMS(axis=2):.2f})")

        ge_.Draw("P")
        ge_.SetMarkerStyle(21)
        ge_.SetMarkerSize(0.6)

        # leg_ratio.SetFillColor(ROOT.kWhite)
        # leg_ratio.AddEntry(g_ratio,'ratio '+ cfg['h1_title'] + '/' + cfg['h2_title'],"P")
        # leg_ratio.AddEntry(ge_sig1,"1 #sigma","F")
        # leg_ratio.AddEntry(ge_sig2,"2 #sigma","F")
        # leg_ratio.AddEntry(ge_sig3,"3 #sigma","F")
        # leg_ratio.Draw()

        # print('draw_pulls_plot', 'Finished drawing pulls plot')

        pass

    def draw_poisson_errors(self, can, cfg):
        error_points = self.get_bin_probabilities()
        graph = self.g_poisson_errors = ROOT.TGraph()

        for i in range(len(error_points)):
            graph.SetPoint(i, error_points[i][0], error_points[i][1])
            
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()

        graph.SetMarkerStyle(21)
        graph.SetMarkerSize(0.6)

        graph.GetXaxis().SetTitle(cfg['xaxis_title'])
        graph.GetYaxis().SetTitle('Poisson probability')
        graph.GetXaxis().SetLabelSize(0.04)
        graph.GetYaxis().SetLabelSize(0.04)
        graph.GetXaxis().SetTitleSize(0.04)
        graph.GetXaxis().SetRangeUser(self.h1_.GetBinLowEdge(1), self.h1_.GetBinLowEdge(self.h1_.GetNbinsX()+1))
        graph.GetYaxis().SetTitleSize(0.04)
        graph.GetXaxis().SetTitleOffset(0.9)
        graph.GetYaxis().SetTitleOffset(0.8)
        graph.Draw("2AP")

    def draw_comparison_plot(self, can1, cfg=None):

        if not cfg: cfg = self.config_

        can1.Divide(1,2)
        can1.cd(1)
        self.draw_hist_comparison(can1, cfg)
        
        can1.cd(2)
        # self.draw_poisson_errors(can1, cfg)
        self.draw_pulls_plot(can1, cfg, 2)
        # self.draw_goodness_plot(can1, cfg)

        if cfg['out_fpath']: 
            #can1.Draw()
            can1.SaveAs(cfg['out_fpath'])

    def get_bin_probabilities(self):
        hmeas = self.h1_
        hfit = self.h2_
        assert hmeas.GetNbinsX() == hfit.GetNbinsX()

        from ROOT.TMath import Poisson
        error_points = [ [hmeas.GetBinCenter(i+1), Poisson(hmeas.GetBinContent(i+1), hfit.GetBinContent(i+1))]
                     for i in range(hmeas.GetNbinsX())]

        return error_points

def plot_tests(can, plot_cfg, tests, data_hist=None):
    name = plot_cfg['name']
    hist_key = plot_cfg['hist_key']

    # can.cd()

    
    if plot_cfg.get('log_yaxis_hist', None): ROOT.gPad.SetLogy()
            
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetBottomMargin(-0.05)


    hist_stats = plot_cfg['hist_stats'] if 'hist_stats' in plot_cfg else 0
    ROOT.gStyle.SetOptStat(hist_stats)

    leg = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    leg.SetFillColor(ROOT.kWhite)

    # colors = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kBlack]
    fillstyles = [3002, 3003, 3005, 3006]

    cdrawn = False
    for i in range(len(tests)):
        test = tests[i]
        hist = test[hist_key]

        leg.AddEntry(hist, test['name'], "L")

        # hist.SetLineColor(colors[i])
        # hist.SetFillColor(colors[i])
        hist.SetFillStyle(fillstyles[i])

        if cdrawn: hist.Draw('h plc pfc sames')
        else: 
            hist.Draw('h plc pfc')
            cdrawn = True
            
        can.Update()

                
    if data_hist:
        leg.AddEntry(data_hist, 'data', "L")
        data_hist.SetLineColor(ROOT.kBlack)
        data_hist.SetFillStyle(0)
        data_hist.Draw('h sames')

    first_hist = tests[0][hist_key]
    
    # adjust y axis range if it's not log scale
    # if plot_cfg.get('log_yaxis_hist', None):
    #     ymax = max([ test[hist_key].GetMaximum() for test in tests ])
    #     first_hist.GetYaxis().SetRangeUser(0, ymax*1.1)

    first_hist.GetXaxis().SetTitle(plot_cfg.get('xaxis_title', ''))
    first_hist.GetYaxis().SetTitle(plot_cfg.get('yaxis_title', ''))
    first_hist.GetXaxis().SetLabelSize(0.04)
    first_hist.GetYaxis().SetLabelSize(0.04)
    first_hist.GetXaxis().SetTitleSize(0.04)
    first_hist.GetYaxis().SetTitleSize(0.04)
    first_hist.GetYaxis().SetTitleOffset(0.5);
    first_hist.SetTitle(plot_cfg.get('plot_title', 'Comparison'))

    leg.Draw()
    can.Update()

def plot_tests_pad(can, i, plot_cfg, tests, data_hist=None):
    name = plot_cfg['name']
    hist_key = plot_cfg['hist_key']

    # can.cd(i)

    
    if plot_cfg.get('log_yaxis_hist', None): ROOT.gPad.SetLogy()
            
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetBottomMargin(-0.05)


    hist_stats = plot_cfg['hist_stats'] if 'hist_stats' in plot_cfg else 0
    ROOT.gStyle.SetOptStat(hist_stats)

    leg = ROOT.TLegend(0.8, 0.75, 0.9, 0.9)
    leg.SetFillColor(ROOT.kWhite)

    # colors = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kBlack]
    fillstyles = [3351, 3315, 3375, 3357]

    cdrawn = False
    for i in range(len(tests)):
        test = tests[i]
        hist = test[hist_key]

        leg.AddEntry(hist, test['name'], "L")

        # hist.SetLineColor(colors[i])
        # hist.SetFillColor(colors[i])
        hist.SetFillStyle(fillstyles[i])

        if cdrawn: hist.Draw('h plc pfc sames')
        else: 
            hist.Draw('h plc pfc')
            cdrawn = True
            
        ROOT.gPad.Update()

                
    if data_hist:
        leg.AddEntry(data_hist, 'data', "L")
        data_hist.SetLineWidth(2)
        data_hist.SetLineColor(ROOT.kBlack)
        data_hist.SetFillStyle(0)
        data_hist.Draw('h sames')

    first_hist = tests[0][hist_key]
    
    # adjust y axis range if it's not log scale
    # if plot_cfg.get('log_yaxis_hist', None):
    #     ymax = max([ test[hist_key].GetMaximum() for test in tests ])
    #     first_hist.GetYaxis().SetRangeUser(0, ymax*1.1)

    first_hist.GetXaxis().SetTitle(plot_cfg.get('xaxis_title', ''))
    first_hist.GetYaxis().SetTitle(plot_cfg.get('yaxis_title', ''))
    first_hist.GetXaxis().SetLabelSize(0.04)
    first_hist.GetYaxis().SetLabelSize(0.04)
    first_hist.GetXaxis().SetTitleSize(0.04)
    first_hist.GetYaxis().SetTitleSize(0.04)
    first_hist.GetYaxis().SetTitleOffset(0.5);
    first_hist.SetTitle(plot_cfg.get('plot_title', 'Comparison'))

    leg.DrawClone('same')
    # leg.Draw()
    ROOT.gPad.Update()
