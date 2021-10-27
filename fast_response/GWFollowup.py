from .FastResponseAnalysis import PriorFollowup

class GWFollowup(PriorFollowup):
    _dataset = 'GFUOnline_v001_p02'
    _fix_index = False
    _float_index = True
    _index = 2.0
    _pixel_scan_nsigma = 3.0
    
    def run_background_trials(self, ntrials):
        pass

    def write_circular(self):
        pass

    def inject_scan(self):
        pass

    def find_coincident_events(self):
        exp_theta = 0.5*np.pi - self.llh.exp['dec']
        exp_phi   = self.llh.exp['ra']
        exp_pix   = hp.ang2pix(self.nside, exp_theta, exp_phi)
        overlap   = np.isin(exp_pix, self.ipix_90)

        t_mask=(self.llh.exp['time'] <= self.stop) & (self.llh.exp['time'] >= self.start)
        events = self.llh.exp[t_mask]

        events = append_fields(
            events, names=['in_contour', 'ts', 'ns', 'gamma', 'B'],
            data=np.empty((5, events['ra'].size)),
            usemask=False)

        for i in range(events['ra'].size):
            events['in_contour'][i]=overlap[i]
            events['B'][i] = self.llh.llh_model.background(events[i])

        #THIS ISN"T DONE YET

    def upper_limit(self):
        pass

    def ps_sens_range(self):
        pass    

    def per_event_pvalue(self):
        if self.p < 0.05:
            pass
        pass

    def make_dec_pdf(self):
        pass