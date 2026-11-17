class SampleSkymaps():
    def __init__(self):
        import fast_response
        skymap_path = fast_response.__file__
        skymap_path = skymap_path[:skymap_path.find('__init__')]
        self.skymap_path = skymap_path + 'sample_skymaps/'
        self.sample_gw = self.skymap_path + 'MS181101ab-1-Preliminary.xml'
        self.sample_gw_initial = self.skymap_path + 'S191216ap_initial.xml'
        self.sample_gw_update = self.skymap_path + 'S191216ap_update.xml'
        self.sample_cascade = self.skymap_path + 'sample_cascade.txt'
        self.sample_track = self.skymap_path + 'sample_astrotrack_alert.xml'
        self.set_sample_dict()

    def get_sample_gw(self):
        return self.sample_gw

    def get_sample_initial_gw(self):
        return self.sample_gw_initial

    def get_sample_update_gw(self):
        return self.sample_gw_update

    def get_sample_cascade(self):
        return self.sample_cascade

    def get_sample_track(self):
        return self.sample_track

    def set_sample_dict(self):
        sample_dict = {'GW': self.get_sample_gw(),
                    'GW_initial': self.get_sample_initial_gw(),
                    'GW_update': self.get_sample_update_gw(),
                    'cascade': self.get_sample_cascade(),
                    'track': self.get_sample_track()}
        self.sample_dict = sample_dict