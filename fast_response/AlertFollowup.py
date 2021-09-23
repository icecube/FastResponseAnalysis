from .FastResponseAnalysis import PriorFollowup

class AlertFollowup(PriorFollowup):
    pass

class TrackFollowup(AlertFollowup):
    pass

    def _scale_2d_gauss(self, arr, sigma_arr, new_sigma):
        tmp = arr**(sigma_arr**2. / new_sigma**2.)/(np.sqrt(2.*np.pi)*new_sigma)* \
                        np.power(np.sqrt(2.*np.pi)*sigma_arr, (sigma_arr**2. / new_sigma**2.)) 
        return tmp / np.sum(tmp)

class CascadeFollowup(AlertFollowup):
    pass