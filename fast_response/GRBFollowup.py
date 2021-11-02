from .FastResponseAnalysis import PriorFollowup

class GRBFollowup(PriorFollowup):
    _dataset = 'GFUOnline_v001_p03'
    _fix_index = True
    _index = 2.0

    def run_background_trials(self, ntrials):
        bg_path = 'data/user/efried/IceCube/GRB_coincidence/allsky_scans/'
        pass

    def write_grb_tables(ts_list, ns_list, pval_list, UL_list):
        # """
        # Print latex for a GRB table
        
        # ts_list: unblinded TS values for 10 time windows
        # ns_list: fit ns values for 10 time windows
        # pval_list: unblinded p-values for 10 time windows (corrected for looking in 10 TWs)
        # """
    
        # print '\\begin{table}\label{}'
        # print '    \centering'
        # print '    \caption{The results from all ten time windows.  The time window is centered on the reported $T_{90}$ of the GRB.  P-values have been corrected for searching 10 time windows.}'
        # print '    \\begin{tabular}{| *{1}{p{3.5cm}} *{1}{p{3.5cm}} *{1}{p{2.5cm}} *{1}{p{2.5cm}} *{1}{p{4.0cm}}|}'
        # print '        \hline'
        # print r'       \centering Time Window & \centering Test Statistic & \centering $\hat{n}_s$ & \centering P$_{\mathrm{post}}$ & \ \ Upper Limit (GeV cm$^2$) \\'
        # print '        \hline'
        # for TW_index in range(10):
        
        #     print '        \centering {} & \centering {:.2f} & \centering {:.2f} & \centering {:.2e} & {} {:.2e} \\\ '.format(
        #                                                             tw_labels[TW_index],
        #                                                             ts_list[TW_index], 
        #                                                             ns_list[TW_index], 
        #                                                             pval_list[TW_index], 
        #                                                             '\hspace{1cm}', 
        #                                                             UL_list[TW_index])
        #     print '        \hline'
        # print '    \end{tabular}'
        # print '\end{table}'
        pass

    # index fixed to 2
    # Dataset: with zenith_smoothing
    # Different precomputed trials