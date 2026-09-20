
'''
Logging utilities for FRA
'''
import logging, os

class LogFileWriter(logging.FileHandler):
    '''
    Make a new logging handler that writes to a file
    and flushes all messges to the file as they are written.
    '''
    def emit(self, record):
        super().emit(record)
        self.flush()

class FRA_Logger(object):
    '''
    Define a logging object for use in FRA listeners.
    '''

    def __init__(self, file='', **fmt_kwargs):
        '''
        Initialize a logger object for FRA.

        Parameters:
        -----------
        file (str): 
            Path or file name to log to. If not passed, logs only to stout
        fmt_kwargs (optional):
            Args passed to self.logger_format
        '''
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        self.logger_format(**fmt_kwargs)
        logger.setFormatter(self.logger_fmt)
        if file and os.path.dirname(file):
            if not os.path.exists(os.path.dirname(file)):
                raise Exception('Unable to find parent directory for log file!')
            
            filelogger = LogFileWriter(file, mode='a+')
            filelogger.setFormatter(self.logger_fmt)
            
            # adds the logfile as an additional logger. this will also log to stout
            logger.addHandler(filelogger)

        self.logger = logger

    def logger_format(self, **fmt_kwargs):
        '''
        Format the loggers used in FRA

        Optional arguments:
        -------------------
        fmt (str):
            Full format string to pass to logger.setFormatter. 
            If this is passed, other format options are ignored and only this format is used
            with datefmt if asctime is included in fmt.
        datefmt (str):
            Format for dates to use. Default %Y/%m/%d %H:%M:%S.
        use_time (bool):
            include the time in the format (default True)
        incl_level (bool):
            include the logger level (info, warning, error) (default True)
        '''

        fmt = fmt_kwargs.pop('fmt', None)
        datefmt = fmt_kwargs.pop('datefmt', '%Y/%m/%d %H:%M:%S')
        if fmt is not None:
            self.logger_fmt=logging.Formatter(fmt=fmt, datefmt=datefmt)
            return

        fmt = ''
        if fmt_kwargs.pop('use_time', True):
            fmt = fmt + '[%(asctime)s] '
        if fmt_kwargs.pop('incl_level', True):
            fmt = fmt + '%(levelname)s '
        if fmt: # if not empty, add tab break
            fmt = fmt + '\t '
        fmt = fmt + '%(message)s'

        self.logger_fmt = logging.Formatter(fmt=fmt, datefmt=datefmt)
