import logging

def get_log(name):
    """
    Configures logger for Tauren MD.
    """
    
    logging.getLogger(name)
    log.setLevel(logging.DEBUG)
    
    # create a file handler
    debug_ = logging.FileHandler('tauren-md.log')
    debug_.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    
    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(filename)s:%(name)s:%(funcName)s:%(lineno)d - %(message)s')
    debug_.setFormatter(formatter)
    
    # add the handlers to the logger
    log.addHandler(debug_)
    log.addHandler(ch)
    
    return log
