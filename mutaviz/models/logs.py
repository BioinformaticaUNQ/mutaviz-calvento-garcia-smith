import logging


def logger():
    mutaviz_logger = logging.getLogger('mutaviz')
    mutaviz_logger.setLevel(logging.INFO)
    fh = logging.FileHandler('./logs.log')
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    mutaviz_logger.addHandler(fh)
    mutaviz_logger.addHandler(ch)
    return mutaviz_logger
