
# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - get_gpu_name
# - get_gpu_memory
# - setup_logger
# - logger
# - write_logger
# - log_predictions
# - log_error_sequences

import logging

# I think logger utils should check for write log directives and cache/save those

def get_gpu_name():
    import subprocess
    cmd = "cat $COBALT_NODEFILE"
    #cmd = "cat $PBS_NODEFILE"
    a = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')
    # cmd = "nvidia-smi -L" #"nvidia-smi --query-gpu=gpu_bus_id --format=csv"
    # b = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')

    # return(f"{a}-{b}")
    return(a)

def get_gpu_memory():
    import subprocess
    command = "nvidia-smi --query-gpu=memory.free,memory.used --format=csv"
    memory_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[1].split(',')
    memory_values = [int(x.strip().split(' ')[0]) for i, x in enumerate(memory_info)]
    return memory_values

# def count_number_of_gpus():
#     import subprocess
#     command = "nvidia-smi --list-gpus | wc -l"
#     return()

def setup_logger(logger_name, log_file, level=logging.INFO):
    log_setup = logging.getLogger(logger_name)
    formatter = logging.Formatter('[%(levelname)s: %(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fileHandler = logging.FileHandler(log_file, mode='a')
    fileHandler.setFormatter(formatter)
    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)
    log_setup.setLevel(level)
    log_setup.addHandler(fileHandler)
    log_setup.addHandler(streamHandler)

    return None

def logger(msg, level, logfile):
    if logfile == 'memory'   : log = logging.getLogger('memory_log')
    if logfile == 'cache'   : log = logging.getLogger('cache_log') 
    if logfile == 'run_error'    : log = logging.getLogger('error_log')
    if logfile == 'time'    : log = logging.getLogger('time_log')
    if logfile == 'run_summary' : log = logging.getLogger('summary_log')

    if level == 'info'    : log.info(msg) 
    if level == 'warning' : log.warning(msg)
    if level == 'error'   : log.error(msg)
    if level == 'time'    : log.info(msg)
    #if level == 'summary' : log.info(msg)

    return None

def write_logger(log_msg_type, logfile, message):
    if log_msg_type == 'memory': setup_logger('memory_log', logfile) ; logger(message, 'info', 'memory')
    if log_msg_type == 'cache': setup_logger('cache_log', logfile) ; logger(message, 'info', 'cache')
    if log_msg_type == 'error': setup_logger('error_log', logfile) ; logger(message, 'error', 'run_error')
    if log_msg_type == 'time' : setup_logger('time_log', logfile) ; logger(message, 'info', 'time')
    if log_msg_type == 'info' : setup_logger('summary_log', logfile) ; logger(message, 'info', 'run_summary')
    if log_msg_type == 'warning' : setup_logger('summary_log', logfile) ; logger(message, 'warning', 'run_summary')

def log_predictions(predictions_log_file, what_to_write):
    #print('Writing log file')

    import os, csv
    #logfile_csv = f'{predictions_log_dir}/{id}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(predictions_log_file) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(predictions_log_file, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['region', 'sample', 'status', 'sequence_source', 'predict_time', 'retrieve_time']) # 
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

    return(0)

def log_error_sequences(error_file, what_to_write):
    import os, csv

    open_mode = 'a' if os.path.isfile(error_file) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(error_file, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file, delimiter='\t')
        if open_mode == 'w':
            logwriter.writerow(['locus', 'reason']) # 
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

    return(0)