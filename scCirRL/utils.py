import sys, os
import time
from threading import Thread
# import psutil
import datetime

def open_err_log(log_fn):
    if log_fn:
        try:
            log_fp = open(log_fn, 'a')
            fps = [sys.stderr, log_fp]
        except:
            sys.stderr.write('Error: Can not open log file: {}\n'.format(log_fn))
            fps = [sys.stderr]
    else:
        fps = [sys.stderr]
    return fps

def close_fp(fp):
    if fp is sys.stderr or fp is sys.stdout:
        fp.flush()
    else:
        fp.close()


def start_time():
    global ut_st_time
    ut_st_time = time.time()

def err_log_end_time(log_fn):
    global ut_st_time
    fps = open_err_log(log_fn)
    # process = psutil.Process(os.getpid())
    for fp in fps:
        fp.write('Elapsed time: {}\n'.format(datetime.timedelta(seconds=int(time.time()-ut_st_time))))
        # fp.write('Peak memory : {:.2f} GB\n\n'.format(process.memory_info().rss / 1024 / 1024 / 1024))
        close_fp(fp)

def err_log_progress_bar(log_fn):
    fps = open_err_log(log_fn)
    # 0% to 100%
    for fp in fps:
        fp.write('\n')
        for i in range(11):
            if i == 0 or i == 10:
                fp.write('{}%'.format(i*10))
            else:
                fp.write('{}'.format(i*10))
            fp.write(' ' * 3)
        fp.write('\n')
        # print |----|
        for i in range(10):
            fp.write('|----')
        fp.write('|\n')
        fp.write('*')
        close_fp(fp)

def err_log_progress_star(log_fn, n_total, n_finished, n_existing_stars):
    fps = open_err_log(log_fn)
    n_to_output_stars = int((n_finished / n_total) * 50)
    if n_to_output_stars > n_existing_stars:
        # append stars
        for fp in fps:
            fp.write('*' * (n_to_output_stars - n_existing_stars))
            if n_to_output_stars == 50:
                fp.write('\n\n')
            close_fp(fp)
        n_existing_stars = n_to_output_stars
    return n_existing_stars

def format_time(fp, header='INFO', str=''):
    fp.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')

def out_format_time(header='INFO', str=''):
    sys.stdout.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')

def err_log_format_time(log_fn, header='INFO', str=''):
    fps = open_err_log(log_fn)
    for fp in fps:
        fp.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')
        close_fp(fp)
    if header == 'ERROR':
        sys.exit(1)

def err_format_time(header='INFO', str=''):
    sys.stderr.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')


def fatal_format_time(fp, header='ERROR', str=''):
    format_time(fp, header, str)
    sys.exit(1)


def err_fatal_format_time(header='ERROR', str=''):
    err_format_time(header, str)
    sys.exit(1)


def exec_cmd(fp, header, cmd):
    format_time(fp, header, cmd)
    ret = os.system(cmd)
    # sys.stderr.write('ret: {}'.format(ret) + '\n')
    if ret != 0:
        fatal_format_time(header, 'Error: ' + cmd)


def exec_cmd_parall(fp_list, header_list, cmd_list):
    T = []
    for fp, header, cmd in zip(fp_list, header_list, cmd_list):
        t = Thread(target=exec_cmd, args=(fp, header, cmd))
        t.start()
        T.append(t)
    for t in T:
        t.join()
