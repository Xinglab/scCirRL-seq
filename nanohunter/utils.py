import sys, os
import time
from threading import Thread


def format_time(fp, header, str):
    fp.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')


def out_format_time(header, str):
    sys.stdout.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')

def err_log_format_time(log_fn, header, str):
    sys.stderr.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')
    if log_fn:
        try:
            log_fp = open(log_fn, 'a')
            log_fp.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')
            log_fp.close()
        except:
            pass

def err_format_time(header, str):
    sys.stderr.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str + '\n')


def fatal_format_time(fp, header, str):
    format_time(fp, header, str)
    sys.exit(1)


def err_fatal_format_time(header, str):
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
