import argparse
import os
import re
import subprocess
import sys
import time

import pexpect

def setup_conda_env(shell, conda_env_name='AmberTools'):
    shell.sendline(f'source ~/anaconda3/etc/profile.d/conda.sh')
    shell.expect(r'\$ ')

    shell.sendline(f"conda activate {conda_env_name}")
    shell.expect(r'\$ ')

    # return to simple shell with $
    shell.sendline('PS1="$ "') # setup a simple promptshell.expect('\$ ')  # Wait for the new prompt to appear
    shell.expect('\$ ')  # Wait for the new prompt to appear

    return shell


def compile_nab(shell, nab_path):
    print("C")
    shell.sendline(f'./bin/nab {nab_path}')
    shell.expect(r'\$ ', timeout=5)

    print("D")
    returncode = check_success(shell)

    if returncode != 0:
        raise ValueError(f"Error code: {returncode}")

    shell.sendline("./a.out")
    shell.expect(r'\$ ')

def check_success(shell, timeout=10):
    # convoluted STATUS message to filter accurately
    shell.sendline('echo __STATUS__$?__')
    shell.expect(r'__STATUS__(\d+)__', timeout=timeout)

    success_status = int(shell.match.group(1)) # Match the number in the regex
    return success_status

def setup_amber_classic(shell, amber_classic_path):
    shell.sendline(f'source {os.path.join(amber_classic_path, "AmberClassic.sh")}')
    shell.expect(r'\$ ')

    success_status = check_success(shell)

    print(shell.before.strip())
    if success_status != 0:
        raise ValueError(shell.stderr.strip())

    # if shell.exitstatus != 0:
    #     # print("Output:", shell.stdout)
    #     raise ValueError(shell.stderr)

def setup_shell():
    # --noprofile and --norc to remove colors
    shell = pexpect.spawn("/bin/bash --noprofile --norc", encoding='utf-8', echo=False)
    shell.logfile = sys.stdout

    shell.sendline('PS1="$ "') # setup a simple promptshell.expect('\$ ')  # Wait for the new prompt to appear

    shell.expect('\$ ')  # Wait for the new prompt to appear

    shell.sendline('echo -e "\\e[?2004l"')
    shell.expect('\$ ')  # Wait for the new prompt to appear
    return shell


def run_tleap(shell):
    shell.sendline('/home/michael/anaconda3/envs/AmberTools/bin/tleap -s -f demo/usr_leap.in')
    shell.expect('\$ ')  # Wait for the new prompt to appear
    success_status = check_success(shell)
    if success_status != 0:
        raise ValueError("Error")

def run_minimization1(shell,driver):
    shell.sendline(f'{driver} -O -i demo/polyAT_wat_min1.in'
                   ' -o polyAT_wat_min1.out -p polyAT_wat.prmtop -c polyAT_wat.rst7 -r polyAT_wat_min1.ncrst -ref polyAT_wat.rst7')
    timeout = 600
    shell.expect('\$ ', timeout=timeout)

    success_status = check_success(shell, timeout=timeout) # minimization is slow, increase timeout for expect
    if success_status != 0:
        raise ValueError("Error")


def run_minimization2(shell, driver):
    shell.sendline(f'{driver} -O -i demo/polyAT_wat_min2.in -o polyAT_wat_min2.out'
                   ' -p polyAT_wat.prmtop -c polyAT_wat_min1.ncrst -r polyAT_wat_min2.ncrst')
    timeout = 600
    shell.expect('\$ ', timeout=timeout)

    success_status = check_success(shell, timeout=timeout) # minimization is slow, increase timeout for expect
    if success_status != 0:
        raise ValueError(f"Error {success_status}")

def run_minimization3(shell, driver):
    shell.sendline(f'{driver} -O -i demo/polyAT_wat_md1.in -o polyAT_wat_md1.out -p polyAT_wat.prmtop'
                   ' -c polyAT_wat_min2.ncrst -r polyAT_wat_md1.ncrst -x polyAT_wat_md1.nc -ref polyAT_wat_min2.ncrst')
    timeout = 600
    shell.expect('\$ ', timeout=timeout)

    success_status = check_success(shell, timeout=timeout) # minimization is slow, increase timeout for expect
    if success_status != 0:
        raise ValueError(f"Error {success_status}")

def run_minimization4(shell, driver):
    shell.sendline(f'{driver} -O -i demo/polyAT_wat_md2.in -o polyAT_wat_md2.out '
                   '-p polyAT_wat.prmtop -c polyAT_wat_md1.ncrst -r polyAT_wat_md2.ncrst -x polyAT_wat_md2.nc')
    timeout = 3600
    shell.expect('\$ ', timeout=timeout)
    success_status = check_success(shell, timeout=timeout)
    if success_status != 0:
        raise ValueError(f'Error {success_status}')


def main():
    cpu_driver = 'home/michael/anaconda3/envs/AmberTools/bin/sander'
    gpu_driver = '/home/michael/molecular_dynamics/install/pmemd24_src/build/src/pmemd/src/pmemd.cuda_DPFP'

    driver = gpu_driver


    shell = setup_shell()
    shell.logfile = sys.stdout

    AMBER_CLASSIC_PATH = '/home/michael/molecular_dynamics/AmberClassic'

    setup_conda_env(shell, conda_env_name='AmberTools')
    setup_amber_classic(shell, AMBER_CLASSIC_PATH)


    compile_nab(shell, 'demo/demo_nuc.nab')
    run_tleap(shell)
    start = time.time()
    run_minimization1(shell, driver)
    end = time.time()
    print(f"Min1 took: {end - start}")
    start = time.time()
    run_minimization2(shell, driver)
    end = time.time()
    print(f"Min2 took: {end - start}")
    start = time.time()
    run_minimization3(shell, driver)
    end = time.time()
    print(f"Min3 took: {end - start}")
    start = time.time()
    run_minimization4(shell, driver)
    end = time.time()
    print(f"Min4 took: {end - start}")


if __name__ == '__main__':
    main()