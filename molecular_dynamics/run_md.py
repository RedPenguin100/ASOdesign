import argparse
import os
import re
import subprocess
import sys

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

def check_success(shell):
    shell.sendline('echo __STATUS__$?__')
    shell.expect(r'__STATUS__(\d+)__')

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
    shell.sendline('/home/michael/anaconda3/envs/AmberTools/bin/tleap -s -f usr_leap.in')
    shell.expect('\$ ')  # Wait for the new prompt to appear
    success_status = check_success(shell)
    if success_status != 0:
        raise ValueError("Error")

def main():
    shell = setup_shell()

    AMBER_CLASSIC_PATH = '/home/michael/molecular_dynamics/AmberClassic'

    shell.logfile = sys.stdout

    #######################################
    ##### Setup conda environment
    setup_conda_env(shell, conda_env_name='AmberTools')

    #######################################
    #### AmberClassic
    setup_amber_classic(shell, AMBER_CLASSIC_PATH)

    compile_nab(shell, 'demo_nuc.nab')
    run_tleap(shell)


if __name__ == '__main__':
    main()