""" Run bash scripts
"""

import os
import subprocess
import warnings
import stat


def run_script(script_str, run_dir):
    """ run a program from a script
    """

    script_name = 'build.sh'
    with _EnterDirectory(run_dir):
        # Write the submit script to the run directory
        try:
            os.remove('build.sh')
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)
        except:
            with open(script_name, 'w') as script_obj:
                script_obj.write(script_str)

        # Make the script executable
        os.chmod(script_name, mode=os.stat(script_name).st_mode | stat.S_IEXEC)

        # Call the program
        try:
            subprocess.check_call('./{:s}'.format(script_name))
        except subprocess.CalledProcessError:
            # If the program failed, continue with a warning
            warnings.warn("run failed in {}".format(run_dir))


class _EnterDirectory():

    def __init__(self, directory):
        assert os.path.isdir(directory)
        self.directory = directory
        self.working_directory = os.getcwd()

    def __enter__(self):
        os.chdir(self.directory)

    def __exit__(self, _exc_type, _exc_value, _traceback):
        os.chdir(self.working_directory)
