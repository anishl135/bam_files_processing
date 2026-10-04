import subprocess
import streamlit as st

def run_command(command, env_name=None):
    """Runs a bash command, activating a conda environment if specified."""
    if env_name:
        # Check if it's an absolute path or a standard env name
        env_flag = "-p" if "/" in env_name else "-n"
        # The eval hook ensures conda activate works inside a subprocess
        full_cmd = f"eval \"$(conda shell.bash hook)\" && conda activate {env_name} && {command}"
    else:
        full_cmd = command
        
    st.text(f"Executing:\n{full_cmd}")
    
    process = subprocess.Popen(
        full_cmd, 
        shell=True, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        executable='/bin/bash'
    )
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        st.error(f"Error encountered:\n{stderr.decode('utf-8')}")
        return False
    else:
        st.success("Command completed successfully!")
        with st.expander("View Execution Logs"):
            st.text(stdout.decode('utf-8'))
            if stderr:
                st.text("Warnings/Info:\n" + stderr.decode('utf-8'))
        return True