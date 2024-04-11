import os
# THIS FILE DOESN'T WORK YET
# AT SOME POINT, USE THIS TO SPAWN A BUNCH OF JOBS
# AND PUT THE DATA IN A UNIQUE FOLDER
# STORE ONE CONFIG FOR EACH JOB (with the different command for the benchmark, and inputs)

# keep trying to make data1, data2, etc. until one doesn't exist
ctr = 1
path = f"data{ctr}"
while
    ctr += 1
    path = 

from multiprocessing import Pool

def run_command(path):
    command = "spark-submit driver.py -i {}".format(path)
    subprocess.Popen(command, shell=True)

pool = Pool()
pool.map(run_command, paths)