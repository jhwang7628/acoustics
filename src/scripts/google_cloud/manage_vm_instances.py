#!/usr/bin/env python 
from subprocess import call
import os

################################################################################
# Class VM
################################################################################
class VM: 
    def __init__(self): 
        self.name = ''
        self.zone = ''
        self.machine_type = ''
        self.preemptible = False 
        self.external_ip = ''
        self.status = ''
    def __str__(self): 
        s = 'VM instance %s: %s' %(self.name, self.external_ip)
        return s

################################################################################
# Execute_Shell_Cmd
################################################################################
def Execute_Shell_Cmd(cmd): 
    print cmd
    call(cmd, shell=True)

################################################################################
# Function Create_Disk_From_Snapshot
################################################################################
def Create_Disk_From_Snapshot(disk_name, snapshot_name, zone_name='us-central1-c'): 
    cmd = 'gcloud compute disks create %s --source-snapshot %s --type=pd-ssd --zone=%s' %(disk_name, snapshot_name, zone_name)
    Execute_Shell_Cmd(cmd)

################################################################################
# Function Create_Instance_With_Disk
################################################################################
def Create_Instance_With_Disk(instance_name, disk_name, machine_type='n1-standard-16'):
    cmd = 'gcloud compute instances create %s --machine-type %s --disk name=%s,boot=yes' %(instance_name, machine_type, disk_name)
    Execute_Shell_Cmd(cmd)

################################################################################
# Function Parse_Instances_From_File
################################################################################
def Parse_Instances_From_File(filename): 
    VMs = []
    with open(filename, 'r') as stream: 
        lines = stream.readlines()
        for line in lines[1:]:
            tokens = line.split()
            vm = VM()
            vm.name = tokens[0]
            vm.zone = tokens[1]
            vm.machine_type = tokens[2]
            # vm.preemptible = tokens[3]
            vm.external_ip = tokens[-2]
            vm.status = tokens[-1]
            VMs.append(vm)
    return VMs

################################################################################
# Function Query_Instances
################################################################################
def Query_Instances(query_file = None):
    # tmp filename
    count = 0
    while (True): 
        tmp_file = '_tmp_%u.txt' %(count)
        if os.path.isfile(tmp_file): 
            count += 1
        else: 
            break 
    if query_file is not None: 
        tmp_file = query_file

    # write instance state to file
    cmd = 'gcloud compute instances list > %s' %(tmp_file)
    Execute_Shell_Cmd(cmd) 

    # read/parse status from file
    VMs = Parse_Instances_From_File(tmp_file)
    for vm in VMs: 
        print vm

    # clear tmp file 
    if query_file is None: 
        cmd = 'rm %s' %(tmp_file)
        Execute_Shell_Cmd(cmd)

    return VMs

################################################################################
# Function Add_VMs_To_ForwardAgent
################################################################################
def Add_VMs_To_ForwardAgent():
    ##
    config_file = '/home/jui-hsien/.ssh/config'
    start_str = '### GCP automatic START ###'
    end_str = '### GCP automatic END ###'
    ##
    with open(config_file) as stream: 
        lines = stream.readlines()
        start_line_idx = 0
        end_line_idx = 0
        count = 0
        for l in lines: 
            if l.find(start_str)!=-1: start_line_idx = count
            elif l.find(end_str)!=-1: end_line_idx = count
            count += 1

    print '\nEditing ssh config file to allow forwarding for the current set of VMs'
    VMs = Query_Instances() 
    with open(config_file, 'w') as stream: 
        for ii in range(start_line_idx): 
            stream.write(lines[ii])
        stream.write(start_str)
        stream.write('\n')
        for vm in VMs: 
            stream.write('Host %s\n' %(vm.external_ip))
            stream.write('  ForwardAgent yes\n')
        stream.write(end_str)
        stream.write('\n')

################################################################################
# Function Main
################################################################################
if __name__ == '__main__': 

    # create vm instances
    # snapshot_name = 'shell-wavesolver-jui-3'
    # zone_name = 'us-central1-c'
    # machine_type = 'n1-standard-16'
    # 
    # instance_name = 'instance-6'
    # Create_Disk_From_Snapshot(instance_name, snapshot_name, zone_name)
    # Create_Instance_With_Disk(instance_name, instance_name, machine_type)

    # query instance status
    Query_Instances() 

    # Add_VMs_To_ForwardAgent()
