#!/usr/bin/env python 
from subprocess import call
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
# Function Main
################################################################################
if __name__ == '__main__': 
    snapshot_name = 'shell-wavesolver-jui-3'
    zone_name = 'us-central1-c'
    machine_type = 'n1-standard-16'
    
    instance_name = 'instance-6'
    # Create_Disk_From_Snapshot(instance_name, snapshot_name, zone_name)
    # Create_Instance_With_Disk(instance_name, instance_name, machine_type)

