#!/usr/bin/env python
from subprocess import call
from manage_vm_instances import *

################################################################################
# Function Add_VMs_To_ForwardAgent
#   Note: start_str and end_str needs to be at the bottom of the config file
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
    Add_VMs_To_ForwardAgent()
