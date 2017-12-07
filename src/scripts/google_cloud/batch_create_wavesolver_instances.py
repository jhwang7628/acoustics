#!/usr/bin/env python
from manage_vm_instances import *

snapshot_name = 'shell-wavesolver-jui-4'
zone_name = 'us-central1-c'
machine_type = 'n1-standard-32'

# instance_id_list = range(17,26)
instance_id_list = [25,26]
print instance_id_list
for iid in instance_id_list:
    instance_name = 'instance-%u' %(iid)
    Create_Disk_From_Snapshot(instance_name, snapshot_name, zone_name)
    Create_Instance_With_Disk(instance_name, instance_name, machine_type)
