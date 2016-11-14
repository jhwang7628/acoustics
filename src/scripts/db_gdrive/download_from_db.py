#!/usr/bin/env python 
import modeldb
from modeldb import create_service,get_folder_id,get_all_models,download_file
import os

################################################################################
## Private
################################################################################
def GuardedMkdir(directory): 
    if not os.path.exists(directory): 
        os.makedirs(directory)
    return directory

def WriteToFile(service, file_id, local_fd):
    download_file(service, file_id, local_fd)

################################################################################
## Public API
################################################################################
def Download(service, download_dir, download_files, model): 
    children = service.children().list(folderId=model['id']).execute()
    # print 'children=', children
    children = children.get('items', [])
    all_ids = [child['id'] for child in children]
    all_files = [service.files().get(fileId=f).execute() for f in all_ids]

    for m in all_files:
        if m['title'] in download_files:
            dirname = GuardedMkdir('%s/%s' %(download_dir, model['title']))
            filename = '%s/%s' %(dirname, m['title'])
            if not os.path.exists(filename):
                with open(filename, 'w+') as out_file: 
                    print ' - Downloading file %s -> %s' %(m['title'], filename)
                    WriteToFile(service, m['id'], out_file)
            else: 
                print ' - File exists, skipping: %s' %(filename)
        # if m['title'] == 'proj.obj':
            # # checks to see if image already exists for item
            # image_dir = os.path.join(os.getcwd(), 'imgs', res['name']) + '.png'
            # if not os.path.exists(image_dir):
            #     with open('temp/' + res['name'] + '.obj', 'w+') as f:
            #         download_file(service, m['id'], f)
            #     # paraview is weird so you have to call it like this
            #     # just take my word for it
            #     call('python sser.py ' + res['name'], shell=True)
            #     # now delete the obj since we don't need it anymore
            #     os.remove('temp/' + res['name'] + '.obj')

        # # deal with sound-gen.cfg
        # if m['title'] == 'sound-gen.cfg':
        #     contents = get_file_content(service, m['id'])
        #     fields = get_soundgencfg_fields(service, contents)
        #     # adds a dictionary to another
        #     res = dict(res, **fields)

    # # add images
    # res['image'] = 'imgs/' + res['name'] + '.png'
    # return res

if __name__ == '__main__': 
    ## customized fields 
    download_dir = 'downloaded'
    query_models = [
        'p_B',
        'p_O',
        'p_U',
        'p_N',
        'p_C',
        'p_E',
        'p_M',
        'p_A',
        'p_P',
    ]
    download_files = [
        'proj.obj',
        #'proj.tet',
        'proj.modes',
        #'proj.geo.txt'
    ]

    ## 
    service = create_service()
    print('Finding folder modal_models...')
    modal_models = get_folder_id(service, 'modal_models')
    print('Finding all models...')
    models = get_all_models(service, modal_models)
    N_models = len(models)
    N_query_models = len(query_models)
    print('Number of models: %u' %(N_models))

    query_found = [False] * N_query_models
    for j in range(N_query_models): 
        for i in range(N_models):
            if models[i]['title'] == query_models[j]:
                query_found[j] = True
                print 'Found query model %s in the database' %(query_models[j])
                Download(service, download_dir, download_files, models[i])
    print('Done')
    
