# from __future__ import print_function
import httplib2
import os
from subprocess import call

from apiclient import discovery, errors, http
import oauth2client
from oauth2client import client
from oauth2client import tools

try:
    import argparse
    flags = argparse.ArgumentParser(parents=[tools.argparser]).parse_args()
except ImportError:
    flags = None

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/drive-python-quickstart.json
SCOPES = 'https://www.googleapis.com/auth/drive.readonly'
CLIENT_SECRET_FILE = 'client_secret.json'
APPLICATION_NAME = 'MODEL TABLE GENERATOR'

TYPES = ["models",
         "ellipsoids",
         "rocks"
]

# relevant files
# order in the list is the order in the table
HEADERS = ["name",
           "image",
           "proj.obj",
           "proj_extracted.obj",
           "proj.modes",
           "proj_surf.modes",
           "proj.mass.csc",
           "proj.stiff.csc",
           "proj_mass.vector",
           "proj_materialProperties.vector",
           "proj_centerOfMass.3vector",
           "proj_inertia.3matrix",
           "drop.cfg"
]

# relevant fields, grabbing entire drop.cfg file
SOUND_GEN_CFG_FIELDS = ["alpha",
                        "beta",
                        "density"
]

HEADERS_SET = set(HEADERS)
SOUND_GEN_SET = set(SOUND_GEN_CFG_FIELDS)

def get_credentials():
    """
    Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.

    Taken from the Google Drive v2 Python quickstart guide
    """
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir,
        'drive-python-quickstart.json')
    store = oauth2client.file.Storage(credential_path)
    credentials = store.get()
    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else: # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print('Storing credentials to ' + credential_path)
    return credentials


def create_service():
    credentials = get_credentials()
    http = credentials.authorize(httplib2.Http())
    return discovery.build('drive', 'v2', http=http)


def get_folder_id(service, folder_name):
    results = service.files().list(q="title='"+folder_name+"'").execute()
    items = results.get('items', [])

    if len(items) == 0:
        raise Exception, 'modal_models folder not found'
    if len(items) != 1:
        raise Exception, 'more than one folder named modal_models'

    return items[0]['id']


def get_all_models(service, folder_id):
    """
    Gets all models in the given folder id
    """
    try:
        children = service.children().list(folderId=folder_id).execute()
        all_model_ids = [child['id'] for child in children.get('items', [])]
        return [service.files().get(fileId=f).execute() for f in all_model_ids]
    except errors.HttpError, error:
        print('An error occurred: %s' % error)


def generate_table_dict(service, model):
    """
    Generates a dictionary where the key is the table header, and the
    value is the corresponding link/string (in the case of name) of the
    key.

    TODO: deal with ellipsoids, renderData, rocks, and tets folders
    """
    # ignore models that don't start with p
    if model['title'][:2] != 'p_':
        return

    res = {'name': model['title'],
           'image': None,#'image not made'
    }

    # get all children of the folder
    children = service.children().list(folderId=model['id']).execute()
    children = children.get('items', [])
    # convert the list to a list of file objects
    # does the same thing in get_all_models
    all_ids = [child['id'] for child in children]
    all_files = [service.files().get(fileId=f).execute() for f in all_ids]

    # looks for the view links for the files
    for m in all_files:
        if m['title'] in HEADERS_SET:
            res[m['title']] = m['alternateLink']
        if m['title'] == 'proj.obj':
            # checks to see if image already exists for item
            image_dir = os.path.join(os.getcwd(), 'imgs', res['name']) + '.png'
            if not os.path.exists(image_dir):
                with open('temp/' + res['name'] + '.obj', 'w+') as f:
                    download_file(service, m['id'], f)
                # paraview is weird so you have to call it like this
                # just take my word for it
                call('python sser.py ' + res['name'], shell=True)
                # now delete the obj since we don't need it anymore
                os.remove('temp/' + res['name'] + '.obj')

        # deal with sound-gen.cfg
        if m['title'] == 'sound-gen.cfg':
            contents = get_file_content(service, m['id'])
            fields = get_soundgencfg_fields(service, contents)
            # adds a dictionary to another
            res = dict(res, **fields)

    # add images
    res['image'] = 'imgs/' + res['name'] + '.png'
    return res


def get_soundgencfg_fields(service, contents):
    """
    Gets the relevant fields in sound-gen.cfg as specified in SOUND_GEN_CFG_
        FIELDS

    Args:
        service: Drive API Service instance
        contents: contents of sound-gen.cfg

    Returns:
        Dictionary where the key is a relevant field, and the value is the
        value of the field.

    Quick and dirty way of doing it, could probably make it easier
    Could also add comment support
    """
    res = {}
    contents = contents.split('\n')
    for line in contents:
        for field in SOUND_GEN_CFG_FIELDS:
            if field in line:
                # yea... hopefully it's always in a standard format
                data = line[line.find('=')+1:line.find(';')].strip()
                res[field] = data

    return res


def download_file(service, file_id, local_fd):
    """
    Download a Drive file's contents to the local filesystem

    Args:
        service: Drive API Service instance
        file_id: ID of the drive file that will be downloaded
        local_fd: io.Base or file object, the stream that the Drive
            file's contents will be written to.

    Taken from https://developers.google.com/drive/v2/reference/files/get
    """
    request = service.files().get_media(fileId=file_id)
    media_request = http.MediaIoBaseDownload(local_fd, request)

    while True:
        try:
            download_progress, done = media_request.next_chunk()
        except errors.HttpError, error:
            print('**ERROR** An error occurred during download: %s' % error)
            return
        if download_progress:
            print('     Progress: %d%%... ' % int(download_progress.progress() * 100)), 
        if done:
            print('[DONE]')
            return

def print_file_content(service, file_id):
    """
    Print a file's content.

    Args:
        service: Drive API service instance
        file_id: ID of the file

    Returns:
        File's content if successful, None otherwise.

    Liberally taken from
        https://developers.google.com/drive/v2/reference/files/get
    """
    try:
        print service.files().get_media(fileId=file_id).execute()
    except errors.HttpError, error:
        print('An error occurred: %s' % error)


def get_file_content(service, file_id):
    """
    Get a file's content.

    Args:
        service: Drive API service instance
        file_id: ID of the file

    Returns:
        File's content if successful, None otherwise.

    Taken with slight modifications from
        https://developers.google.com/drive/v2/referance/files/get
    """
    try:
        return service.files().get_media(fileId=file_id).execute()
    except errors.HttpError, error:
        print('An error occurred: %s' % error)
        return None


def link_string(link, name='link'):
    """
    Generates a string to represent a link because I'm too lazy to
    type it all out every time.
    """
    return '<a href="' + link + '">' + name + '</a>\n'


def generate_HTML(models):
    """
    Generates the HTML of the table.
    """
    with open('index.html', 'w+') as f:
        f.write('<!DOCTYPE html>\n')
        f.write('<html>\n')
        f.write('<head></head>\n')
        f.write('<body>\n')

        f.write('<table border="1" cellpadding="5" cellspacing="5">\n')
        # write headers (including drop.cfg) and the sound-gen.cfg fields
        f.write('<tr>\n')
        for h in HEADERS + SOUND_GEN_CFG_FIELDS:
            f.write('<th>' + h + '</th>\n')
        f.write('</tr>\n')

        # write models in
        for m in models:
            f.write('<tr>\n')

            f.write('<td>' + m['name'] + '</td>\n')

            # if no image, write "no image", else display image
            if m['image'] is None:
                f.write('<td>no image</td>\n')
            else:
                f.write('<td><img src="' + m['image'] + '">\n')

            # now write in headers and drop.cfg
            for h in HEADERS[2:]:
                try:
                    f.write('<td>' + link_string(m[h], h) + '</td>\n')
                except TypeError: # Directory didn't start with p_
                    f.write('<td>File not found</td>\n')
                except KeyError: # Model doesn't have the file in it
                    f.write('<td>' + m['name'] + ' does not contain file '
                        + h + '</td>\n')

            for h in SOUND_GEN_CFG_FIELDS:
                try:
                    f.write('<td>' + m[h] + '</td>\n')
                except:
                    f.write('<td></td>\n')
                    print(m['name'])

            f.write('</tr>\n')

        f.write('</table>\n')
        f.write('</body>\n')
        f.write('</html>')


if __name__ == '__main__':
    service = create_service()
    print('Finding folder modal_models...')
    modal_models = get_folder_id(service, 'modal_models')
    print('Finding all models...')
    models = get_all_models(service, modal_models)
    print('Number of models: %i' % len(models))

    all_model_dicts = []
    #testing
    for i in range(len(models)):
        if models[i]['title'][:2] == 'p_':
            print('Generating dict for model %s, %i/%i' % (models[i]['title'],
                i+1, len(models)))
            all_model_dicts.append(generate_table_dict(service, models[i]))
    # resize images to be 100x100
    print('Resizing images')
    call('./resize.sh', shell=True)

    #print_file_content(service, all_model_dicts[0]['proj.obj'])
    print('Generating HTML')
    generate_HTML(all_model_dicts)
    print('Done')

