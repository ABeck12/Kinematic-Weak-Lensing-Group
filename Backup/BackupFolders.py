'''
Setup:
Download the google python packages with
pip install --upgrade google-api-python-client google-auth-httplib2 google-auth-oauthlib

Go to https://console.cloud.google.com/welcome/new and create a new project
On the left under the navigation menu go to APIs and Services -> OAuth Consent Screen
    User type should be external
    Fill out App name (Something like python drive backup), User support email, and Developer Contact Information
Go to APIs and Services -> Credentials
    Press the Create Credentials button -> OAuth client ID
    Select Desktop app
    Once the popup appears hit the download json button, rename the file to credentials.json and put it in the same directory as this file
    NOTE: credentials.json SHOULD BE KEPT SECRET! DO NOT SHARE IT WITH OTHERS!
Go to APIs and Services -> Enable APIs & Services
    Press the Enable APIs & Services button
    Search for "Google Drive API"
    Click it and enable it

How to back up files:
There are 3 backup functions
    SingleFolderBackup
    MultiFolderBackup
    MultiFolderBackupFromTxt
All of these functions do a recursive search for all files and subfolders within a directory and backs it up to a drive folder

SingleFolderBackup(rootFilepath, gDriveFolderName)
    rootFilepath - This is the directory that you want to use.
    gDriveFolderName - This is the name of the root folder that it will save things to. All subfolders names will be copied over to the google drive subfolders

MultiFolderBackup(rootFilepathsList, gDriveFolderName)
    rootFilepathsList - This is a list of the the directories that you want to use.
    gDriveFolderName - This is the name of the root folder that it will save things to. All subfolders names will be copied over to the google drive subfolders

MultiFolderBackupFromTxt(txtFileName, gDriveFolderName)
    txtFileName - Filepath to a file with a list of directories that you want to use.
    gDriveFolderName - This is the name of the root folder that it will save things to. All subfolders names will be copied over to the google drive subfolders
'''


import os

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaFileUpload
import time
from tqdm import tqdm


def TimeFunction(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        out = func(*args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {end-start:.3f} seconds to run')
        return out
    return wrapper

SCOPES = ["https://www.googleapis.com/auth/drive"]

def SetupGDriveCreds():
    creds = None

    if os.path.exists("token.json"):
        creds = Credentials.from_authorized_user_file("token.json", SCOPES)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                "credentials.json", SCOPES
            )
            creds = flow.run_local_server(port=0)

    with open("token.json", "w") as token:
        token.write(creds.to_json())
    return creds

def CreateFolder(service, name, parentID=None):
    file_metadata = {
                "name": f'{name}',
                "mimeType": "application/vnd.google-apps.folder"
            }
    if parentID is not None:
        file_metadata["parents"] = [parentID]
    print(f'Making new folder "{name}" with parent {parentID}')
    file = service.files().create(body=file_metadata).execute()
    folder_id = file.get('id')
    return folder_id

def ResumeableFileUpload(service, folderpath, filename, parentID):
    try:
        file_metadata = {"name": filename, "parents": [parentID],}
        media = MediaFileUpload(f'{folderpath}/{filename}',
                            chunksize=1024 * 1024,
                            mimetype='text/plain',
                            resumable=True)
        request = service.files().create(body=file_metadata, media_body=media)
        response = None

        print(f'Uploading {folderpath}/{filename}')
        with tqdm(total=100) as pbar:
            while response is None:
                status, response = request.next_chunk()
                if status:
                    currentPercentage = int(status.progress() * 100)
                    pbar.n = currentPercentage
                    pbar.refresh()
            pbar.n = 100
            pbar.refresh
    except:
        print(f'Could not back up {folderpath}/{filename}')

def RecursiveBackupFolder(service, folder_id, folderpath: str):
    for path in os.listdir(folderpath):
        if os.path.isdir(f'{folderpath}/{path}'):
            subfolderID = CreateFolder(service, path, folder_id)
            RecursiveBackupFolder(service, subfolderID, f'{folderpath}/{path}')
        if os.path.isfile(f'{folderpath}/{path}'):
            ResumeableFileUpload(service, folderpath, path, folder_id)
            # try:
            #     print(f'Starting backup of file: {folderpath}/{path}')
            #     file_metadata = {"name": path, "parents": [folder_id],}
            #     media = MediaFileUpload(f'{folderpath}/{path}')
            #     upload_file = service.files().create(body=file_metadata, media_body=media, fields = "id").execute()
            #     print(f'Backed up file: {folderpath}/{path}')
            # except:
            #     print(f'Could not backup file: {folderpath}/{path}')
    return

@TimeFunction
def SingleFolderBackup(rootFilepath, gDriveFolderName):
    creds = SetupGDriveCreds()

    try:
        service = build("drive", "v3", credentials=creds)
        response = service.files().list(q=f"name = '{gDriveFolderName}' and mimeType = 'application/vnd.google-apps.folder' and trashed = false",spaces='drive').execute()
        if not response['files']:
            folder_id = CreateFolder(service, gDriveFolderName)
        else:
            overwrite = input(f'WARNING: Google Drive folder "{gDriveFolderName}" already exists! Overwrite the existing folder? (y/n)\n')
            if overwrite == 'n': 
                print(f'Backup will not overwrite folder "{gDriveFolderName}". Exiting...')
                return
            elif overwrite == 'y':
                folder_id = response['files'][0]['id']
                body_value = {'Trashed':True}
                response = service.files().delete(fileId=folder_id).execute()
                folder_id = CreateFolder(service, gDriveFolderName)
            else:
                print('Unknown answer, exiting...')
                return

        RecursiveBackupFolder(service, folder_id, rootFilepath)

    except HttpError as error:
        print(f"An error occurred: {error}")

    print('Done')

@TimeFunction
def MultiFolderBackup(rootFilepathsList, gDriveFolderName):
    creds = SetupGDriveCreds()

    try:
        service = build("drive", "v3", credentials=creds)
        response = service.files().list(q=f"name = '{gDriveFolderName}' and mimeType = 'application/vnd.google-apps.folder' and trashed = false",spaces='drive').execute()
        if not response['files']:
            folder_id = CreateFolder(service, gDriveFolderName)
        else:
            overwrite = input(f'WARNING: Google Drive folder "{gDriveFolderName}" already exists! Overwrite the existing folder? (y/n)\n')
            if overwrite == 'n': 
                print(f'Backup will not overwrite folder "{gDriveFolderName}". Exiting...')
                return
            elif overwrite == 'y':
                folder_id = response['files'][0]['id']
                body_value = {'Trashed':True}
                response = service.files().delete(fileId=folder_id).execute()
                folder_id = CreateFolder(service, gDriveFolderName)
            else:
                print('Unknown answer, exiting...')
                return
            
        for rootFilepath in rootFilepathsList:
            RecursiveBackupFolder(service, folder_id, rootFilepath)

    except HttpError as error:
        print(f"An error occurred: {error}")

    print('Done')
    return

@TimeFunction
def MultiFolderBackupFromTxt(txtFilename, gDriveFolderName):
    dirlist = []
    try:
        with open(txtFilename,'r') as file:
            for line in file:
                if '\n' in line:
                    dirlist.append(line.replace("\n",""))
                else:
                    dirlist.append(line)
    except FileNotFoundError:
        print(f'{txtFilename} was not found!')
        return
    
    MultiFolderBackup(dirlist, gDriveFolderName)

def main():
    # SingleFolderBackup(r'C:\Users\Alden\OneDrive\Documents\StonyBrook\Spectroscopy\BackupTesting\largeFiles', 'LargeFilesTest')
    SingleFolderBackup(r'C:\Users\Alden\OneDrive\Documents\StonyBrook\Spectroscopy\Specpro', 'Astro Backup')
    return

if __name__ == "__main__":
    main()