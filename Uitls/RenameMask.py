import os
import shutil

OutputDir = r'C:\Users\Alden\OneDrive\Documents\StonyBrook\Spectroscopy\Specpro\MS0451 dat files\All'
newMask = 'J'

assert (newMask == os.path.basename(os.getcwd()))

files = [file for file in os.listdir(os.getcwd()) if file.endswith('.dat')]


for file in files:
    newFileName = file[0:-17] + newMask + '.dat'
    print(newFileName)
    shutil.copyfile(file,f'{OutputDir}\{newFileName}')
     