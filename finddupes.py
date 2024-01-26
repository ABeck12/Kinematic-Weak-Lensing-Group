import os
import numpy as np


filenamePerObject = 'DuplicatesByObject.txt'
filenamePerSlit = 'DiplicatesBySlit.txt'


class SlitObject:
    def __init__(self, ra, dec, tag):
        self.ra = float(ra)
        self.dec = float(dec)
        self.tag = tag
        self.tags = [tag]
        self.key = (ra, dec)
    
    def UpdateTag(self, newTag):
        self.tags.append(newTag)
        return self
    
    def SqrDistance(self, other):
        radist = abs(self.ra - other.ra)
        decdist = abs(self.dec - other.dec)
        return (radist**2 + decdist**2)
    
    def WithinThreshold(self, other):
        threshold = 0.000138889 **2
        return self.SqrDistance(other) <= threshold
    
    def __str__(self):
        # TODO: Swap 
        # return f'RA: {self.ra}, DEC: {self.dec}, Found in: {self.tags}'
        string = f'{self.ra}, {self.dec}'
        for tag in self.tags:
            string += f', {tag}'
        return string
    
    def __repr__(self):
        return str(self)
    
    def __lt__(self, other):
        return len(self.tags) < len(other.tags)
    
    

    

objectDict = {}
numFiles = len([file for file in os.listdir(os.getcwd()) if file.endswith('.dat')])

for file in os.listdir(os.getcwd()):
    if file.endswith(".dat"):

        datFile = np.loadtxt(file, dtype='str')

        # TODO: Swap 
        currentObject = SlitObject(datFile[1,1], datFile[2,1], file[13:-4] + ':' + datFile[0,1] + ':' + datFile[9,1])
        # currentObject = SlitObject(datFile[1,1], datFile[2,1], file[13:-4] + ':' + datFile[0,1])
        
        if objectDict.get(currentObject.key) != None: # if we find it exactly in the dict already
            objectDict[currentObject.key] = objectDict[currentObject.key].UpdateTag(currentObject.tag)
        elif objectDict.get(currentObject.key) == None: #check in dict to find close match
            addedByDist = False
            for key, item in objectDict.items():
                if item.WithinThreshold(currentObject):
                    item.UpdateTag(currentObject.tag)
                    addedByDist = True
                    break
            if not addedByDist:
                objectDict[currentObject.key] = currentObject

                





# Check that we didnt miss any obejcts/slits
counter = 0
for key, item in objectDict.items():
    counter += len(item.tags)

    object = item

print(f'matched {counter} of {numFiles} objects')
assert counter == numFiles


objectList = list(objectDict.values())
objectList.sort(reverse=True)


newFile = open(filenamePerObject,'w')
for item in objectList:
    newFile.write(str(item))
    newFile.write('\n')
newFile.close()













exit()




#---------------------------------------------------------------------------------------------

class PerSlitObject:
    def __init__(self, tag):
        self.tag = tag
        self.otherSlits = []

    def Unpack(self):
        # return str(*self.otherSlits,)
        return ", ".join(map(str, self.otherSlits))

    def __lt__(self, other):
        if len(self.tag) != len(other.tag):
            return len(self.tag) < len(other.tag)
        else:
            return self.tag[:-4] < other.tag[:-4]
    
    def __str__(self):
        # other = 
        # return f'{self.tag}, {self.otherSlits}'
        return f'{self.tag}, ' + str(self.otherSlits)[1:-1]
    
    def __repr__(self):
        return str(self)

# test = 'A22612015D-001' 
# print(test[:-4])
# print('A22612015D' > 'A22612015C')

# Make a second list for objects by slits---------
perSlitList = []
for file in os.listdir(os.getcwd()):
    if file.endswith(".dat"):

        datFile = np.loadtxt(file, dtype='str')

        currentSlit = PerSlitObject(file[13:-4] + '-' + datFile[0,1])
        # print(currentSlit.tag)
        for key, item in objectDict.items():
            if currentSlit.tag in item.tags:
                tempTag = item.tags.copy()
                tempTag.remove(currentSlit.tag)
                currentSlit.otherSlits = tempTag
                break
        perSlitList.append(currentSlit)

        # print(str(currentSlit.otherSlits)[1:-1])



# print(perSlitList.sort())
perSlitList.sort()

newFile = open(filenamePerSlit,'w')
for item in perSlitList:
    newFile.write(str(item))
    newFile.write('\n')
newFile.close()