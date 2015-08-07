######################################
## Read in Scattering Cross Section ##
######################################
def scatter_read_in(file_path,reverse=0):
    ''' read in scattering data base in format given by http://www.nndc.bln.gov/'''
    ''' Input: File path [type=string], Reverse[set to 0 or 1]; Output: Energy [eV], Cosine, Probability '''
    ''' reverse is needed because the database calles a backscatter +1 and a forward scatter -1 and I use a backscatter = -1 and a forward scatter = +1'''
    f = open(file_path,'r')
    raw = f.read().splitlines()
    cosine , probs , energy= [] , [] , []
    index = 1
    while index < len(raw)-2:
        if raw[index][0] == 'I':
            temp = raw[index].split(':')
            temper = temp[1].split('eV')
            energy.append(float(temper[0]))
            index += 3
            new_line = raw[index].split('   ')
            cos , prob = [] , []
            while raw[index] != raw[0] and index < len(raw)-1:
                new_line = raw[index].split('   ')
                cos.append(float(new_line[0]))
                prob.append(float(new_line[1]))
                index += 1
            if reverse == 1:   
                cos.reverse()
            cosine.append(cos)
            probs.append(prob)
        else:
            index += 1  
    return energy,cosine,probs
###################################
## Read in Elastic Cross Section ##
###################################
def elastic_read_in(file_path):
    ''' Read in database data in format given by http://www.nndc.bnl.gov/ '''
    ''' Input: file_path [type=string]; Output: Energy [eV], Cross Section [b] ''' 
    f = open(file_path,'r')
    f = f.read().splitlines()
    f = f[1:]
    Data = []
    for i in f:
        Data.append(i.split(','))
    Data   = [[float(i)for i in line] for line in Data]
    Energy = [col[0] for col in Data]
    XC     = [col[1] for col in Data]
    return Energy,XC
