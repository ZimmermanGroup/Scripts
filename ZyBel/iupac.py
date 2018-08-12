import requests
from urllib2 import urlopen
import re

def iupac_name(smile,substrates,g,sub_id,group_id,pos):
    try:
        name=urlopen('http://cactus.nci.nih.gov/chemical/structure/' + smile + '/iupac_name').read().decode('utf-8') 
        if (len(name)>100):
            try:
                #print(substrates)
                name1=urlopen('http://cactus.nci.nih.gov/chemical/structure/' + substrates + '/iupac_name').read().decode('utf-8') 
            except:
                name1="substrates"+str(sub_id)
            try:
                #print(g)
                name2="_pos"+str(pos)+"_"++urlopen('http://cactus.nci.nih.gov/chemical/structure/' + g + '/iupac_name').read().decode('utf-8') 
                if len(name2)>100:
                    name2="_pos"+str(pos)+"_group"+str(group_id)
                    group_id+=1
                if name2.find("$l^{1}-")==-1: 
                    name2=name2
                else:
                    name2="_pos"+str(pos)+"_dimethylamine"
            except:
                name2="_pos"+str(pos)+"_group"+str(group_id)
            name=name1+name2
    except:
        try:
            #print(substrates)
            name1=urlopen('http://cactus.nci.nih.gov/chemical/structure/' + substrates + '/iupac_name').read().decode('utf-8') 
        except:
            name1="substrates"+str(sub_id)
        try:
            #print(g)
            name2="_pos"+str(pos)+"_"+urlopen('http://cactus.nci.nih.gov/chemical/structure/' + g + '/iupac_name').read().decode('utf-8') 
            if len(name2)>200:
                name2="_pos"+str(pos)+"_group"+str(group_id)
            if name2.find("$l^{1}-")==-1: 
                name2=name2
            else:
                name2="_pos"+str(pos)+"_"+"dimethylamine"
        except:
            name2="_pos"+str(pos)+"_"+"group"+str(group_id)
        name=name1+name2
    name = re.sub(' ', '', name)
    name = re.sub('\)', '', name)
    name = re.sub('\(', '', name)
    return name
