#!/usr/bin/env python3

'''
Deze programma zal een csi.fa bestand inlezen en doormiddel van restrictie enzymen op patronen en specifieke locaties gaan knippen.
de fragmenten die hiervan overblijven worden meegenomen in een dictionary: " csi_seq ".
van de fragmenten uit deze dictionary zal het moleculaire massa worden berekend en dat zal worden opgeslagen in een andere dictionary: " gui_dict ".
Er is een belangrijke functie gemaakt die geimporteed wordt vanuit de restriction_GUI_v2.py, deze functie wordt performRestriction genoemd die al het werk dan verricht.
de gui_dict wordt verder gegeven aan de restriction_GUI_v2.py script om het te visualiseren.

Doormiddel van deze programma kan achterhaald worden wie de dader is bij de CSI Crime Scene.
'''


import sys
import re
import os 

_Writer = 'Armin Dedic'
_Group = 'BFV1'
_datum = "23 - 03 - 2016"


def Reading(DNAfile):
    '''
    Deze functie zal de csi.fa bestand inlezen en de personen als key
    en hun bijbehorende sequentie als value opslaan in de seqdict.
    '''

    DNA_files = open(DNAfile)
    sequence = []
    seqdict = {}
    header = None  
    
    for line in DNA_files:
        line = line.strip()
        if line:
            if line.startswith(">"):
                header = line
                sequence = []
            else:
                sequence.append(line)
                seqdict[header] = sequence
    return seqdict
            
enzyme_dict = {
'EcorI'    :{'PATROON':'GAATTC',        'SITE':1},  
'BamHI'    :{'PATROON':'GGATCC',        'SITE':1},
'EcoRII'   :{'PATROON':'CC[AT]GG',      'SITE':0},
'HindIII'  :{'PATROON':'AAGCTT',        'SITE':1},
'TaqI'     :{'PATROON':'TCGA',          'SITE':1},
'NotI'     :{'PATROON':'GCGGCCGC',      'SITE':2},
'HinfI'    :{'PATROON':'GA[ATCG]TCA',   'SITE':1},
'Sau3A'    :{'PATROON':'GATC',          'SITE':0},
'HaeIII'   :{'PATROON':'GGCC',          'SITE':2},
'EcoRV'    :{'PATROON':'GATATC',        'SITE':3},
'PstI'     :{'PATROON':'CTGCAG',        'SITE':4},
'XbaI'     :{'PATROON':'TCTAGA',        'SITE':1},
'MboI'     :{'PATROON':'GATC',          'SITE':0},
'mst2'     :{'PATROON':'CCT[ATCG]AGG',   'SITE':2},
}


#~ def changing():
    #~ changed_dict = {}
    #~ for enzyme in enzyme_selected:
        #~ print(enzyme)
        #~ enzyme_pat = enzyme_selected.get(enzyme)['PATROON']
        #~ for change in enzyme_pat:
            

def extracting_fragments(csi, enzymelist):
    '''
    hier worden de sequenties van de personen geknipt door enzymen
    de overgebleven fragmenten worden in de vorm van een lijst 
    toegevoegd aan de dictionary als value
    '''
    csi_seq = csi
    matches = []
    for enzyme in enzymelist:
        for sequentie in csi:
            fragmenten = list()
            for seq_frag in csi_seq[sequentie]:
                enzyme_seq = enzyme_dict.get(enzyme)['PATROON']
                Pat = ""
                for char in enzyme_seq:
                    if char == "N":
                        char = "[CGTA]"
                    elif char == "W":
                        char = "[AT]"
                    Pat += char
                matches = re.finditer(Pat,seq_frag)
                cut_site = enzyme_dict.get(enzyme)['SITE']
                cut_start = 0
                count = 0
                for match in matches:
                    cut_end = int(match.span()[0]) + int(cut_site)
                    seq = seq_frag[cut_start:cut_end]
                    cut_start = cut_end
                    fragmenten.append(seq)
                    count += 1
                fragmenten.append(seq_frag[cut_start:])
                csi_seq[sequentie] = fragmenten
            if count == 0:
                fragmenten = csi_seq[sequentie]
    return csi_seq

def calculation(fragments):
    '''
    deze functie zal van de fragmenten de moleculaire massa berekenen
    en opslaan in de gui_dict
    '''
    
    amino_weights = {"A":135,"C":111,"T":126,"G":151}
    gui_dict = {}
    
    for person in fragments:
        fragmentlijst = []
        for frags in fragments[person]:
            weight = 0
            for item in frags:
                weight += amino_weights[item]
            fragmentlijst.append(int(weight))
        gui_dict[person] = fragmentlijst      
    return gui_dict

def performRestriction(DNAfile, enzymes_selected):
    '''
    deze functie zal al het werk verrichten.
    '''
    
    argv = [DNAfile]
    for enzyme in enzymes_selected:
        argv.append(enzyme)
    result = main(argv)
    return result
    
def main(argv=None):
    if argv == None: 
        argv = sys.argv[1:]
    
    if len(argv) == 0:
        print(__doc__)
        return 1
    
    DNAfile = argv[0]
    enzymelist = argv[1:]
    csi = Reading(DNAfile)
    fragments = extracting_fragments(csi, enzymelist)
    result = calculation(fragments)
    return result

if __name__ == '__main__':
    sys.exit(main())    

    
