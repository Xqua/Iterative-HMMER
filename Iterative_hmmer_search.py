#!/usr/bin/env python

# This scripts takes an alignment to start with then take a sliding window over it and run a HMMER search of the sequence against the selected database.
# It then extracts the number of Eurkaryote hits versus bacterial hits.
# It should then plot the relative frequencies on top of the sequence logo or similar

import pandas as pd
import numpy as np
from optparse import OptionParser
from Bio import AlignIO, SeqIO
import os
import sys
import matplotlib.pyplot as plt
import urllib2
import matplotlib.patches as mpatches
import progressbar


class Iterative_search:

    def __init__(self, alignpath, window, name, db, osk=False, preprocess=True):
        # self.idmap_path = 'idmapping_reduced.tab'
        # self.db = './Data/uniprot_trembl.fasta'
        self.osk = osk
        if preprocess is None:
            self.preprocess = True
        else:
            self.preprocess = False
        # if db is not None:
        self.db = db

        idmap_path = "./Data/uniprot_species.tsv"
        idmap2_path = "./Data/uniprot_ID_taxa.tsv"
        # self.UniprotID_collums = ["UniProtKB-ID", "NCBI-taxon"]
        # self.UniprotID_collums = ["UniProtKB-AC", "UniProtKB-ID", "GeneID (EntrezGene)", "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90", "UniRef50", "UniParc", "PIR", "NCBI-taxon",
        #                           "MIM", "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional PubMed"]
        # self.NCBI_taxonomy_path = 'NCBI_taxID.tsv'
        # print "Loading NCBI taxonomy"
        # NCBI_taxonomy = {}
        # f = open(self.NCBI_taxonomy_path)
        # for line in f.readlines():
        #     s = line.strip().split('\t')
        #     NCBI_taxonomy[s[1]] = s[0]

        print "Loading ID mapping file"
        self.ID2Taxa = {}
        f = open(idmap_path)
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            # print s
            self.ID2Taxa[s[0].strip()] = s[1].strip()
        f.close()
        self.ID2Arthro = {}
        f = open(idmap2_path)
        for line in f.readlines():
            s = line.strip().split('\t')
            ID = s[0]
            taxas = s[2].strip().split(',')
            if len(taxas) >= 4:
                if taxas[3] == "Arthropoda":
                    self.ID2Arthro[ID] = True
                else:
                    self.ID2Arthro[ID] = False
        # idmap = pd.read_csv(idmap_path, header=None, names=self.UniprotID_collums, sep='\t')
        # NCBI_taxonomy = pd.read_csv(NCBI_taxonomy_path, header=None, names=['Kingdom', 'TaxaID', 'NCBI_ID'], sep='\t')
        print "Loading Alignment"
        handle = open(alignpath, 'rU')
        self.alignment = AlignIO.read(handle, "fasta")
        self.AlLen = self.alignment.get_alignment_length()
        print("Alignment length %i" % self.alignment.get_alignment_length())
        self.window = window
        self.name = name
        self.uniprot_dl = []

    def Parse_HMMER_output(self, path):
        f = open(path)
        lines = f.readlines()
        f.close()
        collumns = ['target name', 'Prot_ID', 'Specie_ID', 'accession', 'query name', 'accession', 'Pre E-value', 'Pre score', 'Pre bias', 'E-value', 'score', 'bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description of target']
        res = []
        for line in lines:
            if line[0] != '#':
                s = [i for i in line.split(' ') if i]
                s = [s[0]] + [s[0].split('|')[2]] + [s[0].split('|')[2].split('_')[1]] + s[1:18] + [' '.join(s[18:])]
                res.append(s)
        df = pd.DataFrame(res, columns=collumns)
        return df

    def Run_Search(self, alignment, p):
        AlignIO.write(alignment, open(self.name + '_' + str(p) + ".fasta", 'w'), "fasta")
        print "Building HMM"
        os.system('hmmbuild -o atej -n Iter_run_%s %s %s' % (p, self.name + '_' + str(p) + ".hmm", self.name + '_' + str(p) + ".fasta"))
        print "Searching against database ..."
        os.system('hmmsearch --cpu 7 -o atej --tblout ' + self.name + '_' + str(p) + '.out ' + self.name + '_' + str(p) + '.hmm ' + self.db)
        # result = self.Parse_HMMER_output(self.name + '_' + str(p) + ".out")
        # return result

    def Test_Oskar(self, accID):
        exist = False
        if accID not in self.uniprot_dl:
            if not os.path.isfile('%s.fasta' % accID):
                url = "http://www.uniprot.org/uniprot/%s.fasta" % (accID)
                handle = urllib2.urlopen(url)
                fasta = handle.read()
                f = open('%s.fasta' % accID, 'w')
                f.write(fasta)
                f.close()
                handle = SeqIO.parse('%s.fasta' % accID, 'fasta')
                try:
                    handle.next()
                    exist = True
                    self.uniprot_dl.append(accID)
                except:
                    print accID
            else:
                exist = True
        else:
            exist = True
        if exist:
            os.system('hmmsearch --cpu 7 -o atej --tblout ' + 'test_oskar_lotus.out ' + '/home/lblondel/Documents/Harvard/ExtavourLab/Project_Oskar_HGT/Sequences/HMM/LOTUS-refined.hmm ' + '%s.fasta' % accID)
            os.system('hmmsearch --cpu 7 -o atej --tblout ' + 'test_oskar_sgnh.out ' + '/home/lblondel/Documents/Harvard/ExtavourLab/Project_Oskar_HGT/Sequences/HMM/SGNH-refined.hmm ' + '%s.fasta' % accID)
            lotus = self.Parse_HMMER_output("test_oskar_lotus.out")
            sgnh = self.Parse_HMMER_output("test_oskar_sgnh.out")
            if len(sgnh) > 0 and len(lotus) == 0:
                f = open('toCheck', 'a')
                f.write(accID + '\n')
            if len(lotus) > 0:
                if len(sgnh) > 0:
                    return True
        return False

    def Extract_kingdoms(self, df):
        result = []
        for i in range(len(df)):
            ID = df['Specie_ID'][i].split(' ')[0]
            kg = self.ID2Taxa[ID]
            # print ID
            if kg == 'E':
                if ID in self.ID2Arthro:
                    if self.ID2Arthro[ID]:
                        kg = 'A'
                        if self.osk:
                            if 'osk' in df['description of target'][i].lower():
                                kg = 'O'
                            elif self.Test_Oskar(df['Prot_ID'][i]):
                                kg = 'O'
            result.append(kg)
        return result

    def Get_kingdom(self, ID, description, protID):
        kg = self.ID2Taxa[ID]
        # print ID
        if kg == 'E':
            if ID in self.ID2Arthro:
                if self.ID2Arthro[ID]:
                    kg = 'A'
                    if self.osk:
                        if 'osk' in description.lower():
                            kg = 'O'
                        elif self.Test_Oskar(protID):
                            kg = 'O'
        return kg

    def Preprocess(self):
        trembl = SeqIO.parse(self.db, format='fasta')
        match = {}
        print "Running searches against trembl"
        tot = len(np.arange(0, self.AlLen, 250)) + len(np.arange(125, self.AlLen, 250))
        j = 1
        for i in np.arange(0, self.AlLen, 250):
            print "Search %s/%s" % (j, tot)
            tmp_align = self.alignment[:, i:max(self.AlLen, i + 250)]
            self.Run_Search(tmp_align, "pre")
            df = self.Parse_HMMER_output(self.name + '_pre.out')
            for prot in df['Prot_ID'].values:
                match[prot] = True
            j += 1
        for i in np.arange(125, self.AlLen, 250):
            print "Search %s/%s" % (j, tot)
            tmp_align = self.alignment[:, i:max(self.AlLen, i + 250)]
            self.Run_Search(tmp_align, "pre")
            df = self.Parse_HMMER_output(self.name + '_pre.out')
            for prot in df['Prot_ID'].values:
                match[prot] = True
            j += 1
        print "Extracting subDatabase"
        db = []
        match = set(match.keys())
        for seq in trembl:
            if seq.name.split('|')[-1] in match:
                db.append(seq)
        SeqIO.write(db, self.name + '_database.fasta', format='fasta')
        self.db = self.name + '_database.fasta'

    def Run(self):
        if self.preprocess:
            print "Preprocessing ..."
            self.Preprocess()
        for i in range(self.AlLen - self.window):
            # for i in range(282, 290):
            print "Doing %s out of %s iterations ..." % (i, self.AlLen - self.window)
            # take a slice
            print "Building Alignment"
            tmp_align = self.alignment[:, i:i + self.window]
            # run HMMER and parse results
            self.Run_Search(tmp_align, i)

    def Analyze(self):

        print "Loading results ..."
        bar = progressbar.ProgressBar()
        results = []
        for i in bar(range(self.AlLen - self.window)):
            # print "Loading %s/%s" % (i, self.AlLen - self.window)
            df = self.Parse_HMMER_output('%s_%s.out' % (self.name, i))
            res = self.Extract_kingdoms(df)
            results.append(res)

        print "Compute ratio"
        # save ratio in array
        df = []
        if os.path.isfile(self.name + '_database.fasta'):
            res = {'B': 0, 'E': 0, 'A': 0, 'O': 0, 'Other': 0}
            tot = 0
            seqs = SeqIO.parse(self.name + '_database.fasta', format='fasta')
            for i in seqs:
                specID = i.name.split('|')[2].split('_')[1]
                kg = self.Get_kingdom(specID, i.description, i.name.split('|')[1])
                if kg in res.keys():
                    res[kg] += 1
                else:
                    res['Other'] += 1
                tot += 1
            if self.osk:
                df.append(['db', res['B'], res['E'], res['A'], res['O'], res['Other'], 0])
            else:
                df.append(['db', res['B'], res['E'], res['A'], res['Other'], 0])
        counts = []
        j = 0
        for i in results:
            tmp = []
            tmp.append(i.count('B'))
            tmp.append(i.count('E'))
            tmp.append(i.count('A'))
            if self.osk:
                tmp.append(i.count('O'))
                tmp.append(len(i) - (i.count('B') + i.count('A') + i.count('E') + i.count('O')))
            else:
                tmp.append(len(i) - (i.count('B') + i.count('A') + i.count('E')))
            counts.append(tmp)
            p = 1 - self.alignment[:, j].count('-') / float(len(self.alignment[:, j]))
            df.append([j] + tmp + [p])
            j += 1
        counts = np.array(counts).T
        if self.osk:
            df = pd.DataFrame(df, columns=['Position', 'Bacteria', 'Eukaryotes', 'Arthropods', 'Oskar', 'Other', 'Occupancy'])
        else:
            df = pd.DataFrame(df, columns=['Position', 'Bacteria', 'Eukaryotes', 'Arthropods', 'Other', 'Occupancy'])
        df.to_csv(self.name + '.csv', index=False)
        # return results, counts
        # plot ratio
        # sns.
        # http://matplotlib.org/examples/pylab_examples/stackplot_demo.html
        print "Plotting results"
        fig, ax = plt.subplots()
        if self.osk:
            t = ax.stackplot(range(len(counts[0])), counts[0], counts[1], counts[2], counts[3], counts[4])
            leg = ['Bacteria', 'Eukaryotes', 'Arthropoda', 'Oskar', 'Other']
        else:
            t = ax.stackplot(range(len(counts[0])), counts[0], counts[1], counts[2], counts[3])
            leg = ['Bacteria', 'Eukaryotes', 'Arthropoda', 'Other']
        handles = []
        for i in range(len(t)):
            handles.append(mpatches.Patch(color=t[i].get_facecolor()[0], label=leg[i]))
        ax.legend(handles=handles)
        plt.title(self.name)
        fig.savefig("%s.pdf" % self.name)
        fig.savefig("%s.png" % self.name)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-a", "--alignment", dest="alignpath", default=None,
                      help="[Required] Location of the FASTA alignement file.")
    parser.add_option("-w", "--window", dest="window", default="None",
                      help="[Required] Size of the sliding window (aka minimum nb of columns to be used for HMM generation)")
    parser.add_option("-n", "--name", dest="name", default=None,
                      help="[Required] Name of the analysis)")
    parser.add_option("-l", "--analyze", dest="anal", action="store_true",
                      help="[Optional] Analyze only)")
    parser.add_option("-d", "--db", dest="db", default=None,
                      help="[Required] Path to DB (default to Trembl)")
    parser.add_option("-o", "--oskar", dest="osk", action="store_true",
                      help="Is it the Oskar Gene ?")
    parser.add_option("-p", "--preprocess", dest="preprocess", action="store_true",
                      help="Do NOT preprocess the search. (slower and more prone to false positives)")

    # Parse options into variables
    (options, args) = parser.parse_args()

    alignpath = options.alignpath
    if not alignpath:
        print("You must provide an alignement.")
        sys.exit(1)
    name = options.name
    if not name:
        print("You must provide a name.")
        sys.exit(1)
    anal = options.anal
    db = options.db
    if not db:
        print("You must provide a database.")
        sys.exit(1)
    osk = options.osk
    window = options.window
    preprocess = options.preprocess
    if window:
        try:
            window = int(options.window)
        except:
            print "window must be an integer"
    else:
        print "window must be specified"
        sys.exit(1)

    Iter = Iterative_search(alignpath, window, name, db, osk, preprocess)
    print "Starting !"
    if not anal:
        Iter.Run()
    print "Analyzing !"
    Iter.Analyze()

# slices for alignements : align[line_start:line_end, col_start:col_end]


# alignment = []
# for record in SeqIO.parse(handle, "fasta"):
#     print record.id
# handle.close()
