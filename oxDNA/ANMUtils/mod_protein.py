import models as m
import sys
import Bio.PDB
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
import subprocess as sp
import os


def write_PIR(seqs, outfile, labels=[]):
    o = open(outfile, 'w')
    if type(seqs) is str:
        if len(labels) == 1:
            print('>P1;' + labels[0], file=o)
            print('sequence:' + labels[0] + ':::::::0.00: 0.00', file=o)
            print(seqs, file=o)
        else:
            print('>P1;seq1', file=o)
            print('sequence:seq1:::::::0.00: 0.00', file=o)
            print(seqs, file=o)
    else:
        if len(labels) == len(seqs):
            for i in range(len(seqs)):
                print('>P1;' + labels[i], file=o)
                print('sequence:'+labels[i]+':::::::0.00: 0.00', file=o)
                print(seqs[i], file=o)
        else:
            for i in range(len(seqs)):
                comm = 'seq' + str(i+1)
                print('>P1;seq' + str(i+1), file=o)
                print('sequence:'+comm+':::::::0.00: 0.00', file=o)
                print(seqs[i], file=o)
    o.close()


def clustal_align(infasta, outclustal):
    cline = ClustalOmegaCommandline(infile=infasta, outfile=outclustal, outfmt='clustal', verbose=True, auto=False)
    sp.check_call(str(cline), shell=True)


# def get_pdb_seq(pdbfile):
#     if "/" in pdbfile:
#         pdbid = pdbfile.rsplit('/', 1)[1].split('.')[0]
#     else:
#         pdbid = pdbfile.split('.')[0]
#     s1 = m.get_pdb_info(pdbid, pdbfile, returntype=3)
#     trash, seq = zip(*s1)
#     seq = ''.join(list(seq))
#     return s1


def get_seq_by_chain(pdbfile):
    s1 = m.get_pdb_info(pdbfile, returntype='S')
    #print(s1)
    chainids = list(set([x for x, y in s1]))
    chains = {}
    for i in chainids:
        res = [y for x, y in s1 if x == i]
        chains[i] = ''.join(res)
    return chains

def write_chain_seq(chaindict, pdbid='', fmt='fasta'):
    keys = chaindict.keys()
    for i in keys:
        if fmt=='pir':
            write_PIR(chaindict[i], pdbid+'chain'+str(i)+'.pir', labels=[pdbid + 'chain' + str(i)])
        elif fmt=='fasta':
            filename = pdbid + 'chain' + str(i) + '.fasta'
            o = open(filename, 'w')
            print('>' + str(pdbid) + 'chain' + str(i), file=o)
            print(chaindict[i], file=o)
            o.close()

def fasta_read(fastafile):
    o = open(fastafile)
    titles = []
    data = o.read()
    o.close()
    datas = data.split('\n')
    for lid, line in enumerate(datas):
        if line.startswith('>'):
            titles.append(str(line.rstrip().split('>')[1]))
            datas[lid] = '>'
    strdatas = ''.join(datas)
    seqdatas = strdatas.split('>')
    seqs = [x for x in seqdatas if x != '']
    return seqs, titles

def write_fasta(seqs, titles, out):
    o = open(out, 'w')
    for xid, title in enumerate(titles):
        print('>' + str(title), file=o)
        print(seqs[xid], file=o)
    o.close()

def combine_fastas(outfile, *fastafiles):
    seqs, titles = [], []
    for ff in fastafiles:
        seq, title = fasta_read(ff)
        seqs += seq
        titles += title
    print(seqs)
    write_fasta(seqs, titles, outfile)



def conv_clustal_to_PIR(clustalfile, pdbid, pdbcid, pdblb, pdbup, outfile):
    try:
        from Bio import SeqIO
    except ImportError:
        print('Check that Biopython is installed')
    records = SeqIO.parse(clustalfile, "clustal")
    for record in records:
        #PDB sequence
        a = str(record.seq)
        if '-' in a:
            pdb_d = 'structure:'+pdbid+':'+str(pdblb)+':'+pdbcid+':'+str(pdbup)+':'+pdbcid+'::::'
            pdbrecord = SeqRecord(Seq(a, generic_protein), id=pdbid, name='', description=pdb_d)
        else:
            seq_d = 'seqeunce:'+record.id+'::::::::'
            seqrecord = SeqRecord(Seq(a, generic_protein), id=record.id, name='', description=seq_d)
    nrecords=[pdbrecord, seqrecord]
    SeqIO.write(nrecords, outfile, "pir")

#Prior to this step must download Unitprot sequences and rename to chainX.fasta
#Titles must be manually set in these files as well
def auto_prep_files(pdbfile):
    if "/" in pdbfile:
        pdbid = pdbfile.rsplit('/', 1)[1].split('.')[0]
    else:
        pdbid = pdbfile.split('.')[0]
    chaindict = get_seq_by_chain(pdbfile)
    chains = list(chaindict.keys())
    taskmaster = []
    for chain in chains:
        print('Prepping Chain', chain)
        pdbtitle = 'chain' + chain.upper()
        try:
            seqs, titles = fasta_read(pdbtitle+'.fasta')
        except:
            print('Error Reading', pdbtitle+'.fasta')
            print('Make sure each Fasta Uniprot file conforms to the following naming convention chainX.fasta')
        pdbseq = chaindict[chain]
        cseqs, ctitles = [seqs[0], pdbseq], [titles[0], pdbtitle]
        pafile = 'prealign'+chain.upper()+'.fasta'
        cfile = 'align'+chain.upper()+'.clustal'
        write_fasta(cseqs, ctitles, pafile)
        clustal_align(pafile, cfile)
        tl = (chain, pdbid, titles[0])
        taskmaster.append(tl)
    return taskmaster


def auto_models(tasklist, pdbfile, model_number=5):
    if "/" in pdbfile:
        pdbid = pdbfile.rsplit('/', 1)[1].split('.')[0]
    else:
        pdbid = pdbfile.split('.')[0]
    boundsdict = m.get_pdb_info(pdbfile, returntype='p')
    for task in tasklist:
        cid, kid, sid = task
        cfile = 'align' + cid.upper() + '.clustal'
        pirfile = 'mod' + cid.upper() + '.pir'

        lb, ub = boundsdict[cid]
        conv_clustal_to_PIR(cfile, pdbid, cid, lb, ub, pirfile)

        print('Generating Models for Chain', cid)
        ndir = 'mod' + cid
        create_models(pirfile, kid, sid, pdbfile, model_number=model_number, dir=ndir)


def create_models(alnfile, knownid, sequenceid, pdbfile, model_number=5, dir=''):
    if dir:
        cdir = os.getcwd()
        ndir = cdir+'/'+dir
        os.mkdir(ndir)
        os.chdir(ndir)
        sp.check_call("cp ../" + pdbfile + ' ' + pdbfile, shell=True)
    try:
        from modeller import environ
        from modeller.automodel import automodel
    except ImportError:
        print('Make Sure Python Modeller is installed. Double check License Key is in Modeller config.py file')
        sys.exit()
    env = environ()
    if dir:
        a = automodel(env, alnfile='../'+alnfile, knowns=knownid, sequence=sequenceid)
    else:
        a = automodel(env, alnfile=alnfile, knowns=knownid, sequence=sequenceid)
    a.starting_model = 1
    a.ending_model = model_number
    a.make()
    if dir:
        os.chdir(cdir)


def adj_seq(outfile, modpdbfiles):
    #Get Base Structure, we're adding chains to this one
    if "/" in modpdbfiles[0]:
        pdbid = modpdbfiles[0].rsplit('/', 1)[1].split('.')[0]
    else:
        pdbid = modpdbfiles[0].split('.')[0]
    basestruct = Bio.PDB.PDBParser().get_structure(pdbid, modpdbfiles[0])
    basemodel = basestruct[0]
    #Here's where we will store each chain in a dictionary
    chain_cont = {}

    # Get all of the chains (except in OG PDB file) and store
    for i in range(1, len(modpdbfiles)):
        ref = modpdbfiles[i] #its shorter to type

        #open PDB file
        if "/" in ref:
            pdbid = ref.rsplit('/', 1)[1].split('.')[0]
        else:
            pdbid = ref.split('.')[0]
        # Container for Entirety of PDB file
        struct = Bio.PDB.PDBParser().get_structure(pdbid, ref)
        # Get all chains even between different models
        model = Bio.PDB.Selection.unfold_entities(struct, 'M')
        print(model)
        for chain in model:
            print(str(chain.get_id())) #.upper()
            chain_cont[str(chain.get_id().upper())] = chain

    print(chain_cont.keys())
    print(chain_cont.values())
    ck = list(chain_cont.keys())
    ck.sort()

    # Add Each chain into basemodel
    for c in ck:
        basemodel.add(chain_cont[c])

    io = Bio.PDB.PDBIO()
    io.set_structure(basemodel)
    io.save(outfile)





