from scipy.stats import binom_test
import gzip
import sys

vcffile=sys.argv [1]
window_size=int(sys.argv [2])
overlap=int(sys.argv [3])

inf=gzip.open (vcffile, "rt")
out = open ("homozygosity.txt", "w")
out.write ("startpos\tendpos\thap\tfreq\texpt\tobs\tchi\tLOH\tEOH\n")
hap1=list ()
hap2=list ()
pos=list ()
nvar = -1 
minimum_exp =3
nwin = 0

def get_p_values () :
    #midpos = pos [ int((window_size+1)/2) ]
    startpos =pos [0]
    endpos= pos [-1]
    print (nwin,startpos,endpos)
    homo=dict()
    freq = dict ()
    for myhap1,myhap2 in zip (hap1, hap2):
        freq [myhap1] = freq.get (myhap1,0) +1
        freq [myhap2] = freq.get (myhap2,0) +1
        if myhap1 == myhap2:
            homo [myhap1] = homo.get (myhap1,0) + 1
            
    freq = {hap : round(freq[hap]/nhap,4)  for hap in freq}

    for myhap in freq :
        myfreq = freq.get (myhap,0)
        exp = round(myfreq * myfreq * nani)
        obs = homo.get (myhap,0)
        if exp > minimum_exp :
            chi= (exp - obs)**2 / exp 
            loh=binom_test(obs, nani, myfreq**2, alternative='less')            
            eoh=binom_test(obs, nani, myfreq**2, alternative='greater')            
            out.write (f"{startpos}\t{endpos}\t{myhap}\t{myfreq}\t{exp}\t{obs}\t{chi}\t{loh}\t{eoh}\n")

for line in inf:
    if line [0:2] != "##":
        if line [0:6] == "#CHROM":
            spl=line.rstrip ().split ()
            ids = spl [9:]
            nhap = len (ids) *2
            nani = len (ids)
        else:
            spl=line.rstrip ().split ()
            if "," not in spl[3] and "," not in spl[4]:       
                gts = spl [9:]
                nvar+=1
                pos.append (int(spl [1]))
                if nvar == window_size:
                    nwin +=1
                    get_p_values ()
                    nvar = overlap 
                    hap1 = [hap[-overlap:]  for hap in hap1] #keep the last "n=overlap" genotypes in the haplotype
                    hap2 = [hap[-overlap:]  for hap in hap2]
                    pos  = pos [-overlap:]
                for i, gt in enumerate (gts):
                    if nvar ==0:  #this will run only once ! initializes a haplotype for every id
                        hap1.append(gt [0])
                        hap2.append(gt [2])
                    else:
                        hap1 [i]= hap1 [i] + gt [0] 
                        hap2 [i]= hap2 [i] + gt [2]
inf.close()
out.close()

