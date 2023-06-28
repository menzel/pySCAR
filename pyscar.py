#!/usr/bin/env python3

import pandas as pd
import sys

pd.options.mode.chained_assignment = 'raise'

global name
name =sys.argv[1]
copy_numbers = pd.read_csv(name, sep="\t")
sizelimit = 15e+06

chrominfo_full = pd.read_csv("chrominfo", sep=" ")
chrominfo_hg37 = chrominfo_full[chrominfo_full['genome'] == "hg37"][["chrom","centstart","centend"]]
chrominfo = chrominfo_hg37


#################3

def preprocess(dat):
    unknown_chr = set(dat['Chromosome'].unique()) - set(["chr"+ str(c) for c in list(range(1,23))+['X','Y']])

    if len(unknown_chr) > 0:
        print("Chromsome error: " + str(unknown_chr))
        sys.exit(1)

    dat= dat[~dat.Chromosome.isin(["chr"+i for i in ['X','Y','x','y','23','24']])]

    #if len(dat['SampleID'].unique()) > 1:
    if len(dat['SampleID'].unique()) > 1:
        print("Multiple Patients in one file not implemented (yet!)")
        sys.exit(1)


    def switch(row):
        a,b = row["A_cn"], row["B_cn"]
        row["A_cn"] = max(a,b)
        row["B_cn"] = min(a,b)

        return row

    return dat.apply(switch, axis=1)


#shrink
def shrink_chrom(dat):
    #pd.options.mode.chained_assignment = None

    new_dat = []

    if len(dat) < 2 :
        return dat

    tmp_start = dat.iloc[0,:]['Start_position']
    tmp_max = dat.iloc[0,:]['End_position']
    tmp_sum = dat.iloc[0,:]['total_cn']

    for i in range(1,len(dat)):

        last = dat.iloc[i-1,:]
        row = dat.iloc[i,:]

        # if new segment
        if not (row['A_cn'] == last['A_cn'] and row['B_cn'] == last['B_cn']):
            new_dat.append([last['SampleID'],last['Chromosome'],tmp_start, tmp_max, tmp_sum, last['A_cn'], last['B_cn'], last['ploidy']])

            tmp_start = row['Start_position']
            tmp_max = row['End_position']
            tmp_sum = row['total_cn']

        else: 
            tmp_sum += row['total_cn']
            tmp_max = max(tmp_max, row['End_position'])

    new_dat.append([last['SampleID'],last['Chromosome'],tmp_start, tmp_max, tmp_sum, row['A_cn'], row['B_cn'], row['ploidy']])

    df = pd.DataFrame(new_dat)

    df.columns = ["SampleID","Chromosome","Start_position","End_position","total_cn","A_cn","B_cn","ploidy"]
    return df
    

#shrink wrapper
def shrink(dat):

    if len(dat) == 0:
        return dat

    new = []

    for patient in dat['SampleID'].unique():
        for chrom in dat['Chromosome'].unique():
            shrunk = shrink_chrom(dat[(dat['Chromosome'] == chrom) & (dat['SampleID'] == patient)])
            new.append(shrunk)

    return pd.concat(new)
    


def calc_hrd(dat):

    for patient in dat['SampleID'].unique():
        seqSamp = dat.loc[dat['SampleID'] == patient]

        # remove chromsomes where all B_cn is 0 
        # do not count loss of whole chromosome
        for chrom in seqSamp['Chromosome'].unique():
            if seqSamp.loc[dat['Chromosome'] == chrom]['B_cn'].sum() == 0:
                seqSamp = seqSamp[seqSamp['Chromosome'] != chrom]

        seqSamp = preprocess(seqSamp)
        seqSamp = shrink(seqSamp)
        seqSamp.loc[seqSamp["A_cn"] > 1,"A_cn"] = 1
        seqSamp = shrink(seqSamp)



        loh = seqSamp[(seqSamp["B_cn"] == 0) & (seqSamp["A_cn"] != 0)]
        loh = loh[loh["End_position"] - loh["Start_position"] > sizelimit]

    #[print(name.split("/")[-1] + " " + x["Chromosome"] + " " + str(x["Start_position"]) + " " + str(x["End_position"])) for i,x in loh.iterrows()]
    
    return len(loh)


def calc_tai(dat, minsize = 1e+06):

    dat = preprocess(dat)


    ai = []
    events = {}

    dat = dat[dat["End_position"] - dat["Start_position"] > minsize]

    dat = shrink(dat)


    def calc_ai(row):
        if row['B_cn'] == row['A_cn']:
            return 0
        return 2

    def calc_ai_odd(row):
        if row['B_cn'] + row['A_cn'] == ploidy and row['A_cn'] != ploidy:
            return 0
        return 2

    dat['AI'] = -1

    for chrom in dat['Chromosome'].unique():

        slice_c = dat[dat['Chromosome'] == chrom]

        global ploidy
        ploidy = -1
        lengths = 0

        for cnv in set(slice_c["A_cn"].unique())-set([0]):
            cnv_slice = slice_c[slice_c["A_cn"] == cnv]

            tmp = sum(cnv_slice["End_position"] - cnv_slice["Start_position"])
            if tmp > lengths:
                ploidy = cnv
                lengths = tmp

        if ploidy == 1 or ploidy % 2 == 0:
            dat.loc[dat['Chromosome'] == chrom, 'AI'] = dat[dat['Chromosome'] == chrom].apply(calc_ai, axis=1)
        else:
            dat.loc[dat['Chromosome'] == chrom, 'AI'] = dat[dat['Chromosome'] == chrom].apply(calc_ai_odd, axis=1)
            

        slice_c = dat[dat['Chromosome'] == chrom]


        if len(slice_c) == 1 and slice_c.iloc[0]['AI'] != 0:
                dat.iat[(dat['Chromosome'] == chrom).argmax(),8] = 3
        else:
            if slice_c.iloc[0]['AI'] == 2 and slice_c.iloc[0]['End_position'] < chrominfo[chrominfo['chrom'] == chrom]['centstart'].values[0]:
                #argmax of truth series to select the first entry for that chrom
                dat.iat[(dat['Chromosome'] == chrom).argmax(),8] = 1

            #print(chrominfo)
            #print(chrominfo[chrominfo['chrom'] == chrom]['centend'].values[0])

            if slice_c.iloc[-1]['AI'] == 2 and slice_c.iloc[-1]['Start_position'] > chrominfo[chrominfo['chrom'] == chrom]['centend'].values[0]:
                dat.iat[len(dat) - 1 - ((dat['Chromosome'] == chrom).iloc[::-1].argmax()),8] = 1

    events["Telomeric AI"] = len(dat[dat['AI'] == 1])
    events["Interstitial AI"] = len(dat[dat['AI'] == 2])
    events["Whole chr AI"] = len(dat[dat['AI'] == 3])


    #print("TAI:")
    #print(dat)
    #print(dat[dat['AI'] == 1])
    return events["Telomeric AI"]


def calc_lst(dat):

    lst = 0

    dat = preprocess(dat)

    for chrom in dat['Chromosome'].unique():

        current = dat[dat['Chromosome'] == chrom]
        if len(current) < 2 :
            continue

        # split chromosome in arms 
        parm = []
        qarm = []

        parm = current[current['Start_position'] <= chrominfo[chrominfo['chrom'] == chrom]['centstart'].values[0]]
        qarm = current[current['End_position'] >= chrominfo[chrominfo['chrom'] == chrom]['centend'].values[0]]

        qarm = shrink(qarm)
        # set start of second arm to end of center? 
        qarm.iloc[0,2] = chrominfo[chrominfo['chrom'] == chrom]['centend']

        if len(parm) > 0:
            parm = shrink(parm)
            # set end of first arm to start of center? 
            parm.iloc[-1,3] = chrominfo[chrominfo['chrom'] == chrom]['centstart']



        #remove and shrink all intervals below 3e6
        while len(parm[parm['End_position'] - parm['Start_position'] < 3e6]) > 0:
            # only filter first 
            parm = parm.drop(parm[parm['End_position'] - parm['Start_position'] < 3e6].index[0])
            parm = shrink(parm)

        while len(qarm[qarm['End_position'] - qarm['Start_position'] < 3e6]) > 0:
            #tmp = qarm[~(qarm['End_position'] - qarm['Start_position'] < 3e6)]

            qarm = qarm.drop(qarm[qarm['End_position'] - qarm['Start_position'] < 3e6].index[0])
            qarm = shrink(qarm)


        if len(parm) > 1:
            parm = parm.assign(lst=pd.Series(parm.apply(lambda x: x['End_position'] - x['Start_position'] >= 10e6 ,axis=1)))

            for i in range(1,len(parm)):
                last = parm.iloc[i-1]
                current = parm.iloc[i]

                if current["lst"] and last["lst"] and current["Start_position"] - last["End_position"] < 3e6:
                    lst += 1


        if len(qarm) > 1:
            qarm = qarm.assign(lst=pd.Series(qarm.apply(lambda x: x['End_position'] - x['Start_position'] >= 10e6 ,axis=1)))

            for i in range(1,len(qarm)):
                last = qarm.iloc[i-1]
                current = qarm.iloc[i]

                if current["lst"] and last["lst"] and current["Start_position"] - last["End_position"] < 3e6:
                    lst += 1

        #print(chrom,lst)

    return lst


#####################

# call preprocess first here
tai = calc_tai(copy_numbers)
hrd = calc_hrd(copy_numbers)
lst = calc_lst(copy_numbers)

print("LOH TAI LST")
print(hrd,tai,lst,end=" ")
