import pandas as pd
import numpy as np
import os
import glob
import argparse
import re

def gausssmoothing (rawcoveragematrix,shape=(3,3),sigma=1):
    def gaussfilt2D(shape=shape,sigma=sigma):
        m,n=[(edge-1)/2 for edge in shape]
        y,x = np.ogrid[-m:m+1,-n:n+1]
        array=np.exp(-(x*x+y*y)/(2*sigma*sigma))
        array[array<np.finfo(array.dtype).eps * array.max()] =0
        sumarray=array.sum()
        if sumarray !=0:
            array/=sumarray
        return array
    avmatrix=np.zeros_like(rawcoveragematrix)

    paddedrawcoveragematrix=np.concatenate((np.array([rawcoveragematrix[2,:] for i in range(int((shape[0]-1)/2))]),rawcoveragematrix,np.array([rawcoveragematrix[-1,:] for i in range(int((shape[0]-1)/2))])))
    paddedrawcoveragematrix=np.pad(paddedrawcoveragematrix,((0,0),(int((shape[0]-1)/2),int((shape[0]-1)/2))),'constant',constant_values=0)
    for i in range(int((shape[0]-1)/2),int(len(rawcoveragematrix)+(shape[0]-1)/2)):
        print (i,'i')
        for j in range(int((shape[0]-1)/2),int(len(rawcoveragematrix[0])+(shape[0]-1)/2)):
            box=np.ma.masked_invalid(paddedrawcoveragematrix[int(i-(shape[0]-1)/2):int(i+(shape[0]-1)/2+1),int(j-(shape[0]-1)/2):int(j+(shape[0]-1)/2+1)])
            
            avmatrix[int(i-(shape[0]-1)/2),int(j-(shape[0]-1)/2)] = np.nansum(np.multiply(box,gaussfilt2D()))
    return (avmatrix)


def scalingto100range(input):
    a_scaled=np.zeros_like(input)
    for i in range(0,len(input)):
        for j in range (len(input[i])):
            a_scaled[i][j]=(input[i][j]/np.sum(input[:,j]))*100
    return (a_scaled)


# Be aware of file names

def repliseq_normalization(repliseq_data_dir,G1_data,out):
    # read data
    print("Read data....")
    S=["S" + str(i) for i in range(1,17)]
    file_paths=glob.glob(repliseq_data_dir+"/*S*_bin50000.bdg")
    directory=[]
    for i in S:
        directory.append(list(filter(lambda x: re.search(i+"_bin50000",x),file_paths))[0])
    df_dict=dict(zip(S,directory))
    for i in df_dict:
        df_dict[i]=pd.read_csv(df_dict[i],sep="\t",names=["repliseq_chrom","repliseq_start","repliseq_end","repliseq_value"],header=None)
        df_dict[i]["repliseq_fraction"]=i
        total=sum(df_dict[i]["repliseq_value"])
        df_dict[i]["RPM_repliseq_value"]=df_dict[i]["repliseq_value"]/total*1000000
    data=pd.concat(df_dict,axis=0)
    Chromosomes=["chr"+str(i) for i in range(1,20)]
    Chromosomes.append("chrX")
    Chromosomes.append("chrY")
    data=data[data["repliseq_chrom"].isin(Chromosomes)]
    data.to_csv(out+"repliseq_Sfractions.csv",index=False)
    
    # normalized to G1
    print("Normalized to G1 .....")
    G1=pd.read_csv(G1_data,sep="\t",header=None,names=["repliseq_chrom","repliseq_start","repliseq_end","repliseq_value"])
    G1["RPM_repliseq_value"]=G1["repliseq_value"]/sum(G1["repliseq_value"])*1000000
    G1["repliseq_fraction"]="G1"
    G1.to_csv(out+"/G1_repliseq.csv",index=False)
    S_data=data.sort_values(["repliseq_chrom","repliseq_start","repliseq_value"])
    G1=G1.sort_values(["repliseq_chrom","repliseq_start","repliseq_value"])
    data_all=pd.concat([S_data,G1],ignore_index=True)
    wide_t=data_all.pivot(columns="repliseq_fraction",values="RPM_repliseq_value",index=["repliseq_chrom","repliseq_start","repliseq_end"])
    for i in S:
        wide_t[i]=np.log2(wide_t[i].values/wide_t["G1"].values+1)
    data_m=np.array(wide_t[S].values.T)
    print("Smoothing .....")
    data_smooth=gausssmoothing(rawcoveragematrix=data_m,shape=(3,3),sigma=1)
    data_smooth=data_smooth.astype("float")
    data_scaled=scalingto100range(data_smooth)
    data_scaled[np.isnan(data_scaled)]=0
    
    S_fractions=["S"+ str(i) for i in range(1,17)]
    
    igv_out=pd.DataFrame(data_scaled.T)
    igv_out["repliseq_chrom"]=wide_t.index.get_level_values("repliseq_chrom").values
    igv_out["repliseq_start"]=wide_t.index.get_level_values("repliseq_start").values
    igv_out["repliseq_end"]=wide_t.index.get_level_values("repliseq_end").values
    igv_out["repliseq_dtype"]="repliseq"
    
    mapper=dict(zip(list(igv_out.columns[0:16]),list(S_fractions)))
    igv_out=igv_out.rename(columns=mapper)
    igv_out=igv_out[[ 'repliseq_chrom','repliseq_start','repliseq_end','repliseq_dtype','S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11','S12', 'S13', 'S14', 'S15', 'S16']]
    igv_out=igv_out[igv_out["repliseq_chrom"].isin(Chromosomes)]
    print("Writing IGV .....")
    igv_out.to_csv(out+"/repliseq_normalized.igv",sep="\t",index=False)
    
    for i in S:
        result=igv_out[['repliseq_chrom','repliseq_start','repliseq_end',i]]
        outpath=out+"/"+i+"_scaled.bedGraph"
        result.to_csv(outpath,sep="\t",index=False,header=False)
        cmd= "sed -i '1i track type=bedGraph' "+outpath  
        os.system(cmd)
        print("Write bedGraph: ",i)
    



parser=argparse.ArgumentParser(prog="python repliseq_normalization.py",formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",help="The path of repliseq data forlder")
parser.add_argument("-g1","--g1", help="The path of G1 data")
parser.add_argument("-o","--out", help="The path of your output folder")


args=parser.parse_args()

if __name__=="__main__":
    repliseq_data_dir=args.input
    G1_data=args.g1
    out=args.out
    repliseq_normalization(repliseq_data_dir,G1_data,out)
    
    