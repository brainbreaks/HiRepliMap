import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os 
import argparse

def identify_IZs(clusters):
    IZs = []
    for i in range(1,len(clusters)-1):
        if clusters[i-1] > clusters[i]:
            start=i
            while i < len(clusters)-1 and clusters[i] == clusters[i + 1]:
                i+=1
            if i < len(clusters)-1 and clusters[i+1] > clusters[i]:
                for j in range(start,i+1):
                    IZs.append(j)
                
    return IZs
    
def identify_TZs(clusters):
    TZs = []
    for i in range(1,len(clusters)-1):
        if clusters[i-1] < clusters[i]:
            start=i
            while i < len(clusters)-1 and clusters[i] == clusters[i + 1]:
                i+=1
            if clusters[i+1] < clusters[i]:
                for j in range(start,i+1):
                    TZs.append(j)
                
    return TZs

def identify_TTR(clusters):
    L_TTR = []
    R_TTR = []
    for i in range(1,len(clusters)-1):
        if clusters[i-1] <= clusters[i] <= clusters[i+1]:
            for j in range(i,i+1):
                R_TTR.append(j)
        elif clusters[i-1] >= clusters[i] >= clusters[i+1]:
            for k in range(i,i+1):
                L_TTR.append(k)
    
    unique_L_TTR=[]
    [unique_L_TTR.append(e) for e in L_TTR if e not in unique_L_TTR]

    unique_R_TTR=[]
    [unique_R_TTR.append(e) for e in R_TTR if e not in unique_R_TTR]

    return unique_L_TTR, unique_R_TTR

def identify_SBS(clusters): # Small breakage sites
    SBS = []
    for i in range(1,len(clusters)-1):
        if clusters[i-1] == clusters[i]:
            start=i
            while i < len(clusters)-1 and clusters[i] == clusters[i + 1]:
                i+=1
           
            for j in range(start,i+1):
                SBS.append(j)
                
    return SBS



def identify_features(data):
    
    clusters=np.concatenate([data["Max"]])
    
    IZs = identify_IZs(clusters)
    TZs = identify_TZs(clusters)
    L_TTR, R_TTR =identify_TTR(clusters)
    SBS=identify_SBS(clusters)

    data["manually_annot"]="NA"
    data.loc[L_TTR, "manually_annot"]="L_TTR"
    data.loc[R_TTR, "manually_annot"]="R_TTR"
    data.loc[IZs,"manually_annot"]="IZs"
    data.loc[TZs, "manually_annot"]="TZs"
    data.loc[SBS, "manually_annot"]="SBs"

    return(data)


def manually_check(sub,out):
    
    sub["manually_annot2"]=sub["manually_annot"]
    IZs_index=sub[sub["manually_annot2"]=="IZs"]
    TZs_index=sub[sub["manually_annot2"]=="TZs"]
    R_TTR_index=sub[sub["manually_annot2"]=="R_TTR"]
    L_TTR_index=sub[sub["manually_annot2"]=="L_TTR"]
    
    IZs_annot=[]
    TZs_annot=[]
    CTR_annot=[]
    RTTR_temp=[]
    LTTR_temp=[]

    # Check IZs
    for i in IZs_index.index:
        if sub["manually_annot2"][i]=="IZs":
            end=i
            #print(i)
            while end < len(sub)-1 and sub["manually_annot2"][end+1] == "SBs":
                end+=1
              
            for j in range(i,end+1):
                IZs_annot.append(j) 
    
    print("IZs: done")

    # Check TZs and CTR
    
    for i in TZs_index.index:
        if sub["manually_annot2"][i]=="TZs":
            end=i
            while end < len(sub)-1 and sub["manually_annot2"][end+1] == "SBs":
                end+=1
            
            for j in range(i,end+1):
                TZs_annot.append(j) 

    
    print("TZs: done")
    
    # Check TTRs
    
    for i in R_TTR_index.index:
        start=i
        while start < len(sub)-1 and sub["manually_annot2"][start+1] == "SBs":
            start=start+1
        else:
            t=sub["manually_annot2"][i:start+1].values
            RTTR_temp.append(t)

    RTTR_temp_df=pd.DataFrame({"RTTR_temp":RTTR_temp})
    RTTR_temp_df["match"]="No"
    RTTR_temp_df["n"]=0

    for i in range(len(RTTR_temp_df)):
        input_list=RTTR_temp[i]
        input_str = ', '.join(input_list)
        #print(input_str)
        pattern = r"R_TTR(?:, SBs){0,2}, SBs$"
        match = re.search(pattern, input_str)
        if match:
            #print("Pattern matched!")
            RTTR_temp_df.loc[i, "match"] = "Yes"
            RTTR_temp_df.loc[i, "n"] = len(RTTR_temp[i])
        #else:
            #print("No match found.")

    
    R_TTR_annot=[]
    for i in RTTR_temp_df[RTTR_temp_df["match"]=="Yes"].index:
        k=R_TTR_index.index[i]
        for j in range(k,k+RTTR_temp_df["n"][i]):
            R_TTR_annot.append(j)

    

    for i in L_TTR_index.index:
        start=i
        while start > 0 and sub["manually_annot2"][start-1] == "SBs":
            start=start-1
        else:
            t=sub["manually_annot2"][start:i+1].values
            LTTR_temp.append(t)
    
    LTTR_temp_df=pd.DataFrame({"LTTR_temp":LTTR_temp})
    LTTR_temp_df["match"]="No"
    LTTR_temp_df["n"]=0

    for i in range(len(LTTR_temp_df)):
        input_list=LTTR_temp[i]
        input_str = ', '.join(input_list)
        #print(input_str)
        pattern = r"^SBs(?:, SBs){0,2}, L_TTR$"
        match = re.search(pattern, input_str)
        if match:
            #print("Pattern matched!")
            LTTR_temp_df.loc[i, "match"] = "Yes"
            LTTR_temp_df.loc[i, "n"] = len(LTTR_temp[i])
            
        #else:
            #print("No match found.")

    L_TTR_annot=[]
    
    for i in LTTR_temp_df[LTTR_temp_df["match"]=="Yes"].index:
        k=L_TTR_index.index[i]
        for j in range(k-LTTR_temp_df["n"][i]+1,k+1):
            L_TTR_annot.append(j)

    print("TTRs: done")
    annot_lib=dict(zip(list(["L_TTR","R_TTR","TZs","IZs"]),list([L_TTR_annot,R_TTR_annot,TZs_annot,IZs_annot])))
    


    for k in annot_lib.keys():
        sub_list=annot_lib[k]
        sub.loc[sub_list,"manually_annot2"]=k

    print("Correct 1st Annotations: done")


    # 2nd Correction for CTRs and NA

    TZs_index=sub[sub["manually_annot2"]=="TZs"]
    TZs_index=TZs_index[TZs_index["Max"]>=14]
    
    CTR_temp_df=pd.DataFrame({"CTR_temp":TZs_index.index})
    CTR_temp_df["match"]="No"
    CTR_temp_df["n"]=0

    for i in TZs_index.index:
        end=i
        while end < len(sub)-1 and sub["Max"][end] == sub["Max"][end+1]:
                end+=1
        if abs(end-i)>3:
            for j in range(i,end+1):
                    CTR_annot.append(j)

    sub.loc[np.unique(CTR_annot),"manually_annot2"]="CTRs"

    na_0=sub[sub[["S" + str(i) for i in range(1,17)]].sum(axis=1)==0].index

    sub.loc[na_0,"manually_annot2"]="ND"

    print("Correct 2nd Annotations: done")

    sub.to_csv(out+"RepliFeatures.csv",index=False)

    return sub





parser = argparse.ArgumentParser(prog='python 1_Annot_features.py',formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input",help="Input-file to be analyzed. Format accepted is igv")
parser.add_argument("-o","--output",help="Path to output the result: RepliFeatures.csv")
args = parser.parse_args()

if __name__=="__main__":
    igv = args.input
    out = args.output
    data=pd.read_csv(igv,sep="\t")
    data["MaxS"]=data.loc[:,["S" + str(i) for i in range(1,17)]].idxmax(axis=1)
    data["MaxS"]=data["MaxS"].astype(str)
    data["Max"]=[re.sub("S","",i) for i in data["MaxS"]]
    data["Max"]=data["Max"].astype(float)
    
    res=identify_features(data)
    res_df=manually_check(sub=res,out=out)