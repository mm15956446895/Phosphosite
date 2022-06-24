# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import requests
import pandas as pd

#打开phosphosite_list.txt文件
with open('c:\\Pho\\phosphosite_list.txt') as f:
    wordlist=f.read().split('\n')
print(wordlist)


#uniprot上循环抓取完整蛋白序列
results = []
for pho in wordlist:
    try:
        response = requests.get("https://www.uniprot.org/uniprot/"+pho+".fasta")
        result = response.text
        print(result)
        results.append(result)
    except:
        pass
    continue


#将抓取的蛋白序列保存到phosphosite_complete.txt文件
f2 = open('c:\\Pho\\phosphosite_complete.txt','w')
for line in results:
    f2.write(line+'\n')
f2.close()



#将下载的fasta文件中的序列多行并一行
seq = {}
file=open('c:\\Pho\\phosphosite_complete.txt','r')
fw = open('c:\\Pho\\phosphosite_complete(1).txt','w')
for line in file:
    if line.startswith('>'):   
        name=line.split()[0]
        seq[name]=''
    else:
        seq[name]+=line.replace('\n', '')
file.close() 

for i in seq.keys():
    fw.write(i)
    fw.write('\n')
    fw.write(seq[i])
    fw.write('\n')
fw.close()


#将ID提出来
x = []
fr = pd.read_csv('c:\\Pho\\phosphosite_complete(1).txt',sep = '\n',header = None)
fr_list = list(fr.values)
n = len(fr_list)
for i in range(n):
    if ((i+1)%2):
        x.append(fr_list[i])
Id = []
m = len(x)
for i in range(m):
    Id.append(x[i][0][4:10])
 

#对每一个名字找位点
url = 'https://www.uniprot.org/uniprot/XXX.txt'
Id_site= {}
for t in range(len(Id)):
    uniprot_id=Id[t]
    target = url.replace('XXX',uniprot_id)
    r = requests.get(target)#通过r=request.get（url）构造一个向服务器请求资源的url对象
    playFile = open('s_1.txt','wb')
    for chunk in r.iter_content(50000):## 边下载边存硬盘, chunk_size 可以自由调整为可以更好地适合您的用例的数字
        playFile.write(chunk)
    playFile.close()
    
    df = pd.read_csv('s_1.txt',sep = '\n',header = None) 
    data = df[0]
    cc = []
    for i in range(len(data)):
        if data[i] !='//':
            aa = data[i].split()
            if   aa[0] == 'FT' and aa[1] == 'MOD_RES':
                bb = data[i+1].split()
                if 'Phosphoserine' in bb[1]:
                    if aa[2].isdigit():
                        cc.append(int(aa[2]))
                    else:
                        pass
    Id_site[Id[t]] = cc
    t += 1
    print(t)#t是已输出位点的次数


#找到去同源后的位点   
PS_site ={}
for line in Id:
    PS_site[line] = Id_site[line]


#找所有的S位点
S_site = {}
for i in range(int(n/2)):
    uniprot_Id=Id[i]
    l = fr_list[2*i+1]
    dd=[]
    for j in range(len(l[0])):
        if l[0][j]=='S':
            dd.append(j+1)
    S_site[uniprot_Id] = dd


#负样本位点        
NS_site = {}
ee=[]
for line in Id:
    ee = [q for q in S_site[line] if q not in PS_site[line]]
    NS_site[line]=ee  
   
    
#切片段
lenth = 20    
squence = list(seq.values())   
ee=[]
dd=[]

#正样本位点
for i in range(len(Id)):
    name=list(seq.keys())[i][0:11]
    Max = len(squence[i])                 
    for site in PS_site[Id[i]]:
        Psite = site-1
        if Psite-lenth < 0:
            squence_list = 'O'*lenth +squence[i]
            dd.append(name+'S|'+str(site))
            NPsite = Psite+lenth
            dd.append(squence_list[NPsite-lenth:NPsite+lenth+1])
            continue
        
        if Psite+lenth > Max-1:
            squence_list = squence[i]+'O'*lenth 
            dd.append(name+'S|'+str(site))
            dd.append(squence_list[Psite-lenth:Psite+lenth+1])
            continue
            
        if Psite+lenth > Max-1 and Psite-lenth < 0:
            squence_list = 'O'*lenth +squence[i]+'O'*lenth
            dd.append(name+'S|'+str(site))
            NPsite = Psite+lenth
            dd.append(squence_list[NPsite-lenth:NPsite+lenth+1])
            
        
        else:
            squence_list = squence[i]
            dd.append(name+'S|'+str(site))
            dd.append(squence_list[Psite-lenth:Psite+lenth+1])
     
Plist=[]
for i in range(int(len(dd)/2)):
    if len(dd[i*2+1]) == 41 and dd[i*2+1][20]=='S' :
        Plist.append(dd[i*2])
        Plist.append(dd[i*2+1])

#负样本位点
for i in range(len(Id)):
    name=list(seq.keys())[i][0:11]
    Max = len(squence[i])
    for site in NS_site[Id[i]]:
        Nsite = site-1
        if Nsite-lenth < 0:
            nsquence_list = 'O'*lenth +squence[i]
            ee.append(name+'S|'+str(site))
            NNsite = Nsite+lenth
            ee.append(nsquence_list[NNsite-lenth:NNsite+lenth+1])
            continue
        
        if Nsite+lenth > Max-1:
            nsquence_list = squence[i]+'O'*lenth 
            ee.append(name+'S|'+str(site))
            ee.append(nsquence_list[Nsite-lenth:Nsite+lenth+1])
            continue
        
        if Nsite+lenth > Max-1 and Nsite-lenth < 0:
            nsquence_list = 'O'*lenth +squence[i]+'O'*lenth
            ee.append(name+'S|'+str(site))
            NNsite = Nsite+lenth
            ee.append(nsquence_list[NNsite-lenth:NNsite+lenth+1])
            
        
        else:
            nsquence_list = squence[i]
            ee.append(name+'S|'+str(site))
            ee.append(nsquence_list[Nsite-lenth:Nsite+lenth+1])

Nlist = []        
for i in range(int(len(ee)/2)):
    if len(ee[i*2+1]) == 41 and ee[i*2+1][20]=='S':
        Nlist.append(ee[i*2])
        Nlist.append(ee[i*2+1])
    

#将正负样本位点数据保存到txt文件中
P=open('c:\\Pho\\phosphosite_Positive_sites.txt','a+')#a+ 以附加方式打开可读写的文件。若文件不存在，则会建立该文件;如果文件存在，写入的数据会被加到文件尾后，即文件原先的内容会被保留。
for line in Plist:
    print(line)
    P.writelines(line)
    P.write("\n")
P.close

N=open('c:\\Pho\\phosphosite_Negative_sites.txt','a+')
for line in Nlist:
    print(line)
    N.writelines(line)
    N.write("\n")
N.close()






