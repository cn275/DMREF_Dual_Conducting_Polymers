import os,sys, math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("GTK")


def averager(mat,averager):
    averager=int(averager)
    averaged=[0 for i in range(int(float(len(mat))/float(averager))+1)]
    for i in range(len(mat)):
        averaged[int(i/averager)]=averaged[int(i/averager)]+float(mat[i])/float(averager)
    del averaged[-1]
    return averaged


def mob_fit(x,A):
    return A*x**2


def data_get(file_name,key,keynum,*positional_parameters,**keyword_parameters):

    keycount=0
    file1 = open(str(file_name))
    preamble=[]
    if ('preamble' in keyword_parameters):
        preadd=1
    else:
        preadd=0

    if preadd==1:
        read_line=file1.readline()
        preamble.append(read_line)
    while keycount < keynum:
        read_line=file1.readline()
        split_line=read_line.split()
        if len(split_line)>0:
            if split_line[0]==key:
                keycount+=1

    if preadd==1: del preamble[-1]
    var_name=[i for i in split_line]
    variables=[[] for i in range(len(var_name))]
    Nvar=len(variables)
    while split_line != []:
        read_line=file1.readline()
        split_line=read_line.split()
        if split_line!=[]:
            skip=0
            if len(split_line) < Nvar:
                skip=1
            if skip==0:
                for i in range(len(split_line)):
                    variables[i].append(split_line[i])

    output=[]
    for i in range(len(var_name)):
        output.append(variables[i])
    output.append(var_name)
    if preadd==1: output=[output,preamble]
    file1.close()
    return output






args=sys.argv
data=data_get(args[1],'t',1)

fs_per_step=int(args[2])
field=float(args[3])


print("fs_per_step, field= ", fs_per_step, field)

block_size=1

if len(args)>4:
    block_size=int(args[4])

for i in range(len(data[-1])):
    data[i]=averager(data[i],block_size)  
    



fs_per_s=10**15
ns_per_step=float(fs_per_step)/float(fs_per_s)
m2_in_A2=10**(-20)


print("ns_per_step: ", ns_per_step)
t=[float(i)*ns_per_step for i in data[0]]
print("t=",t)
outfile=open("mobilities_"+args[1],'w')
outfile.write("mobility\n")

mobil_vec=[0 for i in range(len(data[-1])-1)]


num_plots=len(data[-1])-1
colormap = plt.cm.jet
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])



use_vec=[]
tmin=[]
tmax=[]
    

batch_size=5


ncycle=int(math.ceil(float(num_plots/batch_size)))+int((num_plots%batch_size)!=0)

print("ncycle= ",ncycle)
keep_sum=0

for k in range(ncycle):

    added='0'
    temp_use=[]
    temp_min=[]
    temp_max=[] 

    while added=='0':
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, batch_size)])
        for j in range(min([num_plots-k*batch_size,batch_size])):
            msd=[float(i)*m2_in_A2 for i in data[1+j+k*batch_size]]
#            print("fit= ",fit)
#            fit=[k for k in np.polyfit(t,msd, 2)[0]]
            fit=[k for k in curve_fit(mob_fit,t,msd)[0]]
            print("fit= ",fit)
            print(k*batch_size+j+1 ," fit= ",fit)
            outfile.write("%s\n"%((fit[0]**.5)/(field)/10**10))
#            outfile.write("%s\n"%((fit[0]**.5)/(field)/10**10))
            print("mobility= ",(fit[0]**.5)/(field)/10**10)
            plt.plot(t,msd,label=str(j+1),linewidth=2)
            mobil_vec[j+k*batch_size]=(fit[0]**.5)/(field)/10**10


        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, batch_size)])
        for j in range(min([num_plots-k*batch_size,batch_size])):
            msd=[float(i)*m2_in_A2 for i in data[1+j+k*batch_size]]
            fit=[k for k in curve_fit(mob_fit,t,msd)[0]]
            plt.plot(t,[j**2*fit[0] for j in t],'--',label=str(j+1),linewidth=2)
#            plt.plot(t,[j**2*fit[0]+fit[1] for j in t],'--',label=str(j+1),linewidth=2)

        plt.ylabel(r"MSD[m^2]",fontsize=20)
        plt.xlabel(r"t[ns]",fontsize=20)
        plt.legend(loc=2)

        plt.show()

        narg=None
        while narg==None:
            print("How many do we use?")
            raw_in=input(" ")    
            raw_in=raw_in.split()
            if len(raw_in)==1: 
                try:
                    narg=int(raw_in[0])
                except ValueError:
                    print("The input isn't valid")
            else:
                print("we can't take in more than one input")

        keep_sum+=narg


        if narg==0: 
            temp_use=[]
            temp_min=[]
            temp_max=[]
            added=1 
        else:
             
            plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, batch_size)])
            for j in range(min([num_plots-k*batch_size,batch_size])):
                msd=[float(i)*m2_in_A2 for i in data[1+j+k*batch_size]]
#                fit=np.polyfit(t,msd, 2)
                fit=[k for k in curve_fit(mob_fit,t,msd)[0]]
                outfile.write("%s\n"%((fit[0]**.5)/(field)/10**10))
                plt.plot(t,msd,label=str(j+1),linewidth=2)
                mobil_vec[j+k*batch_size]=(fit[0]**.5)/(field)/10**10

            temp_use=[]
            temp_min=[]
            temp_max=[]
            for i in range(narg): 
                input_found=0
                while input_found==0:
                    use_input=input("")
                    use_input=use_input.split()
                    
                    
     
                    if len(use_input)==3:   
                        input_found=1         
                        try:
                            temp_use.append(int(use_input[0])+batch_size*k)
                        except ValueError:
                            input_found=0
                                        
                        try:
                            temp_min.append(float(use_input[1])*10**(-9))
                        except ValueError:
                            input_found=0

                        try:
                            temp_max.append(float(use_input[2])*10**(-9))
                        except ValueError:
                            input_found=0 

                    elif len(use_input)==1:
                        input_found=1
                        try:
                            temp_use.append(int(use_input[0])+batch_size*k)
                        except ValueError:
                            input_found=0

                        temp_min.append(0)
                        temp_max.append(t[-1]**10**(-9)) 
 
                    else: print("retype entry, it's not the right length")                   

                 
                new_t=[float(j) for j in t if j>=temp_min[-1] and j<=temp_max[-1]]
                new_msd=[ float(data[temp_use[-1]][j])*m2_in_A2 for j in range(len(data[0])) if t[j]>=temp_min[-1] and t[j]<=temp_max[-1]]
                fit=[k for k in curve_fit(mob_fit,new_t,new_msd)[0]] 
                print("fit= ",fit, (fit[0]**.5)/(field)/10**10)
                plt.plot(new_t,[j**2*fit[0] for j in new_t],'--', label="fit%s"%(i+1))

            plt.ylabel(r"MSD[m^2]",fontsize=20)
            plt.xlabel(r"t[ns]",fontsize=20)
            plt.legend(loc=2)

            plt.show()
                        
            print("So is the fit satisfactory?")
            added=input("")
        
    use_vec.extend(temp_use)
    tmin.extend(temp_min)
    tmax.extend(temp_max)




        
outfile.close()
mob_sum=0
used_mobil=[]
print("tmin, tmax= ", tmin, tmax)

for i in range(len(use_vec)):
    
    new_t=[float(j) for j in t if j>=tmin[i] and j<=tmax[i]]
    new_msd=[ float(data[use_vec[i]][j])*m2_in_A2 for j in range(len(data[0])) if t[j]>=tmin[i] and t[j]<=tmax[i]]
    fit=[k for k in curve_fit(mob_fit,new_t,new_msd)[0]]

    print(use_vec[i] ," fit= ",fit)
    mobil=((fit[0]**.5)/(field)/10**10)
 
    print(mobil)
    mob_sum+=mobil
    used_mobil.append(mobil_vec[i-1])

print()
print("Avg mobility: ")
print(mob_sum/float(len(use_vec)))
print("n_use= ",len(use_vec))

