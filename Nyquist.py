import matplotlib.pyplot as plt
import numpy as np
import os



def average_multiple_files(folder,data,row, column):
    if data == 'parallel' or 'Parallel':
       M_type = 'CP_RP'
    if data == 'series' or 'Series':
       M_type = 'CS_RS'
    os.chdir(folder)
    List_of_files = os.listdir()
    
    List_of_M_type_files=[];
    for i in range(len(List_of_files)):
        Is_type = M_type in List_of_files[i]
    
        if Is_type == True:
           List_of_M_type_files.append(List_of_files[i])
    
    Frequency =[]
    Capacitance = []
    Resistance = []

    
    for i in range(len(List_of_M_type_files)):
        with open(List_of_M_type_files[i]) as f:
             lines = f.readlines()
             
             for line in lines:
                 split = line.split(" ");
                 Frequency.append(float(split[0]))
                 Capacitance.append(float(split[1]))
                 Resistance.append(float(split[2].replace("\n", " ")))

                                         
    r = row
    c = column
    Frequency = np.array(Frequency)
    Frequency = np.reshape(Frequency,(r,c),order='F')
    Frequency = list(np.average(Frequency, axis=1))            

    Capacitance = np.array(Capacitance)
    Capacitance = np.reshape(Capacitance,(r,c),order='F')
    Capacitance = list(np.average(Capacitance, axis=1))          
                 
    Resistance = np.array(Resistance)
    Resistance = np.reshape(Resistance,(r,c),order='F')
    Resistance = list(np.average(Resistance, axis=1))
    
    filename = open("Average_CP_RP.txt","a")
    for i in range(0,len(Frequency)):
        filename.write("%e %e %e\n"%(Frequency[i],Capacitance[i], Resistance[i]))
    filename.close()
        
                 
    return 0




def convert_avg_origin(folder,filename):
    os.chdir(folder)
    Frequency = [];
    Capacitance = [];
    Resistance = []
    with open(filename) as f:
             lines = f.readlines()
             del (lines[0])
             for line in lines:
                 split = line.split(" ");
                 Frequency.append(float(split[0]))
                 Capacitance.append(float(split[1]))
                 Resistance.append(float(split[3]))
                 
    filename = open("Average_CP_RP.txt","a")
    for i in range(0,len(Frequency)):
        filename.write("%e %e %e\n"%(Frequency[i],Capacitance[i], Resistance[i]))
    filename.close()
    





def get_Z_parallel(folder,filename):
    os.chdir(folder)
    with open(filename) as f:
             lines = f.readlines()
             Frequency =[]
             Z_real = []
             Z_img = []
             Complex = []
            # 0 - frequency, 1- capacitance 2- resistance
             for line in lines:
                 split = line.split(" ");
                 Frequency.append(float(split[0]))
                 Real = float(split[2])/(1 + (2*3.14*float(split[0])*float(split[1])*float(split[2].replace("\n"," ")))**2)
                 Z_real.append(Real)
                 Imaginary = 2 * 3.14 * float(split[0]) * float(split[1]) * float(split[2].replace("\n","")) * Real
                 Z_img.append(Imaginary)
                 Complex.append(complex(Real,-Imaginary))
                 
                 
                 
    return np.array(Frequency), np.array(Complex)


def get_Z_series(folder,filename):
    os.chdir(folder)
    with open(filename) as f:
             lines = f.readlines()
             Frequency =[]
             Z_real = []
             Z_img = []
             Complex = []
            # 0 - frequency, 1- capacitance 2- resistance
             for line in lines:
                 
                 split = line.split(" ");
                 Frequency.append(float(split[0]))
                 Z_real.append(float(split[2].replace("\n"," ")))
                 Imaginary = 1/ (2 * 3.14 * float(split[0]) * float(split[1]))
                 Z_img.append(Imaginary)
                 Complex.append(complex(float(split[2].replace("\n"," ")),Imaginary))
                 
    return np.array(Frequency), np.array(Complex)



def get_Z(folder,filename):
    os.chdir(folder)
    with open(filename) as f:
             lines = f.readlines()
             Frequency =[]
             Complex = []
             Z_real = []
             Z_img = []
            # 0 - frequency, 1- capacitance 2- resistance
             for line in lines:
                 split = line.split("\t");
                 Frequency.append(float(split[0]))
                 Real = float(split[1])
                 Z_real.append(Real)
                 Imaginary = float(split[2].replace("\n",""))
                 Z_img.append(Imaginary)
                 Complex.append(complex(Real,-Imaginary))
    return np.array(Frequency), np.array(Complex)
    

def write(F,Z_1,Z_2,name):
    #name = "SiN_stoi"
    filename = open(name+".txt","a")
    filename.write("Frequency Z-data-real Z-data-imaginary Z-fit-real Z-fit-imaginary \n")
    Frequency = F;
    Z_data = Z_1
    Z_fit = Z_2
    
    for i in range(0,len(Frequency)):
        filename.write("%f %f %f %f %f\n"%(F[i],Z_data.real[i], Z_data.imag[i],Z_fit.real[i], Z_fit.imag[i]))
    filename.close()
    

folder = "C:/Users/JaianthV/switchdrive3/HF_dielectric_chracterization/data/Nyquist/Carlos_data/"
filename = "SiN_non_stoi.txt"
#convert_avg_origin(folder,filename)
#average_multiple_files(folder,'series',22,3)


                                      
frequencies, Z = get_Z(folder,filename)


index = [0,1,2,3,4,5,6]
Z = np.delete(Z,index)

frequencies = np.delete(frequencies,index)

print (Z)
                                         


from impedance.models.circuits import CustomCircuit,Randles





circuit = 'R0-p(R1,C1)-p(W2,C2)'
initial_guess = [450, 10000, None, 1 ,None]
#circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
#initial_guess = [450, 10000, None, 10000, 1, 0.00001 ,None]


'''
initial_guess = [100,0.00001, 10000000,1,1, None]
circuit = 'p(R0,C0)-p(R1,C1)'
initial_guess = [1,0.0000001, 100000, None]

'''


'''
#Al2O3
circuit = 'p(C0,R0,W0)-p(C1,R1,W1)'
initial_guess = [0.00001,104,1,0.000001,990000,10]
#100, 10000, 0.0001, 100, 0.000001, 0.01 ,0.01
'''

'''
#SiN- stoi
circuit = 'p(R0,C0)-p(R1,C1)'
initial_guess = [ 100,0.0001,10000, 0.0000000000001]


circuit = 'p(R1,C1)-p(Wo1,C2)'
initial_guess = [200000,None,1,1, None]

circuit = 'p(R1,C1)-p(W2,C2)'
initial_guess = [20,None,1, None]

'''
'''
#SiN non
circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
initial_guess = [350, 100000, 0.000001, 10000, 1, 0.00001 ,0.0000001]


circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
#initial_guess = [240, 100000, 0.00001, 10000, 10, 0.000001 ,0.0000001]
initial_guess = [540, 100000, 0.000009, 10000, 10, 0.0002, 0.0000001]



circuit = 'R0-p(R1,C1)-p(Wo1,C2)'
initial_guess = [ 1,1000, None, 100, 0.2, None]

circuit = 'R0-p(R1,C1)-p(W2,C2)'
initial_guess = [ 1000000,1000000, None, 1000000, None]


circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
initial_guess = [450, 10000, None, 10000, 1, 0.00001 ,None]



#good one worked with linear scale
#circuit = 'R0-p(R1,C1)-p(R2,W2,C2)'
#initial_guess = [450, 10000, None, 10000, 1 ,None]
circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
initial_guess = [450, 10000, None, 10000, 1, 0.00001 ,None]

'''

circuit = CustomCircuit(circuit, initial_guess=initial_guess, constants={'C1':0.3e-9,'C2':1.1e-9})#, constants={'C1':1e-9}
circuit.fit(frequencies, Z)
print (circuit)
Z_fit = circuit.predict(frequencies)
print ((Z_fit))
'''
randles = Randles(initial_guess=[.001, .0005, .1, .1, 2,1], CPE=True)
randles.fit(frequencies, Z)
Z_fit = randles.predict(frequencies)
print ((Z_fit))

'''
Z_real = Z.real
Z_imaginary = Z.imag

Z_fit_real = Z_fit.real
Z_fit_imaginary = Z_fit.imag


#write(frequencies,Z,Z_fit,"SiN_non_Stoi_fit_fixed_cap_infiniteW_21032021")
plt.plot(Z_real,-Z_imaginary,'o', Z_fit_real,-Z_fit_imaginary, linewidth= 2)
#plt.plot(Z_real,Z_imaginary,'o')
#plt.xscale("log")
#plt.yscale("log")
plt.ylabel('Z Imaginary (Ohms))',fontsize = 12)
plt.xlabel('Z Real (Ohms)',fontsize = 12)
plt.legend(["Data","Fit"], loc='upper left')
plt.title("SiN_non_stoi", fontsize = 15)
plt.show()



import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist

'''
fig, ax = plt.subplots()
plot_nyquist(ax, Z, fmt='o')
plot_nyquist(ax, Z_fit, fmt='-')

#plt.plot(Z.real, Z.imaginary)

plt.legend(['Data', 'Fit'])
plt.show()
'''






