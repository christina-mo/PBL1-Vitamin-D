#Authors: Katherine Bair, Kay Lyon, Christina Mo
#Duke University
#BME 260 PBL #1
#March 3, 2025

#Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import math

#EQUATIONS
#Define variables
D2mw_g = 396.648 #g/mol
D2mw_mcg = 396.648e6 #mcg/mol

D3mw_g = 384.64 #g/mol
D3mw_mcg = 384.64e6 #mcg/mol

Camw_g = 40.078 #g/mol
Camw_mcg = 40.078e6 #mcg/mol

CDmw_g = 400.64 #g/mol
CDmw_mcg = 400.64e6 #mcg/mol

CTmw_g = 416.64 #g/mol
CTmw_mcg = 416.64e6 #mcg/mol

Vol_l = 5 #l of blood (typical concentrations are serum concentrations)
Vol_ml = 5000 #ml of blood (typical concentrations are serum concentrations)

DAbsrate = .80 #MOLAR OR MASS ratio of MOLES OR MASS of vitamin d from diet to mass out of intestine

DtoDserumConvert = 2.17 #(mcg * L)/nmol #value determined from ideal vitamin d flow into liver (on same order of magnitude as literature)
# D --> calcidiol
DSigmaRxn1 = -1
CDSigmaRxn1 = 1
Rrxn1ratio = 0.43

Rrxn2ratio_mass = 2.22 #picograms ct/nanograms cd
Rrnx2ratio_moles = Rrxn2ratio_mass*(CDmw_g/CTmw_g)*(1E-3) #moles ct/moles cd

CaPassiveAbsRate = 0.13

SurfAreaDuo = 18600 #cm^2

SurfAreaWell = 4.67 #cm^2

minPerDay = 1440 #minutes in a day



#Vitamin D and Ca entering intestine

D2_Diet_m = 15 # [mcg/day] Total amount of vitamin D consumed per day
D2_Diet_n = D2_Diet_m/D2mw_mcg #[moles/day]

D3_Sun_m = 60 # [mcg/day]
D3_Sun_n = D3_Sun_m/D3mw_mcg #[moles/day]

Ca_Dose_m = 330000 # [mcg] Total amount of Ca consumed per meal
Ca_Dose_n = Ca_Dose_m/Camw_mcg # [moles] of Ca entering intestine from one meal

Ca_Daily_m = 3*Ca_Dose_m
Ca_Daily_n = 3*Ca_Dose_n


#Model of vitamin D entering and leaving storage
D2_Absorbed_m = DAbsrate * D2_Diet_m #mcg/day --> amount of vit D absorbed by intestine per day
D2_Absorbed_n = D2_Absorbed_m/D2mw_mcg #mol/day

D2_Waste_m = D2_Diet_m - D2_Absorbed_m #mcg/day --> amount of vit D that is immediately excreted from intestine
D2_Waste_n = D2_Waste_m / D2mw_mcg #mols/day

#Flows into Storage box --> Treating storage like mixer
D_Storage_in_n = D2_Absorbed_n + D3_Sun_n #mol/day
D_Storage_in_m = D2_Absorbed_m + D3_Sun_m #mcg/day
DmwWeightAvg_mcg = ((D2_Absorbed_n * D2mw_mcg)+(D3_Sun_n*D3mw_mcg))/D_Storage_in_n #mcg/mol

D_Storage_to_liver_conc = D_Storage_in_m /DtoDserumConvert #nmol/L

D_Storage_to_liver_n = (D_Storage_to_liver_conc * Vol_l)/ (1e9) #mol/day
D_Storage_to_liver_m = D_Storage_to_liver_n * DmwWeightAvg_mcg #mcg/day

D_Storage_out_n = D_Storage_in_n - D_Storage_to_liver_n #mol/day -- WASTE
D_Storage_out_m = D_Storage_in_m - D_Storage_to_liver_m #mcg/day --WASTE


#Model conversion of vitamin D to calcidiol
D_liver_in_n = D_Storage_out_n #mol/day
D_liver_in_m = D_Storage_out_m #mcg/day

Rrxn1 = Rrxn1ratio * D_Storage_out_n #mol/day

D_liver_out_n = D_liver_in_n + Rrxn1* DSigmaRxn1 #mol/day
D_liver_out_m = D_liver_out_n * DmwWeightAvg_mcg #mcg/day

CD_generated_n = CDSigmaRxn1*Rrxn1 #mol/day
CD_generated_m = CD_generated_n * CDmw_mcg #mcg/day

CD_liver_to_kidney_n = CD_generated_n


#Function takes kidney function as an input and returns a factor to multiply by total calcidiol in the body to find calcitriol
def find_CT_generated(input):
  CD_kidney_in = CD_liver_to_kidney_n
  CD_consumed = CD_kidney_in*(input)*((CDmw_g*1e9)/(CTmw_g*1e12))
  CD_waste = CD_kidney_in - CD_consumed
  CT_generated = CD_consumed
  CT_kidney_to_intestine = CT_generated

  return(CT_generated)



#Gives calcitriol as a concentration 
def find_CT_conc(input):
  CT_in = find_CT_generated(input) #moles
  CT_conc_nM = CT_in*1E9/Vol_l #nmoles/L = nM
  return CT_conc_nM


#Returns molar calcium absorption as a function of kidney health
def find_absorbed_calcium(input):
  Ca_absorbed_active = (0.2389*math.log(find_CT_conc(input)) + 1.6)*(SurfAreaDuo/SurfAreaWell)* minPerDay *1E-9 # mol/meal
  print('Ca absorbed active: ' + str(Ca_absorbed_active))
  Ca_absorbed_passive =  CaPassiveAbsRate*Ca_Dose_n
  print('Ca absorbed passive:' + str(Ca_absorbed_passive))
  Ca_absorbed_n_meal = Ca_absorbed_active + Ca_absorbed_passive
  Ca_absorbed_n_day = 3*Ca_absorbed_n_meal
  return Ca_absorbed_n_day


#Returns Calcium absorption perecentage as a function of kidney health
def find_absorption_percent(input):
  Ca_absorption_percent = 100*(find_absorbed_calcium(input)/Ca_Daily_n) # of total intake
  return Ca_absorption_percent


#FIGURES
#Calcitriol concentration vs. calcium absorption
ct_conc = np.arange(0.01, 0.05, 0.001)
ca_abs = np.zeros(len(ct_conc))
ca_abs_pct = np.zeros(len(ct_conc))
for i in range(len(ct_conc)):
  ca_abs[i] = find_absorbed_calcium(ct_conc[i])
  ca_abs_pct[i] = find_absorption_percent(ct_conc[i])

plt.figure(figsize=(12,8))
plt.plot(ct_conc, ca_abs_pct)
plt.xlabel ('Calcitriol Concentration (nM)', fontsize = 22,fontweight='normal')
plt.ylabel ('Ca Intake Absorbed (%)', fontsize = 22, fontweight='normal')
plt.title('Effect of Calcitriol Concentration on Calcium Absorption' , fontsize = 26, fontweight='normal')
plt.rcParams['font.size'] = 20
plt.show()

#Models calcitriol concentration and calcium absorption at four different stages: no disease, severe CKD, moderate CKD, and after kidney transplant
Disease_states = ['No Disease', 'Moderate CKD', 'Severe CKD',  'Kidney Transplant']
Calcitriol_conversion_ratio = [2.22, 1.36, 1.11, 1.77]

Calcitriol_produced = np.zeros(4)
Calcitriol_concentration = np.zeros(4)
Calcium_absorbed = np.zeros(4)
Calcium_absorbed_percent = np.zeros(4)

for i in range(len(Disease_states)):
  Calcitriol_produced[i] = find_CT_generated(Calcitriol_conversion_ratio[i])
  Calcitriol_concentration[i] = find_CT_conc(Calcitriol_conversion_ratio[i])
  Calcium_absorbed[i] = find_absorbed_calcium(Calcitriol_conversion_ratio[i])
  Calcium_absorbed_percent[i] = find_absorption_percent(Calcitriol_conversion_ratio[i])

'''
print(Calcitriol_produced)
print(Calcitriol_concentration)
print(Calcium_absorbed)


'''

print('Calcium absorption fractions:' + str(Calcium_absorbed_percent))


plt.figure(figsize=(12,8))
bars = plt.bar(Disease_states, Calcitriol_produced)
bars[3].set_color('darkblue')
plt.xlabel ('Kidney Condition', fontsize = 22,fontweight='normal')
plt.ylabel ('Calcitriol Produced (mol/day)', fontsize = 22, fontweight='normal')
#plt.title('Predicted Effect of Kidney Disease on Calcitriol Production' , fontsize = 26, fontweight='normal')
plt.rcParams['font.size'] = 20

'''
plt.figure()
plt.bar(Disease_states, Calcitriol_concentration)
plt.xlabel ('Kidney Condition')
plt.ylabel ('Calcitriol concentration [nM]')
plt.title('Predicted Effect of Kidney Disease on Calcitriol concentration')

plt.figure()
plt.bar(Disease_states, Calcium_absorbed)
plt.xlabel ('Kidney Condition')
plt.ylabel ('Moles of calcium absorbed per day')
plt.title('Predicted effect of kidney disease on calcium absorption')
'''

plt.figure(figsize=(12,8))
bars2 = plt.bar(Disease_states, Calcium_absorbed_percent)
bars2[3].set_color('darkblue')
plt.xlabel ('Kidney Condition', fontsize = 22,fontweight='normal')
plt.ylabel ('Ca Intake Absorbed (%)', fontsize = 22, fontweight='normal')
plt.rcParams['font.size'] = 20
#plt.title('Predicted Effect of Kidney Disease on Calcium Absorption' , fontsize = 24, fontweight='normal')

plt.show()



