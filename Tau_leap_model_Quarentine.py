# Tau leap code for SAIRQ  agricultural fair model
# Developed 11/25/24 Dana C. Pittman Ratterree

import numpy as np
import pylab as pl

### Parameters ###
##Movement between compartments
#PIG
R0 = 4 #The reproductive number of influenza in pigs
gammas = 1/5 # recovery rate of swine
beta = R0 * gammas # transmission rate
delta = 0.83 # probability of becoming asymptmatic
chi = .20 # sucess at identifying symptomatic pigs
theta = 1/1 #delay at identifying symptomatic pigs

#HUMAN
Cm = 60 # contact duration per minute members
Pm = 0.00168 # probability of transmission members
Ca = 5 # contact duration per minure attendees
Pa = 0.00095 # probability of transmission attendees
kappa = 1/2 # rate exposed become infected
gammah = 1/5 # recovery rate

ND=MaxTime=8

tau=0.1 # setting the value of tau

##Defining intial conditions
#PIG
Ss0, As0, Is0, Rs0, Qs0 = 203, 5, 0, 0, 0

N=208

#MEMBER
Sm0, Em0, Im0, Rm0 = 90, 0, 0, 10

# ATTENDEE
Sa0, Ea0, Ia0, Ra0 = 10042, 0, 0, 4868

# generating the arrary for the inital population values
INPUT = np.array((Ss0, As0, Is0, Rs0, Qs0, Sm0, Em0, Im0, Rm0, Sa0, Ea0, Ia0, Ra0))


### DEFING RATE EVENTS  ###
def stoch_eqs(INP):
    V=INP
    Rate = np.zeros((11))
    Change = np.zeros((11, 13))
    
    Rate[0]= beta*V[0]*((V[1]+V[2])/N)*delta
    Change[0, :]= ([-1,+1,0,0,0,0,0,0,0,0,0,0,0]);
    
    Rate[1]= beta*V[0]*((V[1]+V[2])/N)*(1-delta)
    Change[1, :]= ([-1,0,+1,0,0,0,0,0,0,0,0,0,0]);
    
    Rate[2]= gammas*V[1]
    Change[2, :]= ([0,-1,0,+1,0,0,0,0,0,0,0,0,0]);
    
    Rate[3]= chi*theta*V[2]
    Change[3, :]= ([0,0,-1,0,+1,0,0,0,0,0,0,0,0]);
    
    Rate[4]= (1-chi)*gammas*V[2] 
    Change[4, :]= ([0,0,-1,+1,0,0,0,0,0,0,0,0,0]);
    
    Rate[5]= Cm*Pm*V[5]*((V[1]+V[2])/N) 
    Change[5, :]= ([0,0,0,0,0,-1,+1,0,0,0,0,0,0]);
    
    Rate[6]= kappa*V[6]
    Change[6, :]= ([0,0,0,0,0,0,-1,+1,0,0,0,0,0]);
    
    Rate[7]= gammah*V[7] 
    Change[7, :]= ([0,0,0,0,0,0,0,-1,+1,0,0,0,0]);
    
    Rate[8]= Ca*Pa*V[9]*((V[1]+V[2])/N)
    Change[8, :]= ([0,0,0,0,0,0,0,0,0,-1,+1,0,0]);
    
    Rate[9]= kappa*V[10]
    Change[9, :]= ([0,0,0,0,0,0,0,0,0,0,-1,+1,0]);
    
    Rate[10]= gammah*V[11]
    Change[10, :]= ([0,0,0,0,0,0,0,0,0,0,0,-1,+1]);
    
    for i in range(11):
        Num=np.random.poisson(Rate[i]*tau);
        Use=min([Num, V[pl.where(Change[i,:]<0)]]);
        V=V+Change[i,:]*Use;
    return V

### SETTING UP THE ITERATIONS ###
def Stoch_Iteration(INPUT):
    lop = 0
    Ss, As, Is, Rs, Qs = [INPUT[0]], [INPUT[1]], [INPUT[2]], [INPUT[3]], [INPUT[4]]
    Sm, Em, Im, Rm = [INPUT[5]], [INPUT[6]], [INPUT[7]], [INPUT[8]]
    Sa, Ea, Ia, Ra = [INPUT[9]], [INPUT[10]], [INPUT[11]], [INPUT[12]]
    
    swine_cumulative_infections = 0  # Initialize cumulative infection counter
    member_cumulative_infections = 0 # Initialize cumulative infection counter
    attendee_cumulative_infections = 0 # Initialize cumulative infection counter
    
    for lop in T:
        res = stoch_eqs(INPUT)
        
        # Count infections as those leaving Sm for As or Is
        swine_cumulative_infections += (INPUT[0] - res[0])  # Count pigs leaving Ss
        member_cumulative_infections += (INPUT[5] - res[5]) # Count members leaving Sm
        attendee_cumulative_infections += (INPUT[9] - res[9]) #Count attendees leaving Sa
        
        # Update the compartment lists
        Ss.append(res[0])
        As.append(res[1])
        Is.append(res[2])
        Rs.append(res[3])
        Qs.append(res[4])
        Sm.append(res[5])
        Em.append(res[6])
        Im.append(res[7])
        Rm.append(res[8])
        Sa.append(res[9])
        Ea.append(res[10])
        Ia.append(res[11])
        Ra.append(res[12])
        
        INPUT = res  # Update INPUT for the next iteration
    
    return [Ss, As, Is, Rs, Qs, Sm, Em, Im, Rm, Sa, Ea, Ia, Ra, 
            swine_cumulative_infections, member_cumulative_infections, attendee_cumulative_infections]


T=np.arange(0.0, ND, tau)
        
[Ss, As, Is, Rs, Qs, Sm, Em, Im, Rm, Sa, Ea, Ia, Ra, 
 swine_cumulative_infections, member_cumulative_infections, attendee_cumulative_infections]= Stoch_Iteration(INPUT)  

# making lists
tT=np.array(T)
tSs=np.array(Ss)[1:]
tAs=np.array(As)[1:]
tIs=np.array(Is)[1:]
tRs=np.array(Rs)[1:]
tQs=np.array(Qs)[1:]
tSm=np.array(Sm)[1:]
tEm=np.array(Em)[1:]
tIm=np.array(Im)[1:]
tRm=np.array(Rm)[1:]
tSa=np.array(Sa)[1:]
tEa=np.array(Ea)[1:]
tIa=np.array(Ia)[1:]
tRa=np.array(Ra)[1:]


### MULTIPLE SIMULATIONS ###
num_iterations = 10000

all_results= [] # store all results
swine_cumulative_infections_list= [] #store swine results
member_cumulative_infections_list= [] # store member results
attendee_cumulative_infections_list= [] # store attendee results

for iteration in range(num_iterations):
    INPUT = np.array((Ss0, As0, Is0, Rs0, Qs0, Sm0, Em0, Im0, Rm0, Sa0, Ea0, Ia0, Ra0))
    
    results = Stoch_Iteration(INPUT) # run sim for this iteration
    
    swine_cumulative_infections = results[-3]  # Extract cumulative infections for swine
    member_cumulative_infections = results[-2]
    attendee_cumulative_infections = results[-1]
    
    all_results.append(results[:-3])  # Append the rest of the results
    swine_cumulative_infections_list.append(swine_cumulative_infections)
    member_cumulative_infections_list.append(member_cumulative_infections)
    attendee_cumulative_infections_list.append(attendee_cumulative_infections)

    

### VIEWING RESULTS ###
iteration_1_results = all_results[0]



### SAVING RESULTS ###
import pandas as pd

## SIMULATIONS ##
# Create a DataFrame for the  simulation results
columns = ['Ss', 'As', 'Is', 'Rs', 'Qs', 'Sm', 'Em', 'Im', 'Rm', 'Sa', 'Ea', 'Ia', 'Ra']
sim_res = pd.DataFrame({col: np.concatenate([res[idx] for res in all_results]) for idx, col in enumerate(columns)})

# Create a time column for the range of days (corresponding to T)
time_column = np.tile(np.arange(0, len(T)+1), num_iterations)  # Repeat T for each iteration

# Add the time column to the DataFrame
sim_res['Time'] = time_column

# Add a simulation iteration identifier
simulation_ids = np.repeat(range(1, num_iterations + 1), len(T)+1)
sim_res['Simulation'] = simulation_ids

# Save the DataFrame 
sim_res.to_csv(r'C:\Users\dcpit\OneDrive\Documents\GRADUATE SCHOOL\PhD-yr1\Project 4\Results\Quarentine\Theta_1\Symptomatic_17\Model_results_with_time_symptom_ID_0.2_17.csv', index=False)

# Display the first few rows of the DataFrame
print(sim_res.head())

## CUMULATIVE ##
#Combing the lists into a single dataframe
cum_res = pd.DataFrame({'Swine':swine_cumulative_infections_list,
                       'Member': member_cumulative_infections_list,
                       'Attendee': attendee_cumulative_infections_list})
# Save the DataFrame 
cum_res.to_csv(r'C:\Users\dcpit\OneDrive\Documents\GRADUATE SCHOOL\PhD-yr1\Project 4\Results\Quarentine\Theta_1\Symptomatic_17\Cumulative_infection_symptom_ID_0.2_17.csv', index=False)
