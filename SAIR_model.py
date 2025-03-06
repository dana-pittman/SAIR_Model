
# SAIR model with Quarentine developed 11/21/24
# Dana C. Pittman Ratterree


# Required libraries 
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint

### Parameters ###
##Movement between compartments
#PIG
R0 = 4 #The reproductive number of influenza in pigs
gammas = 1/5 # recovery rate of swine
beta = R0 * gammas # transmission rate
delta = 0.83 # probability of becoming asymptmatic
chi = 0.5 # sucess at identifying symptomatic pigs
theta = 1/2 #delay at identifying symptomatic pigs

#HUMAN
Cm = 60 # contact duration per minute members
Pm = 0.00168 # probability of transmission members
Ca = 5 # contact duration per minure attendees
Pa = 0.00095 # probability of transmission attendees
kappa = 1/2 # rate exposed become infected
gammah = 1/5 # recovery rate

##Defining intial conditions
#PIG
Ss0, As0, Is0, Rs0, Qs0 = 203, 5, 0, 0, 0

N=Ss0+As0+Is0+Rs0

#MEMBER
Sm0, Em0, Im0, Rm0 = 90, 0, 0, 10

# ATTENDEE
Sa0, Ea0, Ia0, Ra0 = 10042, 0, 0, 4868

### Equations ###
def Model (X, t):
    Ss, As, Is, Rs, Qs, Sm, Em, Im, Rm, Sa, Ea, Ia, Ra, Cpigs, Cmembers, Cattendees = X
    
    #pig equaitons
    dSsdt = -beta *Ss * (Is + As)/N
    dAsdt = beta *Ss * (Is + As)/N * delta - gammas * As
    dIsdt = beta *Ss * (Is + As)/N *(1- delta) - chi* theta * Is -(1 - chi) * gammas *Is
    dRsdt = gammas * As + (1 - chi)  * gammas * Is
    dQsdt = chi * theta * Is
    
    #Member equations
    dSmdt = -Cm * Pm * Sm * (Is + As)/N
    dEmdt = Cm * Pm * Sm * (Is + As)/N - kappa * Em
    dImdt = kappa * Em - gammah * Im
    dRmdt = gammah * Im
    
    #Attendee equations
    dSadt = -Ca * Pa * Sa * (Is + As)/N
    dEadt = Ca * Pa * Sa * (Is + As)/N - kappa * Ea
    dIadt = kappa * Ea - gammah * Ia
    dRadt = gammah * Ia
    
    # Cumulative infection counters
    dCpigsdt = -dSsdt  # New infections in pigs
    dCmembersdt = -dSmdt  # New infections in members
    dCattendeesdt = -dSadt  # New infections in attendees
    
    return np.array([dSsdt, dAsdt, dIsdt, dRsdt, dQsdt, dSmdt, dEmdt, dImdt, dRmdt, dSadt, dEadt, dIadt, dRadt, dCpigsdt, dCmembersdt, dCattendeesdt])

# defining inital condition vector
y0 = [Ss0, As0, Is0, Rs0, Qs0, Sm0, Em0, Im0, Rm0, Sa0, Ea0, Ia0, Ra0, 0, 0, 0]

# time vector
t= np.linspace(0, 9, 9) #simulates 20 days

# solving the model
solution = odeint(Model, y0, t)

# extract results
Ss, As, Is, Rs, Qs, Sm, Em, Im, Rm, Sa, Ea, Ia, Ra, Cpigs, Cmembers, Cattendees = solution.T 


### Plots ###
#Pig population dynamic plot
plt.Figure(figsize=(10,6))

plt.plot(t, Ss, label = "Susceptible Swine", color = "blue")
plt.plot(t, Is, label = "Symptomatic Swine", color = "magenta")
plt.plot(t, As, label = "Asymptomatic Swine", color = "red")
plt.plot(t, Rs, label = "Recovered Swine", color = "black")
plt.plot(t, Qs, label = "Quarentined Swine", color = "green")

plt.xlabel("Time (days)")
plt.ylabel("Count")
plt.legend()
plt.title('Swine population dynamics')
plt.show()

#Member population dynamic plot
plt.Figure(figsize=(10,6))

plt.plot(t, Sm, label = "Susceptible Members", color = "blue")
plt.plot(t, Em, label = "Exposed Members", color = "green")
plt.plot(t, Im, label = "Infected Members", color = "red")
plt.plot(t, Rm, label = "Recovered Members", color = "black")

plt.xlabel("Time (days)")
plt.ylabel("Count")
plt.legend()
plt.title('Member Infection dynamics')
plt.show()

#Attendee population dynamic plot
plt.Figure(figsize=(10,6))

plt.plot(t, Sa, label = "Susceptible Attendees", color = "blue")
plt.plot(t, Ea, label = "Expsoed Attendee", color = "green")
plt.plot(t, Ia, label ="Infected Attendee", color = "red")
plt.plot(t, Ra, label = "Recovered Attendees", color = "black")

plt.xlabel("Time (days)")
plt.ylabel("Count")
plt.legend()
plt.title('Attendee Infection dynamics')
plt.show()

# Cumulative infection plot
plt.plot(t, Cpigs, label="Cumulative Infections (Pigs)")
plt.plot(t, Cmembers, label="Cumulative Infections (Members)")
plt.plot(t, Cattendees, label="Cumulative Infections (Attendees)")

plt.xlabel("Time (days)")
plt.ylabel("Cumulative Infections")
plt.title("Cumulative Infections Over Time by Population")
plt.legend()
plt.show()
