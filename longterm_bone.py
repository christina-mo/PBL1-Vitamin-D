from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

t_span = np.linspace(0,20, 1000)

#Need to manually fill in values from main code
Ca_healthy = #amount of Calcium entering from intestine (Calcium absorbed in intestine)
Ca_diseased = #amount of Calcium entering the intestine if CKD
y0 = 1 #1 * 10^6 mg of Calcium stored in bones for a healthy individual

#Ca_Abs amount of Calcium in mg that enters the bone/day
#kidney_regression = the rate of decrease in GFR in ml/min each year
def Ca_In_Bones(y0, t, Ca_Abs, kidney_regression):
  if t < 10 and kidney_regression == 1:
    dy_dt = ((500-287.1)*365 + (Ca_Abs*365) - 500*365)/1E6
  else:
    dy_dt = ((500-287.1)*365 + Ca_Abs*365*((100-kidney_regression*t)/100) - 500*365)/1E6
  return dy_dt

y_healthy = odeint(Ca_In_Bones, y0, t_span, args=(Ca_healthy,1))
y_diseased = odeint(Ca_In_Bones, y0, t_span, args=(Ca_diseased,2.5)) #args = (Ca_Abs = Ca_diseased, kidney_regression = 2.5)

index = np.argmin(np.abs(t_span - 20))
print(f"Diseased after 10 years: {y_diseased[index]}")

plt.figure(figsize = (10,6))
plt.plot(t_span, y_healthy, label = 'Healthy')
plt.plot(t_span, y_diseased, label = 'CKD')
plt.xlabel("Years After 30", fontsize = 18)
plt.xticks(range(0,21, 2))
plt.yticks(np.arange(0,1.1, 0.1))
plt.ylabel("Calcium in Bones (kg)", fontsize = 18)
plt.title("Calcium Bone Stores in Healthy vs. CKD Individuals Over Time", fontsize = 18)
plt.tick_params(length = 10, labelsize =12)
plt.axvline(10, label = "Age 40", linestyle=':', color='g')
plt.legend(fontsize = 18)
plt.show()
