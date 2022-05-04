
#main program
import numpy as np
import scipy.integrate as integ
import matplotlib as plt
from cmaes import CMA
import optuna

# units in centimeters
########sharpe###########
Ea3 = np.array([7.286722469170652, 10.303713848332231, 14.272126621536309, 18.934320956368634, 23.893216946469305])
E_l_sharpe  = np.array([2,4,6,8,10])
gam_sharpe=1.54
q_hat_sharpe=11.5
##########################

####high gamma###
W16_l_high = np.array([1.3145005934249303, 2.2533265176600508, 3.1473177700386983, 4.110565750507097, 5.1971036254122955, 6.420863583536697])/0.175
E_l  = np.array([2,4,6,8,10,12])
gam_high_gam=1.5
q_hat_high_gam=26.333
#################

####high q###
W16_l_high_q=np.array([2.2283549929063122, 4.436693421408165, 7.072866967085392, 10.276651460191752, 14.12380248850615, 18.569446147381107])/0.175
E_l_high_q = np.array([2,4,6,8,10,12])
gam_high_q=1.333
q_hat_high_q=50
#####################

#########lu et al.###########
W16_l = np.array([2.0099448159652655, 3.8616304616834367, 6.0017247326066645, 8.577999547773649, 11.664230325132511, 15.216815094569808, 19.182325129368397, 23.447555070072138])/0.175
E_l_lu  = np.array([2,4,6,8,10,12,14,16])
gam_lu=1.333
q_hat_lu=26.333
############################

##########q_up#############
W16_l_q_up = np.array([2.0079095789135577, 3.8971280594764726, 6.1017728461465195, 8.771413803600637, 11.975003922447678, 15.682001474450658, 19.811107257088878, 24.248338913680286])/0.175
E_l_q_up  = np.array([2,4,6,8,10,12,14,16])
gam_q_up=1.333
q_hat_q_up=30
###############
####################diffrent values of eff_activation energy##########################3


def objective(trial: optuna.Trial):
    x0 = trial.suggest_float("x0", 0, 5,log=True)
    x1 = trial.suggest_float("x1", 0,5,log=True) 
    x2 = trial.suggest_float("x2", 0,5,log=True)
    
    return 

if __name__ == "__main__":
    sampler = optuna.samplers.CmaEsSampler()#CmaEsSampler(consider_pruned_trials=True)
    study = optuna.create_study(sampler=sampler)
    study.optimize(objective, n_trials=200)
    #fig = optuna.visualization.plot_optimization_history(study)
    #fig.show()
    #print('eff{} :[{},{}], Ea = {}'.format(val,study.best_params,study.best_value,  21.573723742459983*(1/2)*val))
    
    #print('\n')
