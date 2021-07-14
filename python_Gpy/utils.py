import scipy.io
import numpy as np
import matplotlib.pyplot as plt

def get_running_time(dataset,nRep):
    # Provide a rough estimation of running time
    if dataset=='nhp':
        exp_time=np.ceil(22*4.3*nRep*np.array([10, 16])/60/60)
    elif dataset=='rat':
        exp_time=np.ceil(40*0.9*nRep*np.array([10, 16])/60/60)
    print('The procedure can take approximately '+ str(exp_time[0])+ ' to '+str(exp_time[1])+ ' hours on a standard workstation (Intel i7-8700 @ 3.2GHz)')
    
def load_matlab_data(path_to_dataset,dataset,m_i):
    # Load data for each subject separately
    # The order of the variables inside the dict changes between Macaque, Cebus and rats
    if dataset=='nhp': # nhp dataset has 4 subjects
        if m_i==0:
            Cebus1_M1_190221 = scipy.io.loadmat(path_to_dataset+'/Cebus1_M1_190221.mat')
            Cebus1_M1_190221= {'emgs': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][0][0],
           'nChan': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][2][0][0],
           'sorted_isvalid': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][8],
           'sorted_resp': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][9],
           'sorted_respMean': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][10],
           'ch2xy': Cebus1_M1_190221['Cebus1_M1_190221'][0][0][16]}
            SET=Cebus1_M1_190221
        if m_i==1:
            Cebus2_M1_200123 = scipy.io.loadmat(path_to_dataset+'/Cebus2_M1_200123.mat')  
            Cebus2_M1_200123= {'emgs': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][0][0],
           'nChan': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][2][0][0],
           'sorted_isvalid': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][8],
           'sorted_resp': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][9],
           'sorted_respMean': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][10],
           'ch2xy': Cebus2_M1_200123['Cebus2_M1_200123'][0][0][16]}
            SET=Cebus2_M1_200123
        if m_i==2:    
            Macaque1_M1_181212 = scipy.io.loadmat(path_to_dataset+'/Macaque1_M1_181212.mat')
            Macaque1_M1_181212= {'emgs': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][0][0],
           'nChan': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][2][0][0],
           'sorted_isvalid': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][8],
           'sorted_resp': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][9],              
           'sorted_respMean': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][15],
           'ch2xy': Macaque1_M1_181212['Macaque1_M1_181212'][0][0][14]}            
            SET=Macaque1_M1_181212
        if m_i==3:    
            Macaque2_M1_190527 = scipy.io.loadmat(path_to_dataset+'/Macaque2_M1_190527.mat')
            Macaque2_M1_190527= {'emgs': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][0][0],
           'nChan': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][2][0][0],
           'sorted_isvalid': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][8],
           'sorted_resp': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][9],              
           'sorted_respMean': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][15],
           'ch2xy': Macaque2_M1_190527['Macaque2_M1_190527'][0][0][14]}
            SET=Macaque2_M1_190527
    elif dataset=='rat':  # rat dataset has 6 subjects
        if m_i==0:
            rat1_M1_190716 = scipy.io.loadmat(path_to_dataset+'/rat1_M1_190716.mat')
            rat1_M1_190716= {'emgs': rat1_M1_190716['rat1_M1_190716'][0][0][0][0],
           'nChan': rat1_M1_190716['rat1_M1_190716'][0][0][2][0][0],
           'sorted_isvalid': rat1_M1_190716['rat1_M1_190716'][0][0][8],
           'sorted_resp': rat1_M1_190716['rat1_M1_190716'][0][0][9],              
           'sorted_respMean': rat1_M1_190716['rat1_M1_190716'][0][0][15],
           'ch2xy': rat1_M1_190716['rat1_M1_190716'][0][0][14]}            
            SET=rat1_M1_190716
        if m_i==1:
            rat2_M1_190617 = scipy.io.loadmat(path_to_dataset+'/rat2_M1_190617.mat')
            rat2_M1_190617= {'emgs': rat2_M1_190617['rat2_M1_190617'][0][0][0][0],
           'nChan': rat2_M1_190617['rat2_M1_190617'][0][0][2][0][0],
           'sorted_isvalid': rat2_M1_190617['rat2_M1_190617'][0][0][8],
           'sorted_resp': rat2_M1_190617['rat2_M1_190617'][0][0][9],              
           'sorted_respMean': rat2_M1_190617['rat2_M1_190617'][0][0][15],
           'ch2xy': rat2_M1_190617['rat2_M1_190617'][0][0][14]}         
            SET=rat2_M1_190617          
        if m_i==2:
            rat3_M1_190728 = scipy.io.loadmat(path_to_dataset+'/rat3_M1_190728.mat')
            rat3_M1_190728= {'emgs': rat3_M1_190728['rat3_M1_190728'][0][0][0][0],
           'nChan': rat3_M1_190728['rat3_M1_190728'][0][0][2][0][0],
           'sorted_isvalid': rat3_M1_190728['rat3_M1_190728'][0][0][8],
           'sorted_resp': rat3_M1_190728['rat3_M1_190728'][0][0][9],              
           'sorted_respMean': rat3_M1_190728['rat3_M1_190728'][0][0][15],
           'ch2xy': rat3_M1_190728['rat3_M1_190728'][0][0][14]}           
            SET=rat3_M1_190728                       
        if m_i==3:
            rat4_M1_191109 = scipy.io.loadmat(path_to_dataset+'/rat4_M1_191109.mat')
            rat4_M1_191109= {'emgs': rat4_M1_191109['rat4_M1_191109'][0][0][0][0],
           'nChan': rat4_M1_191109['rat4_M1_191109'][0][0][2][0][0],
           'sorted_isvalid': rat4_M1_191109['rat4_M1_191109'][0][0][8],
           'sorted_resp': rat4_M1_191109['rat4_M1_191109'][0][0][9],              
           'sorted_respMean': rat4_M1_191109['rat4_M1_191109'][0][0][15],
           'ch2xy': rat4_M1_191109['rat4_M1_191109'][0][0][14]}            
            SET=rat4_M1_191109                       
        if m_i==4:
            rat5_M1_191112 = scipy.io.loadmat(path_to_dataset+'/rat5_M1_191112.mat')
            rat5_M1_191112= {'emgs': rat5_M1_191112['rat5_M1_191112'][0][0][0][0],
           'nChan': rat5_M1_191112['rat5_M1_191112'][0][0][2][0][0],
           'sorted_isvalid': rat5_M1_191112['rat5_M1_191112'][0][0][8],
           'sorted_resp': rat5_M1_191112['rat5_M1_191112'][0][0][9],              
           'sorted_respMean': rat5_M1_191112['rat5_M1_191112'][0][0][15],
           'ch2xy': rat5_M1_191112['rat5_M1_191112'][0][0][14]}           
            SET=rat5_M1_191112                      
        if m_i==5:
            rat6_M1_200218 = scipy.io.loadmat(path_to_dataset+'/rat6_M1_200218.mat')        
            rat6_M1_200218= {'emgs': rat6_M1_200218['rat6_M1_200218'][0][0][0][0],
           'nChan': rat6_M1_200218['rat6_M1_200218'][0][0][2][0][0],
           'sorted_isvalid': rat6_M1_200218['rat6_M1_200218'][0][0][8],
           'sorted_resp': rat6_M1_200218['rat6_M1_200218'][0][0][9],              
           'sorted_respMean': rat6_M1_200218['rat6_M1_200218'][0][0][15],
           'ch2xy': rat6_M1_200218['rat6_M1_200218'][0][0][14]}          
            SET=rat6_M1_200218               
    else:
        print('Invalid value for dataset variable. Has to be either \'nhp\' or \'rat\'. ')     
        SET=None        
    return SET
        
def param_grid(which_opt, dataset, verbose=True):
    # A selection of values to be tested for the selected hyperparameter
    if which_opt=='nrnd':
        if dataset=='nhp':
            this_opt= np.array([1, 2, 3, 5, 10, 20, 30, 50, 90, 96])
        else:
            this_opt= np.array([1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30, 32])     
    elif which_opt=='rho_low':
        this_opt= np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])     
    elif which_opt=='rho_high':
        this_opt= np.arange(10)+1     
    elif which_opt=='noise_min':
        this_opt= np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5])    
    elif which_opt=='noise_max':
        this_opt= np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10])    
    elif which_opt=='kappa':
        if dataset=='nhp':
            this_opt= np.array([1, 2, 3, 3.5, 3.8, 4.1, 4.4, 4.7, 5, 5.5, 6.5, 7, 8, 9, 10])
        else:
            this_opt= np.array([1, 1.5, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 5, 6, 7, 8, 9, 10])              
    else:
        print('Invalid value for to_optimize. Has to be \'nrnd\', \'rho_low\', \'rho_high\', \'min_noise\', \'max_noise\' or \'kappa\'.')
        this_opt=None   
    if verbose==True:    
        print('Values of ' + which_opt + ' to be tested: '+str(this_opt))   
    return this_opt

def plot_optim_results(Stored_perf_explore,Stored_perf_exploit,which_opt, this_opt, dataset):
    # show performance graphs
    emgs=[]
    for m_i in range(Stored_perf_explore.shape[0]): # total number of subjects
        count=0
        for e_i in range(Stored_perf_explore.shape[1]): # max number of emgs
            if np.all((Stored_perf_explore[m_i,e_i] != 0.0)):
                count+=1
        emgs.append(count)
    n= np.sum(emgs)
           
    FinalMeanPerfExploration= np.zeros((len(this_opt),n))
    FinalMeanPerfExploitation= np.zeros((len(this_opt),n))
    for k_i in range(len(this_opt)):
        jj=0
        for m_i in range(len(emgs)):
            for syn in range(emgs[m_i]):
                if len(this_opt)==1:
                    MeanPerfExploration= np.mean(Stored_perf_explore[m_i,syn], axis=0)
                    MeanPerfExploitation= np.mean(Stored_perf_exploit[m_i,syn], axis=0)
                else:
                    MeanPerfExploration= np.mean(Stored_perf_explore[m_i,syn,k_i], axis=0)
                    MeanPerfExploitation= np.mean(Stored_perf_exploit[m_i,syn,k_i], axis=0)
                FinalMeanPerfExploration[k_i,jj]= MeanPerfExploration[-1] # we will display final performance
                FinalMeanPerfExploitation[k_i,jj]= MeanPerfExploitation[-1]
                jj=jj+1 # replicates are individual muscles
    if which_opt=='rho_low' or which_opt=='noise_min' or which_opt=='noise_max':
        this_opt=np.log10(this_opt) # using log scale
    # using linear scale for all other hyperparameters
    plt.plot(this_opt,np.mean(FinalMeanPerfExploration, axis=1), '-ob', alpha=0.9, label= 'Exploration (knowledge of best channel)')
    plt.plot(this_opt,np.mean(FinalMeanPerfExploitation, axis=1), '-ok', alpha= 0.9, label='Exploitation (stimulation efficacy)')
    plt.fill_between(this_opt,np.mean(FinalMeanPerfExploration, axis=1) - (np.std(FinalMeanPerfExploration, axis=1)/np.sqrt(FinalMeanPerfExploration.shape[1])),
                     np.mean(FinalMeanPerfExploration, axis=1) + (np.std(FinalMeanPerfExploration, axis=1)/np.sqrt(FinalMeanPerfExploration.shape[1])),
                     color='blue', alpha=0.2)
    plt.fill_between(this_opt,np.mean(FinalMeanPerfExploitation, axis=1) - (np.std(FinalMeanPerfExploitation, axis=1)/np.sqrt(FinalMeanPerfExploitation.shape[1])),
                     np.mean(FinalMeanPerfExploitation, axis=1) + (np.std(FinalMeanPerfExploitation, axis=1)/np.sqrt(FinalMeanPerfExploitation.shape[1])),
                     color='black', alpha=0.2)
    print('Exploration: '+ str(np.mean(FinalMeanPerfExploration, axis=1)))
    print('Exploitation: ' + str(np.mean(FinalMeanPerfExploitation, axis=1)))
    plt.ylim([0,1])
    plt.ylabel('Performance')
    if which_opt=='rho_low' or which_opt=='noise_min' or which_opt=='noise_max': # using log scale
        plt.xlabel('Hyperparameter value (in log scale)') 
    else:
        plt.xlabel('Hyperparameter value')
    plt.title('Algorithmic performance for multiple values of '+ which_opt)
    plt.legend(loc='lower center')
