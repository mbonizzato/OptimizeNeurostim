%Fixed-lengthscales (rho) implementation of GP-BO optimization of neurostimulation
%The code below displays algorithmic performance on provided datasets
%For all values of a selected hyperparameter.

%Hyperparameter selection is crucial for GP-BO applications. We recommend
%running this code on own existing or surrogate data to tune at least the
%UCB acquisition function hyperparameter "k" (kappa).

% Select moodality
dataset='nhp';              %selected dataset
which_opt= 'kappa';     %hyperparameter to optimize
nRep=50;                    %number of repetitions

% Load data
if strcmp(dataset,'nhp')
    %nhp dataset has 4 subjects
    load('Macaque1_M1_181212.mat')
    load('Macaque2_M1_190527.mat')
    load('Cebus1_M1_190221.mat')
    load('Cebus2_M1_200123.mat')
    SETS=[Macaque1_M1_181212 Macaque2_M1_190527 Cebus1_M1_190221 Cebus2_M1_200123];
    % Provide a rough estimation of running time
    exp_time=ceil(20*6/50*nRep*[5 10]/60);
    
elseif strcmp(dataset,'rat')  
    %rat dataset has 6 subjects  
    load('rat1_M1_190716.mat')
    load('rat2_M1_190617.mat')
    load('rat3_M1_190728.mat')
    load('rat4_M1_191109.mat')
    load('rat5_M1_191112.mat')
    load('rat6_M1_200218.mat')
    SETS=[rat1_M1_190716 rat2_M1_190617 rat3_M1_190728 rat4_M1_191109 rat5_M1_191112 rat6_M1_200218];
    % Provide a rough estimation of running time
    exp_time=ceil(27*2/50*nRep*[3 5]/60);
end

if exp_time(2)==1
    disp(['The procedure should take less than 1 minute on a standard workstation (Intel i7-6700 @ 3.4GHz).'])
else    
    disp(['The procedure can take approximately ' num2str(exp_time(1)) ' to ' num2str(exp_time(2)) ' minutes on a standard workstation (Intel i7-6700 @ 3.4GHz).'])
end

%rho is the kernel geometrical hyperparameter (lengthscales)
%kappa is the UCB acquisition function hyperparameter
%nrand is the number of random queries performed to initialize the GP-BO
%noisemax is the maximum value for the noise hyperparameter
if strcmp(dataset,'nhp') 
    kappa=4;   
    rho=3; 
    nrnd=1;
    noise_estim= 0.001;         
elseif strcmp(dataset,'rat') 
    kappa=3;   
    rho=5; 
    nrnd=1;
    noise_estim= 0.001;         
end
%kappa and rho need to vary between the nhp and rat implementation

%A selection of values to be tested for the selected hyperparameter
if strcmp(which_opt,'rho')
    this_opt = [1:8]; % rho
elseif strcmp(which_opt,'noise_estim')
    this_opt=[0.001 0.005 0.01 0.05 0.1 0.5 1 5 10]; % noise
elseif strcmp(which_opt,'kappa')
    this_opt=[1:.3:5.5 6:10]; % kappa
elseif strcmp(which_opt,'nrand')
    if strcmp(dataset,'nhp')
        this_opt= [1:5 10:10:90 95 96]; %nrand
    else
        this_opt= [1:5 7 10:5:30 32];
    end
end
    
%prepare results storage
clear PP PT YMU PP_t msr
msr = cell(size(SETS,2),8,numel(this_opt),nRep,SETS(1).nChan); 
PP = cell(numel(SETS),8);
PP_t = cell(numel(SETS),8);

%update timing estimation utilities
tic;
count_perf=0;
tot_emgs=0;
for m_i=1:numel(SETS) 
    tot_emgs=tot_emgs+ length(SETS(m_i).emgs);
end
tot_perf=tot_emgs*length(this_opt);

for m_i=1:numel(SETS)                         %for each subject
    
    subject= SETS(m_i);    
    
    for k_i=1:numel(this_opt)       %for each hyperparameter value
    
        if strcmp(which_opt,'rho')
            rho=this_opt(k_i);
        elseif strcmp(which_opt,'noise_estim')
            noise_estim=this_opt(k_i);
        elseif strcmp(which_opt,'kappa')
            kappa=this_opt(k_i);
        elseif strcmp(which_opt,'nrand')
            nrnd= this_opt(k_i);
        end
        
        % Create the kernel
        DimSearchSpace=subject.nChan;
        dist = zeros(DimSearchSpace);
        distKern = dist;
        for x = 1:DimSearchSpace
            for y = 1:DimSearchSpace         
                dif1 = (subject.ch2xy(x,1)-subject.ch2xy(y,1))^2;
                dif2 = (subject.ch2xy(x,2)-subject.ch2xy(y,2))^2;
                dist(x,y) = sqrt(dif1 + dif2);        
                distKern(x,y) = (1+(sqrt(5)*dist(x,y))/rho + (5*(dist(x,y))^2)/(3*rho^2))*exp(-sqrt(5)*dist(x,y)/rho);
            end
        end

        for e_i=1:numel(subject.emgs) % for each muscle of the given subject

            %display remaining time information
            disp([num2str(count_perf/tot_perf*100) ' % completed'])
            if count_perf>0
                t=toc;                
                hs=floor((t/count_perf)*(tot_perf-count_perf)/60/60);
                mins=floor((t/count_perf)*(tot_perf-count_perf)/60-hs*60);
                disp(['Estimated remaining time: ' num2str(hs) ' hours, ' num2str(mins) ' minutes.'])
            end

            % "Ground truth" map
            MPm=subject.sorted_respMean(:,e_i);
            % Best known channel
            mMPm=max(MPm);

            % Then run the sequential optimization
            perf_explore=[]; %performance
            MaxQueries = DimSearchSpace;
            % Here the maximum number of queries is equal to the maximum
            % number of electrodes, which is the number of queries required
            % for extensive search
            P_test = cell(nRep,1); %storing all queries 
            clear perf_explore perf_exploit perf_rsq

            for rep_i=1:nRep % for each repetition

                MaxSeenResp=0; %maximum response obtained in this round,
                %used to normalize all responses between zero and one.
                q=1; % query number
                clear valid_resp
                P_max=[];   

                order_this=randperm(subject.nChan); %random permutation
                %of each entry of the search space

                while q <= MaxQueries  
                    %We will sample the search space randomly for
                    %exactly nrnd queries
                    if q>nrnd
                        % Find next point (max of acquisition function)
                        kappa_norm=abs(MaxSeenResp)*kappa;
                        %note that in this implementation we chose to
                        %normalize kappa, rather than the collected data.
                        AcquisitionMap = MapPrediction + kappa_norm.*real(sqrt(VarianceMap));
                        NextQuery = find(ismember(AcquisitionMap, max(AcquisitionMap))); 
                        %select next query
                        if length(NextQuery) > 1
                            NextQuery = NextQuery(randi(numel(NextQuery)));
                        end
                        P_test{rep_i}(q,1) = NextQuery; 
                    else 
                        P_test{rep_i}(q,1) = order_this(q); 
                        K_maj=[];
                    end
                    query_elec = P_test{rep_i}(q,1);

                    %This offline optimization code randomly choses one
                    %response among all responses stored in the
                    %selected search space look-up table.
                    valid_resp=subject.sorted_resp{query_elec,e_i}(subject.sorted_isvalid{query_elec,e_i}~=0);
                    r_i=randi(numel(valid_resp));  
                    test_respo= valid_resp(r_i);
                    % done reading response
                    P_test{rep_i}(q,2)=test_respo;
                    %The first element of P_test is the selected search
                    %space point, the second the resulting value

                    if (test_respo>MaxSeenResp) || (MaxSeenResp==0)
                        %updated maximum response obtained in this
                        %round
                        MaxSeenResp=test_respo;
                    end

                    %GP-BO optimization
                    [MapPrediction,VarianceMap,K_maj] = CalcPrediction(P_test{rep_i},noise_estim,distKern,K_maj,DimSearchSpace);  

                    % We only test for gp predictions at electrodes that
                    % we had queried (presumable we only want to return an
                    % electrode that we have already queried).                        
                    Tested=unique(sort(P_test{rep_i}(:,1)));
                    MapPredictionTested=MapPrediction(Tested);
                    BestQuery=Tested(find(ismember(MapPredictionTested, max(MapPredictionTested))));
                    if length(BestQuery) > 1
                        BestQuery = BestQuery(randi(numel(BestQuery)));
                    end

                    %Maximum response at time q
                    P_max(q)= BestQuery; 
                    %store all info
                    msr{m_i,e_i,k_i,rep_i,q} = MaxSeenResp;
                    YMU{m_i,e_i,k_i,rep_i,q}=MapPrediction;
                    q=q+1; 
                end

                %estimate current exploration performance: 
                %knowledge of best stimulation point
                perf_explore(rep_i,:)=MPm(P_max)/mMPm; 
                %estimate current exploitation performance: 
                %knowledge of best stimulation point
                perf_exploit(rep_i,:)= P_test{rep_i,1}(:,1);
                %calculate model fitting of ground truth value map
                mdl=corrcoef(MPm,MapPrediction);        
                perf_rsq(rep_i,:)=mdl(2)^2;

                %store all tests 
                PT{m_i,e_i,k_i,rep_i}=P_test{rep_i};
            end

            %store all performance estimations
            RSQ{m_i,e_i,k_i}=perf_rsq;
            PP{m_i,e_i,k_i}=perf_explore;
            PP_t{m_i,e_i,k_i}= MPm(perf_exploit)/mMPm;
            
            count_perf=count_perf+1;
        end
    end
end

%save all workspace, excluding datasets
dd = datetime('today');
dd1=num2str(year(dd));
dd1=dd1(end-1:end);
dd2=month(dd);
if dd2<10
    dd2=['0' num2str(dd2)];
else
    dd2=[num2str(dd2)];
end
dd3=day(dd);
if dd3<10
    dd3=['0' num2str(dd3)];
else
    dd3=[num2str(dd3)];
end
ddate=[dd1 dd2 dd3];
fn=['basicversion_' dataset '_' which_opt '_' num2str(nRep) '_' ddate]; %file name
%
save(fn,'dataset','kappa','MaxQueries','msr','noise_estim','nRep',...
    'nrnd','PP','PP_t','PT','rho','RSQ','this_opt','which_opt','YMU')


%
%show performance graphs
clear perft perftt
for k_i=1:numel(this_opt)
    jj=0;
    for m_i=1:size(PP,1)
        for e_i=1:size(PP,2)
            jj=jj+1; %replicates are individual muscles
            ppm=mean(PP{m_i,e_i,k_i});
            ppt=mean(PP_t{m_i,e_i,k_i});
            perft(k_i,jj)=[ppm(end)]; %we will display final performance
            perftt(k_i,jj)=[ppt(end)];
        end
    end
end

perft=perft(:,~isnan(perft(1,:)));
perftt=perftt(:,~isnan(perftt(1,:)));


if strcmp(which_opt,'rholow') || strcmp(which_opt,'noisemax')

    %using log scale
    figure
    semilogx(this_opt,mean(perft'),'b')
    hold on
    semilogx(this_opt,mean(perftt'),'k')
    semilogx(this_opt,mean(perft')+ std(perft')/sqrt(size(perft,2)),'b')
    semilogx(this_opt,mean(perft')-std(perft')/sqrt(size(perft,2)),'b')
    semilogx(this_opt,mean(perftt')+ std(perftt')/sqrt(size(perftt,2)),'k')
    semilogx(this_opt,mean(perftt')-std(perftt')/sqrt(size(perftt,2)),'k')
    ylim([0 1])
    xlim([this_opt(1), this_opt(end)])

else
    
    %using linear scale for all other hyperparameters
    figure
    plot(this_opt,mean(perft'),'b')
    hold on
    plot(this_opt,mean(perftt'),'k')
    plot(this_opt,mean(perft')+ std(perft')/sqrt(size(perft,2)),'b')
    plot(this_opt,mean(perft')-std(perft')/sqrt(size(perft,2)),'b')
    plot(this_opt,mean(perftt')+ std(perftt')/sqrt(size(perftt,2)),'k')
    plot(this_opt,mean(perftt')-std(perftt')/sqrt(size(perftt,2)),'k')
    ylim([0 1])

end
%data are displayed as mean +/- SEM

title(['Algorithmic performance for multiple values of ' which_opt])
xlabel('Hyperparameter value')
ylabel('Performance') 





%%
function [NEWMEAN,NEWVAR,K_maj] = CalcPrediction(PERFORMANCE,NoiseEstim,distKern,K_maj,NUMELMAP)
% Takes the performance matrix, the prior prediction and the prior
% variance map and finds the new prediction and variance matrices.
    t = size(PERFORMANCE,1);
    AddedRow = distKern(PERFORMANCE(end,1),PERFORMANCE(1:(end),1));
    K_maj(t,1:t) = AddedRow;
    K_maj(1:t,t) = AddedRow;
    K_maj_n = K_maj + eye(t)*NoiseEstim;
    KInv = K_maj_n^-1;
    NEWMEAN = zeros(1,NUMELMAP);
    NEWVAR = NEWMEAN;
    for l = 1:NUMELMAP    
        k_min = distKern(l,PERFORMANCE(:,1));    
        NEWMEAN(l) = k_min*KInv*PERFORMANCE(:,2);
        NEWVAR(l) = 1 - k_min*KInv*(k_min)';
    end
end
