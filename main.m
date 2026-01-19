close all
clear all
clc

addpath(genpath('T1D_VPP'));

nPats = 99;
samples = 288*2;

bolusModulations = 0.5:0.1:1.5;
basalModulations = 0.5:0.1:1.5;
mealModulations = 0:0.2:2.0;

glucose = zeros(samples,nPats);
bolus = zeros(samples,nPats);
basal = zeros(samples,nPats);
meal = zeros(samples,nPats);
weights = zeros(nPats,1);

names = {'bolusMinus50','bolusMinus40','bolusMinus30','bolusMinus20','bolusMinus10',...
        'baseline',...
        'bolusPlus10','bolusPlus20','bolusPlus30','bolusPlus40','bolusPlus50'};

for n = 1:length(names)

    disp(['Creating DatasetPJ_' names{n} '...']);

    for nn = 1:nPats

        [G, CHO, IB, Ib, BW] = runSimulation(nn,bolusModulations(n),1,1);
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;

    end

    saveResults(glucose, bolus, basal, meal, weights, ['DatasetPJ_' names{n}])

end

names = {'basalMinus50','basalMinus40','basalMinus30','basalMinus20','basalMinus10',...
        'baseline',...
        'basalPlus10','basalPlus20','basalPlus30','basalPlus40','basalPlus50'};
        
for n = 1:length(names)
    
    disp(['Creating DatasetPJ_' names{n} '...']);
    
    for nn = 1:nPats

        [G, CHO, IB, Ib, BW] = runSimulation(nn,1,1,basalModulations(n));
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        
    end
    saveResults(glucose, bolus, basal, meal, weights, ['DatasetPJ_' names{n}])
end

names = {'mealMinus100','mealMinus80','mealMinus60','mealMinus40','mealMinus20',...
        'baseline',...
        'mealPlus20','mealPlus40','mealPlus60','mealPlus80','mealPlus100'};
        
for n = 1:length(names)
    
    disp(['Creating DatasetPJ_' names{n} '...']);
    
    for nn = 1:nPats

        [G, CHO, IB, Ib, BW] = runSimulation(nn,1,mealModulations(n),1);
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        
    end
    saveResults(glucose, bolus, basal, meal, weights, ['DatasetPJ_' names{n}])
end

names = {'hypotreatment'};
        
for n = 1:length(names)
    
    disp(['Creating DatasetPJ_' names{n} '...']);
    
    for nn = 1:nPats

        [G, CHO, IB, Ib, BW] = runSimulation(nn,1,1,1);
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        
    end
    saveResults(glucose, bolus, basal, meal, weights, ['DatasetPJ_' names{n}])
end


names = {'correctionBolus'};
        
for n = 1:length(names)
    
    disp(['Creating DatasetPJ_' names{n} '...']);
    
    for nn = 1:nPats

        [G, CHO, IB, Ib, BW] = runSimulation(nn,1,1,1);
        glucose(:,nn) = G(1:samples);
        bolus(:,nn) = IB(1:samples);
        basal(:,nn) = Ib(1:samples);
        meal(:,nn) = CHO(1:samples);
        weights(nn) = BW;
        
    end
    saveResults(glucose, bolus, basal, meal, weights, ['DatasetPJ_' names{n}])
end



function saveResults(glucose, bolus, basal, meal, weights, label)

    saveFolderName = fullfile('results', label);
    if ~exist(saveFolderName, 'dir')
        mkdir(saveFolderName)
    end
    
    % Cut
    startSample = 1;
    endSample = (11 + (24-6)) * 12;
    glucoseAll = glucose(startSample:endSample,:);
    basalAll = basal(startSample:endSample,:);
    bolusAll = bolus(startSample:endSample,:);
    mealAll = meal(startSample:endSample,:);
    
    for p = 1:size(glucose,2)

        t = datetime(2000,1,1,6,0,0):minutes(5):(datetime(2000,1,1,6,0,0)+minutes(5*size(glucoseAll,1)));
        t = t(1:end-1)';

        glucose = glucoseAll(:,p);
        cho = mealAll(:,p);
        bolus = bolusAll(:,p);
        basal = basalAll(:,p);

        cho_label = strings(length(cho),1);
        bolus_label = strings(length(cho),1);
        idxs = [12, 72, 156, 300];
        cho_label(idxs([1, 4])) = "B";
        cho_label(idxs([2])) = "L";
        cho_label(idxs([3])) = "D";
        bolus_label(idxs([1, 4])) = "B";
        bolus_label(idxs([2])) = "L";
        bolus_label(idxs([3])) = "D";

        data = table(t,glucose,bolus,bolus_label,basal,cho,cho_label);

        writetable(data, fullfile(saveFolderName,['patient_' num2str(p) '.csv']));

    end

    patient = (1:99)';
    bw = weights;
    
    data = table(patient,weights);
    writetable(data, fullfile(saveFolderName,['bw.csv']));
    
end