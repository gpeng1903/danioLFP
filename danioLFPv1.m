%% master EEG

% reads LFP files stored in a folder
% process files and output result tables
% Gang Peng @ 2021



%% read files
% tic;
% bandTable
bandTable = cell(6, 2);
bandTable{1,1} = 'delta'; bandTable{1,2} = [1, 4]; 
bandTable{2,1} = 'theta'; bandTable{2,2} = [4, 8];
bandTable{3,1} = 'alpha'; bandTable{3,2} = [8, 12];
bandTable{4,1} = 'beta';  bandTable{4,2} = [12, 30];
bandTable{5,1} = 'gamma1'; bandTable{5,2} = [30, 44];
bandTable{6,1} = 'gamma2'; bandTable{6,2} = [56, 300];

% parameters
spsRate = 5000; % sampling rate

% select 10 min
testLenMin = 10; % in min
testPadSec = 10; % in sec
testLenPoints = testLenMin * spsRate * 60;
testPadPoints = testPadSec * spsRate;

% detrend parameter
detrend_segment = 20; 
detrend_brkLen = detrend_segment * spsRate;

% initiate para to read data
dataFolder = pwd;
dataList = dir(strcat(dataFolder,'\*.txt'));

dataFileNum = size(dataList, 1);
dataStore = cell(dataFileNum, 1); % data read into dataStore

% to get byte size from dataList
fileBytes = [dataList.bytes].'; % first use [] to get all bytes, then use .' to change to vector
minBytes = min(fileBytes);% 
% fprintf('%s%d\n', 'shortest trace in minutes ', (minBytes/8)/(spsRate*60));

% initiate var to store results
dataResultStore = cell(dataFileNum, 4);


%% read data, loop

for i = 1:dataFileNum
    
    fName = dataList(i).name;
    vaName = fName(1:end-4);
    
        % import data
        % import options
        opts = delimitedTextImportOptions("NumVariables", 1);        
        % import range and separator
        opts.DataLines = [3, Inf];
        opts.Delimiter = ",";
        % data type
        opts.VariableNames = "VarName1";
        opts.VariableTypes = "double";
        % file type
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        % read into cell array, start from 10 sec after recording, and read 10 min data
        readData = table2array(readtable(fName, opts));        
        dataStore{i, 1} = readData(testPadPoints+1:testPadPoints+testLenPoints);
        
        
        % reset options
        clear opts;
        
        fprintf('%s%d%s%d\n', 'data loaded into dataStore ', i, ' out of ', dataFileNum);
        
end

%% pre-processing

for data_i = 1:dataFileNum
    LfpD = dataStore{data_i};
    
    % detrend
    LfpD_detrnd = detrend(LfpD, 1, 1:detrend_brkLen:length(LfpD));
    
    % save detrend data
    fName = dataList(data_i).name;
    vaName = fName(1:end-4);  
%     ori_vName = strcat('ori_', vaName);
%     eval([ori_vName '= LfpD_detrnd;']);    
    
    Lfp_prep0 = LfpD_detrnd;  

    % for 5000 sample rate, use universalThreshold
    Lfp_prep1 = wdenoise(Lfp_prep0, 'DenoisingMethod', 'UniversalThreshold'); 
    
    % lowpass, 
    Lfp_prep = lowpass(Lfp_prep1,300,spsRate,'Steepness',0.5,'StopbandAttenuation',60);    
   
    % save result 1, save preprocessed data
    dataResultStore{data_i, 1} = dataList(data_i).name;
    dataResultStore{data_i, 2} = Lfp_prep;    
    
    % bandpower, and percentage bandpower, in matrix
    % initiate bandpower matrix
    bandN = size(bandTable,1); 
    epBandsRaw = zeros(bandN, 1);    
    epBands = zeros(bandN-1, 2);    
    for band_j = 1:bandN
        epBandsRaw(band_j) = bandpower(Lfp_prep, spsRate, bandTable{band_j, 2});
    end
    epBands(1:end-1, 1) = epBandsRaw(1:end-2);
    epBands(end, 1) = epBandsRaw(end-1) + epBandsRaw(end); 
    epBands(:,2) = epBands(:,1) / sum(epBands(:,1));
    
    % save results 2, save epBands
    dataResultStore{data_i, 3} = epBands;   
    
    
    
    fprintf('%s%d%s%d\n', 'data analyzed and saved ', data_i, ' out of ', dataFileNum);    
    
end

%%  re-organize for excel plot
resultMatrix = zeros(dataFileNum, 11);
for data_i = 1:dataFileNum
    
    epBands = dataResultStore{data_i, 3};
    resultMatrix(data_i, 1:5) = epBands(:,1)';
    resultMatrix(data_i, 7:11) = epBands(:,2)';
    
end

%% main processing

%% find peaks, epileptic discharge calling 
testDate = 'fill_data_description_here';


% delete previous results, if not done so yet
delete D:\tmp\epiLFP\lfp*
delete D:\tmp\epiLFP\*.xlsx

% use loop to write sets of parameters
parameterVault = [0.6, 3,  0.18;
                  0.8, 3,  0.18; 
                  0.6, 4,  0.18;
                  0.8, 4,  0.18;
                  1.0, 4,  0.18];

% additional called criteria             
spikeThreshold = 2.4; 
powerThreshold = 0.18; 
checkN_max = 5;    
              
% manually set excel range marker              
excelRange = {'A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'I1', 'J1', 'K1', 'L1'...
                  'M1', 'N1', 'O1', 'P1', 'Q1', 'R1', 'S1', 'T1', 'U1', 'V1'...
                  'W1', 'X1', 'Y1', 'Z1'};              
              
              
              
% loop through sets of parameters
for paraSetI = 1:size(parameterVault, 1)              

    % write to excel, generate file names
    lfpN_fName = strcat('D:\tmp\epiLFP\lfpN_', testDate, num2str(paraSetI), '.xlsx');
    
    lfpDur_fName = strcat('D:\tmp\epiLFP\lfpDur_', testDate, num2str(paraSetI), '.xlsx');
    lfpRMS_fName = strcat('D:\tmp\epiLFP\lfpRMS_', testDate, num2str(paraSetI), '.xlsx');
    lfpPWR_fName = strcat('D:\tmp\epiLFP\lfpPWR_', testDate, num2str(paraSetI), '.xlsx');   
    
    mean_lfpDur_fName = strcat('D:\tmp\epiLFP\mean_lfpDur_', testDate, num2str(paraSetI), '.xlsx');
    mean_lfpRMS_fName = strcat('D:\tmp\epiLFP\mean_lfpRMS_', testDate, num2str(paraSetI), '.xlsx');
    mean_lfpPWR_fName = strcat('D:\tmp\epiLFP\mean_lfpPWR_', testDate, num2str(paraSetI), '.xlsx');
    
    sum_lfpDur_fName = strcat('D:\tmp\epiLFP\sum_lfpDur_', testDate, num2str(paraSetI), '.xlsx');
    sum_lfpRMS_fName = strcat('D:\tmp\epiLFP\sum_lfpRMS_', testDate, num2str(paraSetI), '.xlsx');
    sum_lfpPWR_fName = strcat('D:\tmp\epiLFP\sum_lfpPWR_', testDate, num2str(paraSetI), '.xlsx');
    
    % set parameter
    lfpThreshold = parameterVault(paraSetI, 1); 
    minimal_IPI = spsRate/parameterVault(paraSetI, 2); 

    durThreshold = parameterVault(paraSetI, 3); 
    tailPadding = spsRate/100;  



    lfpPksStore = cell(dataFileNum, 1);
    lfpPksStoreTrim = cell(dataFileNum, 1);

    lfpFlatTableTrim = table;

    for lfpI = 1:dataFileNum

        lfpTest = dataResultStore{lfpI, 2}; 


        % algorithm to count peaks
        tmpPks = find(abs(lfpTest) > lfpThreshold); 

        % some trace doesn't have any peaks above threshold
        % take care of this

        if ~isempty(tmpPks) 

            tmpPksDiff = diff(tmpPks);
            tmpPksCount = find(tmpPksDiff > minimal_IPI); 

            tmpPksCheck = tmpPks(tmpPksCount);

            PksNum = size(tmpPksCount, 1);
            PksLocation = zeros(PksNum+1, 2);

            pksPos1 = tmpPks(1); 

            
            for i = 1:PksNum
                PksLocation(i, 1) = pksPos1 - tailPadding; 
                PksLocation(i, 2) = tmpPks(tmpPksCount(i)) + tailPadding;    
                pksPos1 = tmpPks(tmpPksCount(i)+1);  
            end
            
            PksLocation(end, 1) = pksPos1 - tailPadding; 
            PksLocation(end, 2) = tmpPks(end) + tailPadding; 
            
            if PksLocation(end, 2) > length(lfpTest)
                PksLocation(end, 2) = length(lfpTest);
            end

            

            % some statistics
            PksDur = (PksLocation(:,2) - PksLocation(:,1))./spsRate;
            
            PksHeight = zeros(size(PksDur,1), 1);  % max amplitude 
            PksRms = zeros(size(PksDur,1), 1); % RMS of signal
            PksPower = zeros(size(PksDur,1), 1); % power of signal

            bandN = size(bandTable,1); % gamma was calculate from 2 entries    
            PksBandPower = zeros(size(PksDur,1), bandN-1);

            for i = 1:size(PksRms,1)
                peakI = lfpTest(PksLocation(i,1):PksLocation(i,2));
                
                % PksHeight(i) = max(abs(peakI));                
                PksHeight(i) = min(maxk(abs(peakI), checkN_max)); % for robustness, check top 3 values                
                PksRms(i) = rms(peakI);
                PksPower(i) = bandpower(peakI);


                % bandpower, and percentage bandpower, in matrix
                % initiate bandpower matrix
            %     bandN = size(bandTable,1); % gamma was calculate from 2 entries    
                epBandsRaw = zeros(bandN, 1);    
                epBands = zeros(bandN-1, 1); % -1, due to gamma split  
                for band_j = 1:bandN
                    epBandsRaw(band_j) = bandpower(peakI, spsRate, bandTable{band_j, 2});
                end
                epBands(1:end-1, 1) = epBandsRaw(1:end-2);
                epBands(end, 1) = epBandsRaw(end-1) + epBandsRaw(end); % gamma band, add two splited range [30, 40], [60, 100];

                PksBandPower(i,:) = epBands';
            end

            % initial peaks
            lfpPksTable = table(PksLocation, PksDur, PksHeight, PksRms, PksPower, PksBandPower);
            lfpPksStore{lfpI, 1} = lfpPksTable;

            
            calledPeaks = (lfpPksTable.PksDur > durThreshold) & ...
                          (lfpPksTable.PksHeight > spikeThreshold) & ...
                          (lfpPksTable.PksPower > powerThreshold);
            lfpPksTable_trim = lfpPksTable(calledPeaks, :); 
            
            
            lfpPksStoreTrim{lfpI, 1} = lfpPksTable_trim;

            % cat table
            if size(lfpPksTable_trim,1) > 0
                tableIndex = lfpI .* ones(size(lfpPksTable_trim,1),1);
                addIdxT2 = addvars(lfpPksTable_trim, tableIndex, 'before', 'PksLocation');    
                lfpFlatTableTrim = cat(1, lfpFlatTableTrim, addIdxT2);
            end


            
            lfpFinal = lfpTest;
            lfpRmTable = lfpPksTable(lfpPksTable.PksDur <= durThreshold, :);

            lfpRmIdx = zeros(size(lfpTest,1), 1);
            for checkI = 1: size(lfpRmTable,1)    
                lfpRmIdx(lfpRmTable.PksLocation(checkI,1): lfpRmTable.PksLocation(checkI,2)) = 1;
            end
            lfpRmIdx = logical(lfpRmIdx);
            lfpFinal(lfpRmIdx) = [];

            dataResultStore{lfpI, 4} = lfpFinal;


        else
            dataResultStore{lfpI, 4} = lfpTest;

        end



        disp(lfpI);

    end



    % size count
    entryN = length(lfpPksStoreTrim);

    lfpNcount = zeros(entryN, 1);

    for i = 1:entryN
        lfpNcount(i) = size(lfpPksStoreTrim{i},1);
    end


    % make geneList, using the first 5 characters
    geneList = cell(dataFileNum, 1);

    for i = 1:dataFileNum    
        geneList{i} = dataList(i).name(1:5);
    end

    
    [uniqGene, uia, uic] = unique(geneList); 
    uiaPad = [uia; dataFileNum+1]; 




    

    % take lfpN
    numOfGene = size(uniqGene, 1);


    % write lfpN excel
    for i = 1:numOfGene    
        gStart = uiaPad(i);
        gEnd = uiaPad(i+1)-1; % -1

        gSymbol = uniqGene{i};
        gLfpN_Val = lfpNcount(gStart:gEnd);

        gTable = table(gLfpN_Val);
        gTable.Properties.VariableNames = {gSymbol};

        writePos = excelRange{i};
        writetable(gTable, lfpN_fName, 'Range', writePos, 'WriteMode', 'inplace');
    end

    
    % write lfpDur, lfpPWR, lfpRMS excel    
    for i = 1:numOfGene    
        gStart = uiaPad(i);
        gEnd = uiaPad(i+1)-1; % -1

        gLfpDur = [];
        gLfpRMS = [];
        gLfpPWR = [];

        gSymbol = uniqGene{i};

        for traceI = gStart:gEnd

            lfpPks_Results = lfpPksStoreTrim{traceI, 1}; 

            if size(lfpPks_Results, 1) > 0
                gLfpDur = cat(1, gLfpDur, lfpPks_Results.PksDur);
                gLfpRMS = cat(1, gLfpRMS, lfpPks_Results.PksRms);
                gLfpPWR = cat(1, gLfpPWR, lfpPks_Results.PksPower);
            end
        end
       

        % write tables 1
        gTable = table(gLfpDur);
        gTable.Properties.VariableNames = {gSymbol};

        writePos = excelRange{i};
        writetable(gTable, lfpDur_fName, 'Range', writePos, 'WriteMode', 'inplace');

        % write tables 2
        gTable = table(gLfpRMS);
        gTable.Properties.VariableNames = {gSymbol};

        writePos = excelRange{i};
        writetable(gTable, lfpRMS_fName, 'Range', writePos, 'WriteMode', 'inplace');


        % write tables 3
        gTable = table(gLfpPWR);
        gTable.Properties.VariableNames = {gSymbol};

        writePos = excelRange{i};
        writetable(gTable, lfpPWR_fName, 'Range', writePos, 'WriteMode', 'inplace');

    end
    
    % write mean and sum of lfpDur, lfpPWR, lfpRMS excel    
    for i = 1:numOfGene    
        gStart = uiaPad(i);
        gEnd = uiaPad(i+1)-1; % -1
        gSymbol = uniqGene{i};
        
        
        % initiate 
        tempSize = gEnd - gStart + 1;        
        mean_gLfpDur = zeros(tempSize, 1);
        mean_gLfpRMS = zeros(tempSize, 1);
        mean_gLfpPWR = zeros(tempSize, 1);
        
        sum_gLfpDur = zeros(tempSize, 1);
        sum_gLfpRMS = zeros(tempSize, 1);
        sum_gLfpPWR = zeros(tempSize, 1);
        
        for traceI = gStart:gEnd

            lfpPks_Results = lfpPksStoreTrim{traceI, 1}; 
            
            temp_k = traceI - gStart + 1; % temp idx
            
            if size(lfpPks_Results, 1) > 0
                % mean
                mean_gLfpDur(temp_k,1) = mean(lfpPks_Results.PksDur);
                mean_gLfpRMS(temp_k,1) = mean(lfpPks_Results.PksRms);
                mean_gLfpPWR(temp_k,1) = mean(lfpPks_Results.PksPower);
                % sum
                sum_gLfpDur(temp_k,1) = sum(lfpPks_Results.PksDur);
                sum_gLfpRMS(temp_k,1) = sum(lfpPks_Results.PksRms);
                sum_gLfpPWR(temp_k,1) = sum(lfpPks_Results.PksPower);
            else
                % mean
                mean_gLfpDur(temp_k,1) = 0;
                mean_gLfpRMS(temp_k,1) = 0;
                mean_gLfpPWR(temp_k,1) = 0;
                % sum
                sum_gLfpDur(temp_k,1) = 0;
                sum_gLfpRMS(temp_k,1) = 0;
                sum_gLfpPWR(temp_k,1) = 0;
            end
            
        end
        
        % write tables 1
        gTable = table(mean_gLfpDur);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, mean_lfpDur_fName, 'Range', writePos, 'WriteMode', 'inplace');
        % write tables 2
        gTable = table(mean_gLfpRMS);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, mean_lfpRMS_fName, 'Range', writePos, 'WriteMode', 'inplace');
        % write tables 3
        gTable = table(mean_gLfpPWR);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, mean_lfpPWR_fName, 'Range', writePos, 'WriteMode', 'inplace');
        
         % write tables 4
        gTable = table(sum_gLfpDur);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, sum_lfpDur_fName, 'Range', writePos, 'WriteMode', 'inplace');
        % write tables 5
        gTable = table(sum_gLfpRMS);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, sum_lfpRMS_fName, 'Range', writePos, 'WriteMode', 'inplace');
        % write tables 6
        gTable = table(sum_gLfpPWR);
        gTable.Properties.VariableNames = {gSymbol};
        writePos = excelRange{i};
        writetable(gTable, sum_lfpPWR_fName, 'Range', writePos, 'WriteMode', 'inplace');

    end
    
    

end






