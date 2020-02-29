clear

trainTxtFile = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\trainset.txt';

mainDir = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\classes\jpg\';

T = readtable(trainTxtFile);

nTrainFiles = size(T,1);

logoDir = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\Logos';
imageDir = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\TestImages';
allFilesDir ='E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\AllImages';

dictionaryPath = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\dictionary.mat';


if ~exist(logoDir)
   mkdir( logoDir);
end

if ~exist(imageDir)
   mkdir( imageDir);
end

if ~exist(allFilesDir)
   mkdir( allFilesDir);
end

logoDictionary = [];

logoList = table2cell(T(:,2));

for i = 1:nTrainFiles
    
    logoName = table2cell(T(i,1));
    logoName = logoName{1};
    fileName = table2cell(T(i,2));
    fileName = fileName{1};
    
    imagePath = [mainDir '\' logoName '\' fileName]; 
    
    if isempty(logoDictionary)  %empty dictionary, add first entry
        
        logoDictionary = {logoName, 1};
            
        idx = 1;
    else
          
        idx = find(strcmp( {logoDictionary{:,1}},logoName));   %find key in dictionary
      
        if idx > 0
            logoDictionary{idx,2} = logoDictionary{idx,2} + 1;  %found key, increment logoNumber
        else 
            idx = size(logoDictionary,1) + 1;               %did not find key, add entry
            
            logoDictionary{idx,1} = logoName;
            logoDictionary{idx,2} = 1;           
            
            
        end
    
    end
    
    logoNumber = logoDictionary{idx,2};
    
    
    newName = [logoName int2str(logoNumber)];
    
    newPath = strrep([logoDir '\' newName '.jpg'],' ','');
    
   copyfile(imagePath,newPath);
    
    
end

imageDictionary = logoDictionary;
Nlogos = size(logoDictionary,1);

for i=1:Nlogos
    imageDictionary{i,2} = 0;    
end

imageDictionaryAll = imageDictionary;

dirList = dir(mainDir);
Ndirs = length(dirList);

isSubDir = [dirList(:).isdir];  

for d =1:Ndirs

directory_name = dirList(d).name;

if(strcmp(directory_name(1),'.'))  %hidden folder
    continue
end

if(strcmp(directory_name,'no-logo'))  %not analyzing no logo folder
    continue
end

if ~isSubDir(d)     %it is a file not a directory
   continue 
end

directory_name = [mainDir '\' directory_name];

[~,logoName] = fileparts(directory_name);

filelist = dir([directory_name '\' '*.jpg']);

NTestFiles = length(filelist);

for i=1:NTestFiles
   
    oldName = filelist(i).name;
    
       idx = find(strcmp( {imageDictionary{:,1}},logoName)); 

   imageDictionaryAll{idx,2} = imageDictionaryAll{idx,2}+1;
    
    if find(strcmp(logoList,oldName)) %test logo file
            
    oldPath = [directory_name '\' oldName];
    newName = [imageDictionaryAll{idx,1} int2str(imageDictionaryAll{idx,2})];
    newPath = [allFilesDir '\' newName '.jpg'];
    copyfile(oldPath,newPath);         
         
        continue 
    end
     
    imageDictionary{idx,2} = imageDictionary{idx,2}+1;
     
    oldPath = [directory_name '\' oldName];
    newName = [imageDictionary{idx,1} int2str(imageDictionary{idx,2})];
    newPath = [imageDir '\' newName '.jpg'];
    copyfile(oldPath,newPath);
    
    newPath = [allFilesDir '\' newName '.jpg'];
    copyfile(oldPath,newPath); 
end

end

save(dictionaryPath, 'logoDictionary');













