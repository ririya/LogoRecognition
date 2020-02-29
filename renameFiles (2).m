mainDir = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\classes\jpg\';

dirList = dir(mainDir);
Ndirs = length(dirList);

isSubDir = [dirList(:).isdir];  

for d =1:Ndirs
    


% directory_name = 'E:\FlickrLogos-32_dataset_v2\FlickrLogos-v2\classes\jpg\adidas';

directory_name = dirList(d).name;

if(strcmp(directory_name(1),'.'))
    continue
end

if ~isSubDir(d)
   continue 
end

directory_name = [mainDir '\' directory_name];

[~,logoName] = fileparts(directory_name);

filelist = dir([directory_name '\' '*.jpg']);

Nfiles = length(filelist);

for i=1:Nfiles
   
    oldName = filelist(i).name;
    oldPath = [directory_name '\' oldName];
    newName = [logoName int2str(i)];
    newPath = [directory_name '\' newName '.jpg'];
    movefile(oldPath,newPath);
end

end


