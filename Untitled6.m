imageFiles = dir([imageDir '\' '*.jpg']);

nImages = length(imageFiles);

logoFiles = dir([logoDir '\' '*.jpg']);

nLogosFiles =  dir([logoDir '\' '*.jpg']);

nLogos = size(logoDictionary,1);
truth = zeros(NTotalFiles,nLogos);

for i=1:nImages

    imageName = imageFiles(i).name;
    
    for j=1:nLogos
   
      if strfind(imageName,logoDictionary{j,1})
         truth(i,j) = 1;         
         break;
         
      end
      
    end
    
end

