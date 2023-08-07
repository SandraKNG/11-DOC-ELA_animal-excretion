% function [X] = combineSampleFiles(fpattern)
% fpattern is the file pattern.  For example, to combine samples in xls 
% files in a given directory, fpattern should be '*.xls'
% 
% Usage: 
% At the console, type bigMatrix = combineSampleFiles('*.xls');
% The data from all the files will be stored in bigHonkinMatrix 
%
% P.S. You don't have to really call it bigHonkinMatrix. You could call it
% X if you'd like. 
%
% Written by Celeste Johanson on a cool summer evening in August 2006
% (64.5F on the roof tonight.) 
function [X] = combineSampleFiles(fpattern)

filenames = dir(fpattern);

%creating a structure to hold the filenames.  It is really just an vector
%of filenames. Don't panic!
sortingMatrix = struct('name', {});   

% sort the file names
for i=1:length(filenames)
    currname = filenames(i).name;               %name of the file for this loop iteration
    endindex = find(currname == '.') - 1;       %the index of the . in the file.  I assume there is only 1 .
    startindex = find(currname == '_') + 1;     %the index of the underscore in the file.  Again, assuming only 1
    
    filenumber = currname(startindex:endindex); %the file number found inbetween the _ and the . in the file name
    %this is how the filenames are sorted.  I put the filename in a matrix in 
    %the position corresponding the the filenumber.  It is OK to have
    %filenumbers that are non sequential, i.e. 1, 2, 4.  This will make a
    %vector of length 4 where the 3rd entry is blank.  Blank names are
    %ingored in the next loop.  
    sortingMatrix(str2num(filenumber)).name = char(currname);   
               
                                                            
end

counter = 1; 

% OK!  Now loop through the sorted filenames
for i=1:length(sortingMatrix)
    currname = sortingMatrix(i).name;
    % check to make sure the filename isn't empty.  
    if(length(currname) > 3)  
            
        %comment this line to supress output of the current file processed
        %in the loop. 
        sprintf('%s', currname)
        
        %load the file and put the contents into the big 3-D matrix.  
        currfile = load(currname);  
        bigMatrix(counter, :,:) = currfile; 
        
        %just counting valid files with this variable for correct indexing
        %in the big Matrix. 
        counter = counter + 1; 
        
        clear currfile; 
    end
end

%comment the following lines (or remove them).  I just put these in for
%your own error checking to make sure I wrote the code right. 
sprintf('%s', 'The dimensions of the big Matrix are: ')
size(bigMatrix)

X = bigMatrix; 