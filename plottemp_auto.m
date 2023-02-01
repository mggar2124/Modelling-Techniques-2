function [timeData,tempData] = plottemp_auto(input_image,time_max,temp_max,temp_min)
% functions takes a photo as an input and extracts the data points of both
% axes in the inputed photo and outputs two vectors (identical in size)
% pretaining to the two axis time and temprature.

%The variable of input_image is the variable for image name to be called
%The variable of time_max is the variable for the maximum time of the x-axis
%The variable of temp_max is the variable for the maximum time of the y-axis

image=imread([input_image '.jpg']); % Loads the image and reads it then stores the extracted image data as the variable image
[y,x] = find(image(:,:,1)>150 & image(:,:,2)<100 & image(:,:,3)<100); %searches and filters for the red pixels within the image

[xu,xindex] = unique(x,'last'); %Selects a unique red pixel within a column that contains red pixels

xy = y(xindex); %indexes the y values by x values into one matrix

tempData = xy;
timeData = xu; 

tempData = tempData*(-1); % flips the data across the horizontal axis
tempData = tempData - min(tempData); %starts at 0
coeff = (temp_max - temp_min)/max(tempData); %variable used to calculate the stretching coefficient
tempData = (tempData*coeff) + temp_min; %strechs for a maximum value that is held between the maximum and minimum temperature values

for x=3:length(tempData)-2 % a moving average to help improve smoothness
    tempData(x) = (tempData(x)+tempData(x+1)+tempData(x-1))/3; % average is done across three data points each time
end

timeData = timeData-timeData(1);
timeData = timeData * (time_max/timeData(end));%this is used help homogenise shape of the time vector to match the inputs

%plot(timeData,tempData) can be uncommented to ensure that the data
%extraction was successful