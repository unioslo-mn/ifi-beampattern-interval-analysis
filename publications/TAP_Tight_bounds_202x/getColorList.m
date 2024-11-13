function colorList = getColorList(ID)

% Link to Color Brewer query:
% https://colorbrewer2.org/index.html?type=testtype&scheme=TestMap&n=3#type=qualitative&scheme=Set1&n=7

    switch ID
        case 0 
            colorList = [166,206,227
                        31,120,180
                        178,223,138
                        51,160,44]/255;
        case 1
            colorList = [228,26,28
                        55,126,184
                        77,175,74
                        152,78,163]/255;    
        case 2
            colorList = [27,158,119
                        217,95,2
                        117,112,179
                        231,41,138]/255;
    end
end
