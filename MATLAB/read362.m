function data = read362(filename)
% filename= 'C:\Users\Daniel\Downloads\viblvb';
Fs = 400;  %hz
endianness = 'l';
dataFormat = 'uint16';
fileID=fopen(filename,'rb','l');
rawdata=fread(fileID,Inf,dataFormat,endianness);

X =[];Y =[];Z =[];T =[];
for ind=1:numel(rawdata)
    bits = bitget(rawdata(ind),1:16);
    sens = bits(15:16); % get sample type
    bits(13:16) =  bits(12);
    % get sign
    if bits(12)== 1
        signN = -1;
        bitsn = ~bits ;
    elseif bits(12) == 0
        signN = 1;
        bitsn = bits;
    else
        error('sign error')
    end
    %cunver
    num = bin2num(bitsn,signN);

    switch num2str(sens)
        case '0  0'
            X(end+1) = num;
        case '1  0'
            Y(end+1) = num;
        case '0  1'
            Z(end+1) = num;
        case '1  1'
            T(end+1) = num;
    end
end

fclose(fileID);
clear fileID

data.t = 0:1/Fs:length(X)/Fs-1/Fs;
data.x = X/1000;
data.y = Y/1000;
data.z = Z/1000;
data.temp = T*0.065;

plot(data.t,data.x,data.t,data.y,data.t,data.z,data.t,data.temp);
legend('X','Y','Z','T'); 
ylabel('G and C'); xlabel('time');title('362');
shg 
end

function num = bin2num(bitsn,s)
temp = 0;
for ind = 1:numel(bitsn)
    temp = bitsn(ind)* 2^(ind-1) + temp;
end
if s == -1
    num = (temp +1) * s ;
else
    num=temp ;
end
end
