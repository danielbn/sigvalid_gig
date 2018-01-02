
function [data,t , data_fft, f, Fs, params] = pars_raw_data(filename,hpfilter,toplot,scaledata)
%hpfilter is the  hipass cutoff frequency in  Hz 
filename = strtrim(filename);
params.filename = filename;

fileID=fopen(filename,'rb','l');
header=fgetl(fileID);
if contains(header,'^^%%') || contains(header,'^^%')% check to see if there is a header to the file
    while contains(class(header), 'char') && ~contains(header,'^^$$')
        %%% read header
        I1 = strfind(header,'</');%
        I2 = strfind(header,'>');%
        for ind=1:numel(I1)
            I3 = I1(ind) - I2;
            [~,I4]  = min( I3(I3>0) );
            if isempty(I4)
                continue
            end
            params.(header(I1(ind)+2:I2(I4+1)-1)) = header(I2(I4)+1:I1(ind)-1); % disp([  header(I1(ind)+2:I2(I4+1)-1) ':' header(I2(I4)+1:I1(ind)-1)])
        end
        
        if exist( 'toplot' , 'var') % if plot flag is set display header to consol too
            disp(header);
        end
        header=fgetl(fileID);
    end
    Fs = str2double(params.Frequency);
else % no header found
    frewind(fileID);
     Fs = 20833; % halo Fs
    disp('coldn''t find header');
end

%%% check header for data format
if isfield(params,'dataFormat')
    dataFormat = lower(params.dataFormat);
else
    dataFormat = 'int16';
end
if isfield(params,'endianness')
    switch params.endianness
        case 'LITTLE'
            endianness = 'l';
        case 'BIG'
            endianness = 'b';
    end
else
    endianness = 'l';
end
position = ftell(fileID);
data=fread(fileID,Inf,dataFormat,endianness);

if ~exist( 'scaledata' , 'var')
    scaledata = 1;
end

%%% scale data
if scaledata == 1
    try
        if contains(params.Sensor,'VIBRATION_1')
            switch Fs
                case 20833% 'halo'
                    data =data*2.5/ 2^15 / 0.04 ;% /sqrt(2)  ; % convert data to g the sqrt(2) is to convert g to g(rms)
                case 49019 % 'auguscope'
                    data = 1*data/2^15 / 0.1 * str2double(params.Gain) ; % convert data to g , auguscop data is in RMS
                otherwise
                    disp('data was not re-scaled')
            end
        elseif contains(params.Sensor,'MAGNETIC_1')
            data = 2.5*data/2^15 / 0.003125;
        elseif contains(params.Sensor,'TEMPERATURE_5')
            B=4485;
            data_mean =mean(data);
            if data_mean>0
                R0=47e3;
                T0=273.15 + 25;
                data = (log( data/R0)/B+1/T0).^-1 - 273.15 ; % scale data to degree c
            end
        end
    catch % no header found
        if Fs == 20833 % halo
            switch  isstrprop(filename(end),'digit') 
                case 0
                    if filename(end) == 'm' || contains(filename,'mag')% if 
                        data = 2.5*data/2^15 / 0.003125;
                        params.Sensor = 'MAGNETIC_1';
                        disp('no header MAGNETIC data')
                        %                         data =[data , data1];
                    elseif filename(end) == 't' % if temp
                        data = temp_scale(fileID,position);
                        params.Sensor = 'TEMPERATURE_5';
                    end
                case 1 % VIBRATION
                    if str2double(filename(end)) < 3
                        data = 2.5*data/2^15 / 0.04 ; % /sqrt(2)  ; % convert data to g the sqrt(2) is to convert g to g(rms)
                    elseif str2double(filename(end)) == 3% MAGNETIC
                        data = 2.5*data/2^15 / 0.003125; %
                    elseif str2double(filename(end)) == 4% TEMPERATURE
                        data = temp_scale(fileID,position);
                    end
                    params.Sensor = 'VIBRATION_1';
                otherwise
                    disp([filename 'data was not re-scaled'])
            end
        end
    end
else
    disp('data was not re-scaled');
end


fclose(fileID);
clear fileID

if exist('hpfilter','var')
        data = hipass_filter(data,Fs,hpfilter );
end

%%% calc fft
[data_fft , f ,~] =  myfft(data,Fs);
t = 0:1/Fs:length(data)/Fs-1/Fs;
params.Frequency = Fs;

% I = strfind( params.filename , strtrim(fullfile(' ',' ')));
% if isempty(I)
% name = params.filename;
% else
%     name = params.filename(I(end)+1:end);
% end
% data_cell.(name).data = data;
% data_cell.(name).time = t;
% data_cell.(name).data_fft = data_fft;
% data_cell.(name).f = f;
% data_cell.(name).params = params;

%%% if to plot data
if exist( 'toplot' , 'var')
    subplot(2,1,1) ;
    plot(t , data); ylabel('g') ;xlabel('time') ;xlim([0,max(t)]) ;hold all ;
    subplot(2,1,2);
    plot(f,data_fft) ; xlabel('Hz') ; hold all ;
    shg;
end

%% scale temp data fuction
function data = temp_scale(fileID,position)
fseek(fileID,position,'bof');
data=fread(fileID,Inf,'uint32','l');
% temp scaling
B=4485;
data_mean =mean(data);
if data_mean>0
    R0=47e3;
    T0=273.15 + 25;
    data = (log( data/R0)/B+1/T0).^-1 - 273.15 ; % scale data to degree c
end

% scale data to degree c
%                 dlTemp=diff(TempK-273.15);

% old read header
%         if strfind(header,'Frequency')
%             ind = strfind(header,'Frequency');
%             Fs = str2double(header(ind(1)+10:ind(2)-3));
%             params.frequency = Fs;
%         elseif strfind(header,'Sensor')
%             ind = strfind(header,'Sensor');
%             params.Sensor = header(ind(1)+7:ind(2)-3);
%         elseif strfind(header,'Date')
%             ind = strfind(header,'Date');
%             params.date = header(ind(1)+5:ind(2)-3);
%         elseif strfind(header,'Time')
%             ind = strfind(header,'Time');
%             params.time = header(ind(1)+5:ind(2)-3);
%         end


%%%  example for the file header
% ^^%%
% <?xml version="1.0" encoding="UTF-8"?>
% <AuguryDat><Date>01/10/2017</Date>
% <Time>05:12:14</Time>
% <hwVersion>C1</hwVersion>
% <Sensor>ULTRASONIC_1</Sensor> VIBRATION_1
% <Robin>100</Robin>
% <FrameLength>65536</FrameLength>
% <clientVersion>0.9.71.2+a42c4f7</clientVersion>
% <augurydatVersion>0.4</augurydatVersion>
% <Frequency>117647</Frequency></AuguryDat>
% % ^^$$
% <?xml version="1.0" encoding="UTF-8"?>
% <AuguryDat>
% 	<augurydatVersion>0.6</augurydatVersion>
% 	<dataFormat>UINT32</dataFormat>
% 	<endianness>LITTLE</endianness>
% 	<timestamp>1510748771</timestamp>
% 	<Frequency>100</Frequency>
% 	<Sensor>TEMPERATURE_5</Sensor>
% 	<SensorVersion>NCP15WB473F03RC</SensorVersion>
% 	<fwVersion>01.01.a074</fwVersion>
% 	<hwVersion>apus_alpha</hwVersion>
% 	<clientVersion>odroidc2_v0.3.9-02.01.1031+release-16bbd81955bf57a3e0cd35f7b30cbda5276c6ae5</clientVersion>
% 	<overflow>false</overflow>
% 	<underrun>false</underrun>
% 	<bist>false</bist>
% </AuguryDat>
% ^^$$
% %%%
% ^^%
% <?xml version="1.0" encoding="UTF-8"?>
% <AuguryDat>
% 	<augurydatVersion>0.6</augurydatVersion>
% 	<dataFormat>INT16</dataFormat>
% 	<endianness>LITTLE</endianness>
% 	<timestamp>1510748771</timestamp>
% 	<Frequency>20833</Frequency>
% 	<Sensor>VIBRATION_1</Sensor>
% 	<SensorVersion>ANALOGDEVICES_ADXL1002-50_1</SensorVersion>
% 	<fwVersion>01.01.a074</fwVersion>
% 	<hwVersion>apus_alpha</hwVersion>
% 	<clientVersion>odroidc2_v0.3.9-02.01.1031+release-16bbd81955bf57a3e0cd35f7b30cbda5276c6ae5</clientVersion>
% 	<overflow>false</overflow>
% 	<underrun>false</underrun>
% 	<bist>false</bist>
% </AuguryDat>
% ^^$$
% ^^%
% <?xml version="1.0" encoding="utf-8"?>
% <AuguryDat>
% 	<Date>03/11/2017</Date>
% 	<Time>09:39:10</Time>
% 	<hwVersion>C1</hwVersion>
% 	<Sensor>VIBRATION_1</Sensor>
% 	<Robin>100</Robin>
% 	<FrameLength>1536</FrameLength>
% 	<clientVersion>1.4.40</clientVersion>
% 	<augurydatVersion>0.5</augurydatVersion>
% 	<Frequency>49019</Frequency>
% 	<Gain>2</Gain>
% </AuguryDat>
% ^^$$

