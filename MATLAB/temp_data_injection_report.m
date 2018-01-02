% ccc; clear results
%%% Fs = 20833+1/3;
%% Main test parameters. And gen test signal
tc_num = 1111 %:4
EP = 1 %[13]%1:10


rec_time = 4;%3.9076 ; % 3.9076 for 1074 , 4.0060 for 1075
sig_time = 4;%3 ; % rec_time ;  %time length of transmited signal in seconds
Fs = 20833;
amp_accuracy = 0.02; % sigg to rec   2%
noise_level = 1e-4; % noise 0.1mg
noise_floor = noise_level * 10;  
hann_win = 1 ; % hann window on = 1 ,  off = 0
h_filter = 5 ; % highpass cutoff frequency in Hz 

FeatherHz = 10 ; % the space from peaks to make the noise window 

% folder_root = 'C:\Users\Daniel\Documents\Raw_data\signal_injection_1074\';

for EP_ID =EP
    for test_ID = tc_num
        switch test_ID
            case 1
                folder_test = ['EP_' num2str(EP_ID) '_testID_1'] ; % 'EP_test_TC43';
%                 folder_test = 'case1';
            case 2
                folder_test = ['EP_' num2str(EP_ID) '_testID_2'] ; % 'EP_test_TC4993';
%                 folder_test = 'case2';
            case 3 % 5
                folder_test = ['EP_' num2str(EP_ID) '_testID_3'] ; % 'EP_test_TCprim';
%                 folder_test = 'case3';
            case 4
                folder_test = ['EP_' num2str(EP_ID) '_testID_4'] ; % 'EP_test_TCprim';
%                 folder_test = 'case4';
            case 5 % 3
                folder_test = ['EP_' num2str(EP_ID) '_testID_5'] ; % 'EP_test_TCprim';
            otherwise 
               folder_test = '';
        end
        
        if ~strcmp(folder_test , '')
            test_name = folder_test;
        else
            test_name = 'EP_testID';
        end
        
        list = ls([ folder_root folder_test '\vib*']);
        if isempty ( list)
            continue
        end
        
        Fss=1e6;
        switch test_ID
            case 1 % Low frequency
                f_source = 43;
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 40e-3; % 10mv = 1/4g
            case 2 % High frequency
                f_source = 4993;
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 40e-3; % 40mv = 1/4g
%             case 3 % Wave packet
%                 fprim = primes(8200); % make only prime frequency this way none of them are a multiple of any other.
%                 f_source = fprim([1:7:floor(end/3)-1 , ceil(end/3):20:end ]); % 7 boom!! :-)
%                 offset=0;
%                 amplitude = ones(1,numel(f_source));
%                 amplitude=amplitude * 1e-3; % 4mv = 1/10g
            case 4 % Wave packet
                f_source = 0;
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 0; % 4mv = 1/10g
            case 3 % Wave packet
                fprim = primes(8200); % make only prime frequency this way none of them are a multiple of any other.
                f_source = fprim(10:100:end); % 7 boom!! :-)
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 1e-3; % 4mv = 1/10g
            case 111 % Wave packet
                f_source = [43 4095]; % 7 boom!! :-)
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 440e-3; % 4mv = 1/10g
            otherwise
                disp('*** bad test_ID ***')
                f_source = 0;
                offset=0;
                amplitude = ones(1,numel(f_source));
                amplitude=amplitude * 0; % 4mv = 1/10g
        end
        
        t_source=(0:1/Fss:sig_time-1/Fss)';
        outputdata = ones(numel(t_source),1);
        outputdata=2.5+offset*outputdata; %2.5 in volt
        
        disp('---> Generation signal');
        w=f_source*2*pi;
        for ind=1:length(w) % genarat the test signal
            outputdata=outputdata+amplitude(ind)*sin(w(ind)*t_source + 2*pi*rand);
        end

        disp('---> Signal ready');
        %%  load data,
        sigg = (outputdata-2.5)/0.04;
        if h_filter
            sigg = hipass_filter(sigg , Fss,h_filter) ;
        end
        % subplot (2,1,1) ; plot (t_source,(outputdata-2.5)/0.04) ; hold all
        [sigg_fft_perF , ff ,~] =  myfft(sigg,Fss,hann_win);
        
        % % <<<<<<<<<<<<<<<< add filter to source sig >>>>>>>>>>>>>>>>> %
        [num,~,~] = xlsread('C:\Users\Daniel\Dropbox (Augury)\RandD\MyMatlabfun\Filters_6K_10K.xlsx');
        end_ind = find (ff> num(end,1),1,'first') ;% remove high frequencies to save MEMORY
        Efilter = interp1([0 ; num(:,1)],[1; num(:,2)],ff(1:end_ind-1),'linear'); % add at f = 0 amp = 1;
        Efilter(isnan(Efilter)) = 0; % remove NaN;
        sigg_fft = Efilter'.*sigg_fft_perF(1:end_ind-1); 
        ff = ff(1:end_ind-1);
        % % <<<<<<<<<<<< end of noise >>>>>>>>>>>>>>> %
        
        % subplot (2,1,2) ; plot (ff,sigg_fft) ; hold all
        f1 = figure;
        subplot (2,1,2) ; plot (ff,sigg_fft) ; hold all
        for ind =1:3
            data{ind} = pars_raw_data([folder_root folder_test '\' strtrim(list(ind,:))],h_filter);
            N(ind) = numel(data{ind});
            %%%%%%%% pad zeros to data %%%%%%%%
%             if  round(N(ind)/Fs,4) < sig_time
%                 data{ind} = [ data{ind} ; zeros(sig_time*Fs -N(ind) -1,1) ]; % add zeros
%             else
%                 data{ind} = data{ind}(1:sig_time*Fs); %remove data
%             end
            Npad = numel(data{ind});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
            t = 0:1/Fs:(Npad-1)/Fs;
            %             axt{ind}= subplot (2,3,ind) ;
            axt{ind} = subplot (2,1,1) ; plot (t,data{ind})  ;hold all
            %                       plot (t_source,sigg,'r') ;
            xlabel('Time [s]');ylabel('g');
            %             [data_fft_abs , f ,~] =  myfft(data{ind},Fs,1);
            %             axf{ind} = subplot (2,3,ind+3) ; %plot (f,data_fft_abs) ; hold all
            axf{ind} = subplot (2,1,2) ;
            [data_fft_abs , f ,~] =  myfft(data{ind},Fs,hann_win);
            %             plot (f,data_fft_abs,'Color',[0 , 113 , 188 ]/255) ; hold al
            plot (f,data_fft_abs) ; hold all
            
            xlabel('Frequency [Hz]');ylabel('g');
            
        end
        
        % plot and save it.
        subplot (2,1,1) ; title(strrep(folder_test,'_',' ')) ; legend('vib0','vib1','vib2'); %xlim([1 1.3])
        subplot (2,1,2) ;  legend('source','vib0','vib1','vib2'); xlim([0 10e3])
        
        linkaxes([axt{1} ,axt{2} ,axt{3} ],'xy')
        linkaxes([axf{1} ,axf{2} ,axf{3} ],'xy')
        
        % % <<<<<<<<<<<<<< save fig >>>>>>>>>>>>>>>>> %
%         savefig(f1,[ folder_root folder_test '_' num2str(hann_win) '.fig']);
%         close(f1)
        % % <<<<<<<<<<<< end save fig >>>>>>>>>>>>>>> %
        
        results.(test_name).vib012_length =  N;
        
        if sum(diff(N)) ~= 0
            warning('^^^^^^^^ data size done match ^^^^^^^^^^^^^')
            for indN = 1:numel(N)
                data{indN} = [  data{indN} ; zeros( max(N)-N(indN) ,1)] ; % pad zeros to make data same size 
                N(indN) = max(N);
            end
            N = numel(data{1});
        end
        
        
        %% <read data from HW>
        
        for ind = 1:3
            data_acc = data{ind};
            [data_acc_fft , f] = myfft(data_acc,Fs,hann_win);
            
            % <test case 1: data length, make sure that no sample points were lost >
%             if abs(rec_time*Fs - N(ind))/(rec_time*Fs) < 0.0001 % diff is less then 0.01%  % rec_time == round(N(ind)/Fs,4)  % rounds to 100 microseconds
            if abs(rec_time*Fs - N(ind)) < 256/2 % the diff is less then the pading that can be applyed 
                results.(test_name).(['vib' num2str(ind)]).test_case_1_data_length =  'Pass'; % pass
            else
                results.(test_name).(['vib' num2str(ind)]).test_case_1_data_length =  'Fail'; % fail
            end
            results.(test_name).(['vib' num2str(ind)]).test_case_1_rec_time = numel(data_acc)/Fs ;
            
            
            % <test case 2: frequency accuracy, make sure frequency is with in 0.25 of intenfed intended value >
            %%% calc fft from data
                        
            %%% find peaks frequencys
            [data_pks,locs] = findpeaks(data_acc_fft,'MinPeakProminence',noise_floor); % need to set val of MinPeakProminence
            [sigg_pks,sigg_locs] = findpeaks(sigg_fft,'MinPeakProminence',noise_floor); % need to set val of MinPeakProminence
            
            try
                % check that peaks are at the the correct amplitude and frequency
                f_res = Fs/numel(data_acc);
                ff_source = f_source(f_source > h_filter);
                for f_ind =1:numel(ff_source)
                    % check frequency accuracy
                    if  f(locs(f_ind))<=(ff_source(f_ind)+f_res) & (ff_source(f_ind)-f_res) <= f(locs(f_ind))
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_2_frequency_accuracy_' num2str(ff_source(f_ind)) ]) = 'Pass';
                    else
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_2_frequency_accuracy_' num2str(ff_source(f_ind)) ]) = 'Fail';
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_2_frequency_accuracy_val' num2str(ff_source(f_ind)) ]) = ...
                            [ 'measured value = ' num2str(abs(f(locs(f_ind))-ff_source(f_ind))) ' needs to be less then ' num2str(f_res) ];
                    end
                    results.(test_name).(['vib' num2str(ind)]).('test_case_2_frequency_accuracy_val')(f_ind,:) = [ ff_source(f_ind) , abs(f(locs(f_ind))-ff_source(f_ind)) ]; 

                    % check amplitude accuracy

                    deltaP = (data_pks(f_ind)./sigg_pks(f_ind));
                    if deltaP>1-amp_accuracy & deltaP<1+amp_accuracy
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_3_amplitude_accuracy_' num2str(ff_source(f_ind)) ]) = 'Pass';
                    else
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_3_amplitude_accuracy_' num2str(ff_source(f_ind)) ]) = 'Fail';
                        results.(test_name).(['vib' num2str(ind)]).(['test_case_3_amplitude_accuracy_present_' num2str(ff_source(f_ind)) ]) = ...
                            [ 'measured diff present = ' num2str(abs(1-deltaP)*100,'%0.2f') '% needs to be less then ' num2str(amp_accuracy*100,'%0.2f') '%' ];
                    end
                    results.(test_name).(['vib' num2str(ind)]).('test_case_3_amplitude_accuracy_present')(f_ind,:) = [ ff_source(f_ind) , abs(1-deltaP)]; 
                end
                
                % check for noise level
                if  numel(locs) == 1
                    n  = round(FeatherHz/f_res); %n index of FeatherHz
                    try % make sure the right window doesn't exceed the length of the vector
                        winR = data_acc_fft( locs+n:locs+n+round(N*0.005));
                    catch
                        winR = data_acc_fft(locs+n:end);
                    end
                    try % make sure the left window shorter than the vector
                        winL = data_acc_fft( locs-n-round(N*0.005):locs-n);
                    catch
                        winL = data_acc_fft(round(6/f_res):locs-n);
                    end
                    results.(test_name).(['vib' num2str(ind)]).noise_mean_stdev_min_max =  ...
                        [ mean([winL; winR]) , std([winL; winR]) , min([winL; winR]) , max([winL; winR])];
                else % check the average noise in between the peaks
                    for ind_peak =1:numel(locs)-1
                        n = locs(ind_peak+1)- locs(ind_peak);%n index between peaks
                        win = data_acc_fft( locs(ind_peak) + round(n*0.1): locs(ind_peak+1) - round(n*0.1) ); 
                        results.(test_name).(['vib' num2str(ind)]).noise_mean_stdev_min_max =  ...
                            [ mean(win) , std(win) , min(win) , max(win)];
                    end
                end
                noiseAVG(ind,:) = results.(test_name).(['vib' num2str(ind)]).noise_mean_stdev_min_max;
                peakind = 1 ; 
            catch
                disp([folder_test ' : failed to check peaks : found-' num2str(numel(locs)) ' expected- ' num2str(numel(ff_source))])
                peakind = 0 ; 
                continue
            end
        end
        % mean noise
        if test_ID == 0 
            test_ID = 1;
        end
        if  peakind == 1
            MnoiseAVG(test_ID,:) = mean(noiseAVG, 1);
            results.(test_name).noise_mean_stdev_min_max  = MnoiseAVG(test_ID,:); % mean noise per test
        end
        %         figure; plot (t, abs(data{1}-data{2}) ,t, abs(medfilt1(data{1}-data{3})) ) ;
        if test_ID ~= 4
            results.(test_name).diff_01_mean_stdev_max = [mean(abs(data{1}-data{2})), std(abs(data{1}-data{2})), max(abs(data{1}-data{2}))];
            results.(test_name).diff_02_mean_stdev_max = [mean(abs(data{1}-data{3})), std(abs(data{1}-data{3})), max(abs(data{1}-data{3}))];
            results.(test_name).diff_12_mean_stdev_max = [mean(abs(data{2}-data{3})), std(abs(data{2}-data{3})), max(abs(data{2}-data{3}))];
            t = [results.(test_name).diff_01_mean_stdev_max ; results.(test_name).diff_02_mean_stdev_max ;...
                results.(test_name).diff_12_mean_stdev_max];
            m(test_ID,:) = mean(t,1);
        end
    end
    results.diff012_mean_stdev_max = mean(m,1);
    results.noise_mean_stdev_min_max = mean(MnoiseAVG,1);

end
% save([ folder_root  'results_' num2str(hann_win) '.mat'],'results');
