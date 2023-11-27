function RemoveIED_Spidetection
% USAGE

% Given a SPI and IED res file for a session, omit IEDs detected as SPIs

SPI_files =dir('*spindles.mat');
IED_files =dir('*_IED.mat');

Rs=1250;
CH_key=dir('*CH_key.mat');
load (CH_key.name)
CH_N=length (CH_key);

SPI_l=length (SPI_files);
IED_l=length (IED_files);

N=max([SPI_l IED_l]);
%%
    for i=1:N
        load (SPI_files(i).name);

        filename=SPI_files(i).name;
        CH_index1 =find(filename=='_', 1, 'last' );
        CH_index2 = strfind(filename,'spindles');

        CH=  str2num(filename (CH_index1+1:CH_index2-1)); 
        
        IED_filename=cat(2,filename(1:CH_index1),'Ch',num2str(CH),'33_IED.mat');
        if ~isempty (dir(IED_filename))
            load (IED_filename);
            [r, c] = size(IED);
            [q, a] = size(spindles);
            res = spindles;

            for ii = 1:q
                for j = 1:r
                    if IED(j, 1) < spindles(ii, 3) &&  IED(j, 1) >   spindles(ii, 1)     %test for overlap with 100ms buffer for detection issues
                        res(ii, :) = [0 0 0 0];
                    end
                end
            end

            spindles_only_index = find(res(:, 1) ~= 0);
            spindles_only = res(spindles_only_index, :);

            spindles_only_file = strcat(filename(1:CH_index1), num2str(CH), '_spionly');   %riponly
            save (spindles_only_file, 'spindles_only') 

            spi_events = strcat(filename(1:CH_index1), num2str(CH), '_spionly','.spe.evt');       % .rik.evt
            SaveRippleEvents(spi_events,spindles_only,1); 

        else
            spindles_only = spindles;

            spindles_only_file = strcat(filename(1:CH_index1), num2str(CH), '_spionly');   %riponly
            save (spindles_only_file, 'spindles_only')


            spi_events = strcat(filename(1:CH_index1), num2str(CH), '_spionly','.spe.evt');       % .rik.evt
            SaveRippleEvents(spi_events,spindles_only,1); 

        end
    end
end