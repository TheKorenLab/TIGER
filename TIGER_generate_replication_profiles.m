%% Extract reads from BAM files

% run on bam files (typically on a server):


%  ls *.bam  and insert list below (for multiple bam files, loop on files)
%  for filename in [bam file list]; do echo $filename; mkdir sam_$filename; for i in {1..22} X Y; do echo $i; samtools view -F 1024 -F 256 -F 128 -q 10 $filename ${i} | cut -f 4 > ./sam_$filename/readsChr${i}.txt; done; done



%  OR, for cases in which the reads are labeled "chrXX", e.g. hg38, mm10 (for mouse, chromosomes are {1..19} X):
%  for filename in [bam file list]; do echo $filename; mkdir sam_$filename; for i in {1..22} X Y; do echo $i; samtools view -F 1024 -F 256 -F 128 -q 10 $filename chr${i} | cut -f 4 > ./sam_$filename/readsChr${i}.txt; done; done

%  OR, for drosophila:
%  for filename in [bam file list]; do echo $filename; mkdir sam_$filename; for i in 2L 2R 3L 3R 4 X Y; do echo $i; samtools view -F 1024 -F 256 -F 128 -q 10 $filename chr${i} | cut -f 4 > ./sam_$filename/readsChr${i}.txt; done; done




%  256: not primary alignment; 1024: PCR duplicate; 128: second in pair




%% Import reads to Matlab and filter for alignability
% code is designed for multiple samples run in a loop as above; for single samples, modify accordingly
% make sure the TIGER folder is in the path

clear;clc

%cd    % ! cd to the directory with the txt file folders ! %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TIGER_folder = '/TIGER';  % ! change as appropriate

genome_build = 'hg19'; % hg19, hg38, mm10, dm6
read_length = 100; % refers to the alignability filter file, not the actual lengths of the reads (read lengths should be same or longer than those used for alignability file). Default = 100
win_size = 10000; % size of the windows in uniquely-alignable bps. Default = 1000


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp/Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp'])
eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp/' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)])


switch genome_build
    case 'mm10'
        Wins=Wins(1:20); % no unique sequences on the Y chromosome
end



% find all sam_* folders (each one is a sample). "sam_" is the notation of the file outputs from the samtools reads command (see section above) 
files = dir;
names = {};
for i = 1:size(files,1)
    file_names = files(i).name;
    if length(file_names)>3 && strcmp(file_names(1:3),'sam')
        names = [names;file_names];
    end
end
clear files file_names i 


Reads = cell(size(Wins,1),length(names));
Windowed_data = cell(size(Wins,1),1);
data_to_segment = cell(size(Wins,1),1);
for samp = 1:length(names)
    eval(['cd ' names{samp}])
    
    for Chr = 1:size(Wins,1)
        disp(['Importing, sample ' num2str(samp) ', Chr ' num2str(Chr)]);
        filename = cell2mat(['readsChr' chrnum(Chr,genome_build) '.txt']);
        fid = fopen(filename);
        Coordinates = textscan(fid,'%f');  % % read coordinate
        ReadsT = cell2mat(Coordinates);
        in = ismember(ReadsT,Coordinates_to_remove{Chr}); % remove reads flagged as not uniquely alignable
        ReadsT(in,:) = [];
        Reads{Chr,samp} = ReadsT;
        fclose('all');
        
        
        % count reads in windows
        if samp==1
            Windowed_data{Chr}(:,1:4) = Wins{Chr};
        end
        win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
        read_count  =  histcounts(Reads{Chr,samp},win_coord);
        Windowed_data{Chr}(:,4+samp) = read_count';
            
    end
    
    % prepare data for segmentation
    mean_coverage = cell2mat(Windowed_data(1:TIGER_last_autosome(genome_build)));   mean_coverage = nanmean(mean_coverage(:,samp+4));
    for Chr = 1:size(Wins,1)
        data_to_segment{Chr}(:,1:4) = Windowed_data{Chr}(:,1:4);
        data_to_segment{Chr}(:,samp+4) = Windowed_data{Chr}(:,samp+4)./mean_coverage.*2; % normalize to CN=2 (read number per window is usually much higher and will fail segment)
    end

    
    names{samp} = names{samp}(5:end-4);
    
    cd ..
    
    clear ReadsT Coordinates filename fid in Chr win_coord read_count mean_coverage 
    
end
params.genome_build = genome_build;
params.read_length = read_length;
params.win_size = win_size;
params.TIGER_folder = TIGER_folder;

% segment
[filt_in, data_in_segs_corrected, segment_params] = TIGER_segment_filt(data_to_segment,[],[],params.genome_build); % ** parfor requires installing the Parallel Computing Toolbox. Go to Home-->Add-Ons
params.segment_params = segment_params;




fold_name = pwd;  a = regexp(fold_name,'/');  fold_name = fold_name(a(end)+1:end);
eval(['save -v7.3 ' fold_name  '_'   params.genome_build '_'   num2str(params.read_length)   'bp_filtered_reads_'  num2str(params.win_size/1000) 'Kb  names  Reads  Windowed_data  params  filt_in'])
clear a fold_name Chr Coordinates_to_remove genome_build read_length samp data_to_segment Wins






%%%% QC figures auto-generated:
%%%% whole-genome plot- check for aneuploidies etc;
%%%% evaluating default segmentation parameters- adjust in section below if needed;
%%%% read counts per 1Kb window- evaluate coverage and adjust window size in section below if needed %%%%



Data_filt_segments = Windowed_data;
for samp = 1:length(names)
    for Chr = 1:size(Windowed_data,1)
        Data_filt_segments{Chr}(filt_in{Chr}{samp},samp+4) = NaN;
    end
end


% read counts per window and segment filtering - whole genome plot continuous
Loc{1,1} = Windowed_data{1}(:,4);
for Chr = 2:size(Windowed_data,1) 
    Loc{Chr,1} = Loc{Chr-1}(end)+Windowed_data{Chr}(:,4);
end

for samp = 1:length(names)
    figure;hold on
    set(gca,'fontsize',20)
    set(gcf,'position',[75         175        2500        1150])
    subplot(4,1,1); hold on; set(gca,'fontsize',16)
    bords = [];
    for Chr = 1:size(Windowed_data,1)
        plot(Loc{Chr,1}./1e6,Windowed_data{Chr}(:,samp+4),'b.')
        plot(Loc{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),1)./1e6,Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),samp+4),'g.')
        %plot(Loc{Chr,1}(:,1)./1e6,data_in_segs_corrected{Chr,1}(:,samp).*al./2,'k.')   % optional- plot the segments
        bords = [bords;Loc{Chr,1}(end)./1e6];
    end
    axis tight
    set(gca,'xtick',[],'xticklabel',[])  
    mean_coverage = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));     mean_coverage = nanmean(mean_coverage(:,samp+4));
    set(gca,'ylim',[0 mean_coverage*4])
    X = get(gca,'xlim');
    Y = get(gca,'ylim');
    plot(bords*ones(1,2),Y,'k-')
    set(gca,'xtick',[],'xticklabel',[])
    bords = [Windowed_data{1,1}(1,4)./1e6;bords];
    bords = bords(1:end-1)+diff(bords)/2;
    text(bords-20,(-Y(2)/20)*ones(size(bords)),chrnum((1:size(Windowed_data,1))'),'fontsize',16) % may need to be adjusted for different genomes
    ylabel('Number of reads');
    title([names{samp} ';  Whole-genome;  segment R2 = ' num2str(params.segment_params(1)) ', window size = ' num2str(params.win_size)])
    
    for Chr = 2:4
        subplot(4,1,Chr); hold on; set(gca,'fontsize',16)
        plot(Windowed_data{Chr}(:,4)./1e6,Windowed_data{Chr}(:,samp+4),'.') % original data
        plot(Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),4)./1e6,Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),samp+4),'.') % data after segment-filtering
        plot(Windowed_data{Chr}(:,4)./1e6,data_in_segs_corrected{Chr}(:,samp).*mean_coverage./2,'k.','markersize',1)  % segments
        axis tight
        set(gca,'ylim',[0 mean_coverage*3])
        if Chr==4
            xlabel('Chromosome coordinate (Mb)')
        end
        ylabel('Number of reads')
        title(['Segmentation, Chromosome ' num2str(Chr) ])
    end

    
end
clear bords Loc X Y mean_coverage 



% distribution of read counts per window
figure;hold on
set(gca,'fontsize',20)
set(gcf,'position',[167         145        2189        1189])
all_reads = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));
numrows = ceil(sqrt(length(names)));
numcols = ceil(length(names)/ceil(sqrt(length(names))));
for samp = 1:length(names)
    subplot(numcols,numrows,samp);hold on; set(gca,'fontsize',16)
    mean_coverage = nanmedian(all_reads(:,samp+4));
    [a b] = histcounts(all_reads(:,samp+4),0:mean_coverage*4);
    bar(b(1:end-1)+1,a/1000)
    set(gca,'ylim',[0 max(a/1000)*1.5])
    text(0.04,0.92,['Median reads per window  =  ' num2str(round(mean_coverage))],'units','normalized','fontsize',16)
    text(0.04,0.78,['Total reads (all chromosomes)  =  ' num2str(nansum(all_reads(:,samp+4))/1e6,3) 'M'],'units','normalized','fontsize',16)
    if ceil(samp/ceil(sqrt(length(names))))==numcols
        xlabel('Number of reads per window')
    end
    if rem(samp,numrows)==1
        ylabel([{'Number of windows'} ; {'(thousands)'}])
    end
    title([names{samp} ', window size = ' num2str(win_size)])
end
clear a b all_reads samp 





% !! then delete sam_ folders with .txt files from computer and from server !!





%% % OPTIONAL: evaluate alternative segmentation parameters and window size
% same as above but without importing the data again; loads filtered reads
clear;clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R2 = []; % optional- change segmentation parameter, e.g. for low-coverage libraries. Higher parameter is longer segments; functions become problematic at high values, e.g. >0.1
seg_std=[ ]; % optional- change standard deviation threshold for filtering segments. Default=1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



files = dir;
file_load = {};
for i = 1:size(files,1)
    file_names = files(i).name;
    if ~isempty(regexp(file_names,'_filtered_reads'))
        file_load = [file_load;file_names];
    end
end
if length(file_load)>1
    disp('More than one file! Manually load the desired one')
else
    eval(['load ' file_load{1} ' names params Reads'])
end
clear files file_names i file_load



win_size=params.win_size; % size of the windows in uniquely-alignable bps

% load window coordinates
eval(['load ' params.TIGER_folder '/Alignability_and_GC_filters/' params.genome_build '_' num2str(params.read_length) 'bp/' params.genome_build '_' num2str(params.read_length) 'bp_wins_'   num2str(win_size)])

if strcmp(params.genome_build,'mm10')
    Wins = Wins(1:20); % no unique sequences on the Y chromosome
end



Windowed_data = cell(size(Wins,1),1);
for samp = 1:length(names)
    % count reads in windows
    for Chr = 1:size(Wins,1)
        if samp==1
            Windowed_data{Chr}(:,1:4) = Wins{Chr};
        end
        win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
        read_count  =  histcounts(Reads{Chr,samp},win_coord);
        Windowed_data{Chr}(:,4+samp) = read_count';
    end
    
    % prepare data for segmentation
    mean_coverage = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));   mean_coverage = nanmean(mean_coverage(:,samp+4));
    for Chr = 1:size(Wins,1)
        data_to_segment{Chr,1}(:,1:4) = Windowed_data{Chr}(:,1:4);
        data_to_segment{Chr}(:,samp+4) = Windowed_data{Chr}(:,samp+4)./mean_coverage.*2; % normalize to CN=2 (read number per window is usually much higher and will fail segment)
    end

end


% segment
[filt_in, data_in_segs_corrected, segment_params] = TIGER_segment_filt(data_to_segment,R2,seg_std,params.genome_build);
params.segment_params = segment_params;

Data_filt_segments = Windowed_data;
for samp = 1:length(names)
    for Chr = 1:size(Windowed_data,1)
        Data_filt_segments{Chr}(filt_in{Chr}{samp},samp+4) = NaN;
    end
end




    
% read counts per window and segment filtering - whole genome plot continuous
Loc{1,1} = Windowed_data{1}(:,4);
for Chr = 2:size(Windowed_data,1) % update
    Loc{Chr,1} = Loc{Chr-1}(end)+Windowed_data{Chr}(:,4);
end

for samp = 1:length(names)
    figure;hold on
    set(gca,'fontsize',20)
    set(gcf,'position',[75         175        2500        1150])
    subplot(4,1,1); hold on; set(gca,'fontsize',16)
    bords = [];
    for Chr = 1:size(Windowed_data,1)
        plot(Loc{Chr,1}./1e6,Windowed_data{Chr}(:,samp+4),'b.')
        plot(Loc{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),1)./1e6,Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),samp+4),'g.')
        %plot(Loc{Chr,1}(:,1)./1e6,data_in_segs_corrected{Chr,1}(:,samp).*al./2,'k.')
        bords = [bords;Loc{Chr,1}(end)./1e6];
    end
    axis tight
    set(gca,'xtick',[],'xticklabel',[])  % ,'ydir','reverse','ylim',[0 100]
    mean_coverage = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));     mean_coverage = nanmean(mean_coverage(:,samp+4));
    set(gca,'ylim',[0 mean_coverage*4])
    X = get(gca,'xlim');
    Y = get(gca,'ylim');
    plot(bords*ones(1,2),Y,'k-')
    set(gca,'xtick',[],'xticklabel',[])
    bords = [Windowed_data{1,1}(1,4)./1e6;bords];
    bords = bords(1:end-1)+diff(bords)/2;
    text(bords-20,(-Y(2)/20)*ones(size(bords)),chrnum((1:size(Windowed_data,1))'),'fontsize',16)
    ylabel('Number of reads');
    title([names{samp} ';  Whole-genome;  segment R2 = ' num2str(params.segment_params(1)) ', window size = ' num2str(params.win_size)])
    
    for Chr = 2:4
        subplot(4,1,Chr); hold on; set(gca,'fontsize',16)
        plot(Windowed_data{Chr}(:,4)./1e6,Windowed_data{Chr}(:,samp+4),'.') % original data
        plot(Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),4)./1e6,Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),samp+4),'.') % data after segment-filtering
        plot(Windowed_data{Chr}(:,4)./1e6,data_in_segs_corrected{Chr}(:,samp).*mean_coverage./2,'k.','markersize',1)  % segments
        axis tight
        set(gca,'ylim',[0 mean_coverage*3])
        if Chr==4
            xlabel('Chromosome coordinate (Mb)')
        end
        ylabel('Number of reads')
        title(['Segmentation, Chromosome ' num2str(Chr) ])
    end

    
end
clear bords Loc X Y mean_coverage 


% plot other chromosomes (optional)
samp = 1;
Chr = 2;
figure;hold on
set(gca,'fontsize',20)
set(gcf,'position',[197         499        2222         667])
mean_coverage = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));     mean_coverage = nanmean(mean_coverage(:,samp+4));
plot(Windowed_data{Chr}(:,4)./1e6,Windowed_data{Chr}(:,samp+4),'.') % original data
plot(Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),4)./1e6,Windowed_data{Chr}(~isnan(Data_filt_segments{Chr}(:,samp+4)),samp+4),'.') % data after segment-filtering
plot(Windowed_data{Chr}(:,4)./1e6,data_in_segs_corrected{Chr}(:,samp).*mean_coverage./2,'k.','markersize',1)  % segments
axis tight
set(gca,'ylim',[0 mean_coverage*3])
xlabel('Chromosome coordinate (Mb)')
ylabel('Number of reads')
title([names{samp} '; Segmentation, Chromosome ' num2str(Chr) ])




    
eval(['save -v7.3 ' fold_name  '_'   params.genome_build '_'   num2str(params.read_length)   'bp_filtered_reads_'  num2str(params.win_size/1000) 'Kb_modified_seg_params  names  Reads  Windowed_data  params  filt_in'])

clear fold_name Chr Coordinates_to_remove genome_build read_length samp data_to_segment Wins





%% GC-normalize
clear;clc


files = dir;
file_load = {};
for i = 1:size(files,1)
    file_names = files(i).name;
    if ~isempty(regexp(file_names,'_filtered_reads', 'once'))
        file_load = [file_load;file_names];
    end
end
if length(file_load)>1
    disp('More than one file! Manually load the desired one')
else
    eval(['load ' file_load{1}])
end
clear files file_names i file_load


eval(['load ' params.TIGER_folder '/Alignability_and_GC_filters/' params.genome_build '_' num2str(params.read_length) 'bp/' params.genome_build '_' num2str(params.read_length) 'bp_' num2str(params.win_size) 'bp_wins_Reads_expected_nominal_401'])   

    

figure;hold on
set(gca,'fontsize',20)
set(gcf,'position',[167         145        2189        1189])
numrows = ceil(sqrt(length(names)));
numcols = ceil(length(names)/ceil(sqrt(length(names))));

Count_expected_data = cell(size(Reads,1),length(names));
All_windows = cell2mat(Windowed_data(1:TIGER_last_autosome(params.genome_build)));
for samp = 1:length(names)
    
    % generate the GC filter- all coordinates in the segment-filtered regions (don't use outlier copy number windows for calculating GC% effects)
    disp([names{samp} ', extracting GC windows to exclude'])
    GC_chr = cell(TIGER_last_autosome(params.genome_build),1);
    for Chr = 1:TIGER_last_autosome(params.genome_build)
        in = filt_in{Chr}{samp};
        filt_bps = cell(length(in),1); 
        for i = 1:length(in)
            filt_bps{i} = (Windowed_data{Chr}(in(i),2):Windowed_data{Chr}(in(i),3))';
        end
        GC_chr{Chr} = cell2mat(filt_bps); % this is the GC filter
    end
    clear in i filt_bps
    
    
    for Chr = 1:TIGER_last_autosome(params.genome_build)
        Reads{Chr,samp} = Reads{Chr,samp}(~ismember(Reads{Chr,samp},GC_chr{Chr}));  % use these reads to calculate GC content effect (apply GC filter)
    end
    
    
    %%%% calculate GC bias factor %%%%%%
    % count reads in each GC content bin
    total_number_of_bps_per_GC_bin = zeros(401,1); % bp belonging to each GC bin, after removing alignability- and GC-filtered bps
    num_reads_per_GC_bin = zeros(401,TIGER_last_autosome(params.genome_build));
    for Chr = 1:TIGER_last_autosome(params.genome_build)
        disp([names{samp}, ', finding reads in GC bins, chromosome ' num2str(Chr)])
        % load coordinates belonging to each GC content bin
        eval(['load ' params.TIGER_folder '/Alignability_and_GC_filters/' params.genome_build '_' num2str(params.read_length) 'bp/' params.genome_build '_' num2str(params.read_length) 'bp_GC_cont_401_filtered_bins_chr' num2str(Chr)])
        
        % filter GC bins 
        GC_cont_win_chr_bins_filtered = cell(401,1);
        GC_bin_lengths = cellfun(@numel,GC_cont_win_chr_bins);
        GC_bin_lengths = [0;cumsum(GC_bin_lengths)];
        GC_all = cell2mat(GC_cont_win_chr_bins);
        [~, in] = intersect(GC_all,GC_chr{Chr});
        GC_all(in) = NaN;
        for i = 1:401
            GC_cont_win_chr_bins_filtered{i} = GC_all(GC_bin_lengths(i)+1:GC_bin_lengths(i+1));
            GC_cont_win_chr_bins_filtered{i}(isnan(GC_cont_win_chr_bins_filtered{i})) = [];
        end
        
        Reads_chr = Reads{Chr,samp};
        % process GC bins by order of their abundance, to save time (code removes data as it goes)
        win_len = cellfun(@length,GC_cont_win_chr_bins);
        [~, win_len_sorted] = sort(win_len,'descend');
        for i = win_len_sorted'  % bit faster than 1:401
            Reads_in_bin = find(ismember(Reads_chr,GC_cont_win_chr_bins_filtered{i}));
            num_reads_per_GC_bin(i,Chr) = length(Reads_in_bin); % number of reads in each bin, used for the GC bias calculation
            Reads_chr(Reads_in_bin) = [];  % saves time to reduce the size of R in each iteration- minimize the search space for a
        end
        
        total_number_of_bps_per_GC_bin = total_number_of_bps_per_GC_bin+cellfun(@numel,GC_cont_win_chr_bins_filtered); % total number of bps in each GC bin, used to normalize the read counts per GC bin
        
    end
    clear GC_cont_win_chr_bins i a R GC_cont_win_chr_bins_filtered f fs  GC_chr in GC_all GC_bin_lengths
    
    subplot(numcols,numrows,samp); hold on; set(gca,'fontsize',16)
    GC_bias = sum(num_reads_per_GC_bin')./total_number_of_bps_per_GC_bin';
    plot(80:320,GC_bias(80:320)) % GC bins below 20% or above 80% are excluded
    title(names{samp})
    if ceil(samp/ceil(sqrt(length(names))))==numcols
        xlabel('GC content')
    end
    if rem(samp,numrows)==1
        ylabel([{'Number of reads'}; {'(normalized by GC bin abundance)'}])
    end
    
    
    gcgenome = total_number_of_bps_per_GC_bin/sum(total_number_of_bps_per_GC_bin); % fraction of bps in each bin
    achr = sum(num_reads_per_GC_bin')'; % add number of reads in each bin for all chromosomes
    GC_fractions = (achr/sum(achr))./gcgenome; % this is the GC bias factor- fraction of read in each GC bin divided by fraction of bps in each bin


    % correction factor- to bring copy number to median of 2
    corr_factor = nanmedian(All_windows(:,samp+4))./params.win_size./2;  %  how many reads expected per bp in a window (median number of reads per windows). ./2: to make the genome diploid 

    for Chr = 1:size(Reads,1)
        % in each read number window, expected number of bps in each GC bin given GC bias factor
        sample_Reads_expected_nominal = NaN(size(Reads_expected_nominal{Chr},1),321);
        for i = 81:321
            sample_Reads_expected_nominal(:,i) = Reads_expected_nominal{Chr}(:,i)*GC_fractions(i);
        end
        sample_Reads_expected_nominal(:,402) = nansum(sample_Reads_expected_nominal(:,81:321)'); % total number of reads expected per window

        Expected = sample_Reads_expected_nominal(:,402).*corr_factor; % actual expected number of reads, considering genome coverage of the library
        Count_expected_data{Chr,samp}(:,1) = Windowed_data{Chr,1}(:,samp+4);
        Count_expected_data{Chr,samp}(:,2) = Expected;
        Windowed_data{Chr}(:,samp+4) = Windowed_data{Chr,1}(:,samp+4)./Expected; % this is normalized read counts (count/expected), used for subsequent analyses

        in_bad_unexpected = Reads_expected_nominal{Chr}(:,402)~= params.win_size; % remove windows with unexpected total read counts of GC bins (mostly found in low-mappability regions)
        Windowed_data{Chr}(in_bad_unexpected,samp+4) = NaN;
        in_zero_expected = find(Expected==0); % remove windows with zero expected reads (result in inf after normalization)
        Windowed_data{Chr}(in_zero_expected,samp+4) = NaN;
    end
    clear  total_read_counts_per_GC_bin Chr i  corr_factor total_number_windows_used total_reads_used  GC_bias co GC_fractions achr gcgenome sample_Reads_expected_nominal Expected in_bad_unexpected num_reads_per_GC_bin total_number_of_bps_per_GC_bin

end



fold_name = pwd;  a = regexp(fold_name,'/');  fold_name = fold_name(a(end)+1:end);
eval(['save -v7.3 ' fold_name  '_'   params.genome_build '_'  num2str(params.read_length) 'bp_GC_corrected_' num2str(params.win_size/1000) 'Kb  Windowed_data Count_expected_data names params '])



%% Filter and smooth 
clear;clc

P_interpolation  =  10^-17; % lower power number is less smooth (i.e. 16 is less smooth than 17). Default: 10^-17


files = dir;
file_load = {};
for i = 1:size(files,1)
    file_names = files(i).name;
    if ~isempty(regexp(file_names,'bp_GC_corrected_', 'once'))
        file_load = [file_load;file_names];
    end
end
if length(file_load)>1
    disp('More than one file! Manually load the desired one')
else
    eval(['load ' file_load{1}])
end
clear files file_names i file_load

params.P_interpolation = P_interpolation;


switch params.genome_build
    case 'hg19'
        load Gap_hg19 Gap
    case 'hg38'
        load Gap_hg38 Gap 
    case 'mm10'
        load Mouse_Gap Gap
    case 'dm6'
        load dm6_Gap Gap
end




% segment-filter 
[filt_in, data_in_segs_corrected, segment_params]=segment_filt_v2(Windowed_data,params.segment_params(1),params.segment_params(2),params.genome_build);

Data_filt_segments = Windowed_data;
for samp = 1:length(names)
    for Chr = 1:size(Windowed_data,1)
        Data_filt_segments{Chr}(filt_in{Chr}{samp},samp+4) = NaN;
    end
end


% plot unfiltered and filtered data
Chr_show = [2 5 8]; % adjust for different genomes
for c = 1:min(20,length(names))
    figure;hold on
    set(gcf,'position',[750   402   941   667 ]); 
    set(gca,'fontsize',20)
    for Chr = Chr_show
        [~, in] = find(Chr==Chr_show);
        subplot(length(Chr_show),1,in); hold on; set(gca,'fontsize',16)
        plot(Windowed_data{Chr}(:,4)./1e6,Windowed_data{Chr}(:,c+4),'.','color',[0.4 0.4 0.6],'markersize',1)
        plot(Windowed_data{Chr}(:,4)./1e6,Data_filt_segments{Chr}(:,c+4),'.','color','g','markersize',1)
        %plot(Windowed_data{Chr}(:,4)./1e6,data_in_segs_corrected{Chr}(:,c),'k.')  % segments
        axis tight
        set(gca,'ylim',[0 8])
        title(cell2mat([names(c) ', chromosome ' num2str(Chr)]))
    end
end



% smooth between gaps >= 50Kb and regions separated by >= 100Kb (ignore gaps <50Kb (only for human and mouse))
Smoothed_data = Data_filt_segments;
num_window_thresh = 20; % only smooth regions with >= 20 windows
gap_thresh = 1e5;
for c = 5:size(Smoothed_data{1},2)
    for Chr = 1:size(Windowed_data,1)
        % don't smooth across windows separated by >= 100Kb
        in_not_nan = find(~isnan(Smoothed_data{Chr}(:,c)));
        in1 = find(diff(Smoothed_data{Chr}(in_not_nan,4))>= gap_thresh);
        data_gaps = [Smoothed_data{Chr}(in_not_nan(in1),3)+1 Smoothed_data{Chr}(in_not_nan(in1+1),2)-1]';
        
        % ignore gaps <50Kb (only human and mouse)
        G = Gap(Gap(:,1)==Chr,2:3)';
        if ismember(params.genome_build,[{'hg19'};{'hg38'};{'mm10'}])
            G(:,diff(G)<50e3) = []; 
        end
        G = sort([ Smoothed_data{Chr}(1,2)-1; G(:);  Smoothed_data{Chr}(end,3); data_gaps(:)]);
        G = [G(1:2:end) G(2:2:end)];
        
        for frag = 1:size(G,1)
            in = find( Smoothed_data{Chr}(:,2)>G(frag,1)  &  Smoothed_data{Chr}(:,3)<= G(frag,2)  &  ~isnan(Smoothed_data{Chr}(:,c)) );
            if ~isempty(in) && length(in)>= num_window_thresh
                F_func   =  csaps(Smoothed_data{Chr}(in,4),Smoothed_data{Chr}(in,c),P_interpolation);
                Smoothed_data{Chr}(in,c) = fnval(Smoothed_data{Chr}(in,4),F_func);
            elseif ~isempty(in)
                Smoothed_data{Chr}(in,c) = NaN;
            end
        end
    end
end


% normalize to autosomal mean = 0 and std = 1 (note: for male samples, consider adding 1 to sex chromosomes before normalizing
M = cell2mat(Smoothed_data(1:TIGER_last_autosome(params.genome_build)));
for Chr = 1:size(Smoothed_data,1)
    for c = 5:size(Smoothed_data{Chr},2)
        Smoothed_data{Chr}(:,c) = (Smoothed_data{Chr}(:,c)-nanmean(M(:,c)))./nanstd(M(:,c));
    end
end



fold_name = pwd;  a = regexp(fold_name,'/');  fold_name = fold_name(a(end)+1:end);
eval(['save -v7.3 ' fold_name  '_'   params.genome_build '_'  num2str(params.read_length) 'bp_processed_' num2str(params.win_size/1000) 'Kb  Data_filt_segments   Smoothed_data   names  params'])  % data_in_segs_corrected Smoothed_data_unnorm


