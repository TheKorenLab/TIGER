%% Generate fastq with all possible sequence reads of a given length per a reference genome

% need a reference genome in Matlab (e.g. use fastaread to import data). Must be the same genome as used to align sequence data to. 
% save as "<reference_name>_sequence" in folder "Alignability_and_GC_filters"
% also save in the same file a variable "chr_length" that contains only the length in bps of each chromosome by order

% need a gap file in Matlab (e.g. convert a UCSC gap file for the same genome build into Matlab, save chromosome, start, end

% code generates a fastq file with all possible sequences ("reads") of a chosen length, not including gaps. Generates arbitrary read metadata


clear;clc

genome_build = 'hg19'; % hg19, hg38, mm10, dm6   % other genomes can also be used- need genome sequence and gap file, and possibly code tweeks depending on chromosome naming convention
read_length = 100; % should be no longer than the lengths of the reads in the sequencing data, but no need to make it longer than 100
TIGER_folder = '/TIGER';  % ! change as appropriate


try
   eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/']) 
catch
    eval(['mkdir ' TIGER_folder '/Alignability_and_GC_filters/']) 
    eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/']) 
end
    
eval(['mkdir ' genome_build '_' num2str(read_length) 'bp']) 
eval(['cd '  genome_build '_' num2str(read_length) 'bp']) 

% load Gap file; structure of file should be: chromosome, start, end
switch genome_build
    case 'hg19'
        load Gap_hg19 Gap
    case 'hg38'
        load Gap_hg38 Gap 
    case 'mm10'
        load Mouse_Gap Gap
    case 'dm6'
        load dm6_Gap Gap
end

line4{1} = repmat('A',1,read_length); % arbitrary metadata
for Chr = 1:max(Gap(:,1))
    disp(['Chromosome ' num2str(Chr)])
    
    % load genome sequence  (loading the entire genome every time and only keeping one chromosome takes more time but saves computer memory) (file needs to be in the path or same folder)
    eval(['load ' genome_build '_sequence Sequence']) 
    Use_sequence = Sequence{Chr};
    clear Sequence
    
    % remove gaps
    G = Gap(Gap(:,1)==Chr,2:3);
    for i = size(G,1):-1:1
        Use_sequence(G(i,1)+1:G(i,2)) = [];
    end
    clear G i
    

    % write fastq file
    filename = cell2mat(['Chr' chrnum(Chr,genome_build) '_' num2str(read_length) 'bp_R1.fastq']);
    fid = fopen(filename, 'w');
    line1txt = ['@NB551191:189:HTHCVBGX5:1:11101:' num2str(Chr) ':%d']; % arbitrary initial parameters, followed by chromosome and coordinate

    % generate sequence reads [ordered by location]
    lh = waitbar(0,['Chr' chrnum(Chr)]);  % chrnum coverts chromosome index number to naming convention (e.g. 23 to X)
    cnt = 0;
    for i = 1:read_length
        waitbar(i/read_length,lh)
        clear Fast_I
        seq_length_use = length(Use_sequence)-rem((length(Use_sequence)-i+1),read_length);
        seqs_temp = ( reshape( Use_sequence(i:seq_length_use),read_length,floor(seq_length_use/read_length) ) )';
        Seqs = cellstr(seqs_temp);

        seq_texts = compose(line1txt,(cnt+1:cnt+length(seqs_temp))');
        Fast_I(1:4:length(Seqs)*4,1) = seq_texts;
        Fast_I(2:4:length(Seqs)*4,1) = Seqs;
        Fast_I(3:4:length(Seqs)*4,1) = {'+'}; %line3; arbitrary
        Fast_I(4:4:length(Seqs)*4,1) = line4;
        
        fprintf(fid,'%s\n',Fast_I{:});
        
        cnt = cnt+length(seqs_temp);
    end
    close(lh)
    eval(['gzip ' filename]) % gzips the file - input for bwa is gzipped
    eval(['delete ' filename]) % delete the original
    fclose('all'); 
    clear fid Fast_I A line1 cnt seqs_temp seq_length_use seq_texts lh Seqs filename i
    
end


eval(['cd ' TIGER_folder])



%% Align fastq and extract alignment information
% This part typically performed on a computer server


% 1. copy fastq files to server

% 2. use BWA-MEM to align to reference genome with duplicates marked. Should be the exact same reference used to generate the fastq files. 

% 3. After alignment, extract the locations of reads that were uniquely aligned (change chromosomes and read length if needed):
% for i in {1..22} X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam ${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done
%  - or - (reads are labeled "chrXX", e.g. hg38, mm10)
% for i in {1..22} X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam chr${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done
%  - or - (drosophila)
% for i in 2L 2R 3L 3R 4 X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam chr${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done

% (-F 4: read unmapped; -F 16: read reverse strand (only the forward strand was written to fastq); -F 1024: PCR or optical duplicate; -q 1: only keep sequences with mapQ of 1 or more; -f 4: only extract the start position of the sequence)



% * can then delete fastq files from computer and from server *



%% Import alignment samtools output to Matlab, save coordinates to remove 
% this codes generates the alignability filter- list of coordinates that are not uniquely alignable
clear;clc

genome_build = 'hg19'; % hg19, hg38, mm10, dm6
read_length = 100; 
TIGER_folder = '/TIGER';  % ! change as appropriate

eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 


% load file with chromosome lengths
eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence chr_length'])
% remove "extra" chromosomes (e.g. unmapped contigs)
switch genome_build 
    case 'hg19'
        chr_length = chr_length(1:24);
    case 'hg38'
        chr_length = chr_length(1:24);
    case 'mm10'
        chr_length = chr_length(1:21);
end


for Chr = 1:length(chr_length)
    disp(['Chromosome ' num2str(Chr)])
    switch genome_build 
        case 'mm10'
            filename = cell2mat(['Chr' chrnum(Chr,'mm10') '_' num2str(read_length) 'bp_filt.txt']); % chrnum coverts chromosome index number to naming convention (e.g. 23 to X) 
        case 'dm6'
            filename = cell2mat(['Chr' chrnum(Chr,'dm6') '_' num2str(read_length) 'bp_filt.txt']);
        otherwise
            filename = cell2mat(['Chr' chrnum(Chr) '_' num2str(read_length) 'bp_filt.txt']);
    end
    fid = fopen(filename);
    Coordinates = textscan(fid,'%f');  % read coordinate
    fclose('all'); clear fid filename

    Coordinates_to_remove{Chr,1}(:,1) = setdiff(1:chr_length(Chr),Coordinates{1}); % all genomic coordinates that are not uniquely alignable
    clear Coordinates
end


eval(['save -v7.3 Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp Coordinates_to_remove'])


eval(['cd ' TIGER_folder])


%  * can then delete samtools files from server and from computer *



%% Define genomic DNA copy number windows from alignability filter coordinates to remove 
% This code defines the windows in which read numbers will be counted ("read number windows").
% smaller windows can be merged to larger ones (function "merge_windows") but not vice versa. Downside of using smaller windows is larger files, longer processing times and higher memory usage
clear;clc

genome_build = 'hg19'; % hg19, hg38, mm10, dm6
read_length = 100; 
win_size = 10000; % this will be the size of the windows in uniquely-alignable bps
TIGER_folder = '/TIGER';  % ! change as appropriate


eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 
eval(['load Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp']) % load alignability filter


% load file with chromosome lengths
eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence chr_length'])
switch genome_build
    case 'hg19'
        chr_length = chr_length(1:24); 
        load Gap_hg19 Gap
    case 'hg38'
        chr_length = chr_length(1:24); 
        load Gap_hg38 Gap 
    case 'mm10'
        chr_length = chr_length(1:21);
        load Mouse_Gap Gap
    case 'dm6'
        load dm6_Gap Gap        
end


for Chr = 1:length(chr_length)
    disp(['Chromosome ' num2str(Chr)])
    kept_coordinates = setdiff(1:chr_length(Chr),Coordinates_to_remove{Chr});
    Wins{Chr,1}(:,2) = kept_coordinates(1:win_size:end-win_size+1); % start 
    Wins{Chr,1}(:,3) = kept_coordinates(win_size+1:win_size:end)-1; % end
    Wins{Chr,1}(:,1) = Chr;
    Wins{Chr,1}(:,4) = round(mean(Wins{Chr}(:,2:3)')'); % window center
    
    % remove windows that span gaps
    G = Gap(Gap(:,1)==Chr,2:3);
    in = [];
    for i = 1:size(G,1)
        in = [in;find(Wins{Chr,1}(:,2)<= G(i,1) & Wins{Chr,1}(:,3)>= G(i,2))];
    end
    Wins{Chr}(in,:) = [];
    
end


eval(['save ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)  ' Wins'])


eval(['cd ' TIGER_folder])



%% GC content normalization (files for calculating "expected" number of reads) 
% Calculates GC content in 401 bp disregarding Ns
% Saves files ("xbp_GC_cont_401_filtered_bins_chrx") with all the genomic coordinates (after alignability filter) that belong to each of the 401 GC% bins
% Saves "Reads_expected_nominal"- number of bps in each GC bin in each read number window
clear;clc

genome_build = 'hg19'; % hg19, hg38, mm10, dm6
read_length = 100;
win_size = 10000; % only needed for the last part of the code
TIGER_folder = '/TIGER';  % ! change as appropriate


eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 
eval(['load Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp']) % load alignability filter

% load window coordinates
eval(['load ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)])
if strcmp(genome_build,'mm10')
    Wins = Wins(1:20); % no unique sequences on the Y chromosome
end

% load genome sequence
eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence Sequence'])



        
bins = 0:1/400:1;  
bins = round(bins*10000)/10000; % prevents cases in which the value of bins isn't exact
for Chr = 1:size(Wins,1)
    disp(['Chromosome ' num2str(Chr)])

    %%%%%%%   calculate GC content for all bps in the genome [% all alignable bp] %%%%%
    Use_sequence = Sequence{Chr};
    GC_cont_win_chr(:,1) = 201:length(Use_sequence)-200;

    % GC content for each bps
    GC_sequence = zeros(length(Use_sequence),1);
    GC_sequence(upper(Use_sequence)=='G' | upper(Use_sequence)=='C') = 1;
    GC_sequence_sum = cumsum(GC_sequence);
    GC_sequence_sum = [GC_sequence_sum(401); GC_sequence_sum(402:end)-GC_sequence_sum(1:end-401)];
    GC_sequence_sum = GC_sequence_sum-GC_sequence(201:end-200); % exclude the center position (the one actually tested)
    
    % number of non-Ns in each GC window
    Non_N = ones(length(Use_sequence),1);
    Non_N(upper(Use_sequence)=='N') = 0;
    Non_N_sum = cumsum(Non_N);
    Non_N_sum = [Non_N_sum(401);Non_N_sum(402:end)-Non_N_sum(1:end-401)];
    Non_N_sum = Non_N_sum-Non_N(201:end-200);
    
    GC_cont_win_chr(:,2) = NaN;
    in = find(Non_N_sum==400); % ignore GC windows that contain any Ns
    GC_cont_win_chr(in,2) = GC_sequence_sum(in)./Non_N_sum(in);


    %%%% filter for alignabilty  %%%%
    [~, in] = intersect(GC_cont_win_chr(:,1),Coordinates_to_remove{Chr});
    GC_cont_win_chr(in,:) = [];
    in = find(isnan(GC_cont_win_chr(:,2)));
    GC_cont_win_chr(in,:) = [];    
    

    % find the genomic coordinates that belong to each of the 401 GC% bins
    lh = waitbar(0,'wait');
    for i = 1:length(bins)
        waitbar(i/length(bins),lh)
        in = find(GC_cont_win_chr(:,2)==bins(i));
        GC_cont_win_chr_bins{i,1} = GC_cont_win_chr(in,1);
    end
    close(lh)
    
    eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_GC_cont_401_filtered_bins_chr' num2str(Chr) ' GC_cont_win_chr_bins'])
    



    %%%%%% count number of bps in each GC bin in each read number window- creates variable "Reads_expected_nominal"  %%%%%
    win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
    for GC_bin = 1:401
        Reads_expected_nominal{Chr,1}(:,GC_bin)  =  histcounts(GC_cont_win_chr_bins{GC_bin},win_coord); % in each TIGER window, how many bp are found in each GC bin (in the genome)
    end
    Reads_expected_nominal{Chr}(:,402) = nansum(Reads_expected_nominal{Chr}'); % this is for testing purposes only. Gets overwritten later. column 402 should be the same as the win_size for most or all windows
    
    clear GC_cont_win_chr GC_cont_win_chr_bins
end


eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_' num2str(win_size) 'bp_wins_Reads_expected_nominal_401 Reads_expected_nominal'])   


eval(['cd ' TIGER_folder])



%% Optional: calculate Reads_expected_nominal for additional windows sizes
% For processing files with non-default read number windows, can run just this (GC% coordinate files stay the same) 
clear;clc

genome_build = 'hg19'; % hg19, hg38, mm10, dm6
read_length = 100;
win_size = 10000; 
TIGER_folder = '/TIGER';  % ! change as appropriate


eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 

% load window coordinates
eval(['load ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)])
if strcmp(genome_build,'mm10')
    Wins = Wins(1:20); % no unique sequences on the Y chromosome
end


for Chr = 1:size(Wins,1)
    disp(['Chromosome ' num2str(Chr)])
    eval(['load ' genome_build '_' num2str(read_length) 'bp_GC_cont_401_filtered_bins_chr' num2str(Chr)])
    win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
    for GC_bin = 1:401
        Reads_expected_nominal{Chr,1}(:,GC_bin)  =  histcounts(GC_cont_win_chr_bins{GC_bin},win_coord); % in each TIGER window, how many bp are found in each GC bin (in the genome)
    end
    Reads_expected_nominal{Chr}(:,402) = nansum(Reads_expected_nominal{Chr}'); % this is for testing purposes only. Gets overwritten later
    
    clear GC_cont_win_chr GC_cont_win_chr_bins
end

eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_' num2str(win_size) 'bp_wins_Reads_expected_nominal_401 Reads_expected_nominal'])   


eval(['cd ' TIGER_folder])

