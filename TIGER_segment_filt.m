function [filt_in, data_in_segs_corrected, segment_params]  =  TIGER_segment_filt(varargin) % Data, segment_thresh, std_thres, genome build (for gap_file)
% uses the Matlab function segment to filter out CNVs and outliers
% Requires Data: four columns of chr,start,end,center locations, and any number of additional data columns
% optional: segment_thresh (R2). Default = 0.04. Increase for noisy data, but be cautious! segment function tends to mis-call segment borders with higher R2 values
% optional: std_thres. Default = 1.5. threshold beyond which to call a segment as an outlier and remove
% optional: genome build
% outputs: indexes for filtering outlying segments; segment locations and values; parameters used
% parfor requires installing the Parallel Computing Toolbox. Go to Home-->Add-Ons. can replace parfor with for, but will be slower


if nargin<4
    genome_build = 'hg19';
else
    genome_build = varargin{4};
end


if nargin==1 || isempty(varargin{2})
    R2 = 0.04;
else
    R2 = varargin{2};
end

if nargin<3 || isempty(varargin{3})
    seg_std_thresh = 1.5;
else
    seg_std_thresh = varargin{3};
end

segment_params = [R2 seg_std_thresh];



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
Gap=Gap; % solves an initiation problem with parfor



Data = varargin{1};

depth = cell2mat(Data(1:TIGER_last_autosome(genome_build)));
depth(depth>4 | depth==0) = NaN;
depth_mean = nanmedian(depth);
depth_std = nanstd(depth); % autosomal standard deviation (after excluding extremes)

data_in_segs_corrected = cell(size(Data, 1),1);
parfor Chr = 1:size(Data,1)  % ** parfor requires installing the Parallel Computing Toolbox. Go to Home-->Add-Ons
    
    % remove windows more than 10% the minimum window length (these are windows with low mappability)
    window_length = diff(Data{Chr}(:,2:3)')';
    if strcmp(genome_build,'mm10')
        long_window_in=find( window_length > min(window_length) * 2 );
    else
        long_window_in=find( window_length > min(window_length) * 1.1 );
    end
    
    
    G = Gap(Gap(:,1)==Chr,2:3)';
    G = [ Data{Chr}(1,2)-1; G(:);  Data{Chr}(end,3)];
    G = [G(1:2:end) G(2:2:end)];
    
    
    for c = 5:size(Data{Chr},2)
        disp(['Segmenting, Sample ' num2str(c-4) ' out of ' num2str(size(Data{Chr},2)-4)  ', Chromosome ' num2str(Chr)])
        data_in_segs_corrected{Chr}(:,c-4) = NaN(size(Data{Chr},1),1);
        
        for frag = 1:size(G,1)  % segment between gaps
            in = find(Data{Chr}(:,2)>G(frag,1) & Data{Chr}(:,3)<=G(frag,2) & Data{Chr}(:,c)<=20); % <=20: remove extreme data points that will fail the segment algorithm
            if length(in)>20
                
                % segment
                data_in_segs_corrected{Chr}(in,c-4) = seg_fun(Data{Chr}(in,c), R2);
                
            end
        end
        
        
        % apply filters
        % first filters segments based on autosomal median, then further filters based on median of the actual chromosome
        if Chr>last_autosome(genome_build) % don't filter sex chromosomes based on CN compared to autosomes
            in2 = [];
        else
            in2 = find(  abs(data_in_segs_corrected{Chr}(:,c-4) - depth_mean(c) ) > seg_std_thresh*depth_std(c) ); % number of stds. Removes data points in segments that are more than this number of stds from the *genome* median
        end
        in3 = find(  abs(data_in_segs_corrected{Chr}(:,c-4) - nanmedian(data_in_segs_corrected{Chr}(setdiff(1:length(data_in_segs_corrected{Chr}),in2),c-4)) ) > seg_std_thresh*depth_std(c) ); % number of stds. Removes data points in segments that are more than this number of stds from the *chromosome* median
        in4 = find(data_in_segs_corrected{Chr}(:,c-4)==0);
        in5 = find(isnan(data_in_segs_corrected{Chr}(:,c-4)));
        all_filt = [long_window_in; in2 ; in3; in4 ; in5 ];
        
        
        filt_in{Chr}{c-4} = all_filt;   % save filtered indexes only (to save space since this parameter is going to be saved)
        
    end
    
end




function data_out = seg_fun(data_in,R2)

% segment data
thm = segment([data_in,ones(length(data_in),1)],[0 1 1],R2); % function segment requires the System Identification Toolbox
thm = thm([1 1:end-1]); % correction: segment always results in a shift of 1 in the coordinates


% failed segmentation resulting in all NaNs: iteratively remove 5 data points at a time until resolved
nan = length(find(isnan(thm)));
iter = 1;
while nan>length(thm)*0.9  % if most windows per region between gaps = NaN, remove additional 5 data points from beginning of region for each iteration
    idx = 5*iter;
    thm = segment([data_in(idx+1:end),ones(length(data_in)-idx,1)],[0 1 1],R2);
    thm = thm([1 1:end-1]); % correction: segment always results in a shift of 1 in the coordinates
    thm = [NaN(idx,1); thm]; % add NaNs in the beginning
    nan = length(find(isnan(thm))); % recount number of NaNs
    iter = iter+1;
end

% re-call segment values as median of segments (because the segment function sometimes miscalls values)
% find discrete segments
d = diff(thm);
in_seg = find(d~=0);
in_seg = [0;in_seg;length(thm)];
for i = 1:length(in_seg)-1
    thm(in_seg(i)+1:in_seg(i+1))  =  nanmedian( data_in(in_seg(i)+1:in_seg(i+1)) );
end

data_out = thm;





