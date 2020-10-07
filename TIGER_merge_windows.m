function [Loc_combined Dat_combined]  =  TIGER_merge_windows (varargin)
% merges read number windows to larger windows. 
% requires at least the first three of following input variables:
% 1. Data (normalized data, used for the locations (four columns))
% 2. Count_Expected_Data- this is what will actually be merged
% 3. Windows_to_combine- number of windows to merge together
% 4. Optional: genome_build- used for removal of windows spanning gaps. hg19 by default
% Returns window coordinates (Loc_combined) and Count_expected_Data after window merging (Dat_combined)
% Code can be modified to provide sliding windows as well


Data = varargin{1};
Count_Expected_Data = varargin{2};
Windows_to_combine = varargin{3};
if nargin<4
    genome_build = 'hg19';
else
    genome_build = varargin{4};
end


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


% combine- no overlapping windows
for Chr = 1:size(Data,1)
    
    Win_locs = Data{Chr}(:,2:3);
    if mod(length(Win_locs),Windows_to_combine)~= 0
        Win_locs(end-mod(length(Win_locs),Windows_to_combine)+1:end,:) = [];
    end
    Loc_combined{Chr,1}(:,2) = Win_locs(1:Windows_to_combine:end,1);
    Loc_combined{Chr,1}(:,3) = Win_locs(Windows_to_combine:Windows_to_combine:end,2);
    Loc_combined{Chr,1}(:,1) = Chr;
    Loc_combined{Chr,1}(:,4) = round(mean(Loc_combined{Chr,1}(:,2:3)')');
    
    for c = 1:size(Count_Expected_Data,2)
        CN_data = Count_Expected_Data{Chr,c};
        if mod(length(CN_data),Windows_to_combine)~= 0
            CN_data(end-mod(length(CN_data),Windows_to_combine)+1:end,:) = [];
        end
        Dat_combined{Chr,c}(:,1) = CN_data(1:Windows_to_combine:end,1);
        for k = 2:Windows_to_combine
            Dat_combined{Chr,c}(:,1) = Dat_combined{Chr,c}(:,1)+CN_data(k:Windows_to_combine:end,1);
        end
        Dat_combined{Chr,c}(:,2) = CN_data(1:Windows_to_combine:end,2);
        for k = 2:Windows_to_combine
            Dat_combined{Chr,c}(:,2) = Dat_combined{Chr,c}(:,2)+CN_data(k:Windows_to_combine:end,2);
        end
    end
 
    
    % remove windows that span gaps
    G = Gap(Gap(:,1)==Chr,2:3);
    in = [];
    for i = 1:size(G,1)
        in = [in;find(Loc_combined{Chr,1}(:,2)<= G(i,1) & Loc_combined{Chr,1}(:,3)>= G(i,2))];
    end
    Loc_combined{Chr}(in,:) = [];
    for c = 1:size(Count_Expected_Data,2)
        Dat_combined{Chr,c}(in,:) = [];
    end

    
end






% plot original and merged windows
if 0
    figure; hold on
    Chr = 2;
    samp = 1;
    plot(Data{Chr}(:,4)./1e6,Data{Chr}(:,samp+4),'.')
    plot(Loc_combined{Chr}(:,4)./1e6,Dat_combined{Chr,samp}(:,1)./Dat_combined{Chr,samp}(:,2),'.')
    axis tight
    legend('Original','Combined')
end

