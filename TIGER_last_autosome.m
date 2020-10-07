function num = TIGER_last_autosome(genome_build)

switch genome_build
    case 'mm10'
        num = 19;
    case 'dm6'
        num = 5;
    otherwise
        num = 22;
end
