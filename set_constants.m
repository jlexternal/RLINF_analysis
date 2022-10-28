% set_constants
clear all
samplename = 'sample2';
ncnd = 2;
nblk = 4;
ntrl = 73;
nepis = 24;
fnr = .30;

condrgb = [93 74 25; 25 42 68]/100;

save(sprintf('./constants/constants_rlinf_%s',samplename));