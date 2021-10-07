% This contains optimal/best values from each dataset
% The data is called in qrun_tests.m. 
%text1 = '\n Loading this file with optimal/best values from each dataset; \n';

%% optimal values, problem names for small size problems
small_opt = [9552;9742;11156;9896;7990;9504;11098;1534;2192;2298;14142;
    17212548;68;292;160;16;28;0;26;996;14;8;1652;2724;3720;5358;6922;
    578;1014;1150;1610;1240;1732;1930;2570;235528;354210;725522;31410;
    51140;110030;135028;224416;388214;491812;703482];

name_small = ["chr12a";"chr12b";"chr12c";"chr15a";"chr15b";"chr15c";"chr18a"
    "chr18b";"chr20a";"chr20b";"chr20c";"els19";"esc16a";"esc16b";"esc16c"
    "esc16d";"esc16e";"esc16f";"esc16g";"esc16h";"esc16i";"esc16j";"had12"
    "had14";"had16";"had18";"had20";"nug12";"nug14";"nug15";"nug16a"
    "nug16b";"nug17";"nug18";"nug20";"rou12";"rou15";"rou20";"scr12"
    "scr15";"scr20";"tai10a";"tai12a";"tai15a";"tai17a";"tai20a"];


%% optimal values, problem names for medium size problems
medium_opt = [ 6156;6194;3796;130;168;642;200;2;6;438;88900;91420;88700;
    2438;3596;3488;3744;5234;5166;6124;9526;15852;8239110;1167256;
    1818146;2422002;3139370;149936;240516];

name_medium = ["chr22a";"chr22b";"chr25a";"esc32a";"esc32b";"esc32c";
    "esc32d";"esc32e";"esc32g";"esc32h";"kra30a";"kra30b";"kra32";"nug21";
    "nug22";"nug24";"nug25";"nug27";"nug28";"nug30";"ste36a";"ste36b";
    "ste36c";"tai25a";"tai30a";"tai35a";"tai40a";"tho30";"tho40*"];

%% optimal values, problem names for large size problems
large_opt = [116;15812;23386;34458;48498;
    4938796; 7205962;1855928;48816];

name_large = ["esc64a";"sko42";"sko49";"sko56";"sko64";
    "tai50a";"tai60a";"tai64c";
    "wil50"];


%% consolidate optimal values, problem names for all problems
all_opt = [small_opt;medium_opt;large_opt];
name_all = [name_small;name_medium;name_large];
