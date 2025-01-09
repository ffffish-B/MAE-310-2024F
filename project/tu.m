%% Matlab mesh
%% quarter-plate-with-hole-quad, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 120;
msh.POS = [
1 -1 0;
1 1 0;
-1 1 0;
-0.7 -1 0;
-1 -0.7 0;
-0.7878679656440357 -0.7878679656440357 0;
-0.7018863370312716 -0.9664106571757137 0;
-0.7075216263762176 -0.9332437196783165 0;
-0.7168350010054405 -0.9009162811335227 0;
-0.7297093397656539 -0.8698348779815372 0;
-0.7459827403349979 -0.8403903768807064 0;
-0.7654505552586767 -0.8129530594435265 0;
-0.8129530594860284 -0.7654505552247826 0;
-0.8403903770942989 -0.7459827402007889 0;
-0.8698348785271101 -0.7297093395029198 0;
-0.9009162817706299 -0.7168350007825071 0;
-0.9332437202484005 -0.7075216262460997 0;
-0.9664106573628806 -0.7018863370101829 0;
-1 -0.4571428571432274 0;
-1 -0.2142857142867694 0;
-1 0.02857142856968053 0;
-1 0.2714285714267739 0;
-1 0.5142857142846606 0;
-1 0.7571428571423737 0;
-0.7142857142865064 1 0;
-0.4285714285730131 1 0;
-0.1428571428595189 1 0;
0.1428571428547671 1 0;
0.4285714285698448 1 0;
0.7142857142849224 1 0;
1 0.7142857142865064 0;
1 0.4285714285730131 0;
1 0.1428571428595189 0;
1 -0.1428571428547671 0;
1 -0.4285714285698448 0;
1 -0.7142857142849224 0;
0.7571428571432384 -1 0;
0.5142857142867494 -1 0;
0.2714285714303195 -1 0;
0.02857142857322614 -1 0;
-0.2142857142846606 -1 0;
-0.4571428571423737 -1 0;
0.744590290622789 0.744590290622789 0;
0.4891805812457544 0.4891805812457544 0;
0.2337708718687292 0.2337708718687292 0;
-0.02163883750893048 -0.02163883750893048 0;
-0.2770485468872206 -0.2770485468872206 0;
-0.5324582562657829 -0.5324582562657829 0;
-0.7503035632973445 0.7568733804267743 0;
-0.7863214123081576 0.5137467608534991 0;
-0.8223392613189693 0.270620141280077 0;
-0.8583571103298714 0.02749352170724417 0;
-0.8943749593408633 -0.2156330978651264 0;
-0.930392808351894 -0.4587597174375417 0;
-0.5006674702407806 0.7560683391074988 0;
-0.5727635119084981 0.5121366782149857 0;
-0.6448595535762125 0.2682050173223519 0;
-0.7169555952441082 0.0242733564301037 0;
-0.7890516369121836 -0.2196583044618639 0;
-0.8611476785803363 -0.4635899653539134 0;
-0.2511513055608888 0.7547378570310109 0;
-0.3594454682621834 0.5094757140620468 0;
-0.4677396309634736 0.2642135710929887 0;
-0.5760337936650357 0.01895142812411227 0;
-0.6843279563668671 -0.2263107148446654 0;
-0.7926221190688147 -0.4715728578135623 0;
-0.001813145913782385 0.7528986657853792 0;
-0.1464834346822316 0.5057973315708214 0;
-0.291153723450675 0.2586959973561953 0;
-0.4358240122194807 0.01159466314154728 0;
-0.5804943009886463 -0.2355066710731841 0;
-0.7251645897579663 -0.4826080052880727 0;
0.2472911706181869 0.7505738942572543 0;
0.06601091266665469 0.5011477885146089 0;
-0.1152693452848703 0.2517216827719216 0;
-0.2965496032368485 0.002295577029007928 0;
-0.4778298611892764 -0.2471305287141711 0;
-0.6591101191418979 -0.4965566344575454 0;
0.4961087466037936 0.7477927778253973 0;
0.2779317789228154 0.4955855556509327 0;
0.05975481124184598 0.2433783334764522 0;
-0.1584221564396675 -0.008828888698458909 0;
-0.376599124121721 -0.2610361108738179 0;
-0.5947760918040071 -0.5132433330494103 0;
0.747792777820677 0.496108746611222 0;
0.7505738942383251 0.2472911706514143 0;
0.7528986657482124 -0.001813145831772334 0;
0.754737856999653 -0.2511513054658023 0;
0.7560683390895246 -0.5006674701566255 0;
0.7568733804245007 -0.7503035632692494 0;
0.4955855556415423 0.2779317789360877 0;
0.501147788476852 0.06601091272994036 0;
0.5057973314966395 -0.1464834345229639 0;
0.509475713999535 -0.3594454680767625 0;
0.5121366781792925 -0.5727635117433569 0;
0.5137467608492589 -0.7863214122535517 0;
0.2433783334624244 0.0597548112609619 0;
0.2517216827154025 -0.1152693451915269 0;
0.2586959972450975 -0.2911537232141502 0;
0.2642135709994545 -0.4677396306877191 0;
0.2682050172691047 -0.6448595533300856 0;
0.270620141274069 -0.8223392612378525 0;
-0.008828888717328707 -0.1584221564147062 0;
0.002295576953316128 -0.2965496031134445 0;
0.01159466299291581 -0.435824011905696 0;
0.01895142799873035 -0.5760337932989446 0;
0.02427335635826845 -0.7169555949169931 0;
0.02749352169822376 -0.8583571102222427 0;
-0.2610361108977315 -0.3765991240909125 0;
-0.2471305288094404 -0.4778298610358091 0;
-0.2355066712599578 -0.5804943005975981 0;
-0.2263107150027089 -0.6843279559104369 0;
-0.2196583045533078 -0.7890516365040785 0;
-0.2156330978783875 -0.8943749592067218 0;
-0.5132433330783422 -0.5947760917673517 0;
-0.4965566345723413 -0.6591101189583679 0;
-0.4826080055129116 -0.7251645892896555 0;
-0.4715728580041645 -0.7926221185220459 0;
-0.4635899654648368 -0.861147678091241 0;
-0.4587597174548882 -0.9303928081912398 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 4 7 13
 7 8 13
 8 9 13
 9 10 13
 10 11 13
 11 12 13
 12 6 13
 6 13 13
 13 14 13
 14 15 13
 15 16 13
 16 17 13
 17 18 13
 18 5 13
 5 19 10
 19 20 10
 20 21 10
 21 22 10
 22 23 10
 23 24 10
 24 3 10
 3 25 9
 25 26 9
 26 27 9
 27 28 9
 28 29 9
 29 30 9
 30 2 9
 2 31 8
 31 32 8
 32 33 8
 33 34 8
 34 35 8
 35 36 8
 36 1 8
 1 37 11
 37 38 11
 38 39 11
 39 40 11
 40 41 11
 41 42 11
 42 4 11
];
msh.QUADS =[
 3 25 49 24 12
 24 49 50 23 12
 23 50 51 22 12
 22 51 52 21 12
 21 52 53 20 12
 20 53 54 19 12
 19 54 18 5 12
 25 26 55 49 12
 49 55 56 50 12
 50 56 57 51 12
 51 57 58 52 12
 52 58 59 53 12
 53 59 60 54 12
 54 60 17 18 12
 26 27 61 55 12
 55 61 62 56 12
 56 62 63 57 12
 57 63 64 58 12
 58 64 65 59 12
 59 65 66 60 12
 60 66 16 17 12
 27 28 67 61 12
 61 67 68 62 12
 62 68 69 63 12
 63 69 70 64 12
 64 70 71 65 12
 65 71 72 66 12
 66 72 15 16 12
 28 29 73 67 12
 67 73 74 68 12
 68 74 75 69 12
 69 75 76 70 12
 70 76 77 71 12
 71 77 78 72 12
 72 78 14 15 12
 29 30 79 73 12
 73 79 80 74 12
 74 80 81 75 12
 75 81 82 76 12
 76 82 83 77 12
 77 83 84 78 12
 78 84 13 14 12
 30 2 43 79 12
 79 43 44 80 12
 80 44 45 81 12
 81 45 46 82 12
 82 46 47 83 12
 83 47 48 84 12
 84 48 6 13 12
 2 43 85 31 12
 31 85 86 32 12
 32 86 87 33 12
 33 87 88 34 12
 34 88 89 35 12
 35 89 90 36 12
 36 90 37 1 12
 43 44 91 85 12
 85 91 92 86 12
 86 92 93 87 12
 87 93 94 88 12
 88 94 95 89 12
 89 95 96 90 12
 90 96 38 37 12
 44 45 97 91 12
 91 97 98 92 12
 92 98 99 93 12
 93 99 100 94 12
 94 100 101 95 12
 95 101 102 96 12
 96 102 39 38 12
 45 46 103 97 12
 97 103 104 98 12
 98 104 105 99 12
 99 105 106 100 12
 100 106 107 101 12
 101 107 108 102 12
 102 108 40 39 12
 46 47 109 103 12
 103 109 110 104 12
 104 110 111 105 12
 105 111 112 106 12
 106 112 113 107 12
 107 113 114 108 12
 108 114 41 40 12
 47 48 115 109 12
 109 115 116 110 12
 110 116 117 111 12
 111 117 118 112 12
 112 118 119 113 12
 113 119 120 114 12
 114 120 42 41 12
 48 6 12 115 12
 115 12 11 116 12
 116 11 10 117 12
 117 10 9 118 12
 118 9 8 119 12
 119 8 7 120 12
 120 7 4 42 12
];
