
        set terminal png size 600,400 truecolor
        set output "bwa_mem_Stats/bwa_mem_Stats-indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "bwa_mem_Stats.log"
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	39322
2	8214
3	3788
4	2705
5	701
6	589
7	213
8	316
9	157
10	149
11	56
12	83
13	34
14	36
15	23
16	37
17	6
18	16
19	14
20	5
21	10
22	4
23	3
24	1
25	2
26	4
27	0
28	1
29	0
30	1
31	0
32	0
33	0
34	0
35	0
36	0
40	0
42	0
43	0
end
1	32727
2	10343
3	7157
4	3487
5	837
6	880
7	335
8	451
9	278
10	213
11	103
12	238
13	71
14	85
15	82
16	49
17	27
18	46
19	19
20	35
21	20
22	18
23	7
24	27
25	11
26	9
27	5
28	7
29	1
30	6
31	5
32	8
33	1
34	1
35	1
36	1
40	1
42	1
43	2
end
1	1.201516
2	0.794160
3	0.529272
4	0.775738
5	0.837515
6	0.669318
7	0.635821
8	0.700665
9	0.564748
10	0.699531
11	0.543689
12	0.348739
13	0.478873
14	0.423529
15	0.280488
16	0.755102
17	0.222222
18	0.347826
19	0.736842
20	0.142857
21	0.500000
22	0.222222
23	0.428571
24	0.037037
25	0.181818
26	0.444444
27	0.000000
28	0.142857
29	0.000000
30	0.166667
31	0.000000
32	0.000000
33	0.000000
34	0.000000
35	0.000000
36	0.000000
40	0.000000
42	0.000000
43	0.000000
end
