
            set terminal png size 600,500 truecolor
            set output "bwa.samtools.stats.plot/gc-depth.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Mapped depth"
            set xlabel "Percentile of mapped sequence ordered by GC content"
            set x2label "GC Content [%]"
            set title "bwa.samtools.stats" noenhanced
            set x2tics ("30" 70.588,"40" 82.353,"50" 100.000)
            set xtics nomirror
            set xrange [0.1:99.9]

            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#dedede" t '10-90th percentile' , \
                 '-' using 1:2:3 with filledcurve lt 1 lc rgb "#bbdeff" t '25-75th percentile' , \
                 '-' using 1:2 with lines lc rgb "#0084ff" t 'Median'
        11.765	0.000	0.000
17.647	0.007	0.007
23.529	0.007	0.007
29.412	0.007	0.007
35.294	0.007	0.007
41.176	0.007	0.007
47.059	0.007	0.007
58.824	0.007	0.007
70.588	0.007	0.007
76.471	0.007	0.007
82.353	0.007	0.007
88.235	0.007	0.007
100.000	0.007	0.007
end
11.765	0.000	0.000
17.647	0.007	0.007
23.529	0.007	0.007
29.412	0.007	0.007
35.294	0.007	0.007
41.176	0.007	0.007
47.059	0.007	0.007
58.824	0.007	0.007
70.588	0.007	0.007
76.471	0.007	0.007
82.353	0.007	0.007
88.235	0.007	0.007
100.000	0.007	0.007
end
11.765	0.000
17.647	0.007
23.529	0.007
29.412	0.007
35.294	0.007
41.176	0.007
47.059	0.007
58.824	0.007
70.588	0.007
76.471	0.007
82.353	0.007
88.235	0.007
100.000	0.007
end
