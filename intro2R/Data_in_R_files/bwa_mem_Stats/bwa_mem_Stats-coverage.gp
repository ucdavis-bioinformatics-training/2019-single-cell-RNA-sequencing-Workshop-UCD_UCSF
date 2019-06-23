
            set terminal png size 600,400 truecolor
            set output "bwa_mem_Stats/bwa_mem_Stats-coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "bwa_mem_Stats.log"
            set xrange [:1]
            plot '-' with lines notitle
        1	300
end
