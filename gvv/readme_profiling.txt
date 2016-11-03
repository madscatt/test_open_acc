
nvprof -o timeline.pgprof

nvprof --analysis-metrics -o metrics.pgprof ./a.out

nvvp timeline.pgprof metrics.pgprof

