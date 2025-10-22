parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k21 \
    ./benchmark_script {} 21 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz

parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k11 \
    ./benchmark_script {} 11 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz

parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k6 \
    ./benchmark_script {} 6 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz
