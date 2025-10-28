parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k21 \
    ./benchmark_script {} 21 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz

parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k11 \
    ./benchmark_script {} 11 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz

parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/logs/k6 \
    ./benchmark_script {} 6 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz


parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/construction \
    ./construct_indices {} 6 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz


parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/construction/k6 \
    ./construct_indices {} 6 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz


parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/construction/k11 \
    ./construct_indices {} 11 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz

parallel --bar --results /storage/mi/belap01/Bachelor-Arbeit/construction/k21 \
    ./construct_indices {} 21 ::: /storage/mi/belap01/Bachelor-Arbeit/input_files/*.freqskmers.fa.ust.fa.gz


./benchmark_queries kestrel.k31.freqskmers.fa.ust.fa 6 2 100000 \
    > logs/kestrel_k6.log 2>&1 &

./benchmark_queries cod.k31.freqskmers.fa.ust.fa 6 2 100000 \
    > logs/cod_k6.log 2>&1 &

./benchmark_queries human.k31.freqskmers.fa.ust.fa 6 2 100000 \
    > logs/human_k6.log 2>&1 &

wait

./benchmark_queries kestrel.k31.freqskmers.fa.ust.fa 11 2 100000 \
    > logs/kestrel_k11.log 2>&1 &

./benchmark_queries cod.k31.freqskmers.fa.ust.fa 11 2 100000 \
    > logs/cod_k11.log 2>&1 &

./benchmark_queries human.k31.freqskmers.fa.ust.fa 11 2 100000 \
    > logs/human_k11.log 2>&1 &

wait

./benchmark_queries kestrel.k31.freqskmers.fa.ust.fa 21 2 100000 \
    > logs/kestrel_k21.log 2>&1 &

./benchmark_queries cod.k31.freqskmers.fa.ust.fa 21 2 100000 \
    > logs/cod_k21.log 2>&1 &

./benchmark_queries human.k31.freqskmers.fa.ust.fa 21 2 100000 \
    > logs/human_k21.log 2>&1 &

wait