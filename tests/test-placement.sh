rm data/D652/placements_query_10k.fasta.jplace 2> /dev/null

epik.py place -i data/D652/DB_k7_o2.0.ipk data/D652/query_10k.fasta -o data/D652 --threads 2

python3 ../scripts/jplace_diff.py data/D652/reference.jplace data/D652/placements_query_10k.fasta.jplace
expected_output="10000/10000 placements match."

if [ "$output" == "$expected_output" ]; then
  echo "OK"
else
  echo "Unpredicted output: $output"
  exit 1
fi