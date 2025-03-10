for i in {1..6}; do
  # Extract cit files (barcodes, features, matrix for cit)
  for file in GSM682*patient_${i}_cit*.tsv.gz GSM682*patient_${i}_cit*.mtx.gz; do
    # Extract the base filename without the .gz extension
    base_name=$(basename "$file" .gz)
    # Decompress and move it into the corresponding folder
    gunzip -c "$file" > "patient${i}_cit/$base_name"
  done

  # Extract 2hr files (barcodes, features, matrix for 2hr)
  for file in GSM682*patient_${i}_2hr*.tsv.gz GSM682*patient_${i}_2hr*.mtx.gz; do
    # Extract the base filename without the .gz extension
    base_name=$(basename "$file" .gz)
    # Decompress and move it into the corresponding folder
    gunzip -c "$file" > "patient${i}_2hr/$base_name"
  done
done

# Loop through patient directories (e.g., patient1_cit, patient1_2hr, etc.)
for i in {1..6}; do
  # Loop through each directory (cit and 2hr)
  for dir in "patient${i}_cit" "patient${i}_2hr"; do
    # Check if the directory exists
    if [ -d "$dir" ]; then
      # Change to the directory
      cd "$dir"

      # Rename files inside the directory (be more specific with patterns)
      mv GSM6820*patient_${i}_cit_barcodes.tsv barcodes.tsv
      mv GSM6820*patient_${i}_cit_matrix.mtx matrix.mtx
      mv GSM6820*patient_${i}_cit_features.tsv features.tsv

      mv GSM6820*patient_${i}_2hr_barcodes.tsv barcodes.tsv
      mv GSM6820*patient_${i}_2hr_matrix.mtx matrix.mtx
      mv GSM6820*patient_${i}_2hr_features.tsv features.tsv

      # Go back to the parent directory
      cd ..
    else
      echo "Directory $dir not found!"
    fi
  done
done
