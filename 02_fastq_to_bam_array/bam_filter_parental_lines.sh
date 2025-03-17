
## Pull the Mean Insert Size +/- 2 standard deviations
awk 'NR==1 {for (i=1; i<=NF; i++) if ($i == "MEAN_INSERT_SIZE") mean_idx = i; if ($i == "STANDARD_DEVIATION") sd_idx = i} NR>1 && NF > 0 {print $mean_idx + 2 * $sd_idx}' $name.output_metrics.txt > upper_limit
awk 'NR==1 {for (i=1; i<=NF; i++) if ($i == "MEAN_INSERT_SIZE") mean_idx = i; if ($i == "STANDARD_DEVIATION") sd_idx = i} NR>1 && NF > 0 {result = $mean_idx + 2 * $sd_idx; print (result < 0 ? 0 : result)}' $name.output_metrics.txt > lower_limit
