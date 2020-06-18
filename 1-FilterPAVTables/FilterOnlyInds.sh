head -n1 NewPAV_Table_nice_names_intersection_of_individuals.filtered.csv > ~/NewPAV_Table_nice_names_intersection_of_individuals.filtered.onlyInds.csv
grep Brassica_napus NewPAV_Table_nice_names_intersection_of_individuals.filtered.csv >> ~/NewPAV_Table_nice_names_intersection_of_individuals.filtered.onlyInds.csv
mv ~/NewPAV_Table_nice_names_intersection_of_individuals.filtered.onlyInds.csv .

head -n1 Oleracea_PAV_nice_names_intersection_of_individuals.filtered.csv > ~/Oleracea_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv
grep Brassica_oleracea Oleracea_PAV_nice_names_intersection_of_individuals.filtered.csv >> ~/Oleracea_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv
mv ~/Oleracea_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv .

head -n 1 Rapa_PAV_nice_names_intersection_of_individuals.filtered.csv > ~/Rapa_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv
grep Brassica_rapa Rapa_PAV_nice_names_intersection_of_individuals.filtered.csv >> ~/Rapa_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv
mv ~/Rapa_PAV_nice_names_intersection_of_individuals.filtered.onlyInds.csv .


grep -v ../../synthetic_lines.txt NewPAV_Table_nice_names_intersection_of_individuals.filtered.onlyInds.csv > NewPAV_Table_nice_names_intersection_of_individuals.filtered.onlyInds.NoSynths.csv
