import os
import scipy.stats
import json
import time


start_time = time.time()

####Input files for script to run are attached in the dir input_file_for script:

chromHMM_file = "/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/E008_25_imputed12marks_dense.bed"  
background_dmr_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/merged_processed_DMR/analyse_DMR_output_data/merged_processed_DMR/hypo_H9_hESC_background_dmr.bed"

##
dmr_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/merged_processed_DMR/analyse_DMR_output_data/merged_processed_DMR/human_ES_cell_final_dmrs.bed"
dmr_file1="/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_dmrs.txt"
dmr_file2="/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_dmrs.txt"
##

####Output files from the script run (provide the path location to save the file):
#Provide the file path location to save the output of ENCODE_TF_cluster file as JSON format (could be useful for future usage):
JSON_dict_file = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/ENCODE_chromHMM_JSON_1.txt"

##
TF_fishers_enrichment = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_chromHMM_fishers_enrichment_result.txt"
##

out_file = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_significant_dmr_overlaps_with_chromHMM.bed"
out_file_2 = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/hypo_background_dmr_overlaps_with_chromHMM.bed"

master_dict_return = {}
TF_list = []
header = ["chrom", "start", "end", "chromHMM_state", "cell_line"]
head_print = "\t".join(header)

def parse_chromHMM(chromHMM_file, Json_output_file):
	with open(chromHMM_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file.readlines()[1:]:
				#print line
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, ch_state = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(ch_state)]

				if ch_state in master_dict_return.keys():
					if chrom in master_dict_return[ch_state].keys():
						master_dict_return[ch_state][chrom].append((chrom, int(start), int(end), "\t".join(line_info)))
					else:
						master_dict_return[ch_state][chrom] = [(chrom, int(start), int(end), "\t".join(line_info))]
				else:
					master_dict_return[ch_state] = {chrom:[(chrom, int(start), int(end), "\t".join(line_info))]}
					master_dict_return[ch_state].update({"significant_dmr_hits":0})
					master_dict_return[ch_state].update({"background_dmr_hits":0})
					master_dict_return[ch_state].update({"custom_overlap_list":[]})
					master_dict_return[ch_state].update({"As_overlap_list":[]})
					master_dict_return[ch_state].update({"Bs_overlap_list":[]})

			json.dump(master_dict_return, outfile)
				#print master_dict_return
			return(master_dict_return)

master_chromHMM_dict = parse_chromHMM(chromHMM_file, JSON_dict_file)
print "\n\nTime to parse the ENCODE_TF_cluster = ", time.time()-start_time


#IF needed, read JSON file parsed earlier from the bedfile:
#with open(JSON_dict_file) as json_file:
#	master_TF_dict_1 = json.load(json_file)


def genome_overlap_1(a, b):
	if (a[2] < b[1]) or (a[1] > b[2]):
		return False 
	else:
		return True 




#################
#################




start_time = time.time()
print "\n\nStep 1.... significant dmr analysis..."


##Following piece of code prevents the append file mode to append the outputs to be appended
#to the same file:
directory = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/sig_DMR_chromHMM"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
file_name_list = [key+".txt" for key in master_chromHMM_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file) 


custom_overlap_list = [] 
with open(out_file, "w") as out_file:
	with open(dmr_file1, "r") as file:
		for line in file:
			splitted = line.strip().split()
			chrom = splitted[0]
			start = splitted[1]
			end = splitted[2]
			qval = splitted[3]
			meth_perc = splitted[4]

			line_coords = [chrom, start, end, qval, meth_perc]
			A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))


			for key,value in master_chromHMM_dict.iteritems():
				each_tf_dict = master_chromHMM_dict[key]

				#name the output file with its keyword:
				file_name = key.replace("/","-") + ".txt"
				outfile_path = directory + "/" + file_name

				B_tuple_list = each_tf_dict.get(chrom)       
				if B_tuple_list:
					for B_tuple in B_tuple_list:               
						if genome_overlap_1(A_tuple, B_tuple):
								master_chromHMM_dict[key]["significant_dmr_hits"] += 1
								overlap_tuples = "%s\t%s\t%s" %(key, A_tuple[3], B_tuple[3])
								custom_overlap_list.append(overlap_tuples)

								master_chromHMM_dict[key]["custom_overlap_list"].append(overlap_tuples)
								master_chromHMM_dict[key]["As_overlap_list"].append(A_tuple[3])
								master_chromHMM_dict[key]["Bs_overlap_list"].append(B_tuple[3])
								#print overlap_tuples

								with open(outfile_path, "a") as appendfile:
									appendfile.write(overlap_tuples + "\n")

		print "Total sig count for last ch_state in dict constructed ::: ", each_tf_dict["significant_dmr_hits"]
		print "But, total background count for last ch_state in dict ::: ", each_tf_dict["background_dmr_hits"]

print "\nsignificant_dmr_overlap completed!!!!!"
print "Time for significant dmr analysis = ", time.time()-start_time

out_file = directory + "/" + "sig_dmrhit_results.txt"
with open(out_file, "w") as outfile:
	for key,value in master_chromHMM_dict.items():
		result = "%s\t%s" %(key, master_chromHMM_dict[key]["significant_dmr_hits"])
		outfile.write(result + "\n")
		print result


#####################
#####################


print "\n\nStep 2.... background dmr analysis....."
start_time = time.time()

##Following piece of code prevents the append file mode to append the outputs to be appended
#to the same file:
directory = "/home/surya/Desktop/ENCODE_poster/parsed_data_R_output/backg_DMR_chromHMM"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
file_name_list = [key+".txt" for key in master_chromHMM_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file) 


custom_overlap_list_2 = []
with open(out_file_2, "w") as outfile:
	with open(background_dmr_file, "r") as file:
		for line in file:
			splitted = line.strip().split()
			chrom = splitted[0]
			start = splitted[1]
			end = splitted[2]

			line_coords = [chrom, start, end]
			A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))


			for key,value in master_chromHMM_dict.iteritems():
				each_tf_dict = master_chromHMM_dict[key]
				
				#name the output file with its keyword:
				file_name = key.replace("/","-") + ".txt"
				outfile_path = directory + "/" + file_name

				B_tuple_list = each_tf_dict.get(chrom)       
				if B_tuple_list:
					for B_tuple in B_tuple_list:               
						if genome_overlap_1(A_tuple, B_tuple):
								master_chromHMM_dict[key]["background_dmr_hits"] += 1
								overlap_tuples = "%s\t%s\t%s" %(key, A_tuple[3], B_tuple[3])
								custom_overlap_list_2.append(overlap_tuples)

								master_chromHMM_dict[key]["custom_overlap_list"].append(overlap_tuples)
								master_chromHMM_dict[key]["As_overlap_list"].append(A_tuple[3])
								master_chromHMM_dict[key]["Bs_overlap_list"].append(B_tuple[3])
								#print overlap_tuples

								with open(outfile_path, "a") as appendfile:
									appendfile.write(overlap_tuples + "\n")

		print "Total background count for last ch_state in dict constructed ::: ", each_tf_dict["background_dmr_hits"]
		print "Also, total significant count for last ch_state in dictionary::: ", each_tf_dict["significant_dmr_hits"]


print "\nbackground_dmr_overlap completed!!!!!"
print "Time for background dmr analysis = ", time.time()-start_time
print("\n\n")



#Append significant dmr and background dmr hits to list_array and
#calculate fisher exact test & bonferroni correction for multiple test correction:

with open(TF_fishers_enrichment, "w") as out_file:	
	FishTable = []
	for i in range(5):
		FishTable.append([])

	for dict_key,dict_value in master_chromHMM_dict.iteritems():
		FishTable[0].append(dict_key)
		FishTable[1].append(master_chromHMM_dict[dict_key]["significant_dmr_hits"])
		FishTable[2].append(master_chromHMM_dict[dict_key]["background_dmr_hits"])

	# FishTable[3] = [items for i, items in zip(FishTable[0], FishTable[1])]
	# FishTable[4] = [items for i, items in zip(FishTable[0], FishTable[1])]

	#Calculating the enrichment pvalue using scipy fishers exact test:
	Significant_hits = len(custom_overlap_list)
	Background_hits = len(custom_overlap_list_2)
	header = ("ChromHMM_state", "Significant_DMR_hits", "Background_dmr_Hits", "Enrichment_pval", "Enrichment_Bonferonni_Adj")
	out_file.write("\t".join(header) + "\n")
	print "\t".join(header)
	for i in range(len(FishTable[0])):
		fishers_pvalue = scipy.stats.fisher_exact([[FishTable[1][i], Significant_hits - FishTable[1][i]], [FishTable[2][i], Background_hits - FishTable[2][i]]], 'greater')[1]
		FishTable[3].append(fishers_pvalue)
		bonferroni_correction = min(scipy.stats.fisher_exact([[FishTable[1][i], Significant_hits - FishTable[1][i]], [FishTable[2][i], Background_hits - FishTable[2][i]]], 'greater')[1] * len(FishTable[0]), 1)
		FishTable[4].append(bonferroni_correction)
		fishers_test_result = (str(FishTable[0][i]), str(FishTable[1][i]), str(FishTable[2][i]), str(FishTable[3][i]), str(FishTable[4][i]))
		fisher_print = "\t".join(fishers_test_result)
		print fisher_print
		out_file.write(fisher_print + "\n")

print "\nEnd of job!!\n"
