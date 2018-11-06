#!/bin/python

import sys
import random

#get the dictionary sample2muts from the mutation matrix

def get_sample2muts(mut_matrix_file):
	sample2muts = dict()
	file_i = open(mut_matrix_file,'r')
	for line in file_i:
        	v = line.split()
	        sampleID = v[0].strip()
        	if sampleID in sample2muts:
                	print "Duplicated sample: "+str(sampleID)
	        sample2muts[sampleID]=set()
        	for i in range(len(v)-1):
                	sample2muts[sampleID].add(v[i].strip())
	file_i.close()
	return sample2muts

#get the classes from the class file
def get_sample2class(class_file):
	file_i = open(class_file,'r')
	line_c = 0
	sample2class = dict()
	for line in file_i:
        	if line_c > 0:
                	v = line.split()
	                sampleID = v[0].strip()
        	        class_sample = v[1].strip()
                	if sampleID in sample2class:
                        	print "Duplicated sample: "+str(sampleID)
	                sample2class[sampleID]=class_sample
        	line_c += 1
	file_i.close()
	return sample2class

#get the list of all genes and build the dictionary to get in which samples they are mutated
def get_muts2sample(sample2muts, sample2class):
	all_muts = list()
	muts2sample = dict()
	total_samples_class0 = 0.
	total_samples_class1 = 0.
	for sample in samples_with_muts_AND_class:
        	muts = sample2muts[sample]
	        for mut in muts:
        	        if mut not in all_muts:
                	        all_muts.append(mut)
	                if mut not in muts2sample:
        	                muts2sample[mut]=list()
                	if sample not in muts2sample[mut]:
                        	muts2sample[mut].append(sample)
		if sample2class[sample] == "0":
			total_samples_class0 += 1.
		if sample2class[sample] == "1":
                        total_samples_class1 += 1.
	return muts2sample, all_muts, total_samples_class0, total_samples_class1

#get the adjacency matrix of network from the edge lists
def get_adj_matrix_network(network_file):
	mut2neighbors = dict()
	file_i = open(network_file,'r')

	for line in file_i:
        	v = line.split()
       		mut1 = v[0].strip().upper()
        	mut2 = v[1].strip().upper()
		if mut1 not in mut2neighbors:
			mut2neighbors[mut1] = list()
		if mut2 not in mut2neighbors[mut1]:
			mut2neighbors[mut1].append(mut2)
		if mut2 not in mut2neighbors:
                        mut2neighbors[mut2] = list()
                if mut1 not in mut2neighbors[mut2]:
                        mut2neighbors[mut2].append(mut1)
	return mut2neighbors

#given two lists of muts, computes the absolute value 
#of the discrepancy of the muts in the union of the two

def compute_weight(mut2samples, sample2class, list1, list2, tot_class0, tot_class1):
	set1 = set(list1)
	set2 = set(list2)
	all_muts_set = set1.union(set2)
	all_muts = list( all_muts_set )
	all_muts.sort()
	count_class0 = 0.
        count_class1 = 0.
	all_samples = set()
	for mut in all_muts:
		if mut in mut2samples:
			all_samples = all_samples.union(mut2samples[mut])
        for sample in all_samples:
        	if sample2class[sample] == "0":
                	count_class0 += 1.
                if sample2class[sample] == "1":
                        count_class1 += 1.
        f0 = count_class0/tot_class0
        f1 = count_class1/tot_class1
	weight = abs(f0 - f1)
	#if ("TP53" in list1 and "TEC" in list2) or ("TEC" in list1 and "TP53" in list2):
        #        print "computing weight for TP53 and TEC, weight is: "+str(weight)
	return weight
	
def get_mut_white_list(mut2samples, muts_black_list_file):
        all_muts = mut2samples.keys()
        all_muts.sort()
        #get all muts that pass the minimum threshold
	mut_white_list = set()
        for mut in all_muts:
                samples_mut = mut2samples[mut]
                #minimum threhsold requirement check
                count0 = 0.
                count1 = 0.
                for sample in samples_mut:
                        if sample2class[sample] == "0":
                                count0 += 1.
                        if sample2class[sample] == "1":
                                count1 += 1.
                #NB:  the minimum threshold requirement must be satisfied for at least one of the two classes, NOT BOTH!
                if count0 > min_freq or count1 > min_freq:
                        mut_white_list.add( mut )
	infile = open(muts_black_list_file, 'r')
	mut_bl = set()
	for line in infile:
		mut = line.strip()
		mut_bl.add(mut)
	infile.close()
	return mut_white_list.difference(mut_bl)

#get all edges with weight above a given threshold
def get_edges_above_threshold(mut2samples, sample2class, tot_class0, tot_class1, mut2neighbors, mut_white_list):
	all_muts = mut2samples.keys()
	all_muts.sort()
	#now consider all edges
	edges= list()
	for mut in all_muts:
		if mut in mut_white_list and mut in mut2neighbors:
			neighbors = mut2neighbors[mut]
			for neigh_mut in neighbors:
				#if (mut == "TP53" and neigh_mut == "TEC") or (neigh_mut == "TP53" and mut == "TEC"):
				#	print "Now should compute weight for TEC and TP53" 
				#	if "TEC" in mut_white_list:
				#		print "TEC is in mut_white_list"
				#	else:
				#		print "TEC is not in mut_white_list!"
				if neigh_mut > mut and neigh_mut in mut_white_list:
					edge_weight = compute_weight(mut2samples, sample2class, list([mut]), list([neigh_mut]), tot_class0, tot_class1) 
				#	if ("TP53" in list([mut]) and "TEC" in list([neigh_mut])) or ("TEC" in list([mut]) and "TP53" in list([neigh_mut])):
				#		print "computed weight for TP53 and TEC, value is: "+str(edge_weight)
					if edge_weight >= threshold:
						edge = list([mut, neigh_mut])
						edge.sort()
						edges.append(edge)
				#		print "partial edges: "+str(edges)
	return edges

#starting from a current solution, return the best extension by a neighboring gene
def extend_solution(mut2samples, sample2class, tot_class0, tot_class1, mut2neighbors, mut_white_list, curr_sol):
	#try to add each neighbour of the solution to the current solution, and keep the one
	#with the highest weight; constraints to add a mutation to a solution: must increase
	#the weight
	tmp_curr_sol = list()
	for mut in curr_sol:
		tmp_curr_sol.append(mut)
	neighbors = set()
	#print tmp_curr_sol
	for mut in tmp_curr_sol:
		#print mut
		for mut_neigh in mut2neighbors[mut]:
			neighbors.add(mut_neigh)
	neighbors = neighbors.difference(set(tmp_curr_sol))
	to_remove = set()
	for mut in neighbors:
		if mut not in mut_white_list:
			to_remove.add(mut)
	neighbors = neighbors.difference(to_remove)
	#try one neighbor at the time, and keep the one with highest weight if it increases the weight
	curr_weight = compute_weight(mut2samples, sample2class, tmp_curr_sol, list(), tot_class0, tot_class1)
	max_val = curr_weight
	to_add = list()
	for neighbor in neighbors:
		new_weight = compute_weight(mut2samples, sample2class, tmp_curr_sol, list([neighbor]), tot_class0, tot_class1)
		if new_weight > max_val:
			max_val = new_weight
			to_add = list([neighbor])
	if len(to_add)>0:
		tmp_curr_sol.append(to_add[0])
	return tmp_curr_sol, max_val

#return the best solution found from an edge
def best_solution_from_edge(mut2samples, sample2class, tot_class0, tot_class1, mut2neighbors, mut_white_list, starting_edge):
	old_sol = starting_edge
	new_sol, new_val = extend_solution(mut2samples, sample2class, tot_class0, tot_class1, mut2neighbors, mut_white_list, old_sol)
	#print "new_sol: "+str(new_sol)+", new val: "+str(new_val)
	#print "old_sol: "+str(old_sol)
	while len(new_sol) > len(old_sol) and len(new_sol) < max_size:
		#print "HERE!"
		old_sol = new_sol
		new_sol, new_val = extend_solution(mut2samples, sample2class, tot_class0, tot_class1, mut2neighbors, mut_white_list, old_sol)
		#print "new_sol: "+str(new_sol)+", new val: "+str(new_val)
        	#print "old_sol: "+str(old_sol)
	return new_sol, new_val

#generate permuted dataset
def generate_permuted_dataset( mut2samples, sample2class, samples_with_muts_AND_class, mut_white_list):
	#only use white list mutations in generating the permuted dataset
	permuted_mut2samples = dict()
	samples_class0 = set()
	samples_class1 = set()
	for elem in samples_with_muts_AND_class:
		if sample2class[elem] == "0":
			samples_class0.add(elem)
		if sample2class[elem] == "1":
                        samples_class1.add(elem)
	for mut in mut_white_list:
		num_samples_class0 = 0 
		num_samples_class1 = 0
		for elem in mut2samples[mut]:
			if sample2class[elem] == "0":
				num_samples_class0 += 1
			if sample2class[elem] == "1":
                                num_samples_class1 += 1
		#the following block of code is to do a permutation test where the frequency
		#of mutation in each class is preserved
		#BEGINS
		#new_samples_class0 = random.sample(samples_class0, num_samples_class0)
		#new_samples_class1 = random.sample(samples_class1, num_samples_class1)
		#permuted_mut2samples[mut] = list()
		#for sample in new_samples_class0:
		#	permuted_mut2samples[mut].append(sample)
		#for sample in new_samples_class1:
                #        permuted_mut2samples[mut].append(sample)
		#END
		#the following block of code is to do a permutation test where the frequency
		#of mutation in each class is NOT PRESERVED
		#BEGINS
		samples_both_classes = samples_class0.union(samples_class1)
		new_samples_both_classes = random.sample(samples_both_classes, num_samples_class0+num_samples_class1)
		permuted_mut2samples[mut] = list()
                for sample in new_samples_both_classes:
                        permuted_mut2samples[mut].append(sample)
		#END
	return permuted_mut2samples

#generate permuted dataset
def generate_permuted_dataset2( sample2class, samples_with_muts_AND_class):
	samples_list = list(samples_with_muts_AND_class)
	classes_labels = list()
	for i in range(len(samples_list)):
		classes_labels.append(sample2class[samples_list[i]])
	random.shuffle(classes_labels)
	permuted_sample2class = dict()
	for i in range(len(samples_list)):
		permuted_sample2class[samples_list[i]] = classes_labels[i]
	return permuted_sample2class

#INPUT PARAMETERS

#path to file with mutations, format: one line = one sample, in each line: sampleID\tgeneID1\tgeneID2... 
mut_matrix_file = sys.argv[1]
#path to file with class of each sample, format: first line is header ("SAMPLE_ID\tCLASS"), following lines
#have:sampleID\tClass
class_file = sys.argv[2]
#path to file describing the network, format: one line = one edge, in each line:geneID1\tgeneID2 
network_file = sys.argv[3]
#path to file where to write output
output_file = sys.argv[4]
#minimum frequency threshold: a gene is retained in the analysis if it has frequency above min_freq in
#at least one class
min_freq = float(sys.argv[5])
#minimum threshold for sets (see the paper)
threshold = float(sys.argv[6])
#path to file with blacklisted mutations (mutations to be excluded from analysis); format: one line = one gene,
#in each line: geneID
muts_black_list_file = sys.argv[7]
#maximum size of sets to be produced in output
max_size = int(sys.argv[8])
#number of permutations for permutation testing
num_permutations = int(sys.argv[9])

sample2muts = get_sample2muts(mut_matrix_file)

sample2class = get_sample2class(class_file)

#open the output file
outfile = open(output_file,'w')

#count the number of samples with both mutations and class
samples_with_muts = set(sample2muts.keys())
samples_with_class = set(sample2class.keys())
samples_with_muts_AND_class = samples_with_muts.intersection(samples_with_class)
outfile.write("Number of samples with both mutations and class: "+str(len(samples_with_muts_AND_class))+"\n")

#now build a gene2sample dictionary using only the samples with both mutations and class

#get the list of all genes and build the dictionary to get in which samples they are mutated
mut2samples, all_muts, total_samples_class0, total_samples_class1 = get_muts2sample(sample2muts, sample2class)

outfile.write("Number of samples in class 0: "+str(total_samples_class0)+"\n")
outfile.write("Number of samples in class 1: "+str(total_samples_class1)+"\n")

#get the list of adjacencies
mut2neighbors = get_adj_matrix_network(network_file)

#get the white list of mutations
mut_white_list = get_mut_white_list(mut2samples, muts_black_list_file)

#get all the pairs of edges with discrepancy above a given value
edges = get_edges_above_threshold(mut2samples, sample2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list)
#print len(edges)
#print edges
#now keep adding nodes to a starting edge until you cannot improve

best_sol = list()
best_val = 0.0

for starting_edge in edges:
	best_sol_edge, best_val_edge = best_solution_from_edge(mut2samples, sample2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list, starting_edge)
	if best_val_edge > best_val:
		best_val = best_val_edge
		best_sol = best_sol_edge

outfile.write("Best solution: "+str(best_sol)+"\n")
outfile.write("Best weight: "+str(best_val)+"\n")

p_val = 0.

#now perform permutation test number 1
for i in range(num_permutations):
	permuted_dataset_mut2samples = generate_permuted_dataset( mut2samples, sample2class, samples_with_muts_AND_class, mut_white_list)
	#get all the pairs of edges with discrepancy above a given value
	edges = get_edges_above_threshold(permuted_dataset_mut2samples, sample2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list)
	#now keep adding nodes to a starting edge until you cannot improve
	best_val_perm  = 0.0

	for starting_edge in edges:
        	best_sol_edge_perm, best_val_edge_perm = best_solution_from_edge(permuted_dataset_mut2samples, sample2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list, starting_edge)
        	if best_val_edge_perm > best_val_perm:
         		best_val_perm = best_val_edge_perm
	if best_val_perm >= best_val:
		p_val += 1.
p_val += 1.
p_val = p_val/float(num_permutations+1)

outfile.write("p-value 1: "+str(p_val)+"\n")

p_val = 0.
#now perform permutation tes number 2
for i in range(num_permutations):
        permuted_dataset_samples2class = generate_permuted_dataset2( sample2class, samples_with_muts_AND_class)
        #get all the pairs of edges with discrepancy above a given value
        edges = get_edges_above_threshold(mut2samples, permuted_dataset_samples2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list)
        #now keep adding nodes to a starting edge until you cannot improve
        best_val_perm  = 0.0

        for starting_edge in edges:
                best_sol_edge_perm, best_val_edge_perm = best_solution_from_edge(mut2samples, permuted_dataset_samples2class, total_samples_class0, total_samples_class1, mut2neighbors, mut_white_list, starting_edge)
                if best_val_edge_perm > best_val_perm:
                        best_val_perm = best_val_edge_perm
        if best_val_perm >= best_val:
                p_val += 1.
p_val += 1.
p_val = p_val/float(num_permutations+1)

outfile.write("p-value 2: "+str(p_val)+"\n")
outfile.close()
