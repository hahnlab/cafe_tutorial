#!/usr/bin/python
########################################################################################
# The new CAFE Report Analysis script.
#
# Gregg Thomas, Spring 2016
########################################################################################

import sys, os, argparse, cafecore as cafecore

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Analyzes a CAFE report file (.cafe)");

	parser.add_argument("-i", dest="report_input_file", help="A CAFE report file (.cafe).");
	parser.add_argument("-r", dest="rapids_out", help="1: Output file will be list of only rapidly changing families on each node. 0: Output file will be list of all changing families on each node. Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="large_fam_file", help="A CAFE report file for the large families of the same set of species.");
	parser.add_argument("-o", dest="output_prefix", help="A prefix string to put on all the output files generated.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.report_input_file == None or args.output_prefix == None:
			cafecore.errorOut(1, "A CAFE report file must be specified with -i and an output file name with -o");
			optParse(1);

		if args.rapids_out not in [1,0]:
			cafecore.errorOut(2, "-r can only take values of 0 or 1");
			optParse(1);

		return args.report_input_file, args.rapids_out, args.large_fam_file, args.output_prefix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

#######################

def formatLineParse(line):
# This function handles CAFE's weird node pair format for the p-values and node ids.

	if "=" in line:
		line = line.split("=")[1];
	if ":" in line:
		line = line.split(": ")[1].strip();
	line = line.replace("(", "").replace(")", "");
	line = line.split(" ");
	line = [f.split(",") for f in line];

	return line;

#######################

def nodeRelabel(treedict):
# Family trees are read with gene counts on the tip labels. This function removes them.

	tmp = {};
	#print treedict;

	for oldkey in treedict:
		if treedict[oldkey][3] == 'tip':
			newkey = oldkey[:oldkey.index("_")];
			tmp[newkey] = treedict[oldkey];
		else:
			tmp[oldkey] = treedict[oldkey];

	return tmp;

#######################

def nodeMap(cafetd, mytd):
# CAFE has pre-determined mappings in the tree. When I read the tree with my own script the mappings
# are different. This function creates a map from my node ids to CAFE's node ids.
# The dictionary nodemap has the following {key:value} format: {my node id:CAFE's node id}

	nodemap = {};
	# The map dictionary.

	##############
	# for node in cafetd:
	# 	if cafetd[node][3] == 'tip':
	# 		spec = node[:node.index("<")];
	# 		cafeid = node[node.index("<")+1:node.index(">")];
	# 		nodemap[cafeid] = spec;

	# while len(nodemap) != len(cafetd):
	# 	for node in cafetd:
	# 		if cafetd[node][3] == 'root':
	# 			continue;

	# 		orignode = node;
	# 		node = node[node.index("<")+1:node.index(">")];

	# 		if node in nodemap:
	# 			if cafetd[orignode][3] == 'tip':
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd[nodemap[node]][1];
	# 			else:
	# 				curanc = cafetd[orignode][1];
	# 				mapanc = mytd["<" + nodemap[node] + ">"][1];

	# 			nodemap[curanc.replace("<","").replace(">","")] = mapanc.replace("<","").replace(">","");
	##############
	# The above formats nodemap with the reverse {key:value} format: {CAFE's node id:my node id}

	for node in cafetd:
		if cafetd[node][3] == 'tip':
			spec = node[:node.index("<")];
			cafeid = node[node.index("<"):];
			nodemap[spec] = node;
	# First map the tips by their unique species labels.

	while len(nodemap) != len(mytd):
		for node in mytd:
			if mytd[node][3] == 'root':
				continue;

			if node in nodemap:
				curanc = mytd[node][1];
				mapanc = cafetd[nodemap[node]][1];

				nodemap[curanc] = mapanc;
	# Then do a post-order traversal and map the current node's ancestor to it's map's ancestor.

	return nodemap;

#######################
def cra(inlines, results, node_fams, linestart, afilename, s_nodes, v):
	numbars = 0;
	donepercent = [];
	i = 0;
	acount = 0;
	afile = open(afilename, "a");

	for inline in inlines:
	# Each line of the report file is read.
		if v == 1:
			numbars, donepercent = cafecore.loadingBar(i, len(inlines), donepercent, numbars);
		i = i + 1;

		if i <= linestart:
			continue;
		# If the line is a CAFE info line, skip it.

		inline = inline.strip().split("\t");
		famid = inline[0];
		famtree = inline[1];
		nodeformat = inline[3].replace("),(", ") (");
		# Parsing the information for the current family.

		outline = famid + "\t";
		outlist = [0 for n in s_nodes];
		# Prep for the anc states file.

		nodes = formatLineParse(nodeformat);

		tlinfo, newfamtree = cafecore.treeParseNew(famtree, 1);
		# Reading the tree and adding my node labels.

		for tlnode in tlinfo:
			if tlinfo[tlnode][3] == 'root':
				tlinfo[tlnode].append(famtree[famtree.rfind("_")+1:]);
			elif tlinfo[tlnode][3] == 'tip':
				tlinfo[tlnode].append(tlnode[tlnode.index("_")+1:]);
			else:
				tlinfo[tlnode][4] = tlinfo[tlnode][4][1:];
		# Gene counts for each node are read as support values for internal nodes, but must
		# have the underscore removed. Tip and root node counts are added here as well.

		tlinfo = nodeRelabel(tlinfo);
		# Removes the gene counts from the tip node labels.

		if i == (linestart + 1):
			maps = nodeMap(tinfo, tlinfo);
		# If this is the first family, we need to build our maps from my node ids to CAFE's.

		for tlnode in tlinfo:
		# For each node in the current gene family tree, we make our counts.

			if tlinfo[tlnode][3] == 'root':
				continue;
			# No count is made at the root of the tree.

			curanc = tlinfo[tlnode][1];
			curmap = maps[tlnode];
			# Get the ancestor and the map of the current node.

			curcount = int(tlinfo[tlnode][4]);
			anccount = int(tlinfo[curanc][4]);
			# Get the gene counts of the current node and the ancestor.

			outlist[s_nodes.index(curmap)] = str(curcount);
			# Save the count of the current node to be sent to the anc states file.

			diff = curcount - anccount;
			# Calculate the difference in gene count.

			typeflag = 0;
			# typeflag tells whether an expansion or contraction has occurred.

			if curcount > anccount:
				typeflag = 1;
				results[curmap][0] += 1;
				results[curmap][1] += diff;

				if r_opt == 0:
					node_fams[curmap][0].append(famid + "[+" + str(diff) + "]");

					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
			# If the difference in gene count between the current node and the ancestor is positive, an
			# expansion has occurred. This makes the appropriate counts.

			elif curcount < anccount:
				typeflag = 2
				results[curmap][3] += 1;
				results[curmap][4] += abs(diff);

				if curcount == 0 and anccount != 0:
					results[curmap][5] += 1;

				if r_opt == 0:
					node_fams[curmap][1].append(famid + "[" + str(diff) + "]");

					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
			# If the difference in gene count between the current node and the ancestor is negative, a
			# contraction has occurred. This makes the appropriate counts. It also checks for family losses
			# along that branch by seeing if the current node has transitioned to a count of 0.

			elif curcount == anccount:
				results[curmap][2] += 1;
			# Otherwise, the counts at the current node and the ancestor are the same and no change has occurred.

			if float(inline[2]) < 0.01:
			# If the family p-value is below a threshold, the family is rapidly evolving.

				if r_opt == 1:
					if famid not in node_fams['total']:
						node_fams['total'].append(famid);
				# Add the family id to the 'total' key of node_fams. This also parses the nodeformat line which is
				# in the paired CAFE node format.
				
				pairnodeid = curmap[curmap.index("<")+1:curmap.index(">")];
				# Since the paired format does not include the brackets that the other node labels do, I have to
				# remove them to check against that format.

				for j in range(len(nodes)):
					for k in range(len(nodes[j])):
						if formatline[j][k] == pairnodeid and float(nodes[j][k]) < 0.01:
							if typeflag == 1:
								results[curmap][7] += 1;
								if r_opt == 1:
									node_fams[curmap][0].append(famid + "[+" + str(diff) + "*]");
								elif r_opt == 0:
									node_fams[curmap][0].pop();
									node_fams[curmap][0].append(famid + "[+" + str(diff) + "*]");
							elif typeflag == 2:
								results[curmap][8] += 1;
								if r_opt == 1:
									node_fams[curmap][1].append(famid + "[" + str(diff) + "*]");
								elif r_opt == 0:
									node_fams[curmap][1].pop();
									node_fams[curmap][1].append(famid + "[" + str(diff) + "*]");
							results[curmap][9] += 1;
				# Runs through the paired format as a list of lists. If the p-value of that node is less than a threshold
				# that branch is rapidly evolving. Based on typeflag, the appropriate counts are made. The family id is
				# also added to the current node in rapids.

		outline += "\t".join(outlist) + "\n";
		afile.write(outline);
		# Write the states of the current family to the anc states file

	if v == 1:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);
	afile.close();
	return results, node_fams;

############################################
#Main block
############################################

infilename, r_opt, largefilename, outprefix = optParse(0);
#Get the input parameters.

famfilename = outprefix + "_fams.txt";
nodefilename = outprefix + "_node.txt";
pubfilename = outprefix + "_pub.txt";
ancfilename = outprefix + "_anc.txt";

print "=======================================================================";
print "\t\tCAFE Report File Analysis"
print "\t\t" + cafecore.getDateTime();
print "---------";
print "Parsing format information...\n";

infile = open(infilename, "r");
inlines_main = infile.readlines();
infile.close();
# Reads the input report file.

if inlines_main[2].find("Lambda tree:") != -1:
	treeline = inlines_main[3];
	formatline = inlines_main[4];
	avgline = inlines_main[6];
	linestart_main = 11;
else:
	treeline = inlines_main[2]
	formatline = inlines_main[3];
	avgline = inlines_main[5];
	linestart_main = 10;
# If CAFE was run with a lambda tree structure as input, the report file places that on the third
# line. This shifts all the other relevant lines down by 1. This if/else accounts for that.

labeled_tree = treeline[treeline.index(":")+1:].strip();
tinfo, newtree = cafecore.treeParseNew(labeled_tree,2);
# This reads the CAFE tree with its node labels.

formatline = formatLineParse(formatline);
# formatline is CAFE's line with its paired node format with node ids. The formatLineParse function
# reads that format and returns it as a list of lists.

avgline = avgline.split(":\t")[1].strip().replace("\t", " ");
avgline = formatLineParse(avgline);
# The line of average expansions for each node, in the paired node format. Again passed to formatLineParse
# to make it interpretable.

print "---------";
print "Initializing output structures...\n";

node_fams_main = {"total" : []};
results_main = {};

sorted_nodes = [];
ancfile = open(ancfilename, "w");
header = "Family ID\t";

for node in tinfo:
	if tinfo[node][3] == 'root':
		continue;
	node_fams_main[node] = [[],[]];
	results_main[node] = [0,0,0,0,0,0,0,0,0,0];

	sorted_nodes.append(node);
	header += node + "\t";

ancfile.write(header[:-1] + "\n");
ancfile.close();
# [expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]
# node_fams and results are the two main dictionaries to store CAFE's results.
# node_fams {key:value} format: {node:list of two lists containing family ids for (rapid) expansions and (rapid) contractions, respectively}
# results {key:value} format: {node:[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]}
# This loop also does the prepping of the header and sorted nodes for the anc count file

for j in range(len(formatline)):
	for k in range(len(formatline[j])):
		n = "<" + formatline[j][k] + ">";
		for r in results_main:
			if n in r:
				results_main[r][6] = avgline[j][k];
# Setting average expansion in results as read from avgline

print "---------";
print "Counting changes per branch...\n";

results_main, node_fams_main = cra(inlines_main, results_main, node_fams_main, linestart_main, ancfilename, sorted_nodes, 1);

if largefilename != None:
	print "\n\n---------";
	print "Parsing large families...\n";

	lfile = open(largefilename, "r");
	llines_main = lfile.readlines();
	lfile.close();
	# Reads the input report file.

	if llines_main[2].find("Lambda tree:") != -1:
		linestart_main = 11;
	else:
		linestart_main = 10;

	results_main, node_fams_main = cra(llines_main, results_main, node_fams_main, linestart_main, ancfilename, sorted_nodes, 0);

print "\nDone!";
print "=======================================================================";

print "Writing output files...";

## Begin fam output block.
outfile = open(famfilename, "w");
outfile.write("");
# Initialize the output file.
outfile.write("# The labeled CAFE tree:\t" + labeled_tree + "\n");

if r_opt == 0:
	desc_str = " ";
elif r_opt == 1:
	desc_str = " rapid ";

outline = "Overall" + desc_str + ":\t"
for f in node_fams_main['total']:
	outline = outline + f + ",";
outline = outline[:-1] + "\n";
outfile.write(outline);

for spec in node_fams_main:
	if spec == 'total':
		continue;

	# for f in range(len(node_fams_main[spec])):
	# 	if f == 0:
	# 		outline = spec + desc_str + "expansions:\t";
	# 	elif f == 1:
	# 		outline = spec + desc_str + "contractions:\t";

	# 	for rapid_f in node_fams_main[spec][f]:
	# 		outline = outline + rapid_f + ",";
	# 	outline = outline[:-1] + "\n";
	# 	outfile.write(outline);
	# For output on separate lines for expansions and contractions per species.

	outline = spec + ":\t";
	outline += ",".join(node_fams_main[spec][0] + node_fams_main[spec][1]) + "\n";
	outfile.write(outline);
	# For output on a single line per species.

outfile.close();
## End fam output block


## Begin node and pub output block
nodefile = open(nodefilename, "w");
nodefile.write("Node\tExpansions\tContractions\tRapidly evolving families\n");

pubfile = open(pubfilename, "w");
pubfile.write("Species\tExpanded fams\tGenes gained\tgenes/expansion\tContracted fams\tGenes lost\tgenes/contraction\tNo change\tAvg. Expansion\n");

for node in results_main:
	outline = node + "\t" + str(results_main[node][0]) + "\t" + str(results_main[node][3]) + "\t" + str(results_main[node][9]) + "\n";
	nodefile.write(outline);

	if node.replace("<","").replace(">","").isdigit():
		continue;

	spec = node.title()[:node.index("<")];
	exp = results_main[node][0];
	con = results_main[node][3];
	outline = spec + "\t" + str(results_main[node][0]) + " (" + str(results_main[node][7]) + ")\t" + str(results_main[node][1]) + "\t";
	if exp != 0:
		outline += str(round(float(results_main[node][1])/float(exp),2)) + "\t";
	else:
		outline += '0' + "\t";
	outline += str(results_main[node][3]) + " (" + str(results_main[node][8]) + ")\t" + str(results_main[node][4]) + "\t"; 
	if con != 0:
		outline += str(round(float(results_main[node][4])/float(con),2)) + "\t";
	else:
		outline += '0' + "\t";
	outline += str(results_main[node][2]) + "\t" + str(results_main[node][6]) + "\n";
	pubfile.write(outline);

nodefile.close();
pubfile.close();
## End node and oub output block

## Begin plot block
# print "Generating plots...";

# x_nodes = [];
# y_rapids = [];
# y_changes = [];
# y_exp = [];
# y_rexp = [];
# y_con = [];
# y_rcon = [];
# for node in results_main:
# 	y_rapids.append(results_main[node][9]);
# 	y_changes.append((results_main[node][0]+results_main[node][3])-results_main[node][9]);

# 	y_rexp.append(results_main[node][7]);
# 	y_exp.append(results_main[node][0]-results_main[node][7]);

# 	y_rcon.append(results_main[node][8]);
# 	y_con.append(results_main[node][3]-results_main[node][8]);

# 	if node[0] != "<":
# 		node = node.title()[:node.index("<")];

# 	x_nodes.append(node);

# #barcols = ['#ffef52','#5b5bd7'];
# barcols = ['#e5653c', '#2aa064'];

# y_data = [y_rapids, y_changes];
# y_names = ['total rapids', 'changes'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# changing families","# of changing families",outprefix+"_change.html",barcols,w=1200);
# # Total plot

# y_data = [y_rexp, y_exp];
# y_names = ['rapid expansions', 'expansions'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# expanding families","# of expanding families",outprefix+"_expand.html",barcols,w=1200);
# # Expansion plot

# y_data = [y_rcon, y_con];
# y_names = ['rapid contractions', 'contractions'];
# crplot.barPlotStack(x_nodes,y_data,y_names,"","# contracting families","# of contracting families",outprefix+"_contract.html",barcols,w=1200);
# # Contraction plot
## End plot block

# node_fams {key:value} format: {node:list of two lists containing family ids for (rapid) expansions and (rapid) contractions, respectively}
# results {key:value} format: {node:[expand,gene expand,equal,contract,gene contract,families lost,avg expansion,sigexpand,sigcontract,total sig changes]}
print "RESULTS TABLE -- tab delimted for easy copy/pasting into your favorite spreadsheet program"
print "\tExpansions\tGenes Gained\tEqual\tContractions\tGenes Lost\tFamilies Lost\tAverage Expansion\tSig Expansions\tSig Contractions\tTotal Sig Changes";
for species in results_main:
	outline = species + "\t";
	for col in results_main[species]:
		outline = outline + str(col) + "\t";
	print outline;

print
print "CAFE labeled tree:\t" + labeled_tree;
# This block simply prints the information stored in results to the screen.
print "=======================================================================";


