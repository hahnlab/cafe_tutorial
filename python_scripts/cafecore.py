#############################################################################
#CAFECORE functions
#Gregg Thomas
#August 2013-present, forked from CORE on 04.13.16
#############################################################################

import sys, datetime, subprocess

##########################################################################################################################################################
# IO functions

def loadingRotator(counter, rotate, divisor):
# Provides a loading rotator for loops. The following line must be used to initialize the function
# before the loop in the main code:
# rotator = 0;

	rotation = ['|', '/', '-', '\\'];

	if counter % divisor == 0:
		sys.stderr.write('\b' + rotation[rotate]);
		rotate = rotate + 1;
		if rotate >= len(rotation):
			rotate = 0;

	return rotate;

#############################################################################

def errorOut(errnum, errmsg):
# Formatting for error messages.
	fullmsg = "|**Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	print "\n" + border + "\n" + fullmsg + "\n" + border + "\n";

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

#############################################################################

def loadingBar(counter, length, done, bars):
# This function serves as a text loading bar for long scripts with counters. The following
# lines must be added within the script to initialize and terminate the script:
# Initilization:
# numlines = core.getFileLen(alnfilename);
# numbars = 0;
# donepercent = [];
# i = 0;
# Termination:
#	pstring = "100.0% complete.";
#	sys.stderr.write('\b' * len(pstring) + pstring);
#	print "\nDone!";
#
# If length is lines in a file use the core.getFileLen function to get the number of lines in the file

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "";
		loading = "[";
		j = 0;
		while j <= bars:
			loading = loading + "*";
			j = j + 1;
		while j < 50:
			loading = loading + "-";
			j = j + 1;
		loading = loading + "]";

		loading = loading + "                 ";
		sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;

	sys.stderr.write('\b' * len(pstring) + pstring);

	return bars, done;

#############################################################################

def getFileLen(i_name):
#Calls 'wc -l' to get the number of lines in the file.
	p = subprocess.Popen(['wc', '-l', i_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	result, err = p.communicate();
	if p.returncode != 0:
		raise IOError(err);
	return int(result.strip().split()[0]);

##########################################################################################################################################################
# Tree parsing functions.

def getBranchLength(bltree, spec_label):
# Returns the branch length of a species given a newick formatted tree. Used by treeParse.

	d = 0;
	startind = 0;

	while d < (len(bltree)-1):
		if bltree[d] == ":":
			current_node = bltree[max(bltree.rfind("(",startind,d),bltree.rfind(")",startind,d),bltree.rfind(",",startind,d))+1:d];
			if current_node == spec_label:

				opind = bltree.find("(",d);
				cpind = bltree.find(")",d);
				coind = bltree.find(",",d);

				indcheck = [opind,cpind,coind];

				for a in xrange(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################

def getDesc(d_spec, d_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds the direct descendants of the species.
	d_list = [];
	for node in d_treedict:
		if d_treedict[node][1] == d_spec:
			d_list.append(node);

	if d_list == []:
		return [d_spec];
	else:
		return d_list;

#############################################################################

def treeParseNew(tree, tree_type):
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# New (current) format:
# node:[branch length, ancestral node, ancestral branch length, node type, node labels (if present)]
#
# Old format
# node:[branch length, ancestral node, ancestral branch length, sister node, sister branch length, descendent 1, descendent 1 branch length, descendent 2, descendent 2 branch length, node type]
#
#
# Tree type 1: tree has branch lengths.
# Tree type 2: tree is just topology.

	tree = tree.replace("\n","");
	if tree[len(tree)-1:] != ";":
		tree = tree + ";";
	##Some string handling

	new_tree = "";
	z = 0;
	numnodes = 1;
	supports = {};

	while z < (len(tree)-1):
		if tree_type == 1:
			if tree[z] == ":" and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
			elif tree[z] == ":":
				tmp_str = tree[:z];
				if tmp_str.rfind(")") > tmp_str.rfind(","):
					new_node = "<" + str(numnodes) + ">";

					supports[new_node] = tmp_str[tmp_str.rfind(")")+1:];

					new_tree = new_tree[:new_tree.rfind(")")+1] + new_node;
					numnodes = numnodes + 1;

		if tree_type == 2:
			if (tree[z] == "," or tree[z] == ")") and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;

	if new_tree[-1] not in [")",">"]:
		if new_tree.rfind(")") > new_tree.rfind(">"):
			last_char = ")";
		else:
			last_char = ">";

		new_tree = new_tree[:new_tree.rfind(last_char)+1];

	if new_tree[-1] == ")":
		rootnode = "<" + str(numnodes) + ">"
		new_tree = new_tree + rootnode;
	else:
		rootnode = new_tree[new_tree.rfind(")")+1:];

	##This first block labels all internal nodes with the format <#>
	
	#print tree;
	#print new_tree;
	#print supports;
	#print "-----------------------------------";

	ancs = {};
	nofo = {};

	z = 0;
	startind = 0;
	while z < (len(new_tree)-1):
	##Here, the ancestral nodes of each node are found

		if tree_type == 1:
		##The major difference between trees with branch lengths (type 1) and without (type 2) is seen here. Finding the ancestral nodes requires
		##almost a completely different set of logic statements.
			if new_tree[z] == ":":
				curnode = new_tree[max(new_tree.rfind("(",startind,z),new_tree.rfind(")",startind,z),new_tree.rfind(",",startind,z))+1:z];
				numcpneeded = 1
				numcp = 0;
				nofo[curnode] = [];

				a = z;

				while a < (len(new_tree)-1):
					if new_tree[a] == "(":
						numcpneeded = numcpneeded + 1;
					if new_tree[a] == ")" and numcpneeded != numcp:
						numcp = numcp + 1;
					if new_tree[a] == ")" and numcpneeded == numcp:
						#if a == (len(new_tree)-5):
						#	curanc = new_tree[a+1:];
						if new_tree[a+1:].find(":") == -1:
							#curanc = new_tree[len(new_tree)-5:];
							curanc = new_tree[new_tree.rfind(")")+1:]
						else:
							curanc = new_tree[a+1:new_tree.index(":", a)];
						a = 100000000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		if tree_type == 2:
			if new_tree[z] == "," or new_tree[z] == ")":
				curnode = new_tree[max(new_tree.rfind("(",startind,z),new_tree.rfind(")",startind,z),new_tree.rfind(",",startind,z))+1:z];
				numcpneeded = 1
				numcp = 0;
				nofo[curnode] = [];

				a = z;

				while a < (len(new_tree)-1):
					if new_tree[a] == "(":
						numcpneeded = numcpneeded + 1;
					if new_tree[a] == ")" and numcpneeded != numcp:
						numcp = numcp + 1;
					if new_tree[a] == ")" and numcpneeded == numcp:
						if a == (len(new_tree)-4):
							curanc = new_tree[a+1:];
						else:
							mindex = 999999999;
							for c in ["(",")",","]:
								cind = new_tree.find(c,a+1);
								if cind < mindex and cind != -1:
									mindex = cind;
									minchar = c;
							curanc = new_tree[a+1:mindex];
						a = 10000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		z = z + 1;
	##End ancestral node block
#	print curanc;
	#for key in ancs:
	#	print key + ":", ancs[key]
	#print "---------";
	#sys.exit()

	##The next block gets all the other info for each node: sister and decendent nodes and branch lengths (if type 1)
	##and node type (tip, internal, root). This is easy now that the ancestral nodes are stored.
	nofo[rootnode] = [];
	for node in nofo:
		if tree_type == 1:
			cur_bl = getBranchLength(new_tree,node);
		elif tree_type == 2:
			if node == rootnode:
				cur_bl = None;
			else:
				cur_bl = "NA";
		nofo[node].append(cur_bl);

		if node != rootnode:
			cur_anc = ancs[node];
			nofo[node].append(cur_anc);
			if tree_type == 1:
				cur_anc_bl = getBranchLength(new_tree,cur_anc);
			elif tree_type == 2:
				cur_anc_bl = "NA";
			nofo[node].append(cur_anc_bl);
		else:
			j = 0;
			while j < 2:
				nofo[node].append("");
				j = j + 1;

	for node in nofo:
		if node == rootnode:
			nofo[node].append("root");
		elif getDesc(node,nofo) == [node]:
			nofo[node].append("tip");
		else:
			nofo[node].append("internal");

		if nofo[node][3] != 'tip' and supports != {} and node in supports:
			nofo[node].append(supports[node]);

	##End info retrieval block.

#	for key in nofo:
#		print key + ":" + str(nofo[key]);


	return nofo, new_tree;

#############################################################################