#---------------------------------------------------
# The function GetText takes a name of a FASTA file with an extention and extracts DNA string from it.
# It outputs a DNA string without metadata.
def GetText(finName):
	fin = open(finName, 'r');
	text = '';
	for line in fin:
		if line[0] != '>':
			line = line.strip()
			text = text + line;
	fin.close();

	return text;	


#---------------------------------------------------
# The function CreateDict takes the following inputs:
#	text	DNA string
#	m 	length of a window
#	gtype	can be either "c" for a circular genome or "l" for a linear genome (a genome is treated as circular by default)
#It creates a dictionary of unique strings of length m for a text as an output.

def CreateDict(text, mm, gtype = "c"):
	
	if mm > len(text):
		#print "N is bigger than genome size!!!";
		return {};
	d_g = dict()
	nn = len(text);
	
	if gtype == "c":
		text = text + text[0:(mm+1)];
		lastpos = nn;
	elif gtype == "l":
			lastpos = nn-mm+1;
	else:
		print "Is this genome linear or circular?";
		return d_g;
		
	for ii in range (lastpos):
		bl = text[ii:(ii+mm)]; 
		bl = bl.lower()
		if bl in d_g :
			d_g[bl] = d_g[bl] + 1;
		else:
			d_g[bl] = 1;

	return d_g;	

	
#-----------------------------------------------
# The function FindIntersection takes as inputs two dictionaries, d_g11 and d_g22, and computes an intersection between them	
# It outputs a list of unique strings in intersection of the two dictionaries.
# It is assumed that the length of unique strings in both dictionaries is the same.	

def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;
	
	
	
#------------------------------------------------
# This function finds numbers of nucleotides within a given sequence.
def CountACGT(seq):
	seq = seq.lower();
	n_a = seq.count("a");
	n_c = seq.count("c");
	n_g = seq.count("g");
	n_t = seq.count("t");
	n_other = len(seq) - (n_a + n_c + n_g + n_t);
	if n_other > 0:
		print "Sequence has symbols other than a,c,g,t!!!"
	return [n_a, n_c, n_g, n_t, n_other];
	
# This function generates a random sequence of length n using init_d as initial settings for the random number generating
def GetRandomSeq(init_d, n, range = 100):
	p = 0;
	t = list();
	t_array = list();
	lim = 0;
	for s in init_d:
		t.append(s);
		p = p + init_d[s];
		st = lim;
		lim = lim+ init_d[s]*range;
		while st < lim:
			#print st, lim
			t_array.append(s);
			#print t_array;
			st = st+1;
	num = len(t_array);
	#print t_array;	
	
	if (1-p) > 0.001:
		print "Check frequences for alphabet"
	
	seq = "";
	j = 0;
	while j < n:
		r = randint(0, num -1);
		seq = seq + t_array[r];
		j = j+1;
	

	#print t;

	return seq;



# This function generates a number (equal to mm) of random files of length nn using uniform distribution of nucleotides
# It outputs a name of a file that contains the list of names of generated random files
 	
def GetRandomFiles_Uniform(nn,mm = 10,label = ""):
	fOutLname = "(Unif)"+label+"list.txt";
	fOutLname1 = "(Unif)"+label+"list(info).txt";
	fOutL = open(fOutLname, "w");
	fOutL1 = open(fOutLname1, "w");
	
	scale = 100;
	
	fOutL1.write("Name"+ "\t" + "a"+ "\t" + "c"+ "\t" + "g"+ "\t" + "t"+ "\t"+ "Length"+ "\n");
	fOutL1.write("Ratio"+ "\t" + "25.0"+ "\t" + "25.0"+ "\t" + "25.0"+ "\t" + "25.0"+ "\t"+ str(nn)+ "\n");

	for j in range(mm):
		#ffName = 'ran'+str(j+1)+'_' + str(nn) + '.txt';
		ffName = label +'_ran(unf)'+str(j+1)+ '.txt';
		ffOut = open(ffName, 'w');
		
		n_a = 0;
		n_c = 0;
		n_g = 0;
		n_t = 0;
		

			
		for i in range(nn):
			num = randint(0,3)
			if num == 0 :
				ffOut.write('A');
				n_a  = n_a + 1;
			elif num == 1:
				ffOut.write('C');
				n_c  = n_c + 1;
			elif num == 2:
				ffOut.write('G');
				n_g  = n_g + 1;
			elif num == 3:
				ffOut.write('T');
				n_t  = n_t + 1;
			else:
				print 'Error! No such symbol'
		ffOut.close()
		fOutL.write(ffName+"\n");
		
		nn_a = round(1.0*n_a/nn*scale,3);
		nn_c = round(1.0*n_c/nn*scale,3);
		nn_g = round(1.0*n_g/nn*scale,3);
		nn_t = round(1.0*n_t/nn*scale,3);
		
		#fOutL1.write(ffName+"\t"+str(n_a)+"\t"+str(n_c)+"\t" + str(n_g)+"\t"+ str(n_t) + "\t"+ str(nn) + "\n");
		fOutL1.write(ffName+"\t"+str(nn_a)+"\t"+str(nn_c)+"\t" + str(nn_g)+"\t"+ str(nn_t) + "\t"+ str(nn) + "\n");
	
	fOutL.close();
	fOutL1.close();
	return fOutLname;	
	

# This function converts frequencies of nucleotides into d_ini (initial settings for random generator)	
def GetD_ini(t_freq):
	d_ini = dict();
	if len(t_freq) == 4:
		total = t_freq[0]+t_freq[1]+t_freq[2]+t_freq[3];
		d_ini["a"] = round(1.0*t_freq[0]/total,4);
		d_ini["c"] = round(1.0*t_freq[1]/total,4);
		d_ini["g"] = round(1.0*t_freq[2]/total,4);
		d_ini["t"] = round(1.0*t_freq[3]/total,4);
	else:
		print "Check frequences t_ACGTfreq, more than 4 items!";
	return d_ini;



# This function generates a number of random files specified by nruns with uniform or skew distribution (distr = "unif" or distr = "skew") of nucleotides.
# The length of each random file and skewness of distribution match to the input file specified by finName	
# It outputs a file that contains a list of names of randomly generated files
	
def GetRandomFiles(distr = "unif", finName = "Phage_type10_E.coli_O157.fasta", nruns = 100, label = "", range = 100):
	
	text = GetText(finName);
	
	nn = len(text);
	
	scale = 100;
	
	if distr == "unif":
		fListName = GetRandomFiles_Uniform(nn,nruns,label);
		return fListName;
	elif distr == "skew":
		t_ACGTfreq = CountACGT(text);
		total = t_ACGTfreq[0]+t_ACGTfreq[1]+t_ACGTfreq[2]+t_ACGTfreq[3];
		if total != nn:
			print "Additional symbols are skipped while generating the distribution!!!", nn-total, t[4], "of", nn;
		
		print t_ACGTfreq[0:4];
		d_ini = GetD_ini(t_ACGTfreq[0:4]);
		
		fOutLname = "(Skew)"+label+"list.txt";
		fOutLname1 = "(Skew)"+label+"list(info).txt";
		
		fOutL = open(fOutLname, "w");
		fOutL1 = open(fOutLname1, "w");
		j = 0;
		
		fOutL1.write("Name"+ "\t" + "a"+ "\t" + "c"+ "\t" + "g"+ "\t" + "t"+ "\t"+ "Length"+ "\n");
		fOutL1.write(finName + "\t" + str(d_ini["a"]*scale)+ "\t" + str(d_ini["c"]*scale) + "\t" + str(d_ini["g"]*scale) + "\t" + str(d_ini["t"]

*scale)+ "\t"+ str(nn) + "\n");
		
		
		while j < nruns:
			ffName = label +'_ran(skw)'+str(j+1)+ '.txt';
			ffOut = open(ffName, 'w');
		
			seq = GetRandomSeq(d_ini, nn, range);
			
			tt_acgt = CountACGT(seq);
					
			ffOut.write(seq);
			
			#fOutL1.write(ffName+"\t"+str(tt_acgt[0])+"\t"+str(tt_acgt[1])+"\t" + str(tt_acgt[2])+"\t"+ str(tt_acgt[3]) + "\t"+ str(nn) + "\n")

			fOutL1.write(ffName+"\t"+str(round(1.0*scale*tt_acgt[0]/nn,3))+"\t"+str(round(1.0*scale*tt_acgt[1]/nn,3))+"\t" + str(round

(1.0*scale*tt_acgt[2]/nn,3))+"\t"+ str(round(1.0*scale*tt_acgt[3]/nn,3)) + "\t"+ str(nn) + "\n")

			
			
			ffOut.close()
			fOutL.write(ffName+"\n");
			j = j + 1;
		fOutL.close();
		fOutL1.close();
		return fOutLname;


		

# It generates a number of random files specified by nruns with uniform or skew distribution (distr = "unif" or distr = "skew") of nucleotides.
# The length of each random file equals to nn and skewness of distribution is specified by d_ini to the input file specified by finName
# It outputs a file that contains a list of names of randomly generated files

def GetRandomFilesByLength(distr = "unif", nruns = 1, label = "", range = 100, nn = 0, d_ini = {}):
	
	#text = GetText(finName);
	
	#nn = len(text);
	
	scale = 100;
	
	if distr == "unif":
		fListName = GetRandomFiles_Uniform(nn,nruns,label);
		return fListName;
	elif distr == "skew":

		if d_ini == {}:
			print "Error!!! Initial parameters for a random generator are missing";
			return;
		
		fOutLname = "(Skew)"+label+"list.txt";
		fOutLname1 = "(Skew)"+label+"list(info).txt";
		
		fOutL = open(fOutLname, "w");
		fOutL1 = open(fOutLname1, "w");
		j = 0;
		
		fOutL1.write("Name"+ "\t" + "a"+ "\t" + "c"+ "\t" + "g"+ "\t" + "t"+ "\t"+ "Length"+ "\n");
		fOutL1.write(finName + "\t" + str(d_ini["a"]*scale)+ "\t" + str(d_ini["c"]*scale) + "\t" + str(d_ini["g"]*scale) + "\t" + str(d_ini["t"]

*scale)+ "\t"+ str(nn) + "\n");
		
		
		while j < nruns:
			ffName = label +'_ran(skw)'+str(j+1)+ '.txt';
			ffOut = open(ffName, 'w');
		
			seq = GetRandomSeq(d_ini, nn, range);
			
			tt_acgt = CountACGT(seq);
					
			ffOut.write(seq);
			
			#fOutL1.write(ffName+"\t"+str(tt_acgt[0])+"\t"+str(tt_acgt[1])+"\t" + str(tt_acgt[2])+"\t"+ str(tt_acgt[3]) + "\t"+ str(nn) + "\n")

			fOutL1.write(ffName+"\t"+str(round(1.0*scale*tt_acgt[0]/nn,3))+"\t"+str(round(1.0*scale*tt_acgt[1]/nn,3))+"\t" + str(round

(1.0*scale*tt_acgt[2]/nn,3))+"\t"+ str(round(1.0*scale*tt_acgt[3]/nn,3)) + "\t"+ str(nn) + "\n")

			
			
			ffOut.close()
			fOutL.write(ffName+"\n");
			j = j + 1;
		fOutL.close();
		fOutL1.close();
		return fOutLname;		
		





