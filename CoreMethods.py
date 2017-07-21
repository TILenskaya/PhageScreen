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
#	m 		length of a window
#	gtype	can be either "c" for a circular genome or "l" for a linear genome (a genome is treated as circular by default)
#It creates a dictionary of unique strings of length m for a text as an output.


def CreateDict(text, m, gtype = "c"):
	d_g = dict();
	n = len(text);
	
	if gtype == "c":
		text = text + text;
		lastpos = n;
	elif gtype == "l":
			lastpos = n-m+1;
	else:
		print "Is this genome linear or circular?";
		return d_g;
		
	for i in range (lastpos):
		bl = text[i:(i+m)]; 
		bl = bl.lower()
		if bl in d_g :
			tt = d_g[bl];
			tt[0] = tt[0] + 1;
			d_g[bl] = tt;
		else:
			tt = list();
			tt.append(1);
			d_g[bl] = tt;
			
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


	