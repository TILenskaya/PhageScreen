import CoreMethods as cm;

m = 40; # a string length for creating a dictionary

finName1 = 'E.coli_O157_H7.fasta'; # E.coli O157:H7
finName2 = 'KP869108.fasta'; # E.coli O157 typing phage 10

text1 = cm.GetText(finName1);
text2 = cm.GetText(finName2);

d1 = cm.CreateDict(text1,m); # Host dictionary
d2 = cm.CreateDict(text2,m); # Parasite dictionary

print "Number of unique strings in the intersection:", len(cm.FindIntersection(d1,d2)); # a number of unique strings of length m in the intersection between host and parasite dictionaries
