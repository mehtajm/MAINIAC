// Mafft converter for sequence insertion displays
// Jai Mehta
// 03/26/2025
// Given a .map file, pre-mafft, and post-mafft output, this program will add the necessary
// insertions and will also offset the gaps for other sequences.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <vector>
using namespace std;

struct mapInserts {
	int preMafftResStart, preMafftResEnd, postMafftIns, postMafftEnd;
};
struct seqParts {
	string seq;
	int location;
	bool operator==(const int other) const{
		return location == other;
	}
};

bool lessSeqParts(const seqParts& s1, const seqParts& s2) {
	return s1.location < s2.location;
}

struct postMafftInsertions {
	vector<seqParts> insertions;
};
mapInserts Parse(string line) { // parses the map file to grab the insertion points
	string lower, upper;
	istringstream s3(line);
	char delim = '\n';
	for (int i = 0; i < line.size(); i++) {
		if (line[i] >= 'A' && line[i] <= 'Z') {
			delim = line[i];
			break;
		}
	}
	mapInserts a;
	getline(s3, lower, delim); // gets first range;
	a.preMafftResStart = stoi(lower);
	getline(s3, lower, '-');
	s3 >> a.preMafftResEnd;
	getline(s3, upper, '>');
	s3 >> a.postMafftIns;
	getline(s3, upper, 'v');
	s3 >> a.postMafftEnd; // don't think this is necessary
	return a;
}

int main(int argc, char **argv) {
	if (argc < 5 || argc > 6) {
		cerr << "Usage: ./mafftConvert [Pre-mafft fasta] [.map output] [mafft output fasta (non-existent for j seq)] [output file] ['j' for j seq (optional)]\n";
		return 1;
	}
	bool jSeq = false;
	if (argc == 6 && tolower(argv[5][0]) == 'j') jSeq = true;
	ifstream preMafftFin(argv[1]);
	if (preMafftFin.fail()) {
		cerr << "Pre-mafft fasta did not open properly. Try again.\n";
		return 1;
	}
	ifstream mapFin(argv[2]);
	if (mapFin.fail()) {
		cerr << "map output did not open properly. Try again.\n";
		return 1;
	}
	ifstream postMafftFin(argv[3]);
	if (jSeq) {
		int failCount = 0;
		while (postMafftFin.fail()) {
			if (failCount >= 2) {
				cerr << "post Mafft output fasta failed to be written to and opened twice. Run the program again with a different pre-mafft fasta.\n";
				return 1;
			}
			cerr << "Mafft output fasta file did not open properly. Writing an empty one for assumed J sequence.\n";
			cerr << "If this is not what you wanted, hit ctrl + c and try again.\n";
			ofstream tempOut(argv[3]);
			if (tempOut.fail()) {
				cerr << "The non-existent file had a bad directory. Make sure the directory is correct if you would like to grab the j sequence.\n";
				return 1;
			}
			string tempLine;
			while (getline(mapFin, tempLine)) {
				if (tempLine[0] == '>') {
					tempOut << tempLine << endl << endl;
				}
			}
			mapFin.clear();
			mapFin.seekg(0); // go back to start of file
			tempOut.close();
			postMafftFin.clear();
			postMafftFin.ignore();
			postMafftFin.seekg(0);
			failCount++;
			postMafftFin.open(argv[3]);
		}
	}
	else {
		if (postMafftFin.fail()) {
			cerr << "Mafft output fasta file did not open properly. Try again.\n";
			return 1;
		}
	}
	ofstream fout(argv[4]); // output csv/fasta file
	if (fout.fail()) {
		cerr << "Output file failed to open properly. Try again.\n";
		return 1;
	}
	ofstream jfout;
	string output = (string)argv[4];
	if (jSeq) {
        if (output.size() < 3 || output.substr(output.size() - 3, 3) != ".fa") {
            cerr << "WARNING: Your output file \"" << output << "\" does not contain \".fa\" file type. ";
            cerr << "The cleaved J sequences will be added to \"" << output << "_jSeqs.fa\", while your output file will remain the same.\n";
        }
        else {
            output.resize(output.size() - 3);
        }
        output += "_jSeqs.fa";
        jfout.open(output);
    }
	// first parse through map csv
	unordered_map<string, vector<mapInserts> > changes;
	unordered_map<int, int> masterset; // contains all input locations and their quantities.
	unordered_map<string, postMafftInsertions> postInserts; // final changes
	string line;
	for (int i = 0; i < 2; i++) getline(mapFin, line);  // extracts id range line
	string ranges;
	while (getline(mapFin, ranges)) {
		if (!ranges.empty() && ranges[0] == '>') {
			line = ranges;
			getline(mapFin, ranges);
		}
		string qid = line.substr(1, line.size() - 1);
		while (qid.back() == '\t' || qid.back() == ' ') qid.pop_back();
		mapInserts m = Parse(ranges);
		int diff = m.preMafftResEnd - m.preMafftResStart + 1;
		if (masterset.count(m.postMafftIns) > 0) {
			if (masterset[m.postMafftIns] < diff) masterset[m.postMafftIns] = diff; // updates highest quantity given postMafft ins. location
		}
		else masterset[m.postMafftIns] = diff;
		changes[qid].push_back(m);
	} // takes in all of the changes for each insertion for an id on the map file
	  // now go through pre-mafft ids, and take in the sequence parts 
	string id;
	int seqLength = -1;
	getline(preMafftFin, line);
	while (getline(preMafftFin, line)) {
		if (line[0] == '>') {
			preMafftFin.clear();
			preMafftFin.seekg(0);
			break;
		}
		seqLength++;
	}
	unordered_map<string, string> jSequences;
	while (getline(preMafftFin, line)) {
		id = line.substr(1, line.size() - 1);
		while (id.back() == '\t' || id.back() == ' ') id.pop_back();
		getline(preMafftFin, line); // gets sequence
		string extraLines;
		for (int i = 0; i < seqLength; i++) {
			getline(preMafftFin, extraLines);
			line += extraLines;
		}
		if (changes.find(id) != changes.end()) {
			// id was found in the hash table. Go through vector of insertions and add to the postMafftInsertions
			vector<mapInserts> vec = changes[id];
			for (int i = 0; i < vec.size(); i++) {
				mapInserts p = vec[i];
				string ins = line.substr(p.preMafftResStart - 1, p.preMafftResEnd - p.preMafftResStart + 1); // residues to add to postMafft
				if (ins.size() < masterset[p.postMafftIns] && !jSeq) ins.resize(masterset[p.postMafftIns], '-'); // pads with '-' if length is less
				seqParts s;
				s.location = p.postMafftIns;
				s.seq = ins;
				postInserts[id].insertions.push_back(s);
				if (jSeq && p.postMafftIns == 0) {
					jSequences[id] = (line.substr(p.preMafftResEnd, line.size() - p.preMafftResEnd));
				}
			}
			for (unordered_map<int, int>::const_iterator it = masterset.begin(); it != masterset.end(); it++) {
				if (find(postInserts[id].insertions.begin(), postInserts[id].insertions.end(), it->first) == postInserts[id].insertions.end()) {
					// means that this insertion point doesn't exist. This will be only gaps.
					string fill(it->second, '-');
					seqParts s;
					s.location = it->first;
					s.seq = fill;
					postInserts[id].insertions.push_back(s);
				}
			} // all insertions should be pre-imposed
			sort(postInserts[id].insertions.begin(), postInserts[id].insertions.end(), lessSeqParts); // sorts the vector by location
		}
	}
	// Now that insertions have been parsed, need to add to output
	vector<string> outputs;
	int longestSeq = 0;
	vector<seqParts> mastervec;
	for (unordered_map<int, int>::const_iterator it = masterset.begin(); it != masterset.end(); it++) {
		string fill(it->second, '-');
		seqParts s;
		s.location = it->first;
		s.seq = fill;
		mastervec.push_back(s);
	}
	sort(mastervec.begin(), mastervec.end(), lessSeqParts); // sorted vector to add fills for non-found sequences
	string oldIdCopy;
	seqLength = -1;
	getline(postMafftFin, line);
	while (getline(postMafftFin, line)) {
		if (line[0] == '>') {
			postMafftFin.clear();
			postMafftFin.seekg(0);
			break;
		}
		seqLength++;
	}
	while(getline(postMafftFin, line)) {
		istringstream sin(line);
		string id;
		getline(sin, oldIdCopy, '|');
		getline(sin, id); // actually gets real id
		if (id.substr(0, 3) == "IGH" || id.substr(0, 3) == "TCR") id = oldIdCopy.substr(1, oldIdCopy.size() - 1);
		if (line.size() <= 12) id = line;
		while (id.back() == '\t' || id.back() == ' ') id.pop_back();
		if (id[0] == '>') id = id.substr(1, id.size() - 1);
		string secondHalf = "";
		getline(postMafftFin, line);
		for (int i = 0; i < seqLength; i++) {
			getline(postMafftFin, secondHalf);
			line += secondHalf;
		} // gets the full sequence if needed
		longestSeq = max(longestSeq, (int)line.size());
		if (postInserts.count(id) > 0) {
			if (jSeq) {
				for (int i = 0; i < postInserts[id].insertions.size(); i++) {
					seqParts s = postInserts[id].insertions[i];
					if (s.location != 0) continue;
					line.insert(s.location, s.seq);
				}  
			}
			else {
				for (int i = (int)postInserts[id].insertions.size() - 1; i >= 0; i--) {
					seqParts s = postInserts[id].insertions[i];
					line.insert(s.location, s.seq);
				} // adds backwards to avoid adding offsets
			}
		}
		else { // not found in postMafft, just print seq with dashes
			for (int i = (int)mastervec.size() - 1; i >= 0; i--) {
				line.insert(mastervec[i].location, mastervec[i].seq);
			}
		}
		string outSeq;
		if (!jSeq) {
			outSeq = id + ',' + line[0];
			for (int i = 0; i < line.size(); i++) {
				outSeq += ',';
				outSeq += line[i];
			}
		}
		else {
			outSeq = '>' + id + '\n' + line;
			jfout << '>' << id << endl << jSequences[id] << endl; // add cleaved j sequences to another file.
		}
		outputs.push_back(outSeq);
	}
	if (!jSeq) {
		fout << "Id";
		for (int i = 0; i <= longestSeq; i++) {
			fout << "," << i;
			if (masterset.count(i) > 0) {
				string outChar = "";
				for (int j = 0; j < masterset[i]; j++) {
					if (j != 0 && j % 26 == 0) {
						int ind = outChar.size() - 1;
						if (outChar.empty()) {
							outChar = "A";
						} else {
							while (ind >= 0) {
								if (outChar[ind] == 'Z') {
									outChar[ind] = 'A'; // Carry over
									ind--; // Move left
								} else {
									outChar[ind]++; // Simple increment
									break;
								}
							}
							if (ind < 0) {
								outChar = "A" + outChar;
							}
						}
					}
					fout << "," << i << outChar << (char) ((j % 26) + 'A');
				}
			}
		} // print out colnames
		fout << endl;
	}
	for (int i = 0; i < outputs.size(); i++) fout << outputs[i] << endl;
	preMafftFin.close();
	postMafftFin.close();
	mapFin.close();
	fout.close();
	jfout.close();
	return 0;
}
