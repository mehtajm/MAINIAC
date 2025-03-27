// Jai Mehta
// 03/26/2025
// Code for concatenating V and J numbering schemes
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

using namespace std;

string parse (int &number, string input) { // to separate numbers and chars
    string val = "";
    string insertChars = "";
    for (size_t i = 0; i < input.size(); i++) {
        if (!isalpha(input[i])) val += input[i];
        else insertChars += input[i];
    }
    number = stoi(val);
    return insertChars;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        cerr << "Usage: ./concatenate [v.csv] [j.csv] [output.csv]\n";
        return 1;
    }
    ifstream vFin(argv[1]);
    if (vFin.fail()) {
        cerr << "V numbering csv did not open properly. Try again.\n";
		return 1;
    }
    ifstream jFin(argv[2]);
    if (jFin.fail()) {
        cerr << "J numbering csv did not open properly. Try again.\n";
		return 1;
    }
    ofstream fout(argv[3]);
    if (fout.fail()) {
        cerr << "Output file failed to open. Try again.\n";
        return 1;
    }
    string colNames, lastColChars = "", lastColNumber = "";
    getline(vFin, colNames); // gets colnames 
    for (int i = colNames.size() - 1; i >= 0; i--) {
        if (colNames[i] == ',') break;
        if (isalpha(colNames[i])) {
            //if (lastColChars.empty()) lastColChars.push_back(colNames[i]);
            lastColChars.insert(lastColChars.begin() + 0, colNames[i]); }
        else lastColNumber.insert(lastColNumber.begin() + 0, colNames[i]);
    }
    int colNum = stoi(lastColNumber);
    string jColNames, colId;
    getline(jFin, jColNames);
    istringstream sin(jColNames);
    getline(sin, colId, ',');  // remove Id from 
    while (getline(sin, colId, ',')) {
        int val = -1;
        string insertChars = parse(val, colId);
        if (val == 0) {
            insertChars = lastColChars;
            int ind = insertChars.size() - 1;
            while (ind >= 0) {
                if (insertChars[ind] == 'Z') {
                    insertChars[ind] = 'A'; // Carry over
                    ind--; // Move left
                    } 
                    else {
                        insertChars[ind]++; // Simple increment
                        break;
                    }
                    if (ind < 0) {
						insertChars = "A" + insertChars;
					}
            }
            lastColChars = insertChars;
        }
        val += colNum;
        colNames += ',' + to_string(val) + insertChars;
    }
    fout << colNames << endl;
    int seqSize = 0, largestSize = 0;
    istringstream s0(colNames);
    while (getline(s0, colId, ',')) seqSize++;
    string seq;
    unordered_map<string, string> outputs;
    while (getline(vFin, seq)) {
        istringstream s1(seq);
        string id;
        getline(s1, id, ','); // grab id
        getline(s1, seq);
        if (id[0] == 'q') outputs[id] += seq + ',';
    }
    while (getline(jFin, seq)) {
        istringstream s1(seq);
        string id;
        getline(s1, id, ',');
        getline(s1, seq);
        if (id[0] == 'q') outputs[id] += seq;
    }
    for (unordered_map<string, string>::const_iterator it = outputs.begin(); it != outputs.end(); it++) {
       if (it->first[0] != 'q') continue; 
       fout << it->first << ',' << it->second << endl;
    }
    vFin.close();
    jFin.close();
    fout.close();
    return 0;
}
