#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <bits/stdc++.h>
#include <map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>


using namespace std;

int min_number_of_keyword_matches = 5;
int max_reads_to_report_for_vntr = 10000;

map<int, int> keyword_to_vntr;
vector<int> vntr_ids;

const int MAXS = 103582 * 21; // Sum of lengths of keywords
const int MAXC = 5; // ACGT and N

// OUTPUT FUNCTION IS IMPLEMENTED USING out[] 
// Bit i in this mask is one if the word with index i 
// appears when the machine enters this state. 
map<int, set<int> > out_map;

// FAILURE FUNCTION IS IMPLEMENTED USING f[] 
int f[MAXS]; 

// GOTO FUNCTION (OR TRIE) IS IMPLEMENTED USING g[][] 
int g[MAXS][MAXC]; 

struct count_t
{
	int value;
	count_t(): value(0){} // default value 0
};

int char_to_num(char ch)
{
    if (ch == 'A')
        return 0;
    if (ch == 'C')
        return 1;
    if (ch == 'G')
        return 2;
    if (ch == 'T')
        return 3;
    return 4;
}

int buildMatchingMachine(vector<string> &arr, int k)
{
    // Initialize all values in goto function as -1.
    memset(g, -1, sizeof g);

    // Initially, we just have the 0 state
    int states = 1;

    // Construct values for goto function, i.e., fill g[][]
    // This is same as building a Trie for arr[]
    for (int i = 0; i < k; ++i)
    {
        const string &word = arr[i];
        int currentState = 0;

        // Insert all characters of current word in arr[]
        for (int j = 0; j < word.size(); ++j)
        {
            int ch = char_to_num(word[j]);

            // Allocate a new node (create a new state) if a
            // node for ch doesn't exist.
            if (g[currentState][ch] == -1)
                g[currentState][ch] = states++;

            currentState = g[currentState][ch];
        }

        // Add current word in output function
        out_map[currentState].insert(i);
    }

    // For all characters which don't have an edge from
    // root (or state 0) in Trie, add a goto edge to state
    // 0 itself
    for (int ch = 0; ch < MAXC; ++ch)
        if (g[0][ch] == -1)
            g[0][ch] = 0;

    // Now, let's build the failure function
    // Initialize values in fail function
    memset(f, -1, sizeof f);

    // Failure function is computed in breadth first order
    // using a queue
    queue<int> q;

     // Iterate over every possible input
    for (int ch = 0; ch < MAXC; ++ch)
    {
        // All nodes of depth 1 have failure function value
        // as 0. For example, in above diagram we move to 0
        // from states 1 and 3.
        if (g[0][ch] != 0)
        {
            f[g[0][ch]] = 0; 
            q.push(g[0][ch]);
        }
    }

    // Now queue has states 1 and 3
    while (q.size())
    {
        // Remove the front state from queue
        int state = q.front() ;
        q.pop();

        // For the removed state, find failure function for
        // all those characters for which goto function is
        // not defined.
        for (int ch = 0; ch < MAXC; ++ch)
        {
            // If goto function is defined for character 'ch'
            // and 'state'
            if (g[state][ch] != -1)
            {
                // Find failure state of removed state
                int failure = f[state];

                // Find the deepest node labeled by proper
                // suffix of string from root to current
                // state.
                while (g[failure][ch] == -1)
                {
                    failure = f[failure];
                }

                failure = g[failure][ch];
                f[g[state][ch]] = failure;

                // Merge output values
                for (std::set<int>::iterator it = out_map[failure].begin(); it != out_map[failure].end(); it++)
                	out_map[g[state][ch]].insert(*it);

                // Insert the next level node (of Trie) in Queue
                q.push(g[state][ch]);
            }
        }
    }

    return states;
}

// Returns the next state the machine will transition to using goto
// and failure functions.
// currentState - The current state of the machine. Must be between
//                0 and the number of states - 1, inclusive.
// nextInput - The next character that enters into the machine.
int findNextState(int currentState, char nextInput)
{
    int answer = currentState;
    int ch = char_to_num(nextInput);

    // If goto is not defined, use failure function
    while (g[answer][ch] == -1)
        answer = f[answer];

    return g[answer][ch];
}

map<int, int> get_keywords(vector<string> &arr)
{
    ifstream in("/dev/stdin", ios::in);

    map<int, int> m;
    string line;
    int vntr_id = 0;
    int counter = 0;
    while(true)
    {
        if(!getline(in, line,'\n')) break;
        istringstream iss(line);
        vector<string> tokens;
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(tokens));
        vntr_id = atoi(tokens[0].c_str());
        vntr_ids.push_back(vntr_id);
        for (int i = 1; i < tokens.size(); i++)
        {
            m[arr.size()] = vntr_id;
            arr.push_back(tokens[i]);
        }
//        if (++counter > 20)
//        	break;
    }
/*    vector<string>::iterator it = std::unique( arr.begin(), arr.end() );
    if (it == arr.end())
    	cerr << "keywords are unique" << endl;
    else
    	cerr << "KEYWORDS ARE NOT UNIQUE!" << endl;
    cerr << "returning the keywords map with " << arr.size() << " keywords" << endl;
*/
    return m;
}

int usage(int status=0)
{
	cerr << "Program: adVNTR-Filtering" << endl;
	cerr << "usage: adVNTR-Filtering sequences.fa < keywords.txt > output.txt" << endl;
	return status;
}

int main(int argc,char **argv)
{
	if (argc < 2)
		return usage(1);
	if (strcmp(argv[1], "--help") == 0)
		return usage(0);
    vector<string> arr;
    keyword_to_vntr = get_keywords(arr);
    int k = arr.size();
    buildMatchingMachine(arr, k);
//    cerr << "Matching Machine has been built" << endl;

    map<int, map<string, count_t> > vntr_read_list;
    map<string, string> read_sequences;
    string line;
    string name, seq;
    ifstream in2(argv[1], ios::in);
//    int max_reads = 83083774 / 1000000;
//    max_reads = 83083774;
    int count = 0;
    while(true)
    {
        if(!getline(in2, name,'\n')) break;
        if(!getline(in2, seq,'\n')) break;
        name = name.substr(1);

        int current_state = 0;
        map<int, count_t> vntr_matches_count;
        for (int i = 0; i < seq.length(); i++)
        {
            current_state = findNextState(current_state, seq[i]);
            if (out_map[current_state].size() == 0)
            	continue;

            // Match found, iterating all matching words of arr[]
            for (std::set<int>::iterator it = out_map[current_state].begin(); it != out_map[current_state].end(); it++)
            {
            	int word_index = *it;
            	int vntr_id = keyword_to_vntr[word_index];
                int value = vntr_matches_count[vntr_id].value;
                vntr_matches_count[vntr_id].value = value + 1;
            }
        }
        // PUT RESULTS IN PERMANENT DATA STRUCTURE
        typedef map<int, count_t>::iterator it_type;
        for (it_type it = vntr_matches_count.begin(); it != vntr_matches_count.end(); it++)
        {
            int vntr_id = it->first;
            int occurrence = it->second.value;
            if (occurrence >= min_number_of_keyword_matches)
            {
                vntr_read_list[vntr_id][name].value = occurrence;
                read_sequences[name] = seq;
            }
        }

//        if (++count > max_reads)
//        	break;
    }
    in2.close();

//    cerr << "Done with Aho-Corasick." << endl;

    set<string> filtered_reads;
    map<int, vector<pair<int, string> > > vntr_filtered_reads;
    for (int i = 0; i < vntr_ids.size(); i++)
    {
    	int vntr_id = vntr_ids[i];
        for (map<string, count_t>::iterator it = vntr_read_list[vntr_id].begin(); it != vntr_read_list[vntr_id].end(); it++)
        {
            string read_name = it->first;
            int occurrence = it->second.value;
            vntr_filtered_reads[vntr_id].push_back(make_pair(occurrence, read_name));
    	}

        int result_size = vntr_filtered_reads[vntr_id].size();
        if (result_size > max_reads_to_report_for_vntr)
            result_size = max_reads_to_report_for_vntr;
		cout << vntr_id << " " << result_size;
    	if (vntr_filtered_reads[vntr_id].size() > 0)
    	{
			std::sort(vntr_filtered_reads[vntr_id].rbegin(), vntr_filtered_reads[vntr_id].rend());
			for (int j = 0; j < vntr_filtered_reads[vntr_id].size(); j++)
			{
				filtered_reads.insert(vntr_filtered_reads[vntr_id][j].second);
				cout << " " << vntr_filtered_reads[vntr_id][j].second;
				if (j >= max_reads_to_report_for_vntr)
					break;
			}
    	}
		cout << endl;
    }

//    cerr << "Printing read sequences" << endl;
    for (std::set<string>::iterator it = filtered_reads.begin(); it != filtered_reads.end(); it++)
    {
    	string name = *it;
    	cout << name << " " << read_sequences[name] << endl;
    }
    return 0;

}
