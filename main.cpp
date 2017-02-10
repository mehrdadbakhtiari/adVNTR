#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

map<char, int> get_nucleotide_map(const string& s)
{
	map<char, int> m;
	m['A'] = m['C'] = m['G'] = m['T'] = 0;
	for (int i = 0; i < s.size(); i++)
	{
		m[s[i]] ++;
	}
	return m;
}

int nucleotide_dist(const map<char, int> &m1, const map<char, int> &m2)
{
	string nucleotides = "ACGT";
	int res = 0;
	for (int i = 0; i < nucleotides.size(); i++)
	{
		res += abs(m1.at(nucleotides[i]) - m2.at(nucleotides[i]));
	}
	return res;
}

bool match_read_by_sliding_window(const string &query, const map<char, int> &query_map, const string &s)
{
	map<char, int> current_map;
	current_map['A'] = current_map['C'] = current_map['G'] = current_map['T'] = 0;
	for (int i = 0; i < s.size(); i++)
	{
		current_map[s[i]]++;
		if (i >= query.size())
		{
			current_map[s[i - query.size()]]--;
		}
		if (nucleotide_dist(current_map, query_map) < 3)
			return true;
	}
	return false;
}

int main()
{
	std::ifstream file("/Users/mehrdad/workspace/VeNTeR/paired_dat1.fasta");
	std::string str;
	vector<string> result;
	string query = "TAGAACAGAAGGACAAGGCCCTGGAACCAAAAGATAAAGACT";
	map<char, int> query_map = get_nucleotide_map(query);
	while (getline(file, str))
	{
		getline(file, str);
		if (match_read_by_sliding_window(query, query_map, str))
		{
			result.push_back(str);
		}
	}
	cout << result.size() << endl;
/*	for (int i = 0; i < result.size(); i++)
	{
		cout << result[i] << endl;
	}
*/
}

