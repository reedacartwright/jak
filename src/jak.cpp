// jak.cpp

/* TODO: Create input file format
reps= 10						
branches= ((A:0.1,B:0.1):0.7,C:0.5);			//quoted labels and keys or relabel w/ perl
thetas= ((A:0.5,B:0.5):0.4,C:0.6):0.8;
A = 10
B = 5
C = 5
seed = 1776
*/

#ifdef _MSC_VER
#	include <process.h>
#	define getpid _getpid
#else
#	include <unistd.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <cfloat>
#include <fstream>

#include "xorshift64.h"
#include "rexp.h"
#include "aliastable.h"

using namespace std;

// Conditional probably of mutation
double mutation_matrix[4][4] = {
	{0.25, 0.25, 0.25, 0.25},
	{0.25, 0.25, 0.25, 0.25},
	{0.25, 0.25, 0.25, 0.25},
	{0.25, 0.25, 0.25, 0.25}
};
alias_table mutation[4];

// data structure to hold node information
struct nodestruct {
    int child_1;
    int child_2;
    int parent;
    char label; // nucleotide
    double time;
    int type; // species/population identifier
};

void coaltree(vector<int>& activelist, double theta, double time, char type,
	          vector<nodestruct>& nodeVector, xorshift64& myrand1);

// Function to convert int type to string
// Cache labels for reuse
string id_to_string(int x) {   
	static vector<string> v;
	static char buffer[64];
	for(int xx = (int)v.size(); xx <= x; ++xx) {
		sprintf(buffer, "%d", xx);
		v.push_back(buffer);
	}
    return v[x];
}

// TODO: Cache this
string speciesLabel(int type)			//Function to convert species number to letter format for tree output
{
    int x=1;
    string ans="";
	int value=1;
	vector<int> index;
	if (type > 26) {
		do {
			value = type%26;
			type = type/26;
			index.push_back(value);
		} while(type!=0);
		reverse(index.begin(),index.end());                           //reverses contents of vector
    }
	else 
		index.push_back(type);
	   
    for(int i=0; i<index.size(); i++)
    {
		if ((index[i] - 1) != -1) 
			ans += 'A' + (index[i] - 1);
		else 
			ans += 'A' - 20;
    }
    return ans;
}
//-----------------------------------------------------------------------------//
//create newick tree from node data
string tree_to_string(const vector<nodestruct>& v) {
    vector<string> node_str(v.size(),"");
    char buffer[16];
    for(int i=0; i<v.size(); i++) {
        string temp = "";

        if(v[i].child_1 != -1 && v[i].child_2 != -1) {
            string convert_node1=id_to_string(v[i].child_1);
            string convert_node2=id_to_string(v[i].child_2);
			temp += "(" + node_str[v[i].child_1] + speciesLabel(v[v[i].child_1].type) + convert_node1 + ":";
            sprintf(buffer, "%0.6f", v[v[i].child_1].time);
            temp += string(buffer) + "," + node_str[v[i].child_2] + speciesLabel(v[v[i].child_2].type) + convert_node2+ ":";
            sprintf(buffer, "%0.6f", v[v[i].child_2].time);
            temp += string(buffer) + ")";
        }
        node_str[i] = temp;
    }
	string temp = speciesLabel(v.back().type) + id_to_string(v.size() - 1);
    return node_str.back() + temp + ";";
}
//-----------------------------------------------------------------------------//
string mutationLabels(const vector<nodestruct>& t)			//Function to construct vector of mutation labels for tree, in numerical order
{
	string temp = "";
	for(int i=0; i<t.size(); i++) {
		temp += t[i].label;
	}
	return "[" + temp + "]";
}
//-----------------------------------------------------------------------------//
// random seed generator
inline unsigned int create_random_seed() {
    unsigned int v = static_cast<unsigned int>(getpid());
    v += ((v << 15) + (v >> 3)) + 0x6ba658b3; // Spread 5-decimal PID over 32-bit number
    v^=(v<<17);
    v^=(v>>13);
    v^=(v<<5);
    v += static_cast<unsigned int>(time(NULL));
    v^=(v<<17);
    v^=(v>>13);
    v^=(v<<5);
    return (v == 0) ? 0x6a27d958 : (v & 0x7FFFFFFF); // return at most a 31-bit seed
}

//-----------------------------------------------------------------------------//
void set_mutations(xorshift64 &myrand1, char &G, double time, int& counter)
{   
    double m = rand_exp(myrand1); //m = total distance travelled along branch length
    while (m <= time) { //if m < branch length --> mutate
        ++counter;  //muation counter
		// use the alias tables to effeciently sample the result of the mutation
		G = static_cast<char>(mutation[G](myrand1.get_uint64()));
        m += rand_exp(myrand1);
    }
}

//-------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    int N1, N2, n, N, trees;
    double theta1, theta2, theta3, t1, t2;

	//TODO: fix input validation
    if (argc == 9) {
        trees=atoi(argv[1]);
			while(trees<1)
			{	cout<<"Error: Enter number of trees: ";
				cin>>trees;
			}
        N1=atoi(argv[2]);
			while(N1<1)
			{	cout<<"Error: Enter number of tips for first species: ";
				cin>>N1;
			}
        N2=atoi(argv[3]);
			while(N2<1)
			{	cout<<"Error: Enter number of tips for second species: ";
				cin>>N1;
			}
        theta1=strtod(argv[4], NULL);
			while(theta1<0.0)
			{	cout<<"Error: Enter theta1: ";
				cin>>theta1;
			}
        theta2=strtod(argv[5], NULL);
			while(theta2<0.0)
			{	cout<<"Error: Enter theta2: ";
				cin>>theta2;
			}
        theta3=strtod(argv[6], NULL);
			while(theta3<0.0)
				{	cout<<"Error: Enter theta3: ";
					cin>>theta3;
				}
        t1=strtod(argv[7], NULL);
			while(t1<0.0)
			{	cout<<"Error: Enter t1: ";
				cin>>t1;
			}
        t2=strtod(argv[8], NULL);
			while(t2<0.0)
			{	cout<<"Error: Enter t2: ";
				cin>>t2;
			}
    }

    else if (argc == 1) {
        cout << "Enter number of trees: ";
        cin >> trees;
        cout << endl << "Enter the number of tips for first species: ";
        cin >> N1;                                                                             //receive # of starting nodes (species 1)
        cout << endl << "Enter number of tips for second species: ";
        cin >> N2;                                                                             //receive # of starting nodes (species 2)
        cout << endl << "Enter theta 1: ";
        cin >> theta1;
        cout << endl << "Enter theta 2: ";
        cin >> theta2;
        cout << endl << "Enter theta 3: ";
        cin >> theta3;
        cout << endl << "Enter t1: ";
        cin >> t1;
        cout << endl << "Enter t2: ";
        cin >> t2;
    }

    else {
        cerr << "error, must have format: prg name, trees, # of tips for first"
                "species, # of tips for second species, theta1, theta2,"
                "theta3, t1, t2" << endl;
		cerr << "Press ENTER to quit." << flush;
		cin.clear();
		cin.sync();
		cin.ignore( numeric_limits<streamsize>::max(), '\n' );
        return EXIT_FAILURE;
    }

    n = N1+N2; //n=total original tips from both species

	//use xorshift64 class for random number generator
    xorshift64 myrand;
    myrand.seed(create_random_seed());

	// construct alias tables for mutation simulation
	for(int i=0;i<4;++i) {
		mutation[i].create(&mutation_matrix[i][0],&mutation_matrix[i][4]);
	}
	
    for(int repeat=0; repeat<trees; repeat++) { //loops once for each tree
        vector<nodestruct> nodevector(n); //create nodevector (vector of structs)
		
		vector<int> active1(N1); //initialize active list for species 1
		for(int i=0; i<N1; i++)
		{
			active1[i]=i;
			nodevector[i].child_1=-1;
			nodevector[i].child_2=-1;
			nodevector[i].parent=-1;
			nodevector[i].label='N';
			nodevector[i].time=0;
			nodevector[i].type=1;
		}


		vector<int> active2(N2); //initialize active list for species 2
		for(int j=N1; j<N1+N2; j++)
		{
			active2[j-N1]=j;
			nodevector[j].child_1=-1;
			nodevector[j].child_2=-1;
			nodevector[j].parent=-1;
			nodevector[j].label='N';
			nodevector[j].time=0;
			nodevector[j].type=2;
		}

//----------------------------------------------------------------------------//
		coaltree(active1, theta1, t1, 1, nodevector, myrand);
		coaltree(active2, theta2, t2, 2, nodevector, myrand);

        active1.insert(active1.end(), active2.begin(), active2.end());
        coaltree(active1, theta3, DBL_MAX, 3, nodevector, myrand);

//----------------------------------------------------------------------------//
//mutations
        nodevector.back().label = 0;
		// start mutations for loop
        for (int i = (int)nodevector.size() - 2; i >= 0; --i) {                                               
            int counter = 0;
            nodevector[i].label = nodevector[nodevector[i].parent].label;
            set_mutations(myrand, nodevector[i].label, nodevector[i].time, counter);
        }
		// Label nodes
        for (int i=0; i < nodevector.size(); i++){
            char s[] = "ATCG";
            nodevector[i].label = s[nodevector[i].label];
        }

//----------------------------------------------------------------------------//

		cout << mutationLabels(nodevector) << "\t";
        cout << tree_to_string(nodevector) << endl; //print newick tree to console
    } //end # of trees loop

#ifndef NDEBUG
    cerr<<"Random seed used: "<<create_random_seed()<<endl;
    cerr << "Press ENTER to quit." << flush;
	cin.clear();
	cin.sync();
    cin.ignore( numeric_limits<streamsize>::max(), '\n' );
#endif
    return EXIT_SUCCESS;
}

void coaltree(vector<int>& activelist, double theta, double time, char type,
	          vector<nodestruct>& nodeVector, xorshift64& myrand1)
{
    double T = 0.0;
    int i = 0;
	int random1, random2;

	size_t size = activelist.size();

	while(size>1)
	{
        // Draw waiting time until next coalescent event
		double mean = (2.0/(size*(size-1.0)))*(theta/2.0);
		double U = rand_exp(myrand1, mean);
		if(T+U>time)
			break;
		T+=U;

		// pick a random pair of nodes
		// TODO: use alias table to optimize this
		random1 = (myrand1.get_uint32() % size);
		do {
			random2 = (myrand1.get_uint32() % size);
		} while(random1==random2);

		//orders two nodes minimum to maximum
		if (random1>random2) 
			swap(random1,random2);

		int newparent = (int)nodeVector.size();
		nodeVector.push_back(nodestruct());
		
		nodeVector[newparent].type = type;

		//update parent node
		nodeVector[newparent].child_1 = activelist[random1];
		nodeVector[newparent].child_2 = activelist[random2];
		nodeVector[newparent].time = T;
		
		//update child nodes
		nodeVector[activelist[random1]].parent = newparent;
		nodeVector[activelist[random2]].parent = newparent;
		nodeVector[activelist[random1]].time = T - nodeVector[activelist[random1]].time;
		nodeVector[activelist[random2]].time = T - nodeVector[activelist[random2]].time;

		//update active vector
		activelist[random1] = newparent;
		activelist.erase (activelist.begin() + random2);

		size--;
	}

	for(int i=0; i<size && time!=DBL_MAX; i++)
		nodeVector[activelist[i]].time = nodeVector[activelist[i]].time - time;
}
