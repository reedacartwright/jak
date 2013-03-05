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
#include <iomanip>

#include <boost/program_options.hpp>

#include "xorshift64.h"
#include "rexp.h"
#include "aliastable.h"

// Utility functions
#define CERRORR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl), false);
#define CERROR(err_msg) (std::cerr << "ERROR: " << err_msg << std::endl);

using namespace std;

// Conditional probably of mutation
double mutation_matrix[4][4] = {
/*A*/ {0.25, 0.25, 0.25, 0.25},
/*C*/ {0.25, 0.25, 0.25, 0.25},
/*G*/ {0.25, 0.25, 0.25, 0.25},
/*T*/ {0.25, 0.25, 0.25, 0.25}
};
alias_table mutation[4];

// data structure to hold node information
struct nodestruct {
    int child_1;   // ID of right child
    int child_2;   // ID of left child
    int parent;    // ID of parent
    double time;   // "distance" to parent
    int species;   // species/population identifier
    char label;    // nucleotide
	
	nodestruct() : child_1(-1), child_2(-1), parent(-1),
	               time(0.0), species(0), label(0) { }
};

void coaltree(xorshift64& myrand1, vector<int>& activelist, double theta, double time, char species,
	          vector<nodestruct>& nodeVector);
int set_mutations(xorshift64 &myrand1, char &G, double time);

string id_to_string(int x);
string species_label(int type);
string tree_to_string(const vector<nodestruct>& v);
string mutation_string(const vector<nodestruct>& t);
uint64_t key_create(const vector<nodestruct>& temp_nodes);

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

//----------------------------------------------------------------------------//
namespace po = boost::program_options;

namespace boost { namespace program_options {
template<>
typed_value<bool>* value(bool* v) {
	return bool_switch(v);
	//typed_value<bool>* r = new typed_value<bool>(v);
    //r->default_value(0, "off");
    //r->implicit_value(1, "on");
	//return r;
}
}}

struct args_t
{			
	// use X-Macros to specify argument variables
#	define XM(lname, sname, desc, type, def) type XV(lname) ;
#	include "jakargs.xmh"
#	undef XM
	po::options_description desc;
	string runname;
};

bool process_args(int argc, char* argv[], args_t &a){
	po::variables_map vm;
	a.runname = argv[0];
	a.desc.add_options()
		#define XM(lname, sname, desc, type, def) ( \
			XS(lname) IFD(sname, "," BOOST_PP_STRINGIZE sname), \
			po::value< type >(&a.XV(lname))->default_value(def), \
			desc )				
		#include "jakargs.xmh"
		#undef XM
	;
	try {
		po::store(po::command_line_parser(argc, argv).options(a.desc).run(), vm);
		po::notify(vm);
	} catch (std::exception &e) {
		CERROR(e.what());
		return CERRORR("unable to process command line");
	}
	return true;
}

//----------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
	int N1, N2, n, trees;
	double theta1, theta2, theta3, t1, t2;
	typedef map<uint64_t,size_t> counts_t;
	counts_t counts;
	args_t args;
	if(!process_args(argc, argv, args))
		return EXIT_FAILURE;
	if(args.help) {
		cerr << "Usage:\n  "
		     << args.runname << " [options]"
			 << endl << endl;
		cerr << args.desc << endl;
		return EXIT_FAILURE;
	}
	N1 = args.tips1;
	N2 = args.tips2;
	trees = args.reps;
	theta1 = args.theta1;
	theta2 = args.theta2;
	theta3 = args.theta3;
	t1 = args.t1;
	t2 = args.t2;


	n = N1+N2; //n=total original tips from both species

	//use xorshift64 class for random number generator
    xorshift64 myrand;
    unsigned int seed = args.seed;
    if(seed == 0)
    	seed = create_random_seed();
    myrand.seed(seed);
	
	cerr << "# Using PRNG seed: " << seed << endl;
	
	// construct alias tables for mutation simulation
	for(int i=0;i<4;++i) {
		mutation[i].create(&mutation_matrix[i][0],&mutation_matrix[i][4]);
	}
	
	// loop once for each tree
    for(int repeat=0; repeat<trees; repeat++)
	{ 
        // create nodevector (vector of structs)
		vector<nodestruct> nodevector(n);
		
		// initialize active list for species 1
		vector<int> active1(N1); 
		for(int i=0; i<N1; i++)
		{
			active1[i]=i;
			nodevector[i].species=1;
		}

		// initialize active list for species 2
		vector<int> active2(N2);
		for(int j=N1; j<N1+N2; j++)
		{
			active2[j-N1]=j;
			nodevector[j].species=2;
		}
		
		// coalesce population 1
		coaltree(myrand, active1, theta1, t1, 1, nodevector);
		
		// coalesce population 2
		coaltree(myrand, active2, theta2, t2, 2, nodevector);

		// coalesce population 3
        active1.insert(active1.end(), active2.begin(), active2.end());
        coaltree(myrand, active1, theta3, DBL_MAX, 3, nodevector);

		// Assume root has nucleotide 'A'
        nodevector.back().label = 0;
		// start mutations for loop
        for (int i = (int)nodevector.size() - 2; i >= 0; --i)
		{                                               
            nodevector[i].label = nodevector[nodevector[i].parent].label;
            set_mutations(myrand, nodevector[i].label, nodevector[i].time);
        }
        if(args.count_tips)
        {
        	uint64_t u = key_create(nodevector);
        	// cout << setw(16) << setfill('0') << hex << u << endl;
        	++counts[u];
        	
        } else {
			// Label nodes
		    for (size_t i=0; i < nodevector.size(); i++){
		        char s[] = "ACGT";
		        nodevector[i].label = s[(size_t)nodevector[i].label];
		    }
		    // print newick tree to console
			cout << mutation_string(nodevector) << "\t";
		    cout << tree_to_string(nodevector) << endl;      
        }
        
    } //end # of trees loop
    if(args.count_tips) {
    	cout << "A1\tC1\tG1\tT1\tA2\tC2\tG2\tT2\tCount\n";
    	for(counts_t::const_iterator it = counts.begin();
    		it != counts.end(); ++it) {
    		uint64_t u = it->first;
    		for(int i=0;i<8;u >>= 8,++i)
    			cout << (u & 0xFF) << "\t";
    		cout << it->second << endl;
    	}
    }

#ifndef NDEBUG
//    cerr << "Press ENTER to quit." << flush;
//	cin.clear();
//	cin.sync();
//    cin.ignore( numeric_limits<streamsize>::max(), '\n' );
#endif
    return EXIT_SUCCESS;
}

// key_create must be used before relabeling
uint64_t key_create(const vector<nodestruct>& temp_nodes){
    union key_union{
        uint64_t key;
        unsigned char count[8];

    } k;

    k.key = 0;

    for (size_t i = 0; i < temp_nodes.size(); ++i)
    {
        assert(temp_nodes[i].label < 4);
        //ensure node is a tip
        if(temp_nodes[i].child_1 != -1 || temp_nodes[i].child_2 != -1)
        	break;
       	switch(temp_nodes[i].species) {
       	case 1:
       		k.count[(size_t)temp_nodes[i].label]++;
       		break;
       	case 2:
       		k.count[4+(size_t)temp_nodes[i].label]++;
       		break;
       	default:
       		// You should never get here.
       		assert(false);
       		break;
        }
    }
	return k.key;
} 

void coaltree(xorshift64& myrand1, vector<int>& activelist, double theta, double time, char species,
	          vector<nodestruct>& nodeVector)
{
    double T = 0.0;
	int random1, random2;

	int size = (int)activelist.size();

	while(size>1)
	{
        // Draw waiting time until next coalescent event
		double mean = (2.0/(size*(size-1.0)))*(theta/2.0);
		T += rand_exp(myrand1, mean);
		if(T > time)
			break;

		// pick a random pair of nodes
		random1 = static_cast<int>(size*myrand1.get_double52());
		random2 = static_cast<int>((size-1)*myrand1.get_double52());
		random2 = (random1+random2+1) % size;
		if(random1 > random2)
			swap(random1,random2);
		
		int newparent = (int)nodeVector.size();
		nodeVector.push_back(nodestruct());
		
		nodeVector[newparent].species = species;

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

	for(int i=0; i<size; i++)
		nodeVector[activelist[i]].time = nodeVector[activelist[i]].time - time;
}

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
//Function to convert species number to letter format for tree output
string species_label(int species)
{
	// Species numbers start at 1.
	assert(species >= 1);
	static vector<string> v;
	species -= 1;
	if(species < 0)
		return "%";
	string ans;
	int zz = 0;
	int aa = 0;
	
	for(int xx=(int)v.size(); xx <= species; ++xx) {
		ans.clear();
		zz = xx;
		
		if (zz <= 25)
			ans += 'A' + zz;
			
		else {
			while(zz >= 26){
			zz = zz - ans.size();
			ans += 'A' + zz%26;
			zz /= 26;
			aa = 1;
			if (zz == 26)
				break;
			} 
			if (aa == 1)
				zz = zz - 1;
			ans += 'A' + zz%26;
		}
		
		reverse(ans.begin(),ans.end());
		v.push_back(ans);
	}
    return v[species];
}
//-----------------------------------------------------------------------------//
//create newick tree from node data
string tree_to_string(const vector<nodestruct>& v) {
    vector<string> node_str(v.size(),"");
    char buffer[16];
    for(size_t i=0; i<v.size(); i++) {
        string temp = "";

        if(v[i].child_1 != -1 && v[i].child_2 != -1) {
            string convert_node1=id_to_string(v[i].child_1);
            string convert_node2=id_to_string(v[i].child_2);
			temp += "(" + node_str[v[i].child_1] + species_label(v[v[i].child_1].species) + convert_node1 + ":";
            sprintf(buffer, "%0.6g", v[v[i].child_1].time);
			temp += buffer;
            temp += "," + node_str[v[i].child_2] + species_label(v[v[i].child_2].species) + convert_node2+ ":";
            sprintf(buffer, "%0.6g", v[v[i].child_2].time);
            temp += buffer;
			temp += ")";
        }
        node_str[i] = temp;
    }
	string temp = species_label(v.back().species) + id_to_string((int)v.size() - 1);
    return node_str.back() + temp + ";";
}
//-----------------------------------------------------------------------------//
// Function to construct vector of mutation labels for tree, in numerical order
string mutation_string(const vector<nodestruct>& t)
{
	string temp = "[";
	for(size_t i=0; i<t.size(); i++) {
		temp += t[i].label;
	}
	temp += ']';
	return temp;
}

//-----------------------------------------------------------------------------//
int set_mutations(xorshift64 &myrand1, char &G, double time)
{   
    int counter = 0;
	double m = rand_exp(myrand1); //m = total distance travelled along branch lengthf
    while (m <= time) { //if m < branch length --> mutate
        ++counter;  //muation counter
		// use the alias tables to effeciently sample the result of the mutation
		G = static_cast<char>(mutation[(size_t)G](myrand1.get_uint64()));
        m += rand_exp(myrand1);
    }
	return counter;
}
