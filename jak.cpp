// jak.cpp

#include "stdafx.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <float.h>
#include <fstream>
#include "xorshift64.h"
#include <process.h>



using namespace std;
//----------------------------------------------------------------------------//
struct nodestruct {                                                             //create node structure
    int child_1;
    int child_2;
    int parent;
    char label;                                                                //gene
    double time;
    char type;                                                                 //species
};
//----------------------------------------------------------------------------//
string convert(int x)					//Function to convert int type to string
{   string convert_x;
    ostringstream convert1;
    convert1<< x;
    convert_x = convert1.str();
    return(convert_x);
}
//-----------------------------------------------------------------------------//
string tree_to_string(const vector<nodestruct>& v) {                            //create newick tree from node data
    vector<string> node_str(v.size(),"");
    char buffer[16];
    for(int i=0; i<v.size(); i++) {
        string temp = "";

        if(v[i].child_1 != -1 && v[i].child_2 != -1) {
            string convert_node1=convert(v[i].child_1);
            string convert_node2=convert(v[i].child_2);
            temp += "(" + node_str[v[i].child_1] + "_" + v[v[i].child_1].type  + "_"+ convert_node1 + ":";
            sprintf(buffer, "%0.6f", v[i].time-v[v[i].child_1].time);
            temp += string(buffer) + "," + node_str[v[i].child_2] + "_" + v[v[i].child_2].type  + "_" + convert_node2+ ":";
            sprintf(buffer, "%0.6f", v[i].time-v[v[i].child_2].time);
            temp += string(buffer) + ")";
        }
        temp += v[i].label;
        node_str[i] = temp;
    }
    return node_str.back() + ";";
}
//----------------------------------------------------------------------------//
//REED: use gidpid instead of _getpid for portability to other operating systems.
//      on windows use _getpid via a define
#ifdef WIN32
#	define getpid _getpid
#endif
inline unsigned int create_random_seed() {															//random seed generator
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
#ifdef WIN32
#	undef getpid
#endif


//-----------------------------------------------------------------------------//
void set_mutations(xorshift64 myrand1, int G, double& m, double branchA, double branchB, int& m_counter, char& label)
{   double mutations[4][4]={
        {0.25,0.50,0.75,1.0},
        {0.25,0.50,0.75,1.0},
        {0.25,0.50,0.75,1.0},
        {0.25,0.50,0.75,1.0}
    };

    m = m - log(myrand1.get_double52());                                            //m = total distance travelled along branch length
    if (m <= branchA - branchB) {                                                   //if m < branch length --> mutate
        m_counter = m_counter + 1;                                                  //muation counter
    }
    double rand3=myrand1.get_double52();
    if (rand3<=mutations[G][0])
        label='A';
    else if (rand3<=mutations[G][1])
        label='T';
    else if (rand3<=mutations[G][2])
        label='C';
    else if (rand3<=mutations[G][3])
        label='G';
}
//-----------------------------------------------------------------------------------------------//

int max_element(vector<int>& activelist)
{
    int max=0;
    for(int i=0; i<activelist.size(); i++)
    {   if(max<activelist[i])
            max=activelist[i];
    }

    return max;
}
//-----------------------------------------------------------------------------------------------//
void coaltree(vector<int>& activelist, double theta, double time,
	          vector<nodestruct>& nodeVector, xorshift64& myrand1)
{
    double T = 0.0;
	int random1, random2;
    double Z = myrand1.get_double52();
	
	int size = activelist.size();  

	while(size>1)
	{
		double mean = (2.0/(size*(size-1.0)))*(theta/2.0);
		double U = (-log(Z))*mean;
		if(T+U>time)
		{
			break;
		}
		T+=U;

		random1 = (myrand1.get_uint32()% size);
		do {
			random2 = (myrand1.get_uint32()% size);
		} while(random1==random2);
		
		if (random1>random2) 															//orders two nodes minimum to maximum
			swap(random1,random2);

		int newparent = nodeVector.size();
		nodeVector.push_back(nodestruct());

		nodeVector[newparent].child_1 = activelist[random1];                        //update parent node
		nodeVector[newparent].child_2 = activelist[random2];
		nodeVector[newparent].time = T;

		nodeVector[activelist[random1]].parent = newparent;                         //update child nodes
		nodeVector[activelist[random2]].parent = newparent;
		nodeVector[activelist[random1]].time = T - nodeVector[activelist[random1]].time;
		nodeVector[activelist[random2]].time = T - nodeVector[activelist[random2]].time;
		
		activelist[random1] = newparent;													 //update active vector
		activelist.erase (activelist.begin() + random2);
		size--;
	} 
	for(int i=0; i<size; i++)
	{
		nodeVector[activelist[i]].time = nodeVector[activelist[i]].time - time;
	}
}

//-------------------------------------------------------------------------------------------------//

int main(int argc, char *argv[])														 //receive inputs
{
    int N1, N2, n, N, trees;
    double mean, theta1, theta2, theta3, t1, t2, total_tree=0;

//fix input validation
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
        cout << "error, must have format: prg name, trees, # of tips for first species, # of tips for second species, theta1, theta2, theta3, t1, t2" << endl;
//REED: system("PAUSE") does not work outside of windows.  See below for a better solution using cin.ignore
        system("PAUSE");
        return EXIT_FAILURE;
    }

    n = N1+N2;																		//n=total original tips from both species
    N = 2*n-1;                                                                      //N = total # of nodes
    mean = 2.0/(n*(n-1));                                                           //calculate mean
    

    xorshift64 myrand;																//use xorshift64 class for random number generator
    myrand.seed(create_random_seed());

	
	//REED: What is this for?  It is going to fail on all systems but Kailey's
    ofstream myfile;                                                                //file
    myfile.open ("C://Users//Kailey//Documents//MATLAB//newickstruct_data.txt");    //open file

    for(int repeat=0; repeat<trees; repeat++) {                                     //loops once for each tree
        vector<nodestruct> nodevector(N);                                           //create nodevector (vector of structs)

		
		vector<int> active1(N1);			//initialize active list for species 1
		for(int i=0; i<N1; i++)
		{
			active1[i]=i;
			nodevector[i].child_1=-1;
			nodevector[i].child_2=-1;
			nodevector[i].parent=-1;
			nodevector[i].label='N';
			nodevector[i].time=0;
		}


		vector<int> active2(N2);			//initialize active list for species 2
		for(int j=N1; j<N1+N2; j++)
		{   
			active2[j-N1]=j;
			nodevector[j].child_1=-1;
			nodevector[j].child_2=-1;
			nodevector[j].parent=-1;
			nodevector[j].label='N';
			nodevector[j].time=0;
		}
//----------------------------------------------------------------------------//

        double t=0.0, tN1=0.0, tN2=0.0, var=0.0;

        vector<int> nodes(2*n-1); 													//n is number of initial nodes	
       
		coaltree(active1, theta1, t1, nodevector, myrand);
		coaltree(active2, theta2, t2, nodevector, myrand);

        vector<int> active3(active1.size() + active2.size());

        active1.insert(active1.end(), active2.begin(), active2.end());

//REED: No don't do this it is going to produce bad results.		
        double t3;
        if(t1>t2)
            t3=t1;
        else
            t3=t2;


        while(active1.size()>1)
        {
            coaltree(active1, theta3, t3, myrand, nodevector);
            //cout<<"species 3 "<<endl;
        }


//----------------------------------------------------------------------------////mutations
        char L1;
        int d  = 0;                                                                     //mutation counters
        int d1 = 0;

//REED: Once you create the tree, it shouldn't matter what population the nodes came from.		
        for (int i = N - 1; i > B; --i) {                                               //start mutations for loop
            double T1 = nodevector[i].time;                                              //time @ current node
            double T2 = nodevector[nodevector[i].child_1].time;                          //time @ child 1
            double T3 = nodevector[nodevector[i].child_2].time;			                //time @ child 2
            double m = 0.0;                                                              //initialize/reset m
//node i child 1
            while (T1 - m >= T2) {                                                      //branch between current node and child_1
                set_mutations(myrand, G, m, T1, T2, d, L1);
            }                                                                           //close while loop (node i child 1)
// node i child 2
            m = 0.0;                                                                    //reset m = 0
            while (T1 - m >= T3) {                                                      //branch between current node and child_2
                set_mutations(myrand, G, m, T1, T3, d1, L1);
                nodevector[i].label = L1;												//store mutated gene as nodevector[i].label
            }                                                                           //close while loop (child 2)
        }                                                                               //end mutations for loop
//----------------------------------------------------------------------------//
        cout << "Newick tree: " << repeat+1<< endl;
        cout << tree_to_string(nodevector) << endl;                                 //print newick tree to console
        myfile << tree_to_string(nodevector)<< " \n";                               //print newick tree to file
        int dtotal = d + d1;
        cout << "Number of mutations: " << dtotal << endl;

        total_tree=total_tree+t;

    }                                                                               //end # of trees loop
    myfile.close();

//double tree_height_avg=total_tree/trees;

    cout<<"Random seed used: "<<create_random_seed()<<endl;

//REED: Use this instead.	
    cin.ignore( numeric_limits<streamsize>::max(), '\n' );
    cout << "Press ENTER to quit.";
    cin.ignore( numeric_limits<streamsize>::max(), '\n' );

    return EXIT_SUCCESS;
}




