// jak_data.cpp

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[])
{
	std::string line;
	while(std::getline(cin, line)) 
	{
		size_t p = 0;
		for(; p<line.length() && line[p]!='['; ++p)
			/*noop*/;
		++p;

		char major_allele=line[p];
		for(; p<line.length() && line[p]==major_allele; ++p)
			/*noop*/;
	
		cout << ( (line[p] == ']') ? "1\t" : "0\t");

		for(; p<line.length() && line[p]!='('; ++p)
			/*noop*/;

		double branch_total=0.0;

		for(;;)
		{
			for(; p<line.length() && line[p]!=':'; ++p)
				/*noop*/;
			if(p>=line.length())
				break;
			p++;

			char *end;
			const char *begin=line.c_str();
			branch_total+=strtod(begin+p, &end);
			if(end==begin+p)
				break;
			p=end-begin;
		}

		cout << branch_total;
		cout << endl;
	}	

    return EXIT_SUCCESS;

}
