#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <map>
		
struct pin {
	std::string m_name;

};	

struct net {

	std::vector<pin> m_pins;
};

size_t tokenlize(std::string &str, std::vector<std::string> *token) {
	std::stringstream ss(str);
	std::string s;

	while (ss>>s) {
		token->push_back(s);
	}
	return token->size();
}


int main(int argc, char *argv[]) {
	

	std::ifstream MGfile;
	MGfile.open(argv[1], std::ifstream::in);
	std::string in;
	std::vector <net> nets;
	std::vector <pin> pins;
	std::map <std::string, int> mod;
	bool flag = false;
	while (MGfile >> in) {
		if(in == "*SIGNAL*" && flag == 0) {
			flag = 1;
		}
		if(in == "*SIGNAL*") {
			net n;
			while(std::getline(MGfile, in)) {
				if(in == "") break;
				std::vector<std::string> token;
				size_t size = tokenlize(in,&token); 
				if (size == 2) {
					for (size_t i = 0; i < size; ++i) {
						pin p;
						std::stringstream _name(token[i]); 
						std::string s;
						std::getline(_name, s, '.');
						p.m_name = s;
						mod.insert(std::pair<std::string, int>(s,0));
						n.m_pins.push_back(p);
						pins.push_back(p);
					}	
				}
			}
			nets.push_back(n);
		}
	}



	std::ofstream netFile;
	std::string netFileName = "temp.nets"; 
	std::cout << netFileName << " writing...\n";
	netFile.open(netFileName, std::ios::out);
	netFile << "UCLA nets 1.0" << std::endl;
	netFile << std::endl;
	netFile << "#Created    :" << std::endl;
	netFile << "#Created by :" << std::endl;
	netFile << std::endl;
	netFile << "NumNets : " << nets.size()-1 << std::endl;
	netFile << "NumPins : " << pins.size() << std::endl;
	netFile << std::endl;


	for (size_t i = 1; i < nets.size(); ++i) {

		netFile << "NetDegree : " << nets[i].m_pins.size() << std::endl;
		for (size_t j = 0; j < nets[i].m_pins.size(); ++j)
			netFile << "            " << nets[i].m_pins[j].m_name << " B : 0.0 0.0\n";// << nets[i].m_pins[j].m_x << " " << nets[i].m_pins[j].m_y << std::endl;
	}

	netFile.close();


	std::ofstream nodeFile;
	std::string nodeFileName = "temp.nodes"; 
	std::cout << nodeFileName << " writing...\n";
	nodeFile.open(nodeFileName, std::ios::out);
	nodeFile << "UCLA nodes 1.0" << std::endl;
	nodeFile << std::endl;
	nodeFile << "#Created    :" << std::endl;
	nodeFile << "#Created by :" << std::endl;
	nodeFile << std::endl;
	nodeFile << "NumNodes : " << mod.size() << std::endl;
	nodeFile << "NumTerminals : 0 " << std::endl;
	nodeFile << std::endl;

	for (std::map<std::string,int>::iterator it=mod.begin(); it!=mod.end(); ++it) {
	
		nodeFile << it->first << "          10          10\n";
	}

	nodeFile.close();


	std::ofstream plFile;
	std::string plFileName = "temp.pl"; 
	std::cout << plFileName << " writing...\n";
	plFile.open(plFileName, std::ios::out);
	plFile << "UCLA pl 1.0" << std::endl;
	plFile << std::endl;
	plFile << "#Created    :" << std::endl;
	plFile << "#Created by :" << std::endl;
	plFile << std::endl;

	for (std::map<std::string,int>::iterator it=mod.begin(); it!=mod.end(); ++it) {
	
		plFile << it->first << "      0.0     0.0 : N\n";
	}

	plFile.close();

/*

	std::ofstream wtsFile;
	std::string wtsFileName = "temp.wts"; 
	std::cout << wtsFileName << " writing...\n";
	wtsFile.open(wtsFileName, std::ios::out);
	wtsFile << "UCLA wts 1.0" << std::endl;
	wtsFile << std::endl;
	wtsFile << "#Created    :" << std::endl;
	wtsFile << "#Created by :" << std::endl;
	wtsFile << std::endl;

	for (size_t i = 1; i < nets.size(); ++i) {
		wtsFile << nets[i].m_name << " 1" << std::endl;
	}

	wtsFile.close();*/





	return 0;
}
