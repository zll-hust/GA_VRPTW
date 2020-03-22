#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
#include <vector>
#include <cstdlib>

using namespace std; 

#define cin fin
#define MAX 0x7fffffff
#define MIN 0x80000000

int const customerNr = 25; // customer numbers
int const capacity = 200 ; // capacity of the vehicle
double distanceMatriax[customerNr + 1][customerNr + 1];
int const interation = 1000;
int const mutationInteration = 50;
double const mutationRate = 0.1;
double const selectRate = 0.1;
int const populationNr = 100; // the number of population 

ifstream fin("C101.txt");

//************************************************************

// the class we will use
class Customer;
class Route;
class Solution;
class Chromosome;

// input instance from .txt file
void inputInstance();

// strat the Genetic Algorithm 
Chromosome GeneticAlgorithm();

// initialize the population
Chromosome * initialize();

// selecte the better chromosomes to crossover
Chromosome * selection(Chromosome *);

// crossover the parents chromosomes to create new chromosomes
Chromosome crossover(Chromosome, Chromosome);

// mutate the chromosome
Chromosome mutation(Chromosome chromosome);

// neighborhood move to select better mutation
Chromosome neighborhoodMove(Chromosome);

// get best chromosome
Chromosome getBest(Chromosome *chromosomes);

//************************************************************
class Customer {
    public:
    	int number; // number of customer s
    	double X, Y; // candidate of customers
    	double begin, end, service; // the time windows and service time
    	double demand; // demand of customers
} ;
Customer customers[customerNr + 1];
//************************************************************
class Route {
	public:
		double cost;
		vector <int> cusList;
		bool feasible;
		bool check();
		double getCost();
};

bool Route::check() {
	bool flag = true;
	
	// check capacity constriant
	int load = 0;
    for(int i = 0; i < cusList.size(); i++)
        load += customers[cusList[i]].demand;
    if (load > capacity) flag = false;
    
    // check time windows constriant
    double time = 0;
    time += distanceMatriax[0][cusList[0]];
    if(time > customers[cusList[0]].end) flag = false;
    for(int i = 1; i < cusList.size(); i++){
    	time += customers[cusList[i - 1]].service;
    	time += distanceMatriax[cusList[i - 1]][cusList[i]];
        if(time > customers[cusList[i]].end) flag = false;
        time = max(customers[cusList[i]].begin , time);
    }
    return flag;
}

double Route::getCost(){
	cost = 0;
    cost += distanceMatriax[0][cusList[0]];
    cost += distanceMatriax[0][cusList[cusList.size() - 1]];
    if(cusList.size() > 1) 
        for (int i = 1; i < cusList.size(); i++) 
             cost += distanceMatriax[cusList[i - 1]][cusList[i]];
    return cost;
}
//************************************************************
class Solution {
	public:
		vector <Route> routeList;
		double totalCost;
		double getCost();
		void print();
		Chromosome toChromosome();
}; 

double Solution::getCost(){
	totalCost = 0;
	for (int i = 0; i < routeList.size(); i++){
		totalCost += routeList[i].getCost();
	}
	return totalCost;
}

void Solution::print(){
	for (int i = 0; i < routeList.size(); i++){
		cout << "Route " << i << " : 0 -> ";
		for (int j = 0; j < routeList[i].cusList.size(); j++){
			cout << routeList[i].cusList[j] << "-> ";
		} 
		cout << "0 ." << endl;
		cout << "\t check : " << routeList[i].check() << "\t route cost : " << routeList[i].cost << endl;
	}
	cout << "total cost of the solution : " << totalCost << endl; 
}
//************************************************************
class Chromosome {
	public:
		vector <int> geneList;
		double fitness; // we use 1 / total cost to represent fitness
		Chromosome();
		void setFitness();
		Solution toSolution();
};

Chromosome :: Chromosome() {
	// randomly create geneList
	vector <int> temp;
	for (int i = 1; i <= customerNr ; i++)
		temp.push_back(i);
	
	for (int i = customerNr; i > 0; i--) {
        int index = rand() % i;
        geneList.push_back(temp[index]);
        temp.erase(temp.begin() + index);
    }

} 
//************************************************************
Solution Chromosome::toSolution(){
	Solution solution;
	double time = 0;
	int cost = 0;
	bool if_fir = true;
	Route route;
	for (int i = 0; i < customerNr; i++){
		if(if_fir){
			route.cusList.push_back(geneList[i]);
			time = distanceMatriax[0][geneList[i]];
			time += customers[geneList[i]].service;
			time += distanceMatriax[geneList[i]][0];
			cost = customers[geneList[i]].demand;
			if_fir = false;
		}else{
			double nextTime = time - distanceMatriax[geneList[i - 1]][0] + distanceMatriax[geneList[i - 1]][geneList[i]];
			int nextCost = cost + customers[geneList[i]].demand;
			if (nextTime > customers[geneList[i]].end || nextCost > capacity ){
				route.getCost();
				solution.routeList.push_back(route);
				Route newroute;
				route = newroute;
				if_fir = true;
				i--;
		    }else{
		    	route.cusList.push_back(geneList[i]);
		    	time = max(nextTime,customers[geneList[i]].begin);
			    time += customers[geneList[i]].service;
			    time += distanceMatriax[geneList[i]][0];
			    cost += customers[geneList[i]].demand;
	      	}
		}
    }
    if (route.cusList.size() != 0){
		route.getCost();
		solution.routeList.push_back(route);
	}
	solution.getCost();
	return solution;	
}

void Chromosome::setFitness(){
	fitness = 1 / toSolution().getCost();	
}

Chromosome Solution::toChromosome(){
	Chromosome chromosome;
	for (int i = 0; i < routeList.size(); i++){
		for (int j = 0; j < routeList[i].cusList.size(); j++){
			chromosome.geneList.push_back(routeList[i].cusList[j]);
		}
	}
	return chromosome;
}

//************************************************************
int main(){
	srand ( ( unsigned ) time ( NULL ) );
	
	time_t start,finish;
    start = clock(); // strat the clock
	
	cout << "waiting for a while ......" << endl << endl; 
	inputInstance();
	Chromosome best = GeneticAlgorithm();
	Solution solution = best.toSolution();
	solution.print();
	
	finish = clock(); // stop the clock
    double duration = ((double)(finish-start))/CLOCKS_PER_SEC;
    cout<<endl;
    cout<<"time used : "<< duration <<" s ."<<endl;
	
	return 0;
}
//************************************************************
void inputInstance(){
	for ( int i = 0; i <= customerNr; i++ )
        cin >> customers[i].number >> customers[i].X >> customers[i].Y >> customers[i].demand
            >> customers[i].begin >> customers[i].end >> customers[i].service;
     	
    for ( int i = 0; i <= customerNr; i++ ){
    	for ( int j = i; j <= customerNr; j++){
    		if ( i == j ){
    			distanceMatriax[i][j] == 0;
			}else{
				distanceMatriax[i][j] = distanceMatriax[j][i] = 
				sqrt((customers[i].X - customers[j].X) * (customers[i].X - customers[j].X) + 
				(customers[i].Y - customers[j].Y) * (customers[i].Y - customers[j].Y));
			}
		}
	}
}
//************************************************************
Chromosome GeneticAlgorithm(){
	Chromosome best;
	double bestFitness = MIN;
    Chromosome *parents = initialize();
    for(int i = 0; i < interation; i++){
    	// select
        Chromosome *mid = selection(parents);
        Chromosome *childrens = new Chromosome[populationNr];
        // cross
        for(int j = 0; j < populationNr; j++){
            int fir = rand() % populationNr;
            int sec = rand() % populationNr;
            childrens[j] = crossover(mid[fir],mid[sec]);
        }
        // mutation
        for(int p=0; p < populationNr; p++){
        	int rate = rand() % 100;
        	if (rate < 100 * mutationRate){
        		childrens[p] = mutation(childrens[p]);
			}
		}
        // update the population
        for(int q = 0; q < populationNr; q++)
            parents[q] = childrens[q];
        // update the maxmum fitness
        if(bestFitness < getBest(parents).fitness){
            best = getBest(parents);
            bestFitness = best.fitness;
        }
    }
    return best;
}
//************************************************************
Chromosome * initialize(){
	Chromosome *parents = new Chromosome[populationNr];
    for (int i = 0; i < populationNr; i++) 
        parents[i].setFitness();
    return parents;
}
//************************************************************
Chromosome * selection(Chromosome *chromosomes){
	 Chromosome *childrens = new Chromosome[populationNr];
     for(int t = 0; t < populationNr; t++) {
        int i = rand() % populationNr;
        int j = rand() % populationNr;
        if (chromosomes[i].fitness > chromosomes[j].fitness)
            childrens[t] = chromosomes[i];
        else
            childrens[t] = chromosomes[j];
    }
    return  childrens;
}
//************************************************************
Chromosome crossover(Chromosome p_one,Chromosome p_two){
	Chromosome children;
	vector<int> vecEmpty;
	children.geneList.swap(vecEmpty);
	//children.geneList.clear();

    // get random number 0 <= i < j < customerNr
    int j = rand() % customerNr;
    int i = rand() % customerNr;
    if (i > j){
    	int temp = i;
    	i = j;
    	j = temp;
	}
    
    // push the gene from i to j to mid
    for(int tmp = i;tmp <= j; tmp++) 
        children.geneList.push_back(p_one.geneList[tmp]);
        
    // change the order and push gene of p_two into parent
    vector <int> parent;
    for(int tmp = j; tmp < customerNr; tmp++)
        parent.push_back(p_two.geneList[tmp]);
    for(int tmp = 0; tmp < j; tmp++)
        parent.push_back(p_two.geneList[tmp]);
    
    int tmp = j;
    bool if_fir = false; // if tmp is the first gene

    for(int t = 0; t < parent.size(); t++){
    	// if tmp arrive the last gene, strart from the first gene
        if(tmp == customerNr - 1){ 
            tmp = 0;
            if_fir = true;
        }
        // if the list doesn't include the gene, add it
        bool if_include = false;
        for (int n = 0; n < children.geneList.size() ; n++)
        	if (parent[t] == children.geneList[n])
        		if_include == true;
        if (!if_include){
        	if (if_fir)
        		children.geneList.insert(children.geneList.begin(), parent[t]);
			else
				children.geneList.push_back(parent[t]);
			tmp++;
		}
    }
    
    children.setFitness();
    return children;
}
//************************************************************
Chromosome mutation(Chromosome chromosome){
	Chromosome newChromosome = chromosome;
	int interation = 0;
    while(interation < mutationInteration){
        if(chromosome.fitness < neighborhoodMove(newChromosome).fitness)
            return newChromosome;
        else
            newChromosome = chromosome;
        interation++;
    }
    return  newChromosome;
}

Chromosome neighborhoodMove(Chromosome chromosome){
    int i = rand() % customerNr;
    int j = rand() % customerNr;
    int temp = chromosome.geneList[i];
    chromosome.geneList[i] = chromosome.geneList[j];
    chromosome.geneList[j] = temp;
    chromosome.setFitness();
    return chromosome;
}
//************************************************************
Chromosome getBest(Chromosome *chromosomes){
    Chromosome chromosome;
    double bestFitness = MIN;
    for(int i = 0; i < populationNr; i++){
        if(bestFitness < chromosomes[i].fitness){
            bestFitness = chromosomes[i].fitness;
            chromosome = chromosomes[i];
        }
    }
    return chromosome;
}
//************************************************************
