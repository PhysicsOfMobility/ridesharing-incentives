#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <ctime>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <stdio.h>
#include <sstream>
#include "PerfectMatching.h"
#include "GEOM/GeomPerfectMatching.h"

struct PerfectMatching::Options options;
using namespace std;

// GLOBAL VARIABLES
random_device rd;
mt19937 gen(rd());

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 6-Branch, 2 Ring-Topology
int ID_origin = 0;
const int V = 13;

// DISTANCE MATRIX
double dist_matrix[][V] = {{0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},{1.0,0.0,1.0471975511965976,2.0,2.0,2.0,1.0471975511965976,1.0,2.0471975511965974,3.0,3.0,3.0,2.0471975511965974},{1.0,1.0471975511965976,0.0,1.0471975511965976,2.0, 2.0, 2.0, 2.0471975511965974,1.0,2.0471975511965974,3.0, 3.0, 3.0},{1.0, 2.0,1.0471975511965976,0.0,1.0471975511965976,2.0,2.0,3.0,2.0471975511965974,1.0,2.0471975511965974,3.0,3.0},{1.0,2.0,2.0,1.0471975511965976,0.0,1.0471975511965976,2.0,3.0,3.0,2.0471975511965974,1.0,2.0471975511965974,3.0},{1.0,2.0,2.0,2.0,1.0471975511965976,0.0,1.0471975511965976,3.0,3.0,3.0,2.0471975511965974,1.0,2.0471975511965974},{1.0,1.0471975511965976,2.0,2.0,2.0,1.0471975511965976,0.0,2.0471975511965974,3.0,3.0,3.0,2.0471975511965974,1.0},{2.0,1.0,2.0471975511965974,3.0,3.0,3.0,2.0471975511965974,0.0,2.0943951023931953,4.0,4.0,4.0,2.0943951023931953},{2.0,2.0471975511965974,1.0,2.0471975511965974,3.0,3.0,3.0,2.0943951023931953,0.0,2.0943951023931953,4.0,4.0,4.0},{2.0,3.0,2.0471975511965974,1.0,2.0471975511965974,3.0,3.0,4.0,2.0943951023931953,0.0,2.0943951023931953,4.0,4.0},{2.0,3.0,3.0,2.0471975511965974,1.0,2.0471975511965974,3.0,4.0,4.0,2.0943951023931953,0.0,2.0943951023931953,4.0},{2.0,3.0,3.0,3.0,2.0471975511965974,1.0,2.0471975511965974,4.0,4.0,4.0,2.0943951023931953,0.0,2.0943951023931953},{2.0,2.0471975511965974,3.0,3.0,3.0,2.0471975511965974,1.0,2.0943951023931953,4.0,4.0,4.0,2.0943951023931953,0.0}};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// PRE-DECLARATION OF FUNCTIONS
class Request;
class Fitness;
class edge;

int hom_dest(int, int);
double provider_dist_saving(int,int,int);

////////////////////////////////////////////////////
class Request
{
    public:

        // Static variables
        int ID; //Request ID
        int destination; // Destination node

        double xi; // Inconvenience per distance
        double zeta; // Opportunity cost per distance

        // Dynamic variables
        int t; // Time
        bool share; // Sharing decision
        double prob; // Willingness to share
        int Match_ID; // Request ID that ID is paired with

        double dist_single; // Shortest path distance
        double detour; // Detour from sharing the trip with Match_ID
        double dist_together; // Distance ID and Match_ID spend together in the vehicle
        double round_trip_shared; // Distance of shared round trip

        double cost_financial; // Financial cost of trip
        double cost_opportunity; // Opportunity cost of trip
        double cost_inconvenience; // Cost of inconvenience for trip
        double u; //Utility for the trip

        // Functions
        Request(); // Constructor for request
        void init_request(int,int,double,double,double);
        void set_decision(); // Sharing decision
        void set_utility(double,double,double); // Update utility
};

        // Constructor
        Request::Request()
        {
            // Static variables
            ID = 0; //Request ID
            destination = 0; // Destination node

            zeta = 0; // Inconvenience per distance
            xi = 0; // Opportunity cost per distance

            // Dynamic variables
            t = 0; // Time
            prob = 0; // Willingness to share
            set_decision(); // Sharing decision
            Match_ID = 0; // Request ID that ID is paired with

            dist_single = 0; // Shortest path distance
            detour = 0; // Detour from sharing the trip with Match_ID
            dist_together = 0; // Distance ID and Match_ID spend together in the vehicle
            round_trip_shared = 0; // Round trip distance

            cost_financial = 0; // Financial cost of trip
            cost_opportunity = 0; // Opportunity cost of trip
            cost_inconvenience = 0; // Cost of inconvenience for trip

            u = 0; //Utility for the trip
        };

        void Request::init_request(int ID0, int destination0, double zeta0, double xi0, double prob0)
        {
                // Static variables
                ID = ID0; //Request ID
                destination = destination0; // Destination node

                zeta = zeta0; // Inconvenience per distance
                xi = xi0; // Opportunity cost per distance

                // Dynamic variables
                prob = prob0; // Willingness to share
                set_decision(); // Sharing decision
                dist_single = dist_matrix[ID_origin][destination]; // Shortest path distance
        };

        void Request::set_decision()
        {
            share = false;
            bernoulli_distribution d(prob);
            share = d(gen);
        };

        void Request::set_utility(double u_Single, double p, double epsilon)
        {
            cost_financial = share*p*dist_single*epsilon; // Financial cost increment of sharing depends on sharing decision
            cost_opportunity = share*xi*detour; // Opportunity cost are proportional to xi and the detour
            cost_inconvenience = share*zeta*dist_together; // Cost of inconvenience increment depends on distance together

            u = u_Single + (cost_financial - cost_opportunity - cost_inconvenience); // Utility of trip

            t++;    // Update time step
        };

////////////////////////////////////////////////////
class Fitness
{
    public:
        // Static variables
        int destination; //Destination ID

        // Dynamic variables
        int t; // Time
        double prob; // Willingness to share
        int counter; // Counter
        double fitness_single; // Fitness of not sharing
        double fitness_share; // Fitness of sharing
        double detour_share; // Expected detour of sharing
        double dist_together; // Expected distance spent sharing

        // Functions
        Fitness(); // Constructor for fitness class
        void init_fitness(int, int, double, double);
        void update_fitness(double, double, double);
        double mean_fitness();
        double mean_detour();
        double mean_dist_together();
        void replicator_dynamics();
        stringstream write(bool,int,double,double,double,double);
};

        // Constructor
        Fitness::Fitness()
        {
            destination = 0;
            t = 0;
            prob = 0;
            counter = 0;
            fitness_single = 0;
            fitness_share = 0;
            detour_share = 0;
            dist_together = 0;
        };

        // Initialization
        void Fitness::init_fitness(int t0, int destination0, double prob0, double u_single)
        {
            t = t0;
            prob = prob0;
            destination = destination0;
            fitness_single = u_single;
        };

        // Integrate results from previous matching
        void Fitness::update_fitness(double fitness_share_T, double detour_share_T, double dist_together_T)
        {
            fitness_share += fitness_share_T;
            detour_share += detour_share_T;
            dist_together += dist_together_T;
            counter++;
        };

        // Calculate conditional averages
        double Fitness::mean_fitness(){
            return fitness_share/counter;
        };

        double Fitness::mean_detour(){
            return detour_share/counter;
        };

        double Fitness::mean_dist_together(){
            return dist_together/counter;
        };

        // Replicator dynamics
        void Fitness::replicator_dynamics(){
            // Replicator step
            prob*=fitness_share/counter/(fitness_single*(1-prob)+fitness_share/counter*prob);

            t++; // Time step
            counter = 0; // Reset counter
            fitness_share = 0; // Reset fitness
            detour_share = 0; // Reset detour
            dist_together = 0; // Reset dist_together
        };

        stringstream Fitness::write(bool mode,int S,double epsilon,double zeta,double xi,double u_single){
            stringstream output;

            if(mode)
                output << "S,Eps,Zeta,Xi,Time,Destination,Prob,Fitness_Single,Fitness_Share,Detour_Share,Distance_Together\n";
            else
                output << S << "," << epsilon << "," << zeta << "," << xi << "," << t << "," << destination << "," << prob << "," << u_single << "," << fitness_share/counter << "," << detour_share/counter << "," << dist_together/counter << endl;

            return output;
        };

class edge
{
    public:
        int i;
        int j;
        double weight;

        void init(int I,int J,double W){
            i=I;
            j=J;
            weight=W;
        };
};

// Initialize request vector
void init_requests(vector<Request> &requests, int l, double zeta, double xi, vector<Fitness> fitness_T)
{
    // Initialize the different requests
    for(int i = 1; i<=requests.size(); i++)
    {
        if(i==1)
            requests[i-1].init_request(i,l,zeta,xi,true); // Destination node l, and sharing
        else
        {
            int dest = hom_dest(ID_origin, V); // homogeneous OD distribution
            requests[i-1].init_request(i,dest,zeta,xi,fitness_T[dest].prob);
        }
    }
};

// Filter request vector by shared ride requests
vector<Request> filter_sharing(vector<Request> requests)
{
    vector<Request> requests_Sharing;

    for (auto& c : requests)
    if (c.share == true)
        requests_Sharing.push_back(c);

    return requests_Sharing;
};

// Filter request vector by single ride requests
vector<Request> filter_single(vector<Request> requests)
{
    vector<Request> requests_Single;

    for (auto& c : requests)
    if (c.share == false)
        requests_Single.push_back(c);

    return requests_Single;
};

//////////////////////////////////////////////////////////////////////////////////
int hom_dest(int o, int V){
    int dest = o;

    while(dest==o)
        dest = rand()%(V);
    return dest;
};

double provider_dist_saving(int i, int j, int k){
    // i is origin node, j,k are destination node labels

    bool match = false;
    bool order = false;

    // Distances between nodes
    double Dij = dist_matrix[i][j];
    double Dji = dist_matrix[j][i];
    double Dik = dist_matrix[i][k];
    double Dki = dist_matrix[k][i];
    double Djk = dist_matrix[j][k];
    double Dkj = dist_matrix[k][j];

    double round_trip_single = Dij+Dji+Dik+Dki;
    double round_trip_shared_ijk = Dij+Djk+Dki;
    double round_trip_shared_ikj = Dik+Dkj+Dji;
    double round_trip_shared = min(round_trip_shared_ijk,round_trip_shared_ikj);

    // If a shared trip is shorter than two single strips, offer sharing
    if (round_trip_single > round_trip_shared)
        match = true;

    if(match)
        return round_trip_single-round_trip_shared;
    else
        return -1;
};

void rider_shared_distances(int i, vector<Request> &requests, int ID1, int ID2)
{
    // Update matching
    requests[ID1-1].Match_ID = ID2;
    requests[ID2-1].Match_ID = ID1;

    // Get destination information (i denotes joint origin)
    int j = requests[ID1-1].destination;
    int k = requests[ID2-1].destination;

    // Order of sharing (true: i->j->k, false: i->k->j)
    bool order = false;

    // Distances between nodes
    double Dij = dist_matrix[i][j];
    double Dji = dist_matrix[j][i];
    double Dik = dist_matrix[i][k];
    double Dki = dist_matrix[k][i];
    double Djk = dist_matrix[j][k];
    double Dkj = dist_matrix[k][j];

    double round_trip_shared_ijk = Dij+Djk+Dki;
    double round_trip_shared_ikj = Dik+Dkj+Dji;

    double round_trip_single = Dij+Dji+Dik+Dki;
    double round_trip_shared = min(round_trip_shared_ijk,round_trip_shared_ikj);

    // Determine order of sharing
    if(round_trip_shared_ijk > round_trip_shared_ikj)
        order = false;
    else if(round_trip_shared_ijk == round_trip_shared_ikj){
        if(Dij == Dik){
            // Toss a fair coin
            bernoulli_distribution d(0.5);
            order = d(gen);
        }
        else if(Dij < Dik){
            // Drop j first and online then drive to k
            order = true;
        }
        else
            order = false;
    }
    else
        order = true;

    // Inconvenience
    double dist_together = order*Dij+(1-order)*Dik;

    requests[ID1-1].dist_together = dist_together;
    requests[ID2-1].dist_together = dist_together;

    // Detour
    requests[ID1-1].detour = (1-order)*(Dik+Dkj-Dij);
    requests[ID2-1].detour = order*(Dij+Djk-Dik);

    // Distance saved
    requests[ID1-1].round_trip_shared = round_trip_shared; // save round trip shared to request
    requests[ID2-1].round_trip_shared = round_trip_shared; // save round trip shared to request
};

void init_Matching_options(){
    options.fractional_jumpstart = false;
    options.dual_greedy_update_option = 0;
    options.dual_LP_threshold = 0.00;
    options.update_duals_before = false;
    options.update_duals_after = false;
    options.single_tree_threshold = 1.00;
    options.verbose = false;
};

void Matching_Algorithm(int S_Share, vector<Request> &requests, vector<Request> &requests_Sharing)
{
    // Create shared request graph for the S_Share ride requests
    edge e;
    vector<edge> edge_list;
    double weight=0;

    // CREATE REQUEST GRAPH (iterates over requests_Sharing)
    for(int i=0; i<S_Share-1; i++){
        for(int j=i+1; j<S_Share; j++)
        {
            // Savings potential for provider (include as neg. cost)
            weight = -provider_dist_saving(ID_origin, requests_Sharing[i].destination,requests_Sharing[j].destination);

            // If provider has a savings potential, include request in graph
            if(weight<0){
                e.init(i,j,weight); // Init two identical pairs of edges, shifted by S_Share request nodes
                edge_list.push_back(e);
                e.init(i+S_Share,j+S_Share,weight);
                edge_list.push_back(e);
            }
        }
        e.init(i,i+S_Share,0);
        edge_list.push_back(e);
    }

    e.init(S_Share-1,2*S_Share-1,0);
    edge_list.push_back(e);

    // COMPUTE MINIMUM COST PERFECT WEIGHT MATCHING
    PerfectMatching *pm = new PerfectMatching(2*S_Share,edge_list.size());
    for (int k=0; k<edge_list.size(); k++) pm->AddEdge(edge_list[k].i,edge_list[k].j, (int)(100*edge_list[k].weight));
    pm->options = options;
    pm->Solve();

    int j = 0;
    for (int i=0; i<S_Share; i++){
		j = pm->GetMatch(i);
        if(i<j & j<S_Share){
            rider_shared_distances(ID_origin, requests, requests_Sharing[i].ID, requests_Sharing[j].ID);
        }
    }

    delete pm;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Ride-Sharing Game
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    // PARAMETERS

    // SET PARAMETERS FOR MAXIMUM WEIGHT MATCHING ALGORITHM
    init_Matching_options();

    // REQUEST NUMBER
    int S = 10;

    // RIDER UTILITY PARAMETERS
    double u_Single = 10; // Utility of riding alone

    // FINANCIAL INCENTIVE PARAMETERS
    double p = 1; // Price per distance
    double epsilon = 0.2;// Financial discount of sharing

    // OTHER INCENTIVE PARAMETERS
    double zeta = 0.3; // Inconvenience of sharing
    double xi = 0.3; // Detour effect of sharing

    // REPLICATOR DYNAMICS
    double prob0 = 0.01; // Initial willingness to share
    int T = 50000; // Total number of replicator steps
    int n_iter = 100; // Iterations per replicator step

    // POPULATION FITNESS VECTOR PER OD-PAIR
    vector<Fitness> fitness_T(V);
    for(int i=0; i<V; i++)
        fitness_T[i].init_fitness(0,i,prob0,u_Single);

    /////////////////////////////////////////////////
    //////// THIS IS THE ACTUAL PROGRAM /////////////
    /////////////////////////////////////////////////

    // Iterate over time steps t
    for(int t=0; t<=T; t++){

        // Print current configuration of choice probabilities
        for(int k=1; k<V; k++)
            cout << "Pi(d=" << k << "; t=" << t << "):\t" << fitness_T[k].prob << endl;

        // Iterate over destination nodes and zeta values
        for(int l=1; l<V; l++){ // Destinations

            // Perform n_iter iterations per destination node before updating the willingness to share
            for(int n=1; n<=n_iter; n++){

                // ###########################################################################################
                // # 1. Request sample and sharing decision
                // ###########################################################################################

                    // # 1a. Realize request sample (a) conditional on node l is part of the request sample, (b) l shares (calls willingness to share from fitness vector)
                    vector<Request> requests(S); // Init a vector of size S
                    init_requests(requests, l, zeta, xi, fitness_T); // Initialize request objects

                    // # 1b. Filter sharing requests from request sample
                    vector<Request> requests_Sharing;
                    requests_Sharing = filter_sharing(requests);
                    int S_Share = requests_Sharing.size();

                // ###############################################################################################
                // # 2. Match shared ride requests
                // ###############################################################################################

                    // # 2a. Create a request graph for sharing requests (i.e. request ID as nodes, edge weights denote distance savings potential for provider (i.e. determine the route and order of sharing))

                    // # 2a (i) Iteratively build graph object by passing sharing request IDs as nodes
                    // # 2a (ii) Iteratively build graph object by computing edge weights
                    // # 2b. Perform maximum weight matching on request graph
                    // # 2c. Pass matched request ID to request class; update distance information

                    Matching_Algorithm(S_Share,requests,requests_Sharing);

                // ###########################################################################################
                // # 3. Realize utility per request
                // ###########################################################################################

                    // # 3a. Set utility per request
                    for(int i=0; i<S; i++){
                        requests[i].set_utility(u_Single, p, epsilon);
                    }

                    // # 3b. Pass utility for selected node l to fitness function
                    fitness_T[l].update_fitness(requests[0].u, requests[0].detour, requests[0].dist_together);
                }
            }

            // ###################################################################################################                    // # 4. Evolve according to replicator dynamics
            // ###################################################################################################

                // # 4a. Average over the n_iterations to determine the population fitness of sharing at destination node
                // # 4b. Call replicator dynamics from fitness array
                // # 4c. Write fitness vector to file
                // # 4d. Reset fitness array for next replicator time step

                for(int k=1; k<V; k++)
                    fitness_T[k].replicator_dynamics();
    }

    return 0;
}
