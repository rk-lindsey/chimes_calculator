#include <mpi.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm> // For sort
#include <iomanip>

#include "chimesFF.h"    

using namespace std;

int nprocs;
int my_rank;

// In .cpp file
namespace GlobalParams {
    vector<double> rcin_list;
    vector<double> rcout_2b_list;
    vector<double> rcout_3b_list;
    vector<double> rcout_4b_list;
    vector<double> morse_lambda_list;
}

int split_lines(string line, vector<string> & items)
{
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    
    string       contents;
    stringstream sstream;

    // Strip comments beginining with ! or ## and terminal new line

    int pos = line.find('!');
      
    if ( pos != string::npos ) 
        line.erase(pos, line.length() - pos);

    pos = line.find("##");
    if ( pos != string::npos ) 
        line.erase(pos, line.length()-pos);

    pos = line.find('\n');
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    sstream.str(line);
     
    items.clear();

    while ( sstream >> contents ) 
        items.push_back(contents);

    return items.size();
}

template <typename T>
int get_index_if(const vector<T>  & vec, const T  & element, vector<bool> & disqualified)
{

    if (disqualified.size() != vec.size())
    {
        cout << "chimesFF: " << "ERROR: get_index_if(...): Qualification criteria does not match vector length" << endl;
        cout << "chimesFF: " << "vec.size(): " << vec.size() << endl;
        cout << "chimesFF: " << "disqualified.size(): " << disqualified.size() << endl;
        exit(0);
    }

    for(int i=0; i<vec.size(); i++)
    {
        if ((vec[i]==element) && (!disqualified[i]))
        {
            disqualified[i] = true;
            return i;
        }
    }

    cout << "chimesFF: " << "ERROR: Could not find element in vector: " << element << endl;
    
    for(int i=0; i<vec.size(); i++)
        cout << "chimesFF: " << "\t" << vec[i] << " " << disqualified[i] << endl;
    
    exit(0);
}

int chimesFF::split_line(string line, vector<string> & items)
{
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    
    string       contents;
    stringstream sstream;

    // Strip comments beginining with ! or ## and terminal new line

    int pos = line.find('!');
      
    if ( pos != string::npos ) 
        line.erase(pos, line.length() - pos);

    pos = line.find("##");
    if ( pos != string::npos ) 
        line.erase(pos, line.length()-pos);

    pos = line.find('\n');
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    sstream.str(line);
     
    items.clear();

    while ( sstream >> contents ) 
        items.push_back(contents);

    return items.size();
}

string chimesFF::get_next_line(istream& str)
{
    // Read a line and return it, with error checking.
    
    string line;

    getline(str, line);
    
    if ( ! str.good() )
    {
        if (my_rank == 0)
            cout << "chimesFF: " << "Error reading line" << line << endl;
        exit(0);
    } 

    return line;
}
void chimesFF::read_parameters(string paramfile)
{
    ifstream param_file;
    param_file.open(paramfile.data());
    
    if (rank == 0)
        cout << "chimesFF: " << "Reading parameters from file: " << paramfile << endl;
    
    if(!param_file.is_open())
    {
        if (rank == 0)
            cout << "chimesFF: " << "ERROR: Cannot open parameter file: " << paramfile << endl;
        exit(0);
    }
    
    bool           found_end = false;
    string         line;
    string         tmp_str;
    vector<string> tmp_str_items;
    int            tmp_no_items;
    int            no_pairs;
    
    while (!found_end)
    {
        line = get_next_line(param_file);

        if(line.find("ENDFILE") != string::npos)
        {
            if (rank == 0)
            {
                cout << "chimesFF: " << "ERROR: Could not find line containing: \" PAIRTYP: CHEBYSHEV\" " << endl;
                cout << "chimesFF: " << "       ...Is this a ChIMES force field parameter file?" << endl;
            }
            exit(0);
        }
        
        if(line.find("PAIRTYP: CHEBYSHEV") != string::npos)
        {
            break;    
        }
    }
    
    param_file.seekg(0);
    
    found_end = false;
    
    while (!found_end)
    {
        line = get_next_line(param_file);
        
        if(line.find("ENDFILE") != string::npos)
            break;        
    
        if(line.find("ATOM PAIRS:") != string::npos)
        {
            tmp_no_items = split_line(line, tmp_str_items);
        
            no_pairs = stoi(tmp_str_items[2]);
        
            if (rank == 0)
                cout << "chimesFF: " << "Will consider " << no_pairs << " atom pair types" << endl;        
        }    
        
        if(line.find("# PAIRIDX #") != string::npos)
        {
            if(line.find("# USEOVRP #") != string::npos)
                continue;
        
            pair_params_atm_chem_1.resize(no_pairs);
            pair_params_atm_chem_2.resize(no_pairs);
            chimes_2b_cutoff.resize(no_pairs);
            morse_var.resize(no_pairs);

            string tmp_xform_style;
            
            for (int i=0; i<no_pairs; i++)
            {
                line = get_next_line(param_file);
                
                tmp_no_items = split_line(line, tmp_str_items);

                int pair_input_version = 0;
				
                if ( tmp_no_items == 8 )
                {
					if ( rank == 0 && i == 0 ) cout << "chimesFF: Detected version 1 pair specification (with S_DELTA)\n";
					pair_input_version = 1;
                }
                else if ( tmp_no_items == 7 )
                {
					if ( rank == 0 && i == 0 ) cout << "chimesFF: Detected version 2 pair specification (no S_DELTA)\n";
					pair_input_version = 2;
                }
                else
                {
					if ( rank == 0 )
					{
						cout << "Incorrect input in line: " << line << endl;
						cout << "Expect 7 or 8 entries\n";
					}
					exit(0);
                }
            
                pair_params_atm_chem_1[i] = tmp_str_items[1];
                pair_params_atm_chem_2[i] = tmp_str_items[2];
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << i << " " << pair_params_atm_chem_1[i] << " " << pair_params_atm_chem_2[i]<< endl;
                
                chimes_2b_cutoff[i].push_back(stod(tmp_str_items[3])); // Inner cutoff    
                chimes_2b_cutoff[i].push_back(stod(tmp_str_items[4])); // Outer cutoff

                int xform_style_idx, morse_idx;
				
                if ( pair_input_version == 1 )
                {
					xform_style_idx = 6;
					morse_idx = 7;
                }
                else if ( pair_input_version == 2 )
                {
					xform_style_idx = 5;
					morse_idx = 6;
                } 
                else
                {
					if ( rank == 0 ) cout << "Bad pair input version\n";
					exit(0);
                }
                    
                if (i==0)
                {
                    tmp_xform_style = tmp_str_items[xform_style_idx];
                }
                else if ( tmp_str_items[xform_style_idx] != tmp_xform_style)    
                {
					if (rank == 0)
						cout << "chimesFF: " << "Distance transformation style must be the same for all pair types" << endl;
					exit(0);
                }

                if (tmp_xform_style == "MORSE" )
                {
					if ( tmp_no_items > morse_idx )
						morse_var[i] = stod(tmp_str_items[morse_idx]);
					else {
						if ( rank == 0 )
							cout << "chimesFF: Missing morse lambda value in line: \n" << line << endl;
						exit(0);
					}
				}
            }
                
            xform_style = tmp_xform_style;
            
            if (rank == 0)
                cout << "chimesFF: " << "Read the following pair type information:" << endl;
            
            for (int i=0; i<no_pairs; i++)
            {
				if (rank == 0)
					cout << "chimesFF: " << "\t" << pair_params_atm_chem_1[i] << " " << pair_params_atm_chem_2[i] << " r_cut_in: " << fixed << right << setprecision(5) << chimes_2b_cutoff[i][0] << " r_cut_out: " << chimes_2b_cutoff[i][1] << " " <<  xform_style;
                
				if (xform_style == "MORSE")
				{
					if (rank == 0)
						cout << " " << morse_var[i] << endl;
				}
				else
					if (rank == 0)
						cout << endl;
			}
        }
    }
    
    param_file.close();    
}

bool get_next_line(istream& str, string & line)
{
    // Read a line and return it, with error checking.
    
        getline(str, line);

        if(!str)
            return false;
    
    return true;
}

// Function to print the adjacency matrix with atom types
void print_adjacency_matrix(const vector<vector<double>>& adjacency_matrix, const vector<int>& atom_types) {
    int n = adjacency_matrix.size();
    for (int i = 0; i < n; i++) {
        // Print atom type first
        cout << setw(4) << atom_types[i] << " |";
        // Print adjacency information
        for (int j = 0; j < n; j++) {
            if(i == j) cout << setw(10) << "N/A" << " ";  // Diagonal is atom type
            else cout << setw(10) << fixed << setprecision(4) << adjacency_matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

double transform(double rcin, double rcout, double lambda, double rij)
{
    double x_min = exp(-1*rcin/lambda);
    double x_max = exp(-1*rcout/lambda);

    double x_avg   = 0.5 * (x_max + x_min);
    double x_diff  = 0.5 * (x_max - x_min);

    x_diff *= -1.0; // Special for Morse style
    
    return (exp(-1*rij/lambda) - x_avg)/x_diff;
}

// Function to read the clusters and construct adjacency matrices
void read_flat_clusters(string clufile, int npairs_per_cluster, vector<pair<vector<vector<double>>, vector<int>>> &adjacency_data, const int body_cnt, bool print = false) {
    ifstream clustream(clufile);
    if (!clustream.is_open()) {
        cerr << "ERROR: Could not open file " << clufile << endl;
        exit(1);
    }

    string line;
    vector<string> line_contents;
    int n_contents;

    while (get_next_line(clustream, line)) {
        n_contents = split_lines(line, line_contents);

        if (n_contents != npairs_per_cluster + body_cnt) { // 4 additional columns for atom types
            cout << "ERROR: Read the wrong number of clusters!" << endl;
            cout << "n_contents: " << n_contents << endl;
            cout << "Expected: " << npairs_per_cluster + 4 << endl;
            exit(0);
        }

        // Extract edge lengths and atom types
        vector<double> edge_lengths(npairs_per_cluster);
        vector<int> atom_types(body_cnt);

        for (int i = 0; i < npairs_per_cluster; i++) {
            // edge_lengths[i] = transform(rcin, rcout, lambda, stod(line_contents[i]));
            edge_lengths[i] = stod(line_contents[i]);
        }

        for (int i = 0; i < body_cnt; i++) {
            atom_types[i] = stoi(line_contents[npairs_per_cluster + i]);
        }

        // Construct the adjacency matrix
        // Assuming a fully connected graphlet with {body_cnt} atoms
        vector<vector<double>> adjacency_matrix(body_cnt, vector<double>(body_cnt, 0.0));

        // Fill the adjacency matrix
        int edge_index = 0;
        for (int i = 0; i < body_cnt; i++) {
            for (int j = i + 1; j < body_cnt; j++) {
                adjacency_matrix[i][j] = edge_lengths[edge_index];
                adjacency_matrix[j][i] = edge_lengths[edge_index]; // Symmetric matrix
                edge_index++;
            }
        }

        // Print with atom types
        if (print == true){
            cout << "Adjacency Matrix (Atom Types: ";
            for (auto at : atom_types) cout << at << " ";
            cout << "):" << endl;
            print_adjacency_matrix(adjacency_matrix, atom_types);
        }
        
        // Store both matrix and types as a pair
        adjacency_data.push_back(make_pair(adjacency_matrix, atom_types));
    }

    clustream.close();
}

int get_bin(double binw, double maxval, double dist)
{
    int bin = floor(dist/binw);

    if (dist == maxval)
        return bin-1;
    else
        return bin;
}

void divide_task(int & my_rank_start, int & my_rank_end, int tasks) 
{
	int procs_used;

	// Deal with no tasks to perform.
	if ( tasks <= 0 ) 
	{
	  my_rank_start = 1 ;
	  my_rank_end = 0 ;
	  return ;
	}

	// Deal gracefully with more tasks than processors.
	if ( nprocs <= tasks ) 
		procs_used = nprocs;
	else
		procs_used = tasks;

	// Use ceil so the last process always has fewer tasks than the other
	// This improves load balancing.
	my_rank_start = ceil( (double) my_rank * tasks / procs_used);

	if ( my_rank > tasks ) 
	{
		my_rank_start = tasks + 1;
		my_rank_end = tasks - 1;
	} 
	else if ( my_rank == procs_used - 1 ) 
	{
		// End of the list.
		my_rank_end = tasks - 1;
	} 
	else 
	{
		// Next starting value - 1 .
		my_rank_end   = ceil( (double) (my_rank+1) * tasks / procs_used ) - 1;
		if ( my_rank_end > tasks - 1 ) 
			my_rank_end = tasks - 1;
	}
}



void gen_flat_hists(vector<double > & clu1, vector<double > & clu2, int n_cluster_pairs, int nbin, double binw, double maxd, string histfile, double rcout, bool same = false)
{
    int                     bin;
    double                  total_dist;
    double                  dist_structure;
    double                  dist_comp_1;
    double                  dist_comp_2;
    vector<long long int>   my_hist(nbin,0);  
    vector<long long int>   hist(nbin,0);  
    long long int           my_nsamples = 0;
    long long int           nsamples = 0;
    int                     my_rank_start;
    int                     my_rank_end;	
    int                     looptwo_start;
    int                     total_tasks; 
    int                     status;
    int                     maxIntValue = numeric_limits<int>::max();
    
    // Distribute outer loop over processors
    
    divide_task(my_rank_start, my_rank_end, clu1.size()/n_cluster_pairs);	// Divide atoms on a per-processor basis.
    total_tasks = my_rank_end-my_rank_start;
    
    if(my_rank ==0)
        cout << "Dividing " << clu1.size()/n_cluster_pairs << " tasks across " << nprocs << " processors" << endl;

    if (total_tasks>0)
    {

    for (int i=my_rank_start; i<=my_rank_end; i++)
    {   

         
        // Print progress
        
        status = double(i-my_rank_start)/(total_tasks)*100.0;
	

	// This logic needed to avoid div by zero when total_tasks/10 is zero (since they are integer types)
        if (my_rank == 0)
		if ((total_tasks/10) == 0)
            		cout << histfile << " Completion percent: " << status << " " << i << " of " << total_tasks << " assigned" << endl;
		else if(i%(total_tasks/10) == 0)
			cout << histfile << " Completion percent: " << status << " " << i << " of " << total_tasks << " assigned" << endl;
        
        
        // Modify bounds in case this is a self-calculation
        
        if (same)
            looptwo_start = i+1;
        else
            looptwo_start = 0;
        
        // Compute the distances
        // Need to determine the flat index for the item
        // Should be i*n_cluster_pairs
        
        for (int j=looptwo_start; j<clu2.size()/n_cluster_pairs; j++)
        {
            total_dist = 0;
    
            for (int k=0; k<n_cluster_pairs; k++)
                // Distance calculation between two clusters
                total_dist += pow(clu1[i*n_cluster_pairs+k] - clu2[j*n_cluster_pairs+k],2.0);

            bin  = get_bin(binw, maxd, sqrt(total_dist));
            
            if (bin > nbin)
            {
                
                cout << "Rank: " << my_rank << " ERROR: computed bin larger than nbins:" << endl;

                cout << "Rank: " << my_rank << " nbin: " << nbin << endl;
                cout << "Rank: " << my_rank << " bin:  " << bin << endl;
                cout << "Rank: " << my_rank << " binw: " << binw << endl;
                cout << "Rank: " << my_rank << " maxd: " << maxd << endl;
                cout << "Rank: " << my_rank << " dist: " << total_dist << endl;
                cout << "Rank: " << my_rank << " clus: " << endl;
                
                for (int m=0;m<n_cluster_pairs; m++)
                    cout << "Rank: " << my_rank << clu1[i*n_cluster_pairs+m] << " <--> " << clu2[j*n_cluster_pairs +m] << endl;
                cout << "Rank: " << my_rank << endl;
                
                exit(0);
            }

            my_hist[bin] += 1;
            my_nsamples += 1;
        }
    }
    }

    if (my_rank == 0)
    {
        cout << "Loop done, printing results: " << endl;
        cout << "Counted nsamples: " << my_nsamples << endl;   
    }

    MPI_Reduce(my_hist.data(), hist.data(), nbin, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_nsamples, &nsamples, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
 

    // Print results
    
    if(my_rank == 0)
    {
        ofstream cluhist;
        cluhist.open(histfile);
        
        for (int i=0; i<hist.size(); i++)
             cluhist << i*binw+0.5*binw <<  " " << double(hist[i]) /nsamples << endl; // nprocs << endl;
        
        cluhist.close();
        
        cout << "Printed." << endl;
    }
 
}

int main(int argc, char *argv[])
{
    my_rank = 0;
    nprocs  = 1;

    MPI_Init     (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

		
    if (my_rank==0)
        cout << "Code compiled in MPI mode."; 
 
    if (my_rank==0)
        cout <<" Will run on " << nprocs << " processor(s)." << endl;
    
    string f1_idx = argv[1]; // "0050"; // Frame 1 of liquid carbon at 1000 K & 0.5 gcc // .2b_clu-r.txt;
    string f2_idx = argv[2]; //"0075"; // Frame 6 of liquid carbon at 1000 K & 0.5 gcc // .2b_clu-r.txt;
    
    // string style = argv[3]; // "s"; // Calc distances based on rij, not transformed sij
         
    double rcout_2b = 5.0;
    double rcout_3b = 5.0;
    double rcout_4b = 4.5;
    
    int nbin_2b = 100;
    int nbin_3b = 100;
    int nbin_4b = 100;

    /////////////////////////////////////////////
    // Read in 2B clusters -- IN TERMS OF sij - determined by user
    /////////////////////////////////////////////
    
    string f1_2b = f1_idx + ".all-2b-clusters.txt"; 
    string f2_2b = f2_idx + ".all-2b-clusters.txt";

    
    vector<pair<vector<vector<double>>, vector<int>>> f1_2b_flat_clusters;
    vector<pair<vector<vector<double>>, vector<int>>> f2_2b_flat_clusters;
    
    int npairs_2b = 1;
    
    read_flat_clusters(f1_2b, npairs_2b, f1_2b_flat_clusters, 2);
    read_flat_clusters(f2_2b, npairs_2b, f2_2b_flat_clusters, 2);        
    

    /////////////////////////////////////////////
    // Read in 3B clusters -- IN TERMS OF rij **OR** sij - determined by user
    /////////////////////////////////////////////
        
    string f1_3b = f1_idx + ".all-3b-clusters.txt"; 
    string f2_3b = f2_idx + ".all-3b-clusters.txt";
    
    // vector<double> f1_3b_flat_clusters;
    // vector<double> f2_3b_flat_clusters;
    vector<pair<vector<vector<double>>, vector<int>>> f1_3b_flat_clusters;
    vector<pair<vector<vector<double>>, vector<int>>> f2_3b_flat_clusters;
    
    int npairs_3b = 3;
    
    read_flat_clusters(f1_3b, npairs_3b, f1_3b_flat_clusters, 3);
    read_flat_clusters(f2_3b, npairs_3b, f2_3b_flat_clusters, 3);   
  
    /////////////////////////////////////////////
    // Read in 4B clusters -- IN TERMS OF rij **OR** sij - determined by user
    /////////////////////////////////////////////  
    
    string f1_4b = f1_idx + ".all-4b-clusters.txt"; 
    string f2_4b = f2_idx + ".all-4b-clusters.txt";   
    
    vector<pair<vector<vector<double>>, vector<int>>> f1_4b_flat_clusters;
    vector<pair<vector<vector<double>>, vector<int>>> f2_4b_flat_clusters;    
    
    int npairs_4b = 6;
    
    read_flat_clusters(f1_4b, npairs_4b, f1_4b_flat_clusters, 4);
    read_flat_clusters(f2_4b, npairs_4b, f2_4b_flat_clusters, 4);       
 
    /////////////////////////////////////////////
    // Determine the max possible distance between two clusters
    // Added 1.0 to each due to multi-element formulation
    // Removed rij transformation since sij are already normalized
    /////////////////////////////////////////////

    // double maxd_2b = 2.0 + 1.0;
    // double maxd_3b = sqrt( 3.0*pow(2.0,2.0) ) + 1.0;
    // double maxd_4b = sqrt( 6.0*pow(2.0,2.0) ) + 1.0;
    
    // if(my_rank==0)
    //     cout << "Setting maximum histogram values: " << maxd_2b <<  " " << maxd_3b << " " << maxd_4b << endl;
    

    // /////////////////////////////////////////////
    // // generate the cluster distance histogram
    // /////////////////////////////////////////////

    // // set up the histograms
    
    // double binw_2b = maxd_2b/nbin_2b; 
    // double binw_3b = maxd_3b/nbin_3b; 
    // double binw_4b = maxd_4b/nbin_4b; 
    
    // bool same = false; if (f1_idx == f2_idx) same = true; 
    
    // gen_flat_hists(f1_2b_flat_clusters, f2_2b_flat_clusters, npairs_2b, nbin_2b, binw_2b, maxd_2b, f1_idx + "-" + f2_idx + ".2b_clu-" + style + ".hist", rcout_2b, same);
    // gen_flat_hists(f1_3b_flat_clusters, f2_3b_flat_clusters, npairs_3b, nbin_3b, binw_3b, maxd_3b, f1_idx + "-" + f2_idx + ".3b_clu-" + style + ".hist", rcout_3b, same);
    // gen_flat_hists(f1_4b_flat_clusters, f2_4b_flat_clusters, npairs_4b, nbin_4b, binw_4b, maxd_4b, f1_idx + "-" + f2_idx + ".4b_clu-" + style + ".hist", rcout_4b, same);


    // MPI_Finalize();
     
}


