#include <mpi.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm> // For sort
#include <iomanip>
#include <dirent.h>
#include <cstring>

#include "chimesFF.h"    

using namespace std;
using namespace GlobalParams;

int nprocs;
int my_rank;

#include <iomanip>
#include "chimesFF.h"  // Make sure this includes the GlobalParams namespace declaration

void print_global_params() {
    
    // Print 2-body cutoffs
    cout << "\n2-Body Cutoffs:\n";
    for (size_t i = 0; i < rcut_2b_list.size(); ++i) {
        cout << "Pair " << i << ": ";
        cout << "Inner = " << fixed << setprecision(4) << rcut_2b_list[i][0];
        cout << ", Outer = " << rcut_2b_list[i][1] << "\n";
    }

    // Print 3-body cutoffs
    cout << "\n3-Body Cutoffs:\n";
    for (size_t i = 0; i < rcut_3b_list.size(); ++i) {
        cout << "Triplet " << i << ":\n";
        for (size_t j = 0; j < rcut_3b_list[i].size(); ++j) {
            cout << "  Pair " << j << ": ";
            cout << "Inner = " << rcut_3b_list[i][j][0];
            cout << ", Outer = " << rcut_3b_list[i][j][1] << "\n";
        }
    }

    // Print 4-body cutoffs
    cout << "\n4-Body Cutoffs:\n";
    for (size_t i = 0; i < rcut_4b_list.size(); ++i) {
        cout << "Quadruplet " << i << ":\n";
        for (size_t j = 0; j < rcut_4b_list[i].size(); ++j) {
            cout << "  Pair " << j << ": ";
            cout << "Inner = " << rcut_4b_list[i][j][0];
            cout << ", Outer = " << rcut_4b_list[i][j][1] << "\n";
        }
    }

    // Print Morse lambda values
    cout << "\nMorse Lambda Values:\n";
    for (size_t i = 0; i < morse_lambda_list.size(); ++i) {
        cout << "Pair " << i << ": λ = " 
                  << fixed << setprecision(4) 
                  << morse_lambda_list[i] << "\n";
    }
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
void read_flat_clusters(string clufile, int npairs_per_cluster, vector<double > & clusters, vector<int> & atom_types, const int body_cnt) {
    ifstream clustream(clufile);
    if (!clustream.is_open()) {
        cerr << "ERROR: Could not open file " << clufile << endl;
        exit(0);
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
        vector<int> typ_idxs(body_cnt);

        for (int i = 0; i < body_cnt; i++) {
            typ_idxs[i] = stoi(line_contents[npairs_per_cluster + i]);
        }
        atom_types.insert(atom_types.end(), typ_idxs.begin(), typ_idxs.end());

        double cutoff_0, cutoff_00;
        double cutoff_1, cutoff_01;
        double cutoff_2, cutoff_02;
        double cutoff_3, cutoff_03;
        double cutoff_4, cutoff_04;
        double cutoff_5, cutoff_05;
        double cutoff_6, cutoff_06;

        double morse_pair_1, morse_pair_2, morse_pair_3, morse_pair_4, morse_pair_5, morse_pair_6;

        if (body_cnt == 2){
            int pair_idx = atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[1] ];
            cutoff_0 = rcut_2b_list[pair_idx][1];
            cutoff_00 = rcut_2b_list[pair_idx][0];
            morse_pair_1 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[1]]];
            edge_lengths[0] = transform(cutoff_00, cutoff_0, morse_pair_1, stod(line_contents[0]));
        } else if (body_cnt == 3){
            // cout << "Entered body_cnt 3 look" << endl;
            int type_idx =  typ_idxs[0]*atom_typ_cnt*atom_typ_cnt + typ_idxs[1]*atom_typ_cnt + typ_idxs[2];
            int tripidx = atom_int_trip_mapping[type_idx];
            vector<int> & mapped_pair_idx_3b = pair_int_trip_mapping[type_idx];
            // Get cutoffs
            // cout << "Getting cutoffs" << endl;
            // cout << type_idx << endl;
            // cout << mapped_pair_idx_3b[0] << endl;
            cutoff_0  = rcut_3b_list[ tripidx ][1][mapped_pair_idx_3b[0]]; // outer cutoff
            // cout << "got cutoff" << endl;
            cutoff_00 = rcut_3b_list[ tripidx ][0][mapped_pair_idx_3b[0]]; // inner cutoff
            cutoff_1  = rcut_3b_list[ tripidx ][1][mapped_pair_idx_3b[1]];
            cutoff_01 = rcut_3b_list[ tripidx ][0][mapped_pair_idx_3b[1]]; 
            cutoff_2  = rcut_3b_list[ tripidx ][1][mapped_pair_idx_3b[2]];
            cutoff_02 = rcut_3b_list[ tripidx ][0][mapped_pair_idx_3b[2]];
            // Get morse variables
            // cout << "Getting morse" << endl;
            morse_pair_1 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[1]]];
            morse_pair_2 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[2]]];
            morse_pair_3 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[1]*atom_typ_cnt + typ_idxs[2]]];
            // Assign edge lengths
            // cout << "Assigning edges" << endl;
            edge_lengths[0] = transform(cutoff_00, cutoff_0, morse_pair_1, stod(line_contents[0]));
            edge_lengths[1] = transform(cutoff_01, cutoff_1, morse_pair_2, stod(line_contents[1]));
            edge_lengths[2] = transform(cutoff_02, cutoff_2, morse_pair_3, stod(line_contents[2]));
        } else {
            int idx = typ_idxs[0]*atom_typ_cnt*atom_typ_cnt*atom_typ_cnt + typ_idxs[1]*atom_typ_cnt*atom_typ_cnt + typ_idxs[2]*atom_typ_cnt + typ_idxs[3] ;
            int quadidx = atom_int_quad_mapping[idx] ;
            vector<int> & mapped_pair_idx_4b = pair_int_quad_mapping[idx] ;
            // Get cutoffs
            cutoff_0  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[0]];
            cutoff_00 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[0]];
            cutoff_1  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[1]];
            cutoff_01 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[1]]; 
            cutoff_2  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[2]];
            cutoff_02 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[2]];
            cutoff_3  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[3]];
            cutoff_03 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[3]]; 
            cutoff_4  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[4]];
            cutoff_04 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[4]];
            cutoff_5  = rcut_4b_list[ quadidx ][1][mapped_pair_idx_4b[5]];
            cutoff_05 = rcut_4b_list[ quadidx ][0][mapped_pair_idx_4b[5]];
            // Get morse variables
            morse_pair_1 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[1]]];
            morse_pair_2 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[2]]];
            morse_pair_3 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[0]*atom_typ_cnt + typ_idxs[3]]];            
            morse_pair_4 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[1]*atom_typ_cnt + typ_idxs[2]]];
            morse_pair_5 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[1]*atom_typ_cnt + typ_idxs[3]]];
            morse_pair_6 = morse_lambda_list[atom_int_pair_mapping[ typ_idxs[2]*atom_typ_cnt + typ_idxs[3]]];
            // Assign edge lengths
            edge_lengths[0] = transform(cutoff_00, cutoff_0, morse_pair_1, stod(line_contents[0]));
            edge_lengths[1] = transform(cutoff_01, cutoff_1, morse_pair_2, stod(line_contents[1]));
            edge_lengths[2] = transform(cutoff_02, cutoff_2, morse_pair_3, stod(line_contents[2]));
            edge_lengths[3] = transform(cutoff_03, cutoff_3, morse_pair_4, stod(line_contents[3]));
            edge_lengths[4] = transform(cutoff_04, cutoff_4, morse_pair_5, stod(line_contents[4]));
            edge_lengths[5] = transform(cutoff_05, cutoff_5, morse_pair_6, stod(line_contents[5]));
        }
        sort(edge_lengths.begin(), edge_lengths.end());
        clusters.insert(clusters.end(), edge_lengths.begin(), edge_lengths.end());
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

double compute_composition_weight_2b(vector<double> edges, vector<int> atom_types){
    double normalized_dist = (atom_types[0] + atom_types[1])/(2.0*atom_typ_cnt);
    return normalized_dist;
}

double compute_composition_weight_3b(vector<double> edges, vector<int> atom_types){

    vector<pair<int, double>> sum_edge_lengths(3);
    sum_edge_lengths[0] = {0, edges[0] + edges[1]};
    sum_edge_lengths[1] = {1, edges[0] + edges[2]};
    sum_edge_lengths[2] = {2, edges[1] + edges[2]};

    // Sort the vector by the summed values in ascending order
    sort(sum_edge_lengths.begin(), sum_edge_lengths.end(),
        [](const pair<int, double>& a, const pair<int, double>& b) {
            return a.second < b.second;
        });

    double total_dist =   1.0*atom_types[sum_edge_lengths[0].first] 
                        + 2.0*atom_types[sum_edge_lengths[1].first] 
                        + 3.0*atom_types[sum_edge_lengths[2].first];
    
    double normalized_dist = total_dist / (6.0*atom_typ_cnt);
    return normalized_dist;
}

double compute_composition_weight_4b(vector<double> edges, vector<int> atom_types){

    vector<pair<int, double>> sum_edge_lengths(4);
    sum_edge_lengths[0] = {0, edges[0] + edges[1] + edges[2]};
    sum_edge_lengths[1] = {1, edges[0] + edges[3] + edges[4]};
    sum_edge_lengths[2] = {2, edges[1] + edges[3] + edges[5]};
    sum_edge_lengths[3] = {3, edges[2] + edges[4] + edges[5]};

    // Sort the vector by the summed values in ascending order
    sort(sum_edge_lengths.begin(), sum_edge_lengths.end(),
        [](const pair<int, double>& a, const pair<int, double>& b) {
            return a.second < b.second;
        });

    double total_dist =   1.0*atom_types[sum_edge_lengths[0].first] 
                        + 2.0*atom_types[sum_edge_lengths[1].first] 
                        + 3.0*atom_types[sum_edge_lengths[2].first] 
                        + 4.0*atom_types[sum_edge_lengths[3].first]; 
    
    double normalized_dist = total_dist / (10.0*atom_typ_cnt);
    return normalized_dist;
}

void gen_flat_hists(vector<double > & clu1, vector<double > & clu2, vector<int> & clu1_atm_types, vector<int> & clu2_atm_types, int n_cluster_pairs, int nbin, double binw, double maxd, string histfile, bool same = false, int body_cnt = 2)
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
    double                  dist_struct;
    double                  d_comp1;
    double                  d_comp2;
    
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
        // cout << "Entering first loop generation" << endl;
        for (int j=looptwo_start; j<clu2.size()/n_cluster_pairs; j++)
        {
            dist_struct = 0;
            vector<double> edge_length_1(n_cluster_pairs);
            vector<double> edge_length_2(n_cluster_pairs);
            vector<int>    atom_list_1(body_cnt);
            vector<int>    atom_list_2(body_cnt);
    
            for (int k=0; k<n_cluster_pairs; k++){   
                // Distance calculation between two clusters             
                edge_length_1.push_back(clu1[i*n_cluster_pairs+k]);
                edge_length_2.push_back(clu2[j*n_cluster_pairs+k]);
                dist_struct += pow(clu1[i*n_cluster_pairs+k] - clu2[j*n_cluster_pairs+k],2.0);
                }
            
            for (int l=0; l<body_cnt; l++){
                atom_list_1.push_back(clu1_atm_types[i*body_cnt+l]);
                atom_list_2.push_back(clu2_atm_types[j*body_cnt+l]);
            }

            if (body_cnt == 2){
                d_comp1 = compute_composition_weight_2b(edge_length_1, atom_list_1);
                d_comp2 = compute_composition_weight_2b(edge_length_2, atom_list_2);
            } else if (body_cnt == 3){
                d_comp1 = compute_composition_weight_3b(edge_length_1, atom_list_1);
                d_comp2 = compute_composition_weight_3b(edge_length_2, atom_list_2);
            } else if (body_cnt == 4){
                d_comp1 = compute_composition_weight_4b(edge_length_1, atom_list_1);
                d_comp2 = compute_composition_weight_4b(edge_length_2, atom_list_2);
            } else {
                cout << "Improper body count" << endl;
                exit(1);
            }
            total_dist = sqrt(dist_struct) + abs(d_comp1 - d_comp2);
            bin  = get_bin(binw, maxd, total_dist);
            
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
    nprocs = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    cout << "Nprocs = " << nprocs << endl;

    vector<string> all_files;
    int num_files = 0;

    // File discovery and sorting (rank 0 only)
    if (my_rank == 0) {
        DIR *dir = opendir(".");
        struct dirent *entry;
        while ((entry = readdir(dir)) != nullptr) {
            string filename = entry->d_name;
            size_t num_len = 0;
            while (num_len < filename.size() && isdigit(filename[num_len]))
                num_len++;
            
            if (num_len > 0 && filename.substr(num_len) == ".all-2b-clusters.txt")
                all_files.push_back(filename.substr(0, num_len));
        }
        closedir(dir);

        sort(all_files.begin(), all_files.end(), [](const string &a, const string &b) {
            return stoi(a) < stoi(b);
        });

        num_files = all_files.size();
        cout << "Found " << num_files << " files" << endl;
    }

    // Broadcast file list
    MPI_Bcast(&num_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
    all_files.resize(num_files);

    for (int i = 0; i < num_files; ++i) {
        int str_len = 0;
        char buffer[256] = {0};
        
        if (my_rank == 0) {
            str_len = all_files[i].size();
            strncpy(buffer, all_files[i].c_str(), 255);
        }
        
        MPI_Bcast(&str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(buffer, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        all_files[i] = string(buffer, str_len);
    }

    // Initialize ChIMES calculator
    chimesFF ff;
    ff.init(my_rank);
    ff.read_parameters("params.txt.reduced");
    ff.build_pair_int_trip_map();
    ff.build_pair_int_quad_map();
    
    print_global_params();

    // Process individual files
    for (const auto &file_idx : all_files) {
        string f1_idx = file_idx;
        string f2_idx = file_idx;  // Same file for both indices

        if (my_rank == 0)
            cout << "\nProcessing file: " << file_idx << endl;

        // File paths - using same index for both
        string f1_2b = f1_idx + ".all-2b-clusters.txt";
        string f2_2b = f2_idx + ".all-2b-clusters.txt";

        // ... rest of your processing code ...

        // Histogram parameters
        const int nbin_2b = 100, nbin_3b = 100, nbin_4b = 100;
        const double maxd_2b = 3.0, 
                     maxd_3b = sqrt(12.0) + 1.0,
                     maxd_4b = sqrt(24.0) + 1.0;

        double binw_2b = maxd_2b/nbin_2b; 
        double binw_3b = maxd_3b/nbin_3b; 
        double binw_4b = maxd_4b/nbin_4b; 

        vector<double> f1_2b_flat_clusters, f2_2b_flat_clusters;
        vector<int> f1_2b_atom_types, f2_2b_atom_types;
        int npairs_2b = 1;
        read_flat_clusters(f1_2b, npairs_2b, f1_2b_flat_clusters, f1_2b_atom_types, 2);
        read_flat_clusters(f2_2b, npairs_2b, f2_2b_flat_clusters, f2_2b_atom_types, 2);

        // Process 3B clusters
        string f1_3b = f1_idx + ".all-3b-clusters.txt"; 
        string f2_3b = f2_idx + ".all-3b-clusters.txt";
        vector<double> f1_3b_flat_clusters, f2_3b_flat_clusters;
        vector<int> f1_3b_atom_types, f2_3b_atom_types;
        int npairs_3b = 3;
        read_flat_clusters(f1_3b, npairs_3b, f1_3b_flat_clusters, f1_3b_atom_types, 3);
        read_flat_clusters(f2_3b, npairs_3b, f2_3b_flat_clusters, f2_3b_atom_types, 3);

        // Process 4B clusters
        string f1_4b = f1_idx + ".all-4b-clusters.txt"; 
        string f2_4b = f2_idx + ".all-4b-clusters.txt";   
        vector<double> f1_4b_flat_clusters, f2_4b_flat_clusters;   
        vector<int> f1_4b_atom_types, f2_4b_atom_types; 
        int npairs_4b = 6;
        read_flat_clusters(f1_4b, npairs_4b, f1_4b_flat_clusters, f1_4b_atom_types, 4);
        read_flat_clusters(f2_4b, npairs_4b, f2_4b_flat_clusters, f2_4b_atom_types, 4);
        
        if (my_rank == 0) {
            cout << "Setting maximum histogram values: " 
                 << maxd_2b << " " << maxd_3b << " " << maxd_4b << endl;
        }

        bool same = (f1_idx == f2_idx);

        gen_flat_hists(f1_2b_flat_clusters, f2_2b_flat_clusters, f1_2b_atom_types, 
                      f2_2b_atom_types, npairs_2b, nbin_2b, binw_2b, maxd_2b, 
                      f1_idx + "-" + f2_idx + ".2b_clu-s.hist", same, 2);
        gen_flat_hists(f1_3b_flat_clusters, f2_3b_flat_clusters, f1_3b_atom_types, 
                      f2_3b_atom_types, npairs_3b, nbin_3b, binw_3b, maxd_3b, 
                      f1_idx + "-" + f2_idx + ".3b_clu-s.hist", same, 3);
        gen_flat_hists(f1_4b_flat_clusters, f2_4b_flat_clusters, f1_4b_atom_types, 
                      f2_4b_atom_types, npairs_4b, nbin_4b, binw_4b, maxd_4b, 
                      f1_idx + "-" + f2_idx + ".4b_clu-s.hist", same, 4);
    }

    MPI_Finalize();
    return 0;
}