#include <mpi.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm> // For sort

using namespace std;


int nprocs;
int my_rank;

int split_line(string line, vector<string> & items)
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


void read_flat_clusters(string clufile, int npairs_per_cluster, vector<double > & clusters)
{
    ifstream clustream;
    clustream.open(clufile); 

    string                  line;
    vector<string>          line_contents;
    int                     n_contents;
    vector<double>          one_cluster(npairs_per_cluster);

    while (get_next_line(clustream, line))
    {
        n_contents = split_line(line, line_contents);
        
        if (n_contents != npairs_per_cluster)
        {
            cout << "ERROR: Read the wrong number of clusters!" << endl;
            cout << "n_contents: " <<  n_contents << endl;
            cout << "npairs_per_cluster: " <<  npairs_per_cluster << endl;
            exit(0);
        }

        for (int i=0; i<npairs_per_cluster; i++)
            one_cluster[i] = stod(line_contents[i]);
    
        sort(one_cluster.begin(), one_cluster.end());
    
        // Tack one_cluster on the end of clusters

        clusters.insert(clusters.end(), one_cluster.begin(), one_cluster.end());
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
    double                  dist;
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
            dist = 0;
    
            for (int k=0; k<n_cluster_pairs; k++)
                dist += pow(clu1[i*n_cluster_pairs+k] - clu2[j*n_cluster_pairs+k],2.0);

            bin  = get_bin(binw, maxd, sqrt(dist));
            
            if (bin > nbin)
            {
                
                cout << "Rank: " << my_rank << " ERROR: computed bin larger than nbins:" << endl;

                cout << "Rank: " << my_rank << " nbin: " << nbin << endl;
                cout << "Rank: " << my_rank << " bin:  " << bin << endl;
                cout << "Rank: " << my_rank << " binw: " << binw << endl;
                cout << "Rank: " << my_rank << " maxd: " << maxd << endl;
                cout << "Rank: " << my_rank << " dist: " << dist << endl;
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
    
    string style = argv[3]; // "s"; // Calc distances based on rij, not transformed sij
         
    double rcout_2b = 5.0;
    double rcout_3b = 5.0;
    double rcout_4b = 4.5;
    
    int nbin_2b = 100;
    int nbin_3b = 100;
    int nbin_4b = 100;

    /////////////////////////////////////////////
    // Read in 2B clusters -- IN TERMS OF rij **OR** sij - determined by user
    /////////////////////////////////////////////
    
    string f1_2b = f1_idx + ".2b_clu-" + style + ".txt"; 
    string f2_2b = f2_idx + ".2b_clu-" + style + ".txt";

    
    vector<double> f1_2b_flat_clusters;
    vector<double> f2_2b_flat_clusters;
    
    int npairs_2b = 1;
    
    read_flat_clusters(f1_2b, npairs_2b, f1_2b_flat_clusters);
    read_flat_clusters(f2_2b, npairs_2b, f2_2b_flat_clusters);        
    

    /////////////////////////////////////////////
    // Read in 3B clusters -- IN TERMS OF rij **OR** sij - determined by user
    /////////////////////////////////////////////
        
    string f1_3b = f1_idx + ".3b_clu-" + style + ".txt"; 
    string f2_3b = f2_idx + ".3b_clu-" + style + ".txt";
    
    vector<double> f1_3b_flat_clusters;
    vector<double> f2_3b_flat_clusters;
    
    int npairs_3b = 3;
    
    read_flat_clusters(f1_3b, npairs_3b, f1_3b_flat_clusters);
    read_flat_clusters(f2_3b, npairs_3b, f2_3b_flat_clusters);   
  
    /////////////////////////////////////////////
    // Read in 4B clusters -- IN TERMS OF rij **OR** sij - determined by user
    /////////////////////////////////////////////  
    
    string f1_4b = f1_idx + ".4b_clu-" + style + ".txt"; 
    string f2_4b = f2_idx + ".4b_clu-" + style + ".txt";   
    
    vector<double> f1_4b_flat_clusters;
    vector<double> f2_4b_flat_clusters;    
    
    int npairs_4b = 6;
    
    read_flat_clusters(f1_4b, npairs_4b, f1_4b_flat_clusters);
    read_flat_clusters(f2_4b, npairs_4b, f2_4b_flat_clusters);       
 

    /////////////////////////////////////////////
    // Determine the max possible distance between two clusters
    /////////////////////////////////////////////
    
    double maxd_2b = rcout_2b;                      if (style == "s") maxd_2b = 2.0;
    double maxd_3b = sqrt( 3.0*pow(rcout_3b,2.0) ); if (style == "s") maxd_3b = sqrt( 3.0*pow(2.0,2.0) );
    double maxd_4b = sqrt( 6.0*pow(rcout_4b,2.0) ); if (style == "s") maxd_4b = sqrt( 6.0*pow(2.0,2.0) );
    
    if(my_rank==0)
        cout << "Setting maximum histogram values: " << maxd_2b <<  " " << maxd_3b << " " << maxd_4b << endl;
    

    /////////////////////////////////////////////
    // generate the cluster distance histogram
    /////////////////////////////////////////////

    // set up the histograms
    
    double binw_2b = maxd_2b/nbin_2b; 
    double binw_3b = maxd_3b/nbin_3b; 
    double binw_4b = maxd_4b/nbin_4b; 
    
    bool same = false; if (f1_idx == f2_idx) same = true; 
    
    gen_flat_hists(f1_2b_flat_clusters, f2_2b_flat_clusters, npairs_2b, nbin_2b, binw_2b, maxd_2b, f1_idx + "-" + f2_idx + ".2b_clu-" + style + ".hist", rcout_2b, same);
    gen_flat_hists(f1_3b_flat_clusters, f2_3b_flat_clusters, npairs_3b, nbin_3b, binw_3b, maxd_3b, f1_idx + "-" + f2_idx + ".3b_clu-" + style + ".hist", rcout_3b, same);
    gen_flat_hists(f1_4b_flat_clusters, f2_4b_flat_clusters, npairs_4b, nbin_4b, binw_4b, maxd_4b, f1_idx + "-" + f2_idx + ".4b_clu-" + style + ".hist", rcout_4b, same);


    MPI_Finalize();
     
}


