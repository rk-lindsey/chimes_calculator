/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
	Contributing Author:  Rebecca K. Lindsey (2020) 
*/

#include<vector>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<map>

using namespace std;

#include "chimesFF.h"    

static const double pi      = 3.14159265359;        


template <typename T>
int get_index(const vector<T>  & vec, const T  & element)
{
    auto it = find(vec.begin(), vec.end(), element);
 
    if (it != vec.end())
    {
        return distance(vec.begin(), it);
    }
    else
    {
        cout << "chimesFF: " << "ERROR: Could not find element in vector" << endl;
        exit(0);
    }
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

int chimesFF::get_proper_pair(string ty1, string ty2)
{

    for(int i=0; i<pair_params_atm_chem_1.size(); i++)
    {
        if (ty1 == pair_params_atm_chem_1[i])
            if (ty2 == pair_params_atm_chem_2[i])
                return i;
        
        if (ty2 == pair_params_atm_chem_1[i])
            if (ty1 == pair_params_atm_chem_2[i])
                return i;
    }
            
    cout << "chimesFF: " << "ERROR: No proper pair name found for atom types" << ty1 << ", " << ty2 << endl;
    exit(0);
}

chimesFF::chimesFF()
{
    natmtyps = 0;
    penalty_params.resize(2);
    
    // Set defaults
    
    fcut_type = "CUBIC";
    
    penalty_params[0] = 1.0E4;
    penalty_params[1] = 0.01;
}
chimesFF::~chimesFF(){}

void chimesFF::init(int mpi_rank)
{
    rank = mpi_rank;
    print_pretty_stuff();
}

void chimesFF::print_pretty_stuff()
{
    if (rank == 0)
    {
        cout << "chimesFF: " <<  endl;
        cout << "chimesFF: " << "01000011011010001001001010011010100010101010011 0100010101101110110011101101001011011101100101  " << endl;
        cout << "chimesFF: " <<  endl;
        cout << "chimesFF: " << "      _____  _      _____  __  __  ______   _____   ______                _                      " << endl;
        cout << "chimesFF: " << "     / ____|| |    |_   _||  \\/  ||  ____| / ____| |  ____|              (_)                    " << endl;
        cout << "chimesFF: " << "    | |     | |__    | |  | \\  / || |__   | (___   | |__    _ __    __ _  _  _ __    ___        " << endl;
        cout << "chimesFF: " << "    | |     | '_ \\   | |  | |\\/| ||  __|   \\___ \\  |  __|  | '_ \\  / _` || || '_ \\  / _ \\ " << endl;
        cout << "chimesFF: " << "    | |____ | | | | _| |_ | |  | || |____  ____) | | |____ | | | || (_| || || | | ||  __/        " << endl;
        cout << "chimesFF: " << "     \\_____||_| |_||_____||_|  |_||______||_____/  |______||_| |_| \\__, ||_||_| |_| \\___|     " << endl;
        cout << "chimesFF: " << "                                                                    __/ |                        " << endl;
        cout << "chimesFF: " << "                                                                   |___/                         " << endl;  
        cout << "chimesFF: " << endl;
        cout << "chimesFF: " << "                     Copyright (C) 2020 R.K. Lindsey, L.E. Fried, N. Goldman                     " << endl;    
        cout << "chimesFF: " << endl;
        cout << "chimesFF: " << "01000011011010001001001010011010100010101010011 0100010101101110110011101101001011011101100101   " << endl;
        cout << "chimesFF: " << endl;
    }
      
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
        if (rank == 0)
            cout << "chimesFF: " << "Error reading line" << line << endl;
        exit(0);
    }

    return line;
}

void chimesFF::read_parameters(string paramfile)
{
    // Open the parameter file, run sanity checks
    
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
    
    // Declare parsing variables
    
    
    bool           found_end = false;
    string         line;
    string         tmp_str;
    vector<string> tmp_str_items;
    int            tmp_no_items;
    int            tmp_int;
    int            no_pairs;
    
    // Check that this is actually a chebyshev parameter set

    while (!found_end)
    {
        line = get_next_line(param_file);

           // Break out of loop

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
            tmp_no_items = split_line(line, tmp_str_items);
            
            if (tmp_no_items < 3)
            {    
                if (rank == 0)
                    cout << "chimesFF: " << "ERROR: \"PAIRTYP: CHEBYSHEV\" line must at least contain the 2-body order" << endl;
                exit(0);
            }
            
            poly_orders.push_back(stoi(tmp_str_items[2]));
            
            if (tmp_no_items >= 4)
                poly_orders.push_back(stoi(tmp_str_items[3]));

            if (tmp_no_items >= 5)
                poly_orders.push_back(stoi(tmp_str_items[4]));    
            
            while (poly_orders.size() < 3)
                poly_orders.push_back(0);
            
            if (rank == 0)
            {
                cout << "chimesFF: " << "Using respective 2, 3, and 4-body orders of: " << poly_orders[0] << " " << poly_orders[1] << " " << poly_orders[2] << endl;
            
                cout << "chimesFF: " << "Note: Ignoring polynomial domain; assuming [-1,1]" << endl;    
            }
            
            break;    
        }
    }
    
    // If we've made it to here, then this should contain Chebyshev params. Rewind and start looking for general information
        
    param_file.seekg(0);
    
    found_end = false;
    
    while (!found_end)
    {
        line = get_next_line(param_file);
        
           if(line.find("ENDFILE") != string::npos)
            break;        
    
        if(line.find("ATOM TYPES:") != string::npos)
        {
            tmp_no_items = split_line(line, tmp_str_items);
        
            natmtyps = stoi(tmp_str_items[2]);
        
            if (rank == 0)
                cout << "chimesFF: " << "Will consider " << natmtyps << " atom types:" << endl;
                
            energy_offsets.resize(natmtyps);
            
            for(int i=0; i<natmtyps; i++)
                energy_offsets[i] = 0.0;
        }
        
        if(line.find("# TYPEIDX #") != string::npos)
        {
            atmtyps.resize(natmtyps);
            for (int i=0; i<natmtyps; i++)
            {
                line = get_next_line(param_file);
                split_line(line, tmp_str_items);
                atmtyps[i] = tmp_str_items[1];
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << i << " " << atmtyps[i] << endl;
            }
            
        }
            
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
            chimes_2b_cutoff      .resize(no_pairs);
            morse_var             .resize(no_pairs);
            
            ncoeffs_2b            .resize(no_pairs);
            chimes_2b_pows        .resize(no_pairs);
            chimes_2b_params      .resize(no_pairs);
            chimes_2b_cutoff      .resize(no_pairs);

            string tmp_xform_style;
            
            for (int i=0; i<no_pairs; i++)
            {
                line = get_next_line(param_file);
                
                tmp_no_items = split_line(line, tmp_str_items);
                
                pair_params_atm_chem_1[i] = tmp_str_items[1];
                pair_params_atm_chem_2[i] = tmp_str_items[2];
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << i << " " << pair_params_atm_chem_1[i] << " " << pair_params_atm_chem_2[i]<< endl;
                
                chimes_2b_cutoff[i].push_back(stod(tmp_str_items[3])); // Inner cutoff    
                chimes_2b_cutoff[i].push_back(stod(tmp_str_items[4])); // Outer cutoff
                
                if (i==0)
                {
                    tmp_xform_style = tmp_str_items[6];
                }
                else
                {
                    if ( tmp_str_items[6] != tmp_xform_style)    
                    {
                        if (rank == 0)
                            cout << "chimesFF: " << "Distance transfomration style must be the same for all pair types" << endl;
                            
                        exit(0);
                    }
                }
                
                if (tmp_no_items >= 8)
                    morse_var[i] = stod(tmp_str_items[7]);
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
            
        if(line.find("FCUT TYPE:") != string::npos)
        {
            tmp_no_items = split_line(line, tmp_str_items);
        
            fcut_type = tmp_str_items[2];
            
            if (rank == 0)
                cout << "chimesFF: " << "Will use cutoff style " << fcut_type;
            
            if (fcut_type == "TERSOFF")
            {
                fcut_var = stod(tmp_str_items[3]);
                
                if (rank == 0)
                    cout << " " << fcut_var << endl;
            }
            else
                if (rank == 0)
                    cout << endl;
        }
        
        if(line.find("PAIR CHEBYSHEV PENALTY DIST:") != string::npos)
        {    
            tmp_no_items = split_line(line, tmp_str_items);
            
            penalty_params[0] = stod(tmp_str_items[4]);
            
            if (rank == 0)
                cout << "chimesFF: " << "Will use penalty distance: " << penalty_params[0] << endl;
        }
        
        if(line.find("PAIR CHEBYSHEV PENALTY SCALING:") != string::npos)
        {    
            tmp_no_items = split_line(line, tmp_str_items);
            
            penalty_params[1] = stod(tmp_str_items[4]);
            
            if (rank == 0)
                cout << "chimesFF: " << "Will use penalty scaling: " << penalty_params[1] << endl;
        }
        
        if(line.find("NO ENERGY OFFSETS:") != string::npos)
        {
            int tmp_no = split_line(line, tmp_str_items);
                        
            if(stoi(tmp_str_items[tmp_no-1]) != natmtyps)
            {
                cout << "chimesFF: " << "ERROR: Number of energy offsets do not match number of atom types" << endl;
                exit(0);
            }

            // Expects atom offsets in the same order as atom types were provided originally
            
            if (rank == 0)
                cout << "chimesFF: " << "Will use single atom energy offsets: "<< endl;
            
            int tmp_idx;
            
            for (int i=0; i<natmtyps; i++)
            {
                line = get_next_line(param_file);
                split_line(line, tmp_str_items);
                tmp_idx = stoi(tmp_str_items[2]);
                
                energy_offsets[tmp_idx-1] = stod(tmp_str_items[3]);
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << tmp_idx << " " << atmtyps[tmp_idx-1] << " " << energy_offsets[tmp_idx-1] << endl;
            }
            
        }                
    }
    
    // Rewind and read the 2-body Chebyshev pair parameters
    
    param_file.seekg(0);
    
    found_end = false;
    
    while (!found_end)
    {
        line = get_next_line(param_file);

           if(line.find("ENDFILE") != string::npos)
            break;            
        
        if(line.find("PAIRTYPE PARAMS:") != string::npos)
        {
            tmp_no_items = split_line(line, tmp_str_items);
            
            tmp_int = stoi(tmp_str_items[2]);
            
            if (rank == 0)
                cout << "chimesFF: " << "Read 2B parameters for pair: " << tmp_int << " " << tmp_str_items[3] << " " << tmp_str_items[4] << endl;
            
            line = get_next_line(param_file);
            
            split_line(line, tmp_str_items); // Empty line
            
            ncoeffs_2b[tmp_int] = poly_orders[0];
            
            for(int i=0; i<poly_orders[0]; i++)
            {
                line = get_next_line(param_file);
                split_line(line, tmp_str_items);
                
                chimes_2b_pows  [tmp_int].push_back(stoi(tmp_str_items[0]));                
                chimes_2b_params[tmp_int].push_back(stod(tmp_str_items[1]));
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << chimes_2b_pows[tmp_int][i] << " " << chimes_2b_params[tmp_int][i] << endl;
            }
        }
        
        if(line.find("PAIRMAPS:") != string::npos)
        {
            // Read the slow map and build the fast map
            
            tmp_no_items = split_line(line, tmp_str_items);
            
            n_pair_maps = stoi(tmp_str_items[1]);
            
            atom_typ_pair_map.resize(n_pair_maps);
            atom_idx_pair_map.resize(n_pair_maps);
            
            atom_int_prpr_map.resize(n_pair_maps);
            
            if (rank == 0)
                cout << "chimesFF: " << "Built the following 2-body pair \"slow\" map:" << endl;
            
            for(int i=0; i<n_pair_maps; i++)
            {
                line = get_next_line(param_file);
                split_line(line, tmp_str_items);
                
                atom_idx_pair_map[i] = stoi(tmp_str_items[0]);
                atom_typ_pair_map[i] =      tmp_str_items[1];
                
                if (rank == 0)
                    cout << "chimesFF: " << "\t" << atom_idx_pair_map[i] << " " << atom_typ_pair_map[i] << "(i: " << i << ")" << endl;

            }

            if (rank == 0)
                cout << "chimesFF: " << "Built the following 2-body pair \"fast\" map:" << endl;
            
            atom_int_pair_map.resize((natmtyps-1)*natmtyps + (natmtyps-1) + 1); // Maximum possible pair value + a small buffer
            

            for(int i=0; i<natmtyps; i++)
            {
                for (int j=0; j<natmtyps; j++)
                {
                    // Get the pair type name for the set of atoms
                    
                    tmp_str = atmtyps[i] + atmtyps[j];

                    tmp_int = get_index(atom_typ_pair_map, tmp_str);
                    
                    atom_int_pair_map[ i*natmtyps + j ] = atom_idx_pair_map[tmp_int];
                    

                    tmp_int = get_proper_pair(atmtyps[i],atmtyps[j]);
                    
                    atom_int_prpr_map [ i*natmtyps + j ] = pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];

                    
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << tmp_str << ": " << i*natmtyps + j << " " << atom_int_pair_map[ i*natmtyps + j ] << endl;

                }
            }                        
        }
    }
    
    
    
    // Rewind and read the 3-body Chebyshev pair parameters
    
    if (poly_orders[1] > 0)
    {
        int ntrips;
        int tmp_idx;
        
        // Read parameters
        
        param_file.seekg(0);
        
        found_end = false;
    
        while (!found_end)
        {
            line = get_next_line(param_file);
        
               if(line.find("ENDFILE") != string::npos)
                break;    
            
            if(line.find("ATOM PAIR TRIPLETS:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                ntrips = stoi(tmp_str_items[3]);
                                
                ncoeffs_3b      .resize(ntrips);
                chimes_3b_powers.resize(ntrips);                
                chimes_3b_params.resize(ntrips);
                chimes_3b_cutoff.resize(ntrips);
    
                
                trip_params_atm_chems.resize(ntrips);                
                trip_params_pair_typs.resize(ntrips);
            }
            
            if(line.find("TRIPLETTYPE PARAMS:") != string::npos)
            {
                vector<int> tmp_int_vec(3);

                line = get_next_line(param_file);
                
                split_line(line, tmp_str_items);
                
                tmp_int = stoi(tmp_str_items[1]);
                
                trip_params_atm_chems[tmp_int].push_back(tmp_str_items[3]);
                trip_params_atm_chems[tmp_int].push_back(tmp_str_items[4]);
                trip_params_atm_chems[tmp_int].push_back(tmp_str_items[5]);

                if (rank == 0)
                    cout << "chimesFF: " << "Read 3B parameters for triplet: " << tmp_int << " " << trip_params_atm_chems[tmp_int][0] << " " << trip_params_atm_chems[tmp_int][1] << " " << trip_params_atm_chems[tmp_int][2] << endl;
                
                line = get_next_line(param_file);
                
                split_line(line, tmp_str_items);
            
                trip_params_pair_typs[tmp_int].push_back(tmp_str_items[1]);
                trip_params_pair_typs[tmp_int].push_back(tmp_str_items[2]);
                trip_params_pair_typs[tmp_int].push_back(tmp_str_items[3]);
            
                ncoeffs_3b[tmp_int] = stoi(tmp_str_items[7]);    

                get_next_line(param_file);
                get_next_line(param_file);
            
                for(int i=0; i<ncoeffs_3b[tmp_int]; i++)
                {
                    line = get_next_line(param_file);
                    split_line(line, tmp_str_items);
                    
                    tmp_int_vec[0] = stoi(tmp_str_items[1]);
                    tmp_int_vec[1] = stoi(tmp_str_items[2]);
                    tmp_int_vec[2] = stoi(tmp_str_items[3]);
                    
                    chimes_3b_powers[tmp_int].push_back(tmp_int_vec);                    
                    chimes_3b_params[tmp_int].push_back(stod(tmp_str_items[6]));
                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << chimes_3b_powers[tmp_int][i][0] << " " << chimes_3b_powers[tmp_int][i][1] << " " << chimes_3b_powers[tmp_int][i][2] << " " << chimes_3b_params[tmp_int][i] << endl;
                }
            }    
            
            if(line.find("TRIPMAPS:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                n_trip_maps = stoi(tmp_str_items[1]);
                
                atom_idx_trip_map.resize(n_trip_maps);
                atom_typ_trip_map.resize(n_trip_maps);
                
                if (rank == 0)                
                    cout << "chimesFF: " << "Built the following 3-body pair \"slow\" map:" << endl;
            
                for(int i=0; i<n_trip_maps; i++)
                {
                    line = get_next_line(param_file);
                    split_line(line, tmp_str_items);
                
                    atom_idx_trip_map[i] = stoi(tmp_str_items[0]);
                    atom_typ_trip_map[i] =      tmp_str_items[1];
                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << atom_idx_trip_map[i] << " " << atom_typ_trip_map[i] << endl;
                }        
                
                if (rank == 0)
                    cout << "chimesFF: " << "Built the following 3-body pair \"fast\" map:" << endl;

                atom_int_trip_map.resize(natmtyps*natmtyps*natmtyps);
                
                for(int i=0; i<natmtyps; i++)
                {
                    for (int j=0; j<natmtyps; j++)
                    {
                        for(int k=0; k<natmtyps; k++)
                        {
                            // Get the trip type name for the set of atoms
                            
                            tmp_str = "";

                            tmp_int  = get_proper_pair(atmtyps[i], atmtyps[j]);
                            tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    
                
                            tmp_int  = get_proper_pair(atmtyps[i], atmtyps[k]);
                            tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    
                
                            tmp_int  = get_proper_pair(atmtyps[j], atmtyps[k]);
                            tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];            
                                            
                            tmp_int = get_index(atom_typ_trip_map, tmp_str);

                            tmp_idx = i*natmtyps*natmtyps + j*natmtyps + k;

                            atom_int_trip_map[ tmp_idx ] = atom_idx_trip_map[tmp_int];
                                                        
                            if (rank == 0)
                                cout << "chimesFF: " << "\t" << tmp_idx << " " << atom_int_trip_map[ tmp_idx  ]  << endl;
                        }
                    }
                }
            }            
        }
        
        // Set up cutoffs ... First set to match 2-body, then read special if they exist
        
        int atmtyp_1,  atmtyp_2,  atmtyp_3;
        int pairtyp_1, pairtyp_2, pairtyp_3;

        for(int i=0; i<ntrips; i++) 
        {
            // Figure out the atom type index for each atom in the triplet type 
                        
            atmtyp_1 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), trip_params_atm_chems[i][0]));    
            atmtyp_2 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), trip_params_atm_chems[i][1]));    
            atmtyp_3 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), trip_params_atm_chems[i][2]));    
                        
            // Figure out the corresponding 2-body pair type
            
            pairtyp_1 = atom_int_pair_map[ atmtyp_1*natmtyps + atmtyp_2 ];
            pairtyp_2 = atom_int_pair_map[ atmtyp_1*natmtyps + atmtyp_3 ];
            pairtyp_3 = atom_int_pair_map[ atmtyp_2*natmtyps + atmtyp_3 ];
    
            // Set the default inner/outer cutoffs to the corresponding 2-body value

            chimes_3b_cutoff[i].resize(2);

            chimes_3b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_1][0]);
            chimes_3b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_2][0]);
            chimes_3b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_3][0]);
            
            chimes_3b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_1][1]);
            chimes_3b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_2][1]);
            chimes_3b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_3][1]);                    
        }
        
        param_file.seekg(0);
        
        int    nentries;
        double cutval;
        
        found_end = false;
        
        while (!found_end)
        {
            line = get_next_line(param_file);
        
               if(line.find("ENDFILE") != string::npos)
                break;                
            
            if(line.find("SPECIAL 3B S_MAXIM:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                if (rank == 0)
                    cout << "chimesFF: " << "Set the following special 3-body outer cutoffs: " << endl;
                
                if(tmp_str_items[3] == "ALL")
                {
                    cutval = stod(tmp_str_items[4]);
                                        
                    for(int i=0; i<ntrips; i++)
                    {
                        chimes_3b_cutoff[i][1][0] = cutval;
                        chimes_3b_cutoff[i][1][1] = cutval;
                        chimes_3b_cutoff[i][1][2] = cutval;                    
                    }
                }
                else
                {
                    nentries = stoi(tmp_str_items[4]);
                    
                    vector<string> pair_name(3);
                    vector<double> cutoffval(3);

    
                    for(int i=0; i<nentries; i++)
                    {
                        line = get_next_line(param_file);
                        
                        split_line(line, tmp_str_items);
                        
                        tmp_int = atom_idx_trip_map[distance(atom_typ_trip_map.begin(), find(atom_typ_trip_map.begin(), atom_typ_trip_map.end(), tmp_str_items[0]))];

                        pair_name[0] = tmp_str_items[1];
                        pair_name[1] = tmp_str_items[2];
                        pair_name[2] = tmp_str_items[3];
                        
                        cutoffval[0] = stod(tmp_str_items[4]);
                        cutoffval[1] = stod(tmp_str_items[5]);
                        cutoffval[2] = stod(tmp_str_items[6]);
                        
                        vector<bool>   disqualified(3,false);
                        
                        chimes_3b_cutoff[tmp_int][1][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[0], disqualified) ] = cutoffval[0];
                        chimes_3b_cutoff[tmp_int][1][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[1], disqualified) ] = cutoffval[1];
                        chimes_3b_cutoff[tmp_int][1][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[2], disqualified) ] = cutoffval[2];
                                        
                    }
                }
                
                for(int i=0; i<ntrips; i++)
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << i << " " << chimes_3b_cutoff[i][1][0] << " " << chimes_3b_cutoff[i][1][1] << " " << chimes_3b_cutoff[i][1][2] << endl;
                
            }

            if(line.find("SPECIAL 3B S_MINIM:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                if (rank == 0)
                    cout << "chimesFF: " << "Set the following special 3-body inner cutoffs: " << endl;
                
                if(tmp_str_items[3] == "ALL")
                {
                    cutval = stod(tmp_str_items[4]);
                    
                    for(int i=0; i<ntrips; i++)
                    {
                        chimes_3b_cutoff[i][0][0] = cutval;
                        chimes_3b_cutoff[i][0][1] = cutval;
                        chimes_3b_cutoff[i][0][2] = cutval;
                        
                    }
                }
                else
                {
                    nentries = stoi(tmp_str_items[4]);
                    
                    vector<string> pair_name(3);
                    vector<double> cutoffval(3);


                    for(int i=0; i<nentries; i++)
                    {
                        line = get_next_line(param_file);
                        
                        split_line(line, tmp_str_items);
                        
                        tmp_int = atom_idx_trip_map[distance(atom_typ_trip_map.begin(), find(atom_typ_trip_map.begin(), atom_typ_trip_map.end(), tmp_str_items[0]))];
                        
                        pair_name[0] = tmp_str_items[1];
                        pair_name[1] = tmp_str_items[2];
                        pair_name[2] = tmp_str_items[3];
                        
                        cutoffval[0] = stod(tmp_str_items[4]);
                        cutoffval[1] = stod(tmp_str_items[5]);
                        cutoffval[2] = stod(tmp_str_items[6]);
                        
                        vector<bool>   disqualified(3,false);
                        
                        chimes_3b_cutoff[tmp_int][0][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[0], disqualified) ] = cutoffval[0];
                        chimes_3b_cutoff[tmp_int][0][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[1], disqualified) ] = cutoffval[1];
                        chimes_3b_cutoff[tmp_int][0][ get_index_if(trip_params_pair_typs[tmp_int], pair_name[2], disqualified) ] = cutoffval[2];
                
                    }
                }
                
                for(int i=0; i<ntrips; i++)
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << i << " " << chimes_3b_cutoff[i][0][0] << " " << chimes_3b_cutoff[i][0][1] << " " << chimes_3b_cutoff[i][0][2] << endl;
            }            
        }    
    }
    
    // Rewind and read the 4-body Chebyshev pair parameters
    
    if (poly_orders[2] > 0)
    {
        int nquads;
        int tmp_idx;
        
        // Read parameters
        
        param_file.seekg(0);
        
        found_end = false;
    
        while (!found_end)
        {
            line = get_next_line(param_file);
        
               if(line.find("ENDFILE") != string::npos)
                break;    
            
            if(line.find("ATOM PAIR QUADRUPLETS:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                nquads = stoi(tmp_str_items[3]);
                                
                ncoeffs_4b      .resize(nquads);                                 
                chimes_4b_powers.resize(nquads);                                              
                chimes_4b_params.resize(nquads);                       
                chimes_4b_cutoff.resize(nquads);                            
                
                quad_params_atm_chems.resize(nquads);                
                quad_params_pair_typs.resize(nquads);
            }
            
            if(line.find("QUADRUPLETYPE PARAMS:") != string::npos)
            {
                
                line = get_next_line(param_file);
                
                split_line(line, tmp_str_items);
                
                tmp_int = stoi(tmp_str_items[1]);
                
                quad_params_atm_chems[tmp_int].push_back(tmp_str_items[3]);
                quad_params_atm_chems[tmp_int].push_back(tmp_str_items[4]);
                quad_params_atm_chems[tmp_int].push_back(tmp_str_items[5]);
                quad_params_atm_chems[tmp_int].push_back(tmp_str_items[6]);

                if (rank == 0)
                    cout << "chimesFF: " << "Read 4B parameters for quadruplets: " << tmp_int << " " << quad_params_atm_chems[tmp_int][0] << " " << quad_params_atm_chems[tmp_int][1] << " " << quad_params_atm_chems[tmp_int][2] << " " << quad_params_atm_chems[tmp_int][3]<< endl;
                
                line = get_next_line(param_file);
                
                split_line(line, tmp_str_items);
            
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[1]);
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[2]);
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[3]);
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[4]);
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[5]);
                quad_params_pair_typs[tmp_int].push_back(tmp_str_items[6]);                
            
                ncoeffs_4b[tmp_int] = stoi(tmp_str_items[10]);    

                get_next_line(param_file);
                get_next_line(param_file);
            
                vector<int> tmp_int_vec(6);
                
                for(int i=0; i<ncoeffs_4b[tmp_int]; i++)
                {                
                    line = get_next_line(param_file);
                    split_line(line, tmp_str_items);
                    
                    tmp_int_vec[0] = stoi(tmp_str_items[1]);
                    tmp_int_vec[1] = stoi(tmp_str_items[2]);
                    tmp_int_vec[2] = stoi(tmp_str_items[3]);
                    tmp_int_vec[3] = stoi(tmp_str_items[4]);
                    tmp_int_vec[4] = stoi(tmp_str_items[5]);
                    tmp_int_vec[5] = stoi(tmp_str_items[6]);
                    
                    chimes_4b_powers[tmp_int].push_back(tmp_int_vec);                 
                    
                    chimes_4b_params[tmp_int].push_back(stod(tmp_str_items[9]));
                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << 
                        chimes_4b_powers[tmp_int][i][0] << " " << 
                        chimes_4b_powers[tmp_int][i][1] << " " << 
                        chimes_4b_powers[tmp_int][i][2] << " " << 
                        chimes_4b_powers[tmp_int][i][3] << " " << 
                        chimes_4b_powers[tmp_int][i][4] << " " << 
                        chimes_4b_powers[tmp_int][i][5] << " " <<                                
                        chimes_4b_params[tmp_int][i] << endl;
                }
            }    
            
            if(line.find("QUADMAPS:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                n_quad_maps = stoi(tmp_str_items[1]);
                
                atom_idx_quad_map.resize(n_quad_maps);
                atom_typ_quad_map.resize(n_quad_maps);
                    
                if (rank == 0)            
                    cout << "chimesFF: " << "Built the following 4-body pair \"slow\" map:" << endl;
            
                for(int i=0; i<n_quad_maps; i++)
                {
                    line = get_next_line(param_file);
                    split_line(line, tmp_str_items);
                
                    atom_idx_quad_map[i] = stoi(tmp_str_items[0]);
                    atom_typ_quad_map[i] =      tmp_str_items[1];
                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << atom_idx_quad_map[i] << " " << atom_typ_quad_map[i] << endl;
                }        
                
                if (rank == 0)
                    cout << "chimesFF: " << "Built the following 4-body pair \"fast\" map:" << endl;

                atom_int_quad_map.resize(natmtyps*natmtyps*natmtyps*natmtyps);
                
                for(int i=0; i<natmtyps; i++)
                {
                    for (int j=0; j<natmtyps; j++)
                    {
                        for(int k=0; k<natmtyps; k++)
                        {
                            for(int l=0; l<natmtyps; l++)
                            {                            
                                // Get the quad type name for the set of atoms
                            
                                tmp_str = "";
                                
                                
                                tmp_int  = get_proper_pair(atmtyps[i], atmtyps[j]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    

                                tmp_int  = get_proper_pair(atmtyps[i], atmtyps[k]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    
                
                                tmp_int  = get_proper_pair(atmtyps[i], atmtyps[l]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    
                            
                                tmp_int  = get_proper_pair(atmtyps[j], atmtyps[k]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];    
                            
                                tmp_int  = get_proper_pair(atmtyps[j], atmtyps[l]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];
                            
                                tmp_int  = get_proper_pair(atmtyps[k], atmtyps[l]);
                                tmp_str += pair_params_atm_chem_1[tmp_int] + pair_params_atm_chem_2[tmp_int];                                                                                                

                                tmp_int = get_index(atom_typ_quad_map, tmp_str);
                            
                                tmp_idx = i*natmtyps*natmtyps*natmtyps + j*natmtyps*natmtyps + k*natmtyps + l;

                                atom_int_quad_map[ tmp_idx ] = atom_idx_quad_map[tmp_int];

                                if (rank == 0)
                                    cout << "chimesFF: " << "\t" << tmp_idx << " " << atom_int_quad_map[ tmp_idx  ]  << endl;
                            }
                        }
                    }
                }
            }            
        }
        
        // Set up cutoffs ... First set to match 2-body, then read special if they exist
        
        int atmtyp_1,  atmtyp_2,  atmtyp_3,  atmtyp_4;
        int pairtyp_1, pairtyp_2, pairtyp_3, pairtyp_4, pairtyp_5, pairtyp_6;
        
        for(int i=0; i<nquads; i++) 
        {
            // Figure out the atom type index for each atom in the quadruplet type 
                        
            atmtyp_1 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), quad_params_atm_chems[i][0]));    
            atmtyp_2 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), quad_params_atm_chems[i][1]));    
            atmtyp_3 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), quad_params_atm_chems[i][2]));    
            atmtyp_4 = distance(atmtyps.begin(), find(atmtyps.begin(), atmtyps.end(), quad_params_atm_chems[i][3]));    
                        
            // Figure out the corresponding 2-body pair type
            
            pairtyp_1 = atom_int_pair_map[ atmtyp_1*natmtyps + atmtyp_2 ];
            pairtyp_2 = atom_int_pair_map[ atmtyp_1*natmtyps + atmtyp_3 ];
            pairtyp_3 = atom_int_pair_map[ atmtyp_1*natmtyps + atmtyp_4 ];
            pairtyp_4 = atom_int_pair_map[ atmtyp_2*natmtyps + atmtyp_3 ];
            pairtyp_5 = atom_int_pair_map[ atmtyp_2*natmtyps + atmtyp_4 ];
            pairtyp_6 = atom_int_pair_map[ atmtyp_3*natmtyps + atmtyp_4 ];            
    
            // Set the default inner/outer cutoffs to the corresponding 2-body value                    

            chimes_4b_cutoff[i].resize(2);            

            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_1][0]);
            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_2][0]);
            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_3][0]);
            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_4][0]);
            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_5][0]);
            chimes_4b_cutoff[i][0].push_back(chimes_2b_cutoff[pairtyp_6][0]);              
            
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_1][1]);
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_2][1]);
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_3][1]);          
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_4][1]);
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_5][1]);
            chimes_4b_cutoff[i][1].push_back(chimes_2b_cutoff[pairtyp_6][1]);                                              
        }
        
        param_file.seekg(0);
        
        int    nentries;
        double cutval;
        
        found_end = false;
        
        while (!found_end)
        {
            line = get_next_line(param_file);
        
               if(line.find("ENDFILE") != string::npos)
                break;                
            
            if(line.find("SPECIAL 4B S_MAXIM:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                if (rank == 0)
                    cout << "chimesFF: " << "Set the following special 4-body outer cutoffs: " << endl;
                
                if(tmp_str_items[3] == "ALL")
                {
                    cutval = stod(tmp_str_items[4]);
                                        
                    for(int i=0; i<nquads; i++)
                    {                
                        chimes_4b_cutoff[i][1][0] = cutval;
                        chimes_4b_cutoff[i][1][1] = cutval;
                        chimes_4b_cutoff[i][1][2] = cutval;
                        chimes_4b_cutoff[i][1][3] = cutval;
                        chimes_4b_cutoff[i][1][4] = cutval;
                        chimes_4b_cutoff[i][1][5] = cutval;                                                      
                    }
                }
                else
                {
                    nentries = stoi(tmp_str_items[4]);
                    
                    vector<string> pair_name(6);
                    vector<double> cutoffval(6);

                    for(int i=0; i<nentries; i++)
                    {
                        line = get_next_line(param_file);
                        
                        split_line(line, tmp_str_items);
                        
                        tmp_int = atom_idx_quad_map[distance(atom_typ_quad_map.begin(), find(atom_typ_quad_map.begin(), atom_typ_quad_map.end(), tmp_str_items[0]))];

                        pair_name[0] = tmp_str_items[1];
                        pair_name[1] = tmp_str_items[2];
                        pair_name[2] = tmp_str_items[3];
                        pair_name[3] = tmp_str_items[4];
                        pair_name[4] = tmp_str_items[5];
                        pair_name[5] = tmp_str_items[6];
                        
                        cutoffval[0] = stod(tmp_str_items[7 ]);
                        cutoffval[1] = stod(tmp_str_items[8 ]);
                        cutoffval[2] = stod(tmp_str_items[9 ]);
                        cutoffval[3] = stod(tmp_str_items[10]);
                        cutoffval[4] = stod(tmp_str_items[11]);
                        cutoffval[5] = stod(tmp_str_items[12]);
                        
                        vector<bool>   disqualified(6,false);
                        
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[0], disqualified) ] = cutoffval[0];
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[1], disqualified) ] = cutoffval[1];
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[2], disqualified) ] = cutoffval[2];    
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[3], disqualified) ] = cutoffval[3];
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[4], disqualified) ] = cutoffval[4];
                        chimes_4b_cutoff[tmp_int][1][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[5], disqualified) ] = cutoffval[5];                                               }
                }
                
                for(int i=0; i<nquads; i++)
                {                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << i << " " 
                        << chimes_4b_cutoff[i][1][0] << " " 
                        << chimes_4b_cutoff[i][1][1] << " " 
                        << chimes_4b_cutoff[i][1][2] << " " 
                        << chimes_4b_cutoff[i][1][3] << " " 
                        << chimes_4b_cutoff[i][1][4] << " " 
                        << chimes_4b_cutoff[i][1][5] << endl;
                }                
            }

            if(line.find("SPECIAL 4B S_MINIM:") != string::npos)
            {
                split_line(line, tmp_str_items);
                
                if (rank == 0)
                    cout << "chimesFF: " << "Set the following special 4-body inner cutoffs: " << endl;
                
                if(tmp_str_items[3] == "ALL")
                {
                    cutval = stod(tmp_str_items[4]);
                    
                    for(int i=0; i<nquads; i++)
                    {                
                        chimes_4b_cutoff[i][0][0] = cutval;
                        chimes_4b_cutoff[i][0][1] = cutval;
                        chimes_4b_cutoff[i][0][2] = cutval;
                        chimes_4b_cutoff[i][0][3] = cutval;
                        chimes_4b_cutoff[i][0][4] = cutval;
                        chimes_4b_cutoff[i][0][5] = cutval;                         
                    }
                }
                else
                {
                    nentries = stoi(tmp_str_items[4]);
                    
                    vector<string> pair_name(6);
                    vector<double> cutoffval(6);

                                        for(int i=0; i<nquads; i++)
                                        {
                                                chimes_4b_cutoff[i][0][0] = -1.0;
                                                chimes_4b_cutoff[i][0][1] = -1.0;
                                                chimes_4b_cutoff[i][0][2] = -1.0;
                                                chimes_4b_cutoff[i][0][3] = -1.0;
                                                chimes_4b_cutoff[i][0][4] = -1.0;
                                                chimes_4b_cutoff[i][0][5] = -1.0;
                                        }

                    for(int i=0; i<nentries; i++)
                    {
                        line = get_next_line(param_file);
                        
                        split_line(line, tmp_str_items);
                        
                        tmp_int = atom_idx_quad_map[distance(atom_typ_quad_map.begin(), find(atom_typ_quad_map.begin(), atom_typ_quad_map.end(), tmp_str_items[0]))];

                        pair_name[0] = tmp_str_items[1];
                        pair_name[1] = tmp_str_items[2];
                        pair_name[2] = tmp_str_items[3];
                        pair_name[3] = tmp_str_items[4];
                        pair_name[4] = tmp_str_items[5];
                        pair_name[5] = tmp_str_items[6];
                        
                        cutoffval[0] = stod(tmp_str_items[7 ]);
                        cutoffval[1] = stod(tmp_str_items[8 ]);
                        cutoffval[2] = stod(tmp_str_items[9 ]);
                        cutoffval[3] = stod(tmp_str_items[10]);
                        cutoffval[4] = stod(tmp_str_items[11]);
                        cutoffval[5] = stod(tmp_str_items[12]);
                        
                        vector<bool>   disqualified(6,false);
                        
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[0], disqualified) ] = cutoffval[0];
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[1], disqualified) ] = cutoffval[1];
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[2], disqualified) ] = cutoffval[2];    
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[3], disqualified) ] = cutoffval[3];
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[4], disqualified) ] = cutoffval[4];
                        chimes_4b_cutoff[tmp_int][0][ get_index_if(quad_params_pair_typs[tmp_int], pair_name[5], disqualified) ] = cutoffval[5];                     
                    }
                }
                
                for(int i=0; i<nquads; i++)
                {                
                    if (rank == 0)
                        cout << "chimesFF: " << "\t" << i << " " 
                        << chimes_4b_cutoff[i][1][0] << " " 
                        << chimes_4b_cutoff[i][1][1] << " " 
                        << chimes_4b_cutoff[i][1][2] << " " 
                        << chimes_4b_cutoff[i][1][3] << " " 
                        << chimes_4b_cutoff[i][1][4] << " " 
                        << chimes_4b_cutoff[i][1][5] << endl;
                }                
            }            
        }    
    }
    
    param_file.close();    
}

inline void chimesFF::set_cheby_polys(double *Tn, double *Tnd, const double dx, const int pair_idx, const double inner_cutoff, const double outer_cutoff, const int bodiedness_idx)
{
    // Currently assumes a Morse-style transformation has been requested
    
    // Sets the value of the Chebyshev polynomials (Tn) and their derivatives (Tnd).  Tnd is the derivative
    // with respect to the interatomic distance, not the transformed distance (x).
    
    // Do the Morse transformation
    
    double x_min = exp(-1*inner_cutoff/morse_var[pair_idx]);
    double x_max = exp(-1*outer_cutoff/morse_var[pair_idx]);
    
    double x_avg   = 0.5 * (x_max + x_min);
    double x_diff  = 0.5 * (x_max - x_min);
    
    x_diff *= -1.0; // Special for Morse style
    
    double exprlen = exp(-1*dx/morse_var[pair_idx]);
    
    double x  = (exprlen - x_avg)/x_diff;

    if ( x < -1.0)
        x =  -1.0;
    else if ( x > 1.0 )
        x =  1.0;                            

    // Generate Chebyshev polynomials by recursion. 
    // 
    // What we're doing here. Want to fit using Cheby polynomials of the 1st kinD[i]. "T_n(x)."
    // We need to calculate the derivative of these polynomials.
    // Derivatives are defined through use of Cheby polynomials of the 2nd kind "U_n(x)", as:
    //
    // d/dx[ T_n(x) = n * U_n-1(x)] 
    // 
    // So we need to first set up the 1st-kind polynomials ("Tn[]")
    // Then, to compute the derivatives ("Tnd[]"), first set equal to the 2nd-kind, then multiply by n to get the der's
     
    // First two 1st-kind Chebys:
        
    Tn[0] = 1.0;
    Tn[1] = x;
    
    // Start the derivative setup. Set the first two 1st-kind Cheby's equal to the first two of the 2nd-kind

    Tnd[0] = 1.0;
    Tnd[1] = 2.0 * x;
    
    // Use recursion to set up the higher n-value Tn and Tnd's

    for ( int i = 2; i <= poly_orders[bodiedness_idx]; i++ ) 
    {
        Tn[i]  = 2.0 * x *  Tn[i-1] -  Tn[i-2];
        Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
    }
    
    // Now multiply by n to convert Tnd's to actual derivatives of Tn
    
    // The following dx_dr compuation assumes a Morse transformation
    // DERIV_CONST is no longer used. (old way: dx_dr = DERIV_CONST*cheby_var_deriv(x_diff, rlen, ff_2body.LAMBDA, ff_2body.CHEBY_TYPE, exprlen);)

    double dx_dr = (-exprlen/morse_var[pair_idx])/x_diff;

    for ( int i = poly_orders[bodiedness_idx]; i >= 1; i-- ) 
        Tnd[i] = i * dx_dr * Tnd[i-1];

    Tnd[0] = 0.0;

}

inline void chimesFF::get_fcut(const double dx, const double outer_cutoff, double & fcut, double & fcutderiv)
{

    static double THRESH;
    static double fcut0;
    static double fcut0_deriv, fcut_deriv;
    
    if(fcut_type == "CUBIC")
    {        
        fcut0 = (1.0 - dx/outer_cutoff);
        fcut        = pow(fcut0,3.0);
        fcutderiv  = pow(fcut0,2.0);
        fcutderiv *= -1.0 * 3.0 /outer_cutoff;

        return;
    }

    if(fcut_type=="TERSOFF")
    {
    
        THRESH = outer_cutoff-fcut_var*outer_cutoff;    
    
    
        if      (dx < THRESH)        // Case 1: Our pair distance is less than the fcut kick-in distance
        {
            fcut       = 1.0;
            fcutderiv = 0.0;
        }
        else if (dx > outer_cutoff)        // Case 2: Our pair distance is greater than the cutoff
        {
            fcut       = 0.0;
            fcutderiv = 0.0;
        }                    
        else                // Case 3: We'll use our modified sin function
        {
            fcut0       = (dx-THRESH) / (outer_cutoff-THRESH) * pi + pi/2.0;
            fcut0_deriv = pi / (outer_cutoff - THRESH);
            
            fcut        = 0.5 + 0.5 * sin( fcut0 );
            fcutderiv  = 0.5 * cos( fcut0 ) * fcut0_deriv; 
        }
        return;
    }
}

inline void chimesFF::get_penalty(const double dx, const int & pair_idx, double & E_penalty, double & force_scalar)
{
    double r_penalty = 0.0;
    
    E_penalty    = 0.0;
    force_scalar = 1.0;
    
    if (dx - penalty_params[1] < chimes_2b_cutoff[pair_idx][0])
        
        r_penalty = chimes_2b_cutoff[pair_idx][0] + penalty_params[0] - dx;
        
    if ( r_penalty > 0.0 ) 
    {        
        E_penalty    = r_penalty * r_penalty * r_penalty * penalty_params[1];
        force_scalar = 3.0 * r_penalty * r_penalty * penalty_params[1];
        
        if (rank == 0)
        {
            cout << "chimesFF: " << "Adding penalty in 2B Cheby calc, r < rmin+penalty_dist " << fixed 
                 << dx << " " 
                 << chimes_2b_cutoff[pair_idx][0] + penalty_params[0]  
                 << " pair type: " << pair_idx << endl;
            cout << "chimesFF: " << "\t...Penalty potential = "<< E_penalty << endl;
        }
    }   
}

inline void chimesFF::build_atom_and_pair_mappers(const int natoms, const int npairs, const vector<int> typ_idxs, const vector<vector<string> > & clu_params_pair_typs, const int & cluidx, vector<int >  & mapped_pair_idx)
{
    // Generate permutations for atoms... all we are doing is permuting the possible indices for typ_idxs


    // build a copy of the atom type vector for permuting

    static vector<int> tmp_typ_idxs;
    static int         nelements;
    
    nelements = typ_idxs.size();
    tmp_typ_idxs.resize(nelements);
    
    for(int i=0; i<nelements; i++)
        tmp_typ_idxs[i] = i;
        
    // Build a copy of the original pairs for comparison against permuted pairs
    
    
    static vector<vector<int> > tmp_pairs;
    tmp_pairs.resize(npairs,vector<int>(2));
    static vector<vector<int> > runtime_pairs;
    runtime_pairs.resize(npairs,vector<int>(2));
    
    int idx = 0;
    
    for(int i=0; i<natoms; i++) 
    {
        for (int j=i+1; j<natoms; j++)
        {
            tmp_pairs[idx][0] = i;
            tmp_pairs[idx][1] = j;
            
            idx++;
        }
    }
        
    vector<string> runtime_pair_typs(npairs);    
        
    while ( true )
    {
        // Check if the permutation leads to pair types that match the order specified by the force field type
    
        idx = 0;
        
        for(int i=0; i<natoms; i++) // Associate the current atom pairs with a "proper" 2-body force field name
        {
            for (int j=i+1; j<natoms; j++)
            {
                runtime_pair_typs[idx] = atom_int_prpr_map[ typ_idxs[tmp_typ_idxs[i]]*natmtyps + typ_idxs[tmp_typ_idxs[j]] ];
                        
                idx++;
            }
        }
        
        bool match = true;
        
        for(int i=0; i<npairs; i++)
        {
            if (clu_params_pair_typs[cluidx][i] != runtime_pair_typs[i])
            {
                match = false;
                break;
            }
        }
        
        if (match) // Then we've found an appropriate atom ordering... now what?
        {
            idx = 0;
        
            for(int i=0; i<natoms; i++) // Associate the current atom pairs with a "proper" 2-body force field name
            {
                for (int j=i+1; j<natoms; j++)
                {
                    runtime_pairs[idx][0] = tmp_typ_idxs[i];
                    runtime_pairs[idx][1] = tmp_typ_idxs[j];
                    
                    idx++;
                }
            }
        
            break;
        }
		
		if(!next_permutation(tmp_typ_idxs.begin(),tmp_typ_idxs.begin()+nelements))
			break;
    }
    
    // Once we've found a re-ordering of atoms that properly maps to the force field pair types, need to figure out how to convert that to a map between *pairs*
    
    idx = 0;
    
    for(int i=0; i<npairs; i++)
        for(int j=0; j<npairs; j++)
            if (((runtime_pairs[i][0] == tmp_pairs[j][0]) &&  (runtime_pairs[i][1] == tmp_pairs[j][1])) || ((runtime_pairs[i][0] == tmp_pairs[j][1]) &&  (runtime_pairs[i][1] == tmp_pairs[j][0])))
                mapped_pair_idx[j] = i;

}

void chimesFF::compute_1B(const int typ_idx, double & energy )
{
    // Compute 1b (input: a single atom type index... outputs (updates) energy

    energy += energy_offsets[typ_idx];
}

void chimesFF::compute_2B(const double dx, const vector<double> & dr, const vector<int> typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy )
{
    // Compute 2b (input: 2 atoms or distances, corresponding types... outputs (updates) force, acceleration, energy, stress
    //
    // Input parameters:
    //
    // dx: Scalar (pair distance)
    // dr: 1d-Array (pair distance: [x, y, and z-component]) 
    // Force: [natoms in interaction set][x,y, and z-component] *note
    // Stress [sxx, sxy, sxz, syy, syz, szz]  *note
    // Energy: Scalar; energy for interaction set
    // Assumes atom indices start from zero
    // Assumes distances are atom_2 - atom_1
    //
    // *note: force's pointer vector of vectors elements must be dereferenced! 

    static bool    called_before = false; // Used to insure we only declare the Tn's and Tnd's one time
    static double  *Tn, *Tnd;             // The Chebyshev polymonials and thier derivatives
    static int     pair_idx;    
    static double  fcut;
    static double  fcutderiv;
    static double  dx_dr;
    static double  E_penalty;
    static double  force_scalar;
    static double  fpair_total;
    static double  deriv;
    static double  coeff_val;

    if (!called_before)     // Set up 2-body polynomials
    {
        called_before = true;

        Tn   = new double [poly_orders[0]+1];
        Tnd  = new double [poly_orders[0]+1];
    }
    
    pair_idx = atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[1] ];

    if (dx >= chimes_2b_cutoff[pair_idx][1])
        return;    

    set_cheby_polys(Tn, Tnd, dx, pair_idx, chimes_2b_cutoff[pair_idx][0], chimes_2b_cutoff[pair_idx][1], 0);
    
    get_fcut(dx, chimes_2b_cutoff[pair_idx][1], fcut, fcutderiv);
    
    for(int coeffs=0; coeffs<ncoeffs_2b[pair_idx]; coeffs++)
    {
        coeff_val = chimes_2b_params[pair_idx][coeffs];        
        
        energy += coeff_val * fcut * Tn[ chimes_2b_pows[pair_idx][coeffs]+1 ];
                                                
        deriv = fcut * Tnd[ chimes_2b_pows[pair_idx][coeffs]+1 ]  + fcutderiv * Tn[ chimes_2b_pows[pair_idx][coeffs]+1 ];    

        force_scalar = coeff_val * deriv; 

        *force[0][0] += force_scalar * dr[0]/ dx;
        *force[0][1] += force_scalar * dr[1]/ dx;
        *force[0][2] += force_scalar * dr[2]/ dx;
        
        *force[1][0] -= force_scalar * dr[0]/ dx;
        *force[1][1] -= force_scalar * dr[1]/ dx;
        *force[1][2] -= force_scalar * dr[2]/ dx;
        
        // xx xy xz yy yz zz
        // 0  1  2  3  4  5
        
        // xx xy xz yx yy yz zx zy zz
        // 0  1  2  3  4  5  6  7  8
        // *           *           *
        
        *stress[0] -= force_scalar / dx * dr[0] * dr[0]; // xx tensor component
        *stress[1] -= force_scalar / dx * dr[0] * dr[1]; // xy tensor component 
        *stress[2] -= force_scalar / dx * dr[0] * dr[2]; // xz tensor component
        
        *stress[3] -= force_scalar / dx * dr[1] * dr[0]; // yx tensor component
        *stress[4] -= force_scalar / dx * dr[1] * dr[1]; // yy tensor component
        *stress[5] -= force_scalar / dx * dr[1] * dr[2]; // yz tensor component
        
        *stress[6] -= force_scalar / dx * dr[2] * dr[0]; // zx tensor component
        *stress[7] -= force_scalar / dx * dr[2] * dr[1]; // zy tensor component
        *stress[8] -= force_scalar / dx * dr[2] * dr[2]; // zz tensor component
            
    }
    
    get_penalty(dx, pair_idx, E_penalty , force_scalar); // Initialize with E_penalty = 0.0, force_scalar = 1.0
    
    if ( E_penalty > 0.0 ) 
    {
        energy += E_penalty;

        *force[0][0] += force_scalar * dr[0]/ dx;
        *force[0][1] += force_scalar * dr[1]/ dx;
        *force[0][2] += force_scalar * dr[2]/ dx;
        
        *force[1][0] -= force_scalar * dr[0]/ dx;
        *force[1][1] -= force_scalar * dr[1]/ dx;
        *force[1][2] -= force_scalar * dr[2]/ dx;
    }
}

void chimesFF::compute_3B(const vector<double> & dx, const vector<vector<double> > & dr, const vector<int> & typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy )
{
    // Compute 3b (input: 3 atoms or distances, corresponding types... outputs (updates) force, acceleration, energy, stress
    //
    // Input parameters:
    //
    // dx_ij: Scalar (pair distance)
    // dr_ij: 1d-Array (pair distance: [x, y, and z-component])
    // Force: [natoms in interaction set][x,y, and z-component] *note
    // Stress [sxx, sxy, sxz, syy, syz, szz]  *note
    // Energy: Scalar; energy for interaction set
    // Assumes atom indices start from zero
    // Assumes distances are atom_2 - atom_1
    //
    // *note: force's pointer vector of vectors elements must be dereferenced! 
    
    static int natoms = 3;                   // Number of atoms in an interaction set
    static int npairs = natoms*(natoms-1)/2; // Number of pairs in an interaction set
    
        
    static bool    called_before = false;        // Used to insure we only declare the Tn's and Tnd's one time
    static double  *Tn_ij,   *Tn_ik,   *Tn_jk;   // The Chebyshev polymonials
    static double  *Tnd_ij,  *Tnd_ik,  *Tnd_jk;  // The Chebyshev polymonial derivatives
    
    static vector<double> fcut     (npairs);
    static vector<double> fcutderiv(npairs);
    static vector<double> deriv    (npairs);
    
    if (!called_before)    // Set up 3-body polynomials
    {
        called_before = true;

        Tn_ij   = new double [poly_orders[1]+1];
        Tn_ik   = new double [poly_orders[1]+1];
        Tn_jk   = new double [poly_orders[1]+1];

        Tnd_ij  = new double [poly_orders[1]+1];
        Tnd_ik  = new double [poly_orders[1]+1];
        Tnd_jk  = new double [poly_orders[1]+1];       
    }
            
    int tripidx = atom_int_trip_map[typ_idxs[0]*natmtyps*natmtyps + typ_idxs[1]*natmtyps + typ_idxs[2]];

    if(tripidx < 0)    // Skipping an excluded interaction
        return;
    
    // Build maps for associating atoms/pairs with corresponding force field atom/pair
    
    vector<int > mapped_pair_idx(npairs);        // Will store the force field index for the atom type
    
    build_atom_and_pair_mappers(natoms, npairs, typ_idxs, trip_params_pair_typs, tripidx, mapped_pair_idx);    
    
    // Check whether cutoffs are within allowed ranges
        
    if (dx[0] >= chimes_3b_cutoff[ tripidx ][1][mapped_pair_idx[0]])    // ij
        return;    
    if (dx[1] >= chimes_3b_cutoff[ tripidx ][1][mapped_pair_idx[1]])    // ik
        return;    
    if (dx[2] >= chimes_3b_cutoff[ tripidx ][1][mapped_pair_idx[2]])    // jk
        return;    
    
    // At this point, all distances are within allowed ranges. We can now proceed to the force/stress/energy calculation
    
    vector<double> dx_dummy(npairs);

    for (int i=0; i<npairs; i++)
        dx_dummy[i] = dx[i];
        
    if (dx_dummy[0] < chimes_3b_cutoff[tripidx][0][mapped_pair_idx[0]])    // ij
        dx_dummy[0] = chimes_3b_cutoff[tripidx][0][mapped_pair_idx[0]];
    
    if (dx_dummy[1] < chimes_3b_cutoff[tripidx][0][mapped_pair_idx[1]])    // ik
        dx_dummy[1] = chimes_3b_cutoff[tripidx][0][mapped_pair_idx[1]];
    
    if (dx_dummy[2] < chimes_3b_cutoff[tripidx][0][mapped_pair_idx[2]])    // jk
        dx_dummy[2] = chimes_3b_cutoff[tripidx][0][mapped_pair_idx[2]];
    

    // Set up the polynomials

    set_cheby_polys(Tn_ij, Tnd_ij, dx_dummy[0], atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[1] ], chimes_3b_cutoff[tripidx][0][mapped_pair_idx[0]], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[0]], 1);
    set_cheby_polys(Tn_ik, Tnd_ik, dx_dummy[1], atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[2] ], chimes_3b_cutoff[tripidx][0][mapped_pair_idx[1]], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[1]], 1);
    set_cheby_polys(Tn_jk, Tnd_jk, dx_dummy[2], atom_int_pair_map[ typ_idxs[1]*natmtyps + typ_idxs[2] ], chimes_3b_cutoff[tripidx][0][mapped_pair_idx[2]], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[2]], 1);
    
    
    // Set up the smoothing functions
        
    get_fcut(dx[0], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[0]], fcut[0], fcutderiv[0]);
    get_fcut(dx[1], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[1]], fcut[1], fcutderiv[1]);
    get_fcut(dx[2], chimes_3b_cutoff[tripidx][1][mapped_pair_idx[2]], fcut[2], fcutderiv[2]);

    // Start the force/stress/energy calculation
        
    double coeff;
    vector<int> powers(npairs);
    vector<double> force_scalar(npairs);
    
    for(int coeffs=0; coeffs<ncoeffs_3b[tripidx]; coeffs++)
    {
        coeff = chimes_3b_params[tripidx][coeffs];
        
        powers[0] = chimes_3b_powers[tripidx][coeffs][mapped_pair_idx[0]];
        powers[1] = chimes_3b_powers[tripidx][coeffs][mapped_pair_idx[1]];
        powers[2] = chimes_3b_powers[tripidx][coeffs][mapped_pair_idx[2]];
        
        energy += coeff * fcut[0] * fcut[1] * fcut[2] * Tn_ij[ powers[0] ] * Tn_ik[ powers[1] ] * Tn_jk[ powers[2] ];    

        deriv[0] = fcut[0] * Tnd_ij[ powers[0] ] + fcutderiv[0] * Tn_ij[ powers[0] ];
        deriv[1] = fcut[1] * Tnd_ik[ powers[1] ] + fcutderiv[1] * Tn_ik[ powers[1] ];
        deriv[2] = fcut[2] * Tnd_jk[ powers[2] ] + fcutderiv[2] * Tn_jk[ powers[2] ];

        force_scalar[0]  = coeff * deriv[0] * fcut[1] * fcut[2] * Tn_ik[powers[1]]  * Tn_jk[powers[2]];
        force_scalar[1]  = coeff * deriv[1] * fcut[0] * fcut[2] * Tn_ij[powers[0]]  * Tn_jk[powers[2]];
        force_scalar[2]  = coeff * deriv[2] * fcut[0] * fcut[1] * Tn_ij[powers[0]]  * Tn_ik[powers[1]];
        
        // Accumulate forces/stresses on/from the ij pair
        
        *force[0][0] += force_scalar[0] * dr[0][0]/ dx[0];
        *force[0][1] += force_scalar[0] * dr[0][1]/ dx[0];
        *force[0][2] += force_scalar[0] * dr[0][2]/ dx[0];

        *force[1][0] -= force_scalar[0] * dr[0][0]/ dx[0];
        *force[1][1] -= force_scalar[0] * dr[0][1]/ dx[0];
        *force[1][2] -= force_scalar[0] * dr[0][2]/ dx[0];   
        
        *stress[0] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][0]; // xx tensor component
        *stress[4] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][1]; // yy tensor component
        *stress[8] -= force_scalar[0] / dx[0] * dr[0][2] * dr[0][2]; // zz tensor component
        
        *stress[1] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][1]; // xy tensor component
        *stress[2] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][2]; // xz tensor component
        *stress[5] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][2]; // yz tensor component
        
        *stress[3] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][1]; // yx tensor component
        *stress[6] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][2]; // zx tensor component
        *stress[7] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][2]; // zy tensor component
        
        
        // Accumulate forces/stresses on/from the ik pair
        
        *force[0][0] += force_scalar[1] * dr[1][0]/ dx[1];
        *force[0][1] += force_scalar[1] * dr[1][1]/ dx[1];
        *force[0][2] += force_scalar[1] * dr[1][2]/ dx[1];

        *force[2][0] -= force_scalar[1] * dr[1][0]/ dx[1];
        *force[2][1] -= force_scalar[1] * dr[1][1]/ dx[1];
        *force[2][2] -= force_scalar[1] * dr[1][2]/ dx[1];   

        *stress[0] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][0]; // xx tensor component
        *stress[4] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][1]; // yy tensor component
        *stress[8] -= force_scalar[1] / dx[1] * dr[1][2] * dr[1][2]; // zz tensor component
        
        *stress[1] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][1]; // xy tensor component
        *stress[2] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][2]; // xz tensor component
        *stress[5] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][2]; // yz tensor component
        
        *stress[3] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][1]; // yx tensor component
        *stress[6] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][2]; // zx tensor component
        *stress[7] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][2]; // zy tensor component                     
        
        // Accumulate forces/stresses on/from the jk pair
        
        *force[1][0] += force_scalar[2] * dr[2][0]/ dx[2];
        *force[1][1] += force_scalar[2] * dr[2][1]/ dx[2];
        *force[1][2] += force_scalar[2] * dr[2][2]/ dx[2];

        *force[2][0] -= force_scalar[2] * dr[2][0]/ dx[2];
        *force[2][1] -= force_scalar[2] * dr[2][1]/ dx[2];
        *force[2][2] -= force_scalar[2] * dr[2][2]/ dx[2];   
        
        *stress[0] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][0]; // xx tensor component
        *stress[4] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][1]; // yy tensor component
        *stress[8] -= force_scalar[2] / dx[2] * dr[2][2] * dr[2][2]; // zz tensor component
        
        *stress[1] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][1]; // xy tensor component
        *stress[2] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][2]; // xz tensor component
        *stress[5] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][2]; // yz tensor component
        
        *stress[3] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][1]; // yx tensor component
        *stress[6] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][2]; // zx tensor component
        *stress[7] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][2]; // zy tensor component                             
    }

    return;    
}

void chimesFF::compute_4B(const vector<double> & dx, const vector<vector<double> > & dr, const vector<int> & typ_idxs, vector<vector<double* > > force, vector<double*> stress, double & energy )
{
    // Compute 3b (input: 3 atoms or distances, corresponding types... outputs (updates) force, acceleration, energy, stress
    //
    // Input parameters:
    //
    // dx_ij: Scalar (pair distance)
    // dr_ij: 1d-Array (pair distance: [x, y, and z-component])
    // Force: [natoms in interaction set][x,y, and z-component] *note
    // Stress [sxx, sxy, sxz, syy, syz, szz]  *note
    // Energy: Scalar; energy for interaction set
    // Assumes atom indices start from zero
    // Assumes distances are atom_2 - atom_1
    //
    // *note: force's pointer vector of vectors elements must be dereferenced! 
    
    int natoms = 4;                     // Number of atoms in an interaction set
    int npairs = natoms*(natoms-1)/2;    // Number of pairs in an interaction set
        
    static bool    called_before = false;            // Used to insure we only declare the Tn's and Tnd's one time
    static double    *Tn_ij,   *Tn_ik,   *Tn_il,  *Tn_jk,  *Tn_jl,  *Tn_kl;    // The Chebyshev polymonials
    static double    *Tnd_ij,  *Tnd_ik,  *Tnd_il,  *Tnd_jk, *Tnd_jl, *Tnd_kl;    // The Chebyshev polymonial derivatives
    
    vector<double> fcut     (npairs);
    vector<double> fcutderiv(npairs);
    vector<double> deriv    (npairs);

    if (!called_before)    // Set up 4-body polynomials
    {
        called_before = true;

        Tn_ij   = new double [poly_orders[2]+1];
        Tn_ik   = new double [poly_orders[2]+1];
        Tn_il   = new double [poly_orders[2]+1];
        Tn_jk   = new double [poly_orders[2]+1];
        Tn_jl   = new double [poly_orders[2]+1];
        Tn_kl   = new double [poly_orders[2]+1];        
                                          
        Tnd_ij  = new double [poly_orders[2]+1];
        Tnd_ik  = new double [poly_orders[2]+1];
        Tnd_il  = new double [poly_orders[2]+1];  
        Tnd_jk  = new double [poly_orders[2]+1];
        Tnd_jl  = new double [poly_orders[2]+1];
        Tnd_kl  = new double [poly_orders[2]+1];              
    }
    
    int quadidx = atom_int_quad_map[typ_idxs[0]*natmtyps*natmtyps*natmtyps + typ_idxs[1]*natmtyps*natmtyps + typ_idxs[2]*natmtyps + typ_idxs[3]];

    if(quadidx < 0)    // Skipping an excluded interaction
        return;
    

    // Build maps for associating atoms/pairs with corresponding force field atom/pair

    vector<int > mapped_pair_idx(npairs);        // Will store the force field index for the atom type
    
    build_atom_and_pair_mappers(natoms, npairs, typ_idxs, quad_params_pair_typs, quadidx, mapped_pair_idx);
        
    
    // Check whether cutoffs are within allowed ranges

    for(int i=0; i<npairs; i++)
        if (dx[i] >= chimes_4b_cutoff[ quadidx ][1][mapped_pair_idx[i]])
            return;    

    // At this point, all distances are within allowed ranges. We can now proceed to the force/stress/energy calculation
    
    vector<double> dx_dummy(npairs);

    for (int i=0; i<npairs; i++)
        dx_dummy[i] = dx[i];
        
    for (int i=0; i<npairs; i++)    
        if (dx_dummy[i] < chimes_4b_cutoff[quadidx][0][mapped_pair_idx[i]])
            dx_dummy[i] = chimes_4b_cutoff[quadidx][0][mapped_pair_idx[i]];        
            

    // Set up the polynomials
    
    set_cheby_polys(Tn_ij, Tnd_ij, dx_dummy[0], atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[1] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[0]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[0]], 2);
    set_cheby_polys(Tn_ik, Tnd_ik, dx_dummy[1], atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[2] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[1]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[1]], 2);
    set_cheby_polys(Tn_il, Tnd_il, dx_dummy[2], atom_int_pair_map[ typ_idxs[0]*natmtyps + typ_idxs[3] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[2]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[2]], 2);
    set_cheby_polys(Tn_jk, Tnd_jk, dx_dummy[3], atom_int_pair_map[ typ_idxs[1]*natmtyps + typ_idxs[2] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[3]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[3]], 2);
    set_cheby_polys(Tn_jl, Tnd_jl, dx_dummy[4], atom_int_pair_map[ typ_idxs[1]*natmtyps + typ_idxs[3] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[4]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[4]], 2);
    set_cheby_polys(Tn_kl, Tnd_kl, dx_dummy[5], atom_int_pair_map[ typ_idxs[2]*natmtyps + typ_idxs[3] ], chimes_4b_cutoff[quadidx][0][mapped_pair_idx[5]], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[5]], 2);     
    
    
    // Set up the smoothing functions
    
    for (int i=0; i<npairs; i++)    
        get_fcut(dx[i], chimes_4b_cutoff[quadidx][1][mapped_pair_idx[i]], fcut[i], fcutderiv[i]);
            
    // Start the force/stress/energy calculation
        
    double coeff;
    vector<int> powers(npairs);
    vector<double> force_scalar(npairs);

    for(int coeffs=0; coeffs<ncoeffs_4b[quadidx]; coeffs++)
    {
        coeff = chimes_4b_params[quadidx][coeffs];
        
        for (int i=0; i<npairs; i++)
            powers[i] = chimes_4b_powers[quadidx][coeffs][mapped_pair_idx[i]];

        energy += coeff * fcut[0] * fcut[1] * fcut[2] * fcut[3] * fcut[4] * fcut[5] 
                    * Tn_ij[ powers[0] ] * Tn_ik[ powers[1] ] * Tn_il[ powers[2] ] 
                * Tn_jk[ powers[3] ] * Tn_jl[ powers[4] ] * Tn_kl[ powers[5] ];        

        deriv[0] = fcut[0] * Tnd_ij[ powers[0] ] + fcutderiv[0] * Tn_ij[ powers[0] ];
        deriv[1] = fcut[1] * Tnd_ik[ powers[1] ] + fcutderiv[1] * Tn_ik[ powers[1] ];
        deriv[2] = fcut[2] * Tnd_il[ powers[2] ] + fcutderiv[2] * Tn_il[ powers[2] ];
        deriv[3] = fcut[3] * Tnd_jk[ powers[3] ] + fcutderiv[3] * Tn_jk[ powers[3] ];
        deriv[4] = fcut[4] * Tnd_jl[ powers[4] ] + fcutderiv[4] * Tn_jl[ powers[4] ];
        deriv[5] = fcut[5] * Tnd_kl[ powers[5] ] + fcutderiv[5] * Tn_kl[ powers[5] ];        

        force_scalar[0]  = coeff * deriv[0] * fcut[1] * fcut[2] * fcut[3] * fcut[4] * fcut[5] * Tn_ik[powers[1]]  * Tn_il[powers[2]] * Tn_jk[powers[3]]  * Tn_jl[powers[4]] * Tn_kl[powers[5]];
        force_scalar[1]  = coeff * deriv[1] * fcut[0] * fcut[2] * fcut[3] * fcut[4] * fcut[5] * Tn_ij[powers[0]]  * Tn_il[powers[2]] * Tn_jk[powers[3]]  * Tn_jl[powers[4]] * Tn_kl[powers[5]];
        force_scalar[2]  = coeff * deriv[2] * fcut[0] * fcut[1] * fcut[3] * fcut[4] * fcut[5] * Tn_ij[powers[0]]  * Tn_ik[powers[1]] * Tn_jk[powers[3]]  * Tn_jl[powers[4]] * Tn_kl[powers[5]];
        force_scalar[3]  = coeff * deriv[3] * fcut[0] * fcut[1] * fcut[2] * fcut[4] * fcut[5] * Tn_ij[powers[0]]  * Tn_ik[powers[1]] * Tn_il[powers[2]]  * Tn_jl[powers[4]] * Tn_kl[powers[5]];
        force_scalar[4]  = coeff * deriv[4] * fcut[0] * fcut[1] * fcut[2] * fcut[3] * fcut[5] * Tn_ij[powers[0]]  * Tn_ik[powers[1]] * Tn_il[powers[2]]  * Tn_jk[powers[3]] * Tn_kl[powers[5]];
        force_scalar[5]  = coeff * deriv[5] * fcut[0] * fcut[1] * fcut[2] * fcut[3] * fcut[4] * Tn_ij[powers[0]]  * Tn_ik[powers[1]] * Tn_il[powers[2]]  * Tn_jk[powers[3]] * Tn_jl[powers[4]];

        // Accumulate forces/stresses on/from the ij pair
        
        *force[0][0] += force_scalar[0] * dr[0][0]/ dx[0];
        *force[0][1] += force_scalar[0] * dr[0][1]/ dx[0];
        *force[0][2] += force_scalar[0] * dr[0][2]/ dx[0];

        *force[1][0] -= force_scalar[0] * dr[0][0]/ dx[0];
        *force[1][1] -= force_scalar[0] * dr[0][1]/ dx[0];
        *force[1][2] -= force_scalar[0] * dr[0][2]/ dx[0];   

        *stress[0] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][0]; // xx tensor component
        *stress[4] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][1]; // yy tensor component
        *stress[8] -= force_scalar[0] / dx[0] * dr[0][2] * dr[0][2]; // zz tensor component
        
        *stress[1] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][1]; // xy tensor component
        *stress[2] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][2]; // xz tensor component
        *stress[5] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][2]; // yz tensor component
        
        *stress[3] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][1]; // yx tensor component
        *stress[6] -= force_scalar[0] / dx[0] * dr[0][0] * dr[0][2]; // zx tensor component
        *stress[7] -= force_scalar[0] / dx[0] * dr[0][1] * dr[0][2]; // zy tensor component
        
        // Accumulate forces/stresses on/from the ik pair
        
        *force[0][0] += force_scalar[1] * dr[1][0]/ dx[1];
        *force[0][1] += force_scalar[1] * dr[1][1]/ dx[1];
        *force[0][2] += force_scalar[1] * dr[1][2]/ dx[1];

        *force[2][0] -= force_scalar[1] * dr[1][0]/ dx[1];
        *force[2][1] -= force_scalar[1] * dr[1][1]/ dx[1];
        *force[2][2] -= force_scalar[1] * dr[1][2]/ dx[1];   
        
        *stress[0] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][0]; // xx tensor component
        *stress[4] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][1]; // yy tensor component
        *stress[8] -= force_scalar[1] / dx[1] * dr[1][2] * dr[1][2]; // zz tensor component
        
        *stress[1] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][1]; // xy tensor component
        *stress[2] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][2]; // xz tensor component
        *stress[5] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][2]; // yz tensor component
        
        *stress[3] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][1]; // yx tensor component
        *stress[6] -= force_scalar[1] / dx[1] * dr[1][0] * dr[1][2]; // zx tensor component
        *stress[7] -= force_scalar[1] / dx[1] * dr[1][1] * dr[1][2]; // zy tensor component             
        
        // Accumulate forces/stresses on/from the il pair
        
        *force[0][0] += force_scalar[2] * dr[2][0]/ dx[2];
        *force[0][1] += force_scalar[2] * dr[2][1]/ dx[2];
        *force[0][2] += force_scalar[2] * dr[2][2]/ dx[2];

        *force[3][0] -= force_scalar[2] * dr[2][0]/ dx[2];
        *force[3][1] -= force_scalar[2] * dr[2][1]/ dx[2];
        *force[3][2] -= force_scalar[2] * dr[2][2]/ dx[2];   
        
        *stress[0] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][0]; // xx tensor component
        *stress[4] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][1]; // yy tensor component
        *stress[8] -= force_scalar[2] / dx[2] * dr[2][2] * dr[2][2]; // zz tensor component           
        
        *stress[1] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][1]; // xy tensor component
        *stress[2] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][2]; // xz tensor component
        *stress[5] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][2]; // yz tensor component
        
        *stress[3] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][1]; // yx tensor component
        *stress[6] -= force_scalar[2] / dx[2] * dr[2][0] * dr[2][2]; // zx tensor component
        *stress[7] -= force_scalar[2] / dx[2] * dr[2][1] * dr[2][2]; // zy tensor component
        
        // Accumulate forces/stresses on/from the jk pair
        
        *force[1][0] += force_scalar[3] * dr[3][0]/ dx[3];
        *force[1][1] += force_scalar[3] * dr[3][1]/ dx[3];
        *force[1][2] += force_scalar[3] * dr[3][2]/ dx[3];

        *force[2][0] -= force_scalar[3] * dr[3][0]/ dx[3];
        *force[2][1] -= force_scalar[3] * dr[3][1]/ dx[3];
        *force[2][2] -= force_scalar[3] * dr[3][2]/ dx[3];   
        
        *stress[0] -= force_scalar[3] / dx[3] * dr[3][0] * dr[3][0]; // xx tensor component
        *stress[4] -= force_scalar[3] / dx[3] * dr[3][1] * dr[3][1]; // yy tensor component
        *stress[8] -= force_scalar[3] / dx[3] * dr[3][2] * dr[3][2]; // zz tensor component
        
        *stress[1] -= force_scalar[3] / dx[3] * dr[3][0] * dr[3][1]; // xy tensor component
        *stress[2] -= force_scalar[3] / dx[3] * dr[3][0] * dr[3][2]; // xz tensor component
        *stress[5] -= force_scalar[3] / dx[3] * dr[3][1] * dr[3][2]; // yz tensor component
        
        *stress[3] -= force_scalar[3] / dx[3] * dr[3][0] * dr[3][1]; // yx tensor component
        *stress[6] -= force_scalar[3] / dx[3] * dr[3][0] * dr[3][2]; // zx tensor component
        *stress[7] -= force_scalar[3] / dx[3] * dr[3][1] * dr[3][2]; // zy tensor component                           

        // Accumulate forces/stresses on/from the jl pair
        
        *force[1][0] += force_scalar[4] * dr[4][0]/ dx[4];
        *force[1][1] += force_scalar[4] * dr[4][1]/ dx[4];
        *force[1][2] += force_scalar[4] * dr[4][2]/ dx[4];

        *force[3][0] -= force_scalar[4] * dr[4][0]/ dx[4];
        *force[3][1] -= force_scalar[4] * dr[4][1]/ dx[4];
        *force[3][2] -= force_scalar[4] * dr[4][2]/ dx[4];     
        
        *stress[0] -= force_scalar[4] / dx[4] * dr[4][0] * dr[4][0]; // xx tensor component
        *stress[4] -= force_scalar[4] / dx[4] * dr[4][1] * dr[4][1]; // yy tensor component
        *stress[8] -= force_scalar[4] / dx[4] * dr[4][2] * dr[4][2]; // zz tensor component
        
        *stress[1] -= force_scalar[4] / dx[4] * dr[4][0] * dr[4][1]; // xy tensor component
        *stress[2] -= force_scalar[4] / dx[4] * dr[4][0] * dr[4][2]; // xz tensor component
        *stress[5] -= force_scalar[4] / dx[4] * dr[4][1] * dr[4][2]; // yz tensor component
        
        *stress[3] -= force_scalar[4] / dx[4] * dr[4][0] * dr[4][1]; // yx tensor component
        *stress[6] -= force_scalar[4] / dx[4] * dr[4][0] * dr[4][2]; // zx tensor component
        *stress[7] -= force_scalar[4] / dx[4] * dr[4][1] * dr[4][2]; // zy tensor component 
        
        // Accumulate forces/stresses on/from the kl pair
        
        *force[2][0] += force_scalar[5] * dr[5][0]/ dx[5];
        *force[2][1] += force_scalar[5] * dr[5][1]/ dx[5];
        *force[2][2] += force_scalar[5] * dr[5][2]/ dx[5];

        *force[3][0] -= force_scalar[5] * dr[5][0]/ dx[5];
        *force[3][1] -= force_scalar[5] * dr[5][1]/ dx[5];
        *force[3][2] -= force_scalar[5] * dr[5][2]/ dx[5];     
        
        *stress[0] -= force_scalar[5] / dx[5] * dr[5][0] * dr[5][0]; // xx tensor component
        *stress[4] -= force_scalar[5] / dx[5] * dr[5][1] * dr[5][1]; // yy tensor component
        *stress[8] -= force_scalar[5] / dx[5] * dr[5][2] * dr[5][2]; // zz tensor component
        
        *stress[1] -= force_scalar[5] / dx[5] * dr[5][0] * dr[5][1]; // xy tensor component
        *stress[2] -= force_scalar[5] / dx[5] * dr[5][0] * dr[5][2]; // xz tensor component
        *stress[5] -= force_scalar[5] / dx[5] * dr[5][1] * dr[5][2]; // yz tensor component
        
        *stress[3] -= force_scalar[5] / dx[5] * dr[5][0] * dr[5][1]; // yx tensor component
        *stress[6] -= force_scalar[5] / dx[5] * dr[5][0] * dr[5][2]; // zx tensor component
        *stress[7] -= force_scalar[5] / dx[5] * dr[5][1] * dr[5][2]; // zy tensor component            
    }

    return;
}

double chimesFF::max_cutoff(int ntypes, vector<vector<vector<double> > > & cutoff_list)
{
    double max = cutoff_list[0][1][0]; 
    
    for (int i=0; i<ntypes; i++)
        for (int j=0; j<cutoff_list[i][1].size(); j++)
            if (cutoff_list[i][1][j] > max)
                max = cutoff_list[i][1][j];
    
    return max;

}

double chimesFF::max_cutoff_2B(bool silent)
{
    double max = chimes_2b_cutoff[0][1]; 
    
    for (int i=0; i<chimes_2b_cutoff.size(); i++)
        if (chimes_2b_cutoff[i][1] > max)
            max = chimes_2b_cutoff[i][1];
    
    if ((rank == 0)&&(!silent))        
        cout << "chimesFF: " << "\t" << "Setting 2-body max cutoff to: " << max << endl;
    
    return max;    
}

double chimesFF::max_cutoff_3B(bool silent)
{
    
    if (poly_orders[1] == 0)
        return 0.0;
    
    double max = max_cutoff(chimes_3b_cutoff.size(), chimes_3b_cutoff);
    
    if ((rank == 0)&&(!silent))    
        cout << "chimesFF: " << "\t" << "Setting 3-body max cutoff to: " << max << endl;
    
    return max;
    
}

double chimesFF::max_cutoff_4B(bool silent)
{
    if (poly_orders[2] == 0)
        return 0.0;
    
    double max =  max_cutoff(chimes_4b_cutoff.size(), chimes_4b_cutoff);
        
    if ((rank == 0)&&(!silent))    
        cout << "chimesFF: " << "\t" << "Setting 4-body max cutoff to: " << max << endl;
    
    return max;
}

void chimesFF::set_atomtypes(vector<string> & type_list)
{
    type_list.resize(natmtyps);
    
    for(int i=0;i<natmtyps;i++)
        type_list[i] = atmtyps[i];
}
