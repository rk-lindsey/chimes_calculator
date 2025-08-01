/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>
#include<fstream>

using namespace std;

#include "serial_chimes_interface.h"

// Simple linear algebra functions (for arbitrary triclinic cell support)

double mag_a    (const vector<double> & a)
{
    double mag = 0;
    
    for(int i=0; i<a.size(); i++)
        mag += a[i]*a[i];
    
    mag = sqrt(mag);
    
    return mag;
}
void   unit_a   (const vector<double> & a, vector<double> & unit)        
{
    double mag = mag_a(a);

    unit.resize(a.size());
    
    for(int i=0; i<a.size(); i++)
        unit[i] += a[i]/mag;
    
    return;
}
double a_dot_b  (const vector<double> & a, const vector<double> & b)
{
    double dot = 0;
    
    if (a.size() != b.size())
    {
        cout << "ERROR in a_dot_b: Vectors of different length!" << endl;
        exit(0);
    }
    
    for(int i=0; i<a.size(); i++)
        dot += a[i]*b[i];
    
    return dot;
}
double angle_ab (const vector<double> & a, const vector<double> & b)
{
    double ang = a_dot_b(a,b);

    ang /= mag_a(a);
    ang /= mag_a(b);

    return acos(ang);    
}        
void   a_cross_b(const vector<double> & a, const vector<double> & b, vector<double> & cross)
{
    if( a.size() != b.size())
    {
        cout << "ERROR in a_cross_b: Vectors of different length!" << endl;
        exit(0);
    }
    if(a.size() != 3)
    {
        cout << "ERROR in a_cross_b: Vectors should be of length 3!" << endl;
        exit(0);
    }
    
    cross.resize(3);
    cross[0] =    (a[1]*b[2] - a[2]*b[1]);
    cross[1] = -1*(a[0]*b[2] - a[2]*b[0]);
    cross[2] =    (a[0]*b[1] - a[1]*b[0]);
    return;
}    
void set_hmat(const vector<double> & cell_a,const vector<double> & cell_b, const vector<double> & cell_c, vector<double> & hmat, vector<double> & invr_hmat, int replicates)
{
    // Define the h-matrix (stores the cell vectors locally)

    hmat[0] = cell_a[0]*(replicates+1);  hmat[3] = cell_a[1]*(replicates+1); hmat[6] = cell_a[2]*(replicates+1);
    hmat[1] = cell_b[0]*(replicates+1);  hmat[4] = cell_b[1]*(replicates+1); hmat[7] = cell_b[2]*(replicates+1);
    hmat[2] = cell_c[0]*(replicates+1);  hmat[5] = cell_c[1]*(replicates+1); hmat[8] = cell_c[2]*(replicates+1);

    // Determine the h-matrix inverse
    
    double hmat_det = hmat[0] * (hmat[4]*hmat[8] - hmat[5]*hmat[7])
                    - hmat[1] * (hmat[3]*hmat[8] - hmat[5]*hmat[6])
                    + hmat[2] * (hmat[3]*hmat[7] - hmat[4]*hmat[6]);

    vector<double> tmp_vec(9);
    
    tmp_vec[0] =      (hmat[4]*hmat[8] - hmat[5]*hmat[7]); tmp_vec[3] = -1 * (hmat[1]*hmat[8] - hmat[2]*hmat[7]); tmp_vec[6] =      (hmat[1]*hmat[5] - hmat[2]*hmat[4]);
    tmp_vec[1] = -1 * (hmat[3]*hmat[8] - hmat[5]*hmat[6]); tmp_vec[4] =      (hmat[0]*hmat[8] - hmat[2]*hmat[6]); tmp_vec[7] = -1 * (hmat[0]*hmat[5] - hmat[2]*hmat[3]);
    tmp_vec[2] =      (hmat[3]*hmat[7] - hmat[4]*hmat[6]); tmp_vec[5] = -1 * (hmat[0]*hmat[7] - hmat[1]*hmat[6]); tmp_vec[8] =      (hmat[0]*hmat[4] - hmat[1]*hmat[3]);
    
    invr_hmat[0] = tmp_vec[0]; invr_hmat[3] = tmp_vec[1]; invr_hmat[6] = tmp_vec[2];
    invr_hmat[1] = tmp_vec[3]; invr_hmat[4] = tmp_vec[4]; invr_hmat[7] = tmp_vec[5];
    invr_hmat[2] = tmp_vec[6]; invr_hmat[5] = tmp_vec[7]; invr_hmat[8] = tmp_vec[8];

    invr_hmat[0] /= hmat_det; invr_hmat[3] /= hmat_det; invr_hmat[6] /= hmat_det;
    invr_hmat[1] /= hmat_det; invr_hmat[4] /= hmat_det; invr_hmat[7] /= hmat_det;
    invr_hmat[2] /= hmat_det; invr_hmat[5] /= hmat_det; invr_hmat[8] /= hmat_det;
    
}
    
    
// simulation_system member functions

simulation_system::simulation_system()
{
    hmat     .resize(9);
    invr_hmat.resize(9);
}
simulation_system::~simulation_system()
{}
void simulation_system::set_atomtyp_indices(vector<string> & type_list)
{
    sys_atmtyp_indices.resize(n_ghost);

    for(int i=0; i<n_ghost; i++)
    {
        sys_atmtyp_indices[i] = -1;
        
        for (int j=0; j<type_list.size(); j++)
        {
            if ( sys_atmtyps[i] == type_list[j])
            {
                sys_atmtyp_indices[i] = j;
                break;
            }
        }
        
        if (sys_atmtyp_indices[i] == -1)
        {
            cout << "ERROR: Couldn't assign an atom type index for (index/type) " << sys_parent[i] << " " << sys_atmtyps[i] << endl;
            exit(0);
        }
    }
}

void simulation_system::init(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, double max_2b_cut, bool small)
{
    allow_replication = small;
    max_cut = max_2b_cut;
    static bool called_before = false;
    
    //////////////////////////////////////////
    // STEP 1: Copy the system
    //////////////////////////////////////////
    
    n_atoms = x_in.size();
    
    // Sanity checks
    
    if (n_atoms != y_in.size())
    {
        cout << "ERROR: x and y coordinate vector lengths do not match!" << endl;
        exit(0);
    }
    if (n_atoms != z_in.size())
    {
        cout << "ERROR: x and z coordinate vector lengths do not match!" << endl;
        exit(0);
    }
    
    // Copy over the system
    
    n_ghost = n_atoms;
    n_repl  = n_atoms;
    
    sys_atmtyp_indices .resize(0);
    
    sys_x              .resize(0);
    sys_y              .resize(0);
    sys_z              .resize(0);
    sys_parent         .resize(0);
    sys_rep_parent     .resize(0);
    atmtyps            .erase(atmtyps.begin() + n_atoms, atmtyps.end());
    sys_atmtyps        .resize(0);
    
    for (int a=0; a<n_atoms; a++)
    {
        sys_atmtyps.push_back(atmtyps[a]);    

        sys_x.push_back( x_in[a] );
        sys_y.push_back( y_in[a] );
        sys_z.push_back( z_in[a] );
        
        sys_parent.push_back(a); // for ghost
        sys_rep_parent.push_back(a);    // for replicates
    }  
    
    // Determine if system is large enough  

    latcon_a = mag_a({cella_in[0], cella_in[1], cella_in[2]});
    latcon_b = mag_a({cellb_in[0], cellb_in[1], cellb_in[2]});
    latcon_c = mag_a({cellc_in[0], cellc_in[1], cellc_in[2]});

    double min_latcon = latcon_a;

    if (latcon_b < min_latcon)
        min_latcon = latcon_b;
    if (latcon_c < min_latcon)
        min_latcon = latcon_c;
    
    n_replicates = 0;
    
    if (allow_replication)
        n_replicates = ceil(max_cut/min_latcon)-1;

    if (n_replicates > 0)
    {
            if (!called_before)
            {
                called_before = true;

                // Avoid printing the number of replicates on each MD step (LEF).
                cout << "SerialchimesFF: " << "Replicating the system " << n_replicates << " times prior to generating ghost atoms" << endl;    
    
                cout << "SerialchimesFF: " << "\t" << "Warning: At least one cell length is smaller than the ChIMES outer cutoff." << endl;
                cout << "SerialchimesFF: " << "\t" << "System will be replicated prior to ghost atom generation." << endl;
                cout << "SerialchimesFF: " << "\t" << "Results will only be correct for perfectly crystalline cells." << endl;
                cout << "SerialchimesFF: " << "\t" << "For any other case, system size should be increased." << endl;
            }
    }

    set_hmat(cella_in, cellb_in, cellc_in, hmat, invr_hmat, 0);
    
    // Build the replicates

    double tmp_x, tmp_y, tmp_z;


    for (int i=0; i<=n_replicates; i++) // x
    {
        for (int j=0; j<=n_replicates; j++) // y
        {
            for (int k=0; k<=n_replicates; k++) // z
            {                
                if ((i==0)&&(j==0)&&(k==0))
                   continue;
                
                for (int a=0; a<n_atoms; a++)
                {
                    n_ghost++;
                    n_repl++;    
                    
                    atmtyps.push_back(atmtyps[a]);
                    sys_atmtyps.push_back(atmtyps[a]);    
                    
                    sys_x.push_back(0.0); // Holder    
                    sys_y.push_back(0.0);
                    sys_z.push_back(0.0);
                    
                    // Transform into inverse space 
                    
                    tmp_x = invr_hmat[0]*sys_x[a] + invr_hmat[1]*sys_y[a] + invr_hmat[2]*sys_z[a];
                    tmp_y = invr_hmat[3]*sys_x[a] + invr_hmat[4]*sys_y[a] + invr_hmat[5]*sys_z[a];
                    tmp_z = invr_hmat[6]*sys_x[a] + invr_hmat[7]*sys_y[a] + invr_hmat[8]*sys_z[a];
                    
                    tmp_x += i;
                    tmp_y += j;    
                    tmp_z += k;
                    
                    sys_x[n_ghost-1] = hmat[0]*tmp_x + hmat[1]*tmp_y + hmat[2]*tmp_z;
                    sys_y[n_ghost-1] = hmat[3]*tmp_x + hmat[4]*tmp_y + hmat[5]*tmp_z;
                    sys_z[n_ghost-1] = hmat[6]*tmp_x + hmat[7]*tmp_y + hmat[8]*tmp_z;    
                    
                    sys_parent.push_back(n_repl-1); // As far as ghosts are concerned, these are real atoms
                    sys_rep_parent.push_back(a);    // replicates know they have a parent
                }
            }
        }
    }

    // Note on replicates: We're essentially tricking the code into thinking the system is bigger than it is, 
    // before any ghost atoms are built 

    n_atoms = n_repl;

    set_hmat(cella_in, cellb_in, cellc_in, hmat, invr_hmat, n_replicates);

    //////////////////////////////////////////
    // STEP 2: Wrap atoms
    //////////////////////////////////////////
        
    // Wrap atoms by converting to inverse space (orthorhombic, cell lengths = unity)
    // Leave in fractional coordinates since we'll need to transform to the new basis later
        
    double tmp_ax, tmp_ay, tmp_az;

    for(int i=0; i<n_atoms; i++)
    {
        // Convert to fractional coordinates
        
        tmp_ax = invr_hmat[0]*sys_x[i] + invr_hmat[1]*sys_y[i] + invr_hmat[2]*sys_z[i];
        tmp_ay = invr_hmat[3]*sys_x[i] + invr_hmat[4]*sys_y[i] + invr_hmat[5]*sys_z[i];
        tmp_az = invr_hmat[6]*sys_x[i] + invr_hmat[7]*sys_y[i] + invr_hmat[8]*sys_z[i];

        // Wrap

        tmp_ax -= floor(tmp_ax);
        tmp_ay -= floor(tmp_ay);
        tmp_az -= floor(tmp_az);
        
        // Revert to absolute coordinates
        
        sys_x[i] = hmat[0]*tmp_ax + hmat[1]*tmp_ay + hmat[2]*tmp_az;    
        sys_y[i] = hmat[3]*tmp_ax + hmat[4]*tmp_ay + hmat[5]*tmp_az;
        sys_z[i] = hmat[6]*tmp_ax + hmat[7]*tmp_ay + hmat[8]*tmp_az;        
    }
    
    //////////////////////////////////////////
    // STEP 3: Determine cell volume 
    //////////////////////////////////////////    
    
    latcon_a = mag_a({hmat[0], hmat[3], hmat[6]});
    latcon_b = mag_a({hmat[1], hmat[4], hmat[7]});
    latcon_c = mag_a({hmat[2], hmat[5], hmat[8]});

    cell_alpha = angle_ab({hmat[1], hmat[4], hmat[7]}, {hmat[2], hmat[5], hmat[8]});
    cell_beta  = angle_ab({hmat[2], hmat[5], hmat[8]}, {hmat[0], hmat[3], hmat[6]});
    cell_gamma = angle_ab({hmat[0], hmat[3], hmat[6]}, {hmat[1], hmat[4], hmat[7]});
    
    vol  = 1;
    vol += 2*cos(cell_alpha)*cos(cell_beta )*cos(cell_gamma);
    vol -=   cos(cell_alpha)*cos(cell_alpha);
    vol -=   cos(cell_beta )*cos(cell_beta );
    vol -=   cos(cell_gamma)*cos(cell_gamma);

    vol = latcon_a * latcon_b * latcon_c * sqrt(vol);   
}
void simulation_system::reorient()
{
    //////////////////////////////////////////
    //  Manipulate atoms and cell to conform with LAMMPS convention
    //////////////////////////////////////////
    // Assumes the cell has already been wrapped
    // This will make linear-scaling neighbor list construction easier down the line
    // See: https://lammps.sandia.gov/doc/Howto_triclinic.html
    
    // Convert to fractional coordinates from the original basis

    double tmp_ax, tmp_ay, tmp_az;
    
    for(int i=0; i<n_atoms; i++)
    {
        tmp_ax = invr_hmat[0]*sys_x[i] + invr_hmat[1]*sys_y[i] + invr_hmat[2]*sys_z[i];
        tmp_ay = invr_hmat[3]*sys_x[i] + invr_hmat[4]*sys_y[i] + invr_hmat[5]*sys_z[i];
        tmp_az = invr_hmat[6]*sys_x[i] + invr_hmat[7]*sys_y[i] + invr_hmat[8]*sys_z[i];
        
        sys_x[i] = tmp_ax;
        sys_y[i] = tmp_ay;
        sys_z[i] = tmp_az;
    }  
    
    
    // Rotate the cell
    
    vector<double> tmp_cella(3);
    vector<double> tmp_cellb(3);
    vector<double> tmp_cellc(3);
    vector<double> tmp_unit (3);
    vector<double> tmp_cross(3);
    
    vector<double> cella_in = {hmat[0], hmat[3], hmat[6]};
    vector<double> cellb_in = {hmat[1], hmat[4], hmat[7]};
    vector<double> cellc_in = {hmat[2], hmat[5], hmat[8]};
    
    
    unit_a   (cella_in, tmp_unit);
    a_cross_b(tmp_unit, cellb_in, tmp_cross);
    
    // Determine cell vectors for the rotated system
    
    tmp_cella[0] = mag_a(cella_in);
    tmp_cella[1] = 0;
    tmp_cella[2] = 0;

    tmp_cellb[0] = a_dot_b(cellb_in, tmp_unit);
    tmp_cellb[1] = mag_a(tmp_cross);
    tmp_cellb[2] = 0;

    tmp_cellc[0] = a_dot_b(cellc_in,tmp_unit);
    tmp_cellc[1] = (a_dot_b(cellb_in,cellc_in) - tmp_cellb[0]*tmp_cellc[0]) / tmp_cellb[1];
    tmp_cellc[2] = sqrt( mag_a(cellc_in)*mag_a(cellc_in) -  tmp_cellc[0]*tmp_cellc[0] - tmp_cellc[1]*tmp_cellc[1]);
            
    // Determine the new cell h-matrix and its inverse

    set_hmat({tmp_cella[0], tmp_cella[1], tmp_cella[2]}, {tmp_cellb[0], tmp_cellb[1], tmp_cellb[2]}, {tmp_cellc[0], tmp_cellc[1], tmp_cellc[2]}, hmat, invr_hmat, 0);

    // Transform to the new nominally rotated cell  

    for(int i=0; i<n_atoms; i++)
    {
        tmp_ax = sys_x[i];
        tmp_ay = sys_y[i];
        tmp_az = sys_z[i];
        
        sys_x[i] = hmat[0]*tmp_ax + hmat[1]*tmp_ay + hmat[2]*tmp_az;    
        sys_y[i] = hmat[3]*tmp_ax + hmat[4]*tmp_ay + hmat[5]*tmp_az;
        sys_z[i] = hmat[6]*tmp_ax + hmat[7]*tmp_ay + hmat[8]*tmp_az;
    }   

    // Determine cell extents and volume for neighbor list constructions
     
    latcon_a = mag_a({hmat[0], hmat[3], hmat[6]});
    latcon_b = mag_a({hmat[1], hmat[4], hmat[7]});
    latcon_c = mag_a({hmat[2], hmat[5], hmat[8]});

    cell_alpha = angle_ab({hmat[1], hmat[4], hmat[7]}, {hmat[2], hmat[5], hmat[8]});
    cell_beta  = angle_ab({hmat[2], hmat[5], hmat[8]}, {hmat[0], hmat[3], hmat[6]});
    cell_gamma = angle_ab({hmat[0], hmat[3], hmat[6]}, {hmat[1], hmat[4], hmat[7]});

    double cell_lx = latcon_a;
    double xy = latcon_b * cos(cell_gamma);
    if(abs(xy) < 1E-12)
        xy = 0.0;

    double xz = latcon_c * cos(cell_beta );
    if(abs(xz) < 1E-12)
        xz = 0.0;

    double cell_ly = sqrt(latcon_b *latcon_b - xy*xy);

    double yz = (latcon_b*latcon_c * cos(cell_alpha) -xy*xz)/cell_ly;
    if(abs(yz) < 1E-12)
        yz = 0.0;

    double cell_lz = sqrt( latcon_c*latcon_c - xz*xz -yz*yz);
    double tmp;
        
    double xlo = 0.0;
    if (xy < xlo)
        xlo = xy;
    if (xz < xlo)
        xlo = xz;
    if(xy+xz < xlo)
        xlo = xy+xz;
        
    double ylo = 0.0;
    if (yz< ylo)
        ylo = yz;
    double zlo = 0.0;
    
    double xhi = hmat[0];
    tmp = 0.0;
    if (xy > tmp)
        tmp = xy;
    if (xz > tmp)
        tmp = xz;
    if (xy+xz> tmp)
        tmp = xy+xz;
    xhi += tmp;
    
    double yhi = hmat[4];
    if(yz > 0.0)
        yhi += yz;
        
    double zhi = hmat[8];
    
    extent_x = xhi - xlo;
    extent_y = yhi - ylo;
    extent_z = zhi - zlo; 


}
void simulation_system::build_layered_system(vector<string> & atmtyps, vector<int> & poly_orders, double max_2b_cut, double max_3b_cut, double max_4b_cut)
{
    
    // use smallest lattice length to determine number of ghost atom layers (n_layers)
     
    vector<double> latdist = {latcon_a,latcon_b,latcon_c};
    
    double lat_min = *min_element(latdist.begin(),latdist.end());
    
    double eff_length = max_2b_cut*2.0;
    
    // n_layers is set to ensure that max 2b rcut is less than half smallest box length
     
    n_layers = ceil(eff_length/lat_min+1);
    
    double eff_lx = latcon_a * (2*n_layers + 1);
    double eff_ly = latcon_b * (2*n_layers + 1);
    double eff_lz = latcon_c * (2*n_layers + 1);
    
    if ((max_2b_cut>0.5*eff_lx)||(max_2b_cut>0.5*eff_ly)||(max_2b_cut>0.5*eff_lz))
    {
        cout << "ERROR: Maximum 2b cutoff is greater than half at least one box length." << endl;
        cout << "       Increase requested n_layers." << endl;
        cout << "       Max 2b cutoff:            " << max_2b_cut << endl;
        cout << "       Effective cell length(x): " << eff_lx << endl;
        cout << "       Effective cell length(y): " << eff_ly << endl;
        cout << "       Effective cell length(z): " << eff_lz << endl;
        cout << "       nlayers:                  " << n_layers << endl;
        exit(0);
    }        
    if (poly_orders[1] >0)
    {
        if ((max_3b_cut>0.5*eff_lx)||(max_3b_cut>0.5*eff_ly)||(max_3b_cut>0.5*eff_lz))
        {
            cout << "ERROR: Maximum 3b cutoff is greater than half at least one box length." << endl;
            cout << "       Increase requested n_layers." << endl;
            cout << "       Max 3b cutoff:            " << max_3b_cut << endl;
            cout << "       Effective cell length(x): " << eff_lx << endl;
            cout << "       Effective cell length(y): " << eff_ly << endl;
            cout << "       Effective cell length(z): " << eff_lz << endl;
            cout << "       nlayers:                  " << n_layers << endl;
            exit(0);
        }
    }
    if (poly_orders[2] >0)
    {
        if ((max_4b_cut>0.5*eff_lx)||(max_4b_cut>0.5*eff_ly)||(max_4b_cut>0.5*eff_lz))
        {
            cout << "ERROR: Maximum 4b cutoff is greater than half at least one box length." << endl;
            cout << "       Increase requested n_layers." << endl;
            cout << "       Max 4b cutoff:            " << max_4b_cut << endl;
            cout << "       Effective cell length(x): " << eff_lx << endl;
            cout << "       Effective cell length(y): " << eff_ly << endl;
            cout << "       Effective cell length(z): " << eff_lz << endl;
            cout << "       nlayers:                  " << n_layers << endl;
            exit(0);
        }
    }

    // Build the layers. First clear out the associated vectors so you don't run into memory issues on multiple calls:
  
    sys_atmtyps.erase(sys_atmtyps.begin() + n_atoms, sys_atmtyps.end());
    sys_x      .erase(sys_x      .begin() + n_atoms, sys_x      .end());
    sys_y      .erase(sys_y      .begin() + n_atoms, sys_y      .end());
    sys_z      .erase(sys_z      .begin() + n_atoms, sys_z      .end());
    sys_parent .erase(sys_parent .begin() + n_atoms, sys_parent .end());
 
    double tmp_x, tmp_y, tmp_z;

    for (int i=-n_layers; i<=n_layers; i++) // x
    {
        for (int j=-n_layers; j<=n_layers; j++) // y
        {
            for (int k=-n_layers; k<=n_layers; k++) // z
            {                
                if ((i==0)&&(j==0)&&(k==0))
                    continue;
                    
                for (int a=0; a<n_atoms; a++)
                {
                    n_ghost++;    
                    
                    sys_atmtyps.push_back(atmtyps[a]);    
                    
                    sys_x.push_back(0.0); // Holder    
                    sys_y.push_back(0.0);
                    sys_z.push_back(0.0);
                    
                    // Transform into inverse space 

                    tmp_x = invr_hmat[0]*sys_x[a] + invr_hmat[1]*sys_y[a] + invr_hmat[2]*sys_z[a];
                    tmp_y = invr_hmat[3]*sys_x[a] + invr_hmat[4]*sys_y[a] + invr_hmat[5]*sys_z[a];
                    tmp_z = invr_hmat[6]*sys_x[a] + invr_hmat[7]*sys_y[a] + invr_hmat[8]*sys_z[a];
                    
                    tmp_x += i;
                    tmp_y += j;    
                    tmp_z += k;

                    sys_x[n_ghost-1] = hmat[0]*tmp_x + hmat[1]*tmp_y + hmat[2]*tmp_z;
                    sys_y[n_ghost-1] = hmat[3]*tmp_x + hmat[4]*tmp_y + hmat[5]*tmp_z;
                    sys_z[n_ghost-1] = hmat[6]*tmp_x + hmat[7]*tmp_y + hmat[8]*tmp_z;    

                    sys_parent.push_back(a);
                    

                }
            }
        }
    }
    
}
void simulation_system::build_neigh_lists(vector<int> & poly_orders, vector<vector<int> > & neighlist_2b, vector<vector<int> > & neighlist_3b, vector<vector<int> > & neighlist_4b, double max_2b_cut, double max_3b_cut, double max_4b_cut)
{
    vector<double> maxpos(3, -1.0e100) ;
    vector<double> minpos(3, +1.0e100) ;

    // Determine limits on position of all particles.
    for(int i=0; i<n_ghost; i++)
    {
        if ( sys_x[i] > maxpos[0] )
            maxpos[0] = sys_x[i] ;
        if ( sys_y[i] > maxpos[1] )
            maxpos[1] = sys_y[i] ;
        if ( sys_z[i] > maxpos[2] )
            maxpos[2] = sys_z[i] ;

        if ( sys_x[i] < minpos[0] )
            minpos[0] = sys_x[i] ;
        if ( sys_y[i] < minpos[1] )
            minpos[1] = sys_y[i] ;
        if ( sys_z[i] < minpos[2] )
            minpos[2] = sys_z[i] ;

    }

    // Make the 2b neighbor lists
    
    neighlist_2b.resize(n_ghost);
    
    for (int i = 0; i < n_ghost; i++) 
        neighlist_2b[i].resize(0,0);
    
    neighlist_3b.resize(0);
    neighlist_4b.resize(0);

    // Determine search distances

    double search_dist = max_2b_cut;
    
    if (max_3b_cut > search_dist)
        search_dist = max_3b_cut;
    if (max_4b_cut > search_dist)
        search_dist = max_4b_cut;    
    
    // Prepare bins

    int nbins_x = ceil((maxpos[0]-minpos[0])/search_dist) ;
    int nbins_y = ceil((maxpos[1]-minpos[1])/search_dist) ;
    int nbins_z = ceil((maxpos[2]-minpos[2])/search_dist) ;
    if ( (nbins_x < 3) || (nbins_y < 3) || (nbins_z < 3) )
    {
       cout << "Error: require at least 3 neighbor bins in all directions.\n" ;
       cout << "The number of layers is not correct\n" ;
    }
    
    int total_bins = nbins_x * nbins_y * nbins_z;
    
    vector<vector<int> > bin(total_bins);

    for (int i=0; i<total_bins; i++) 
        vector<int>().swap(bin[i]);
    
    int bin_x_idx, bin_y_idx, bin_z_idx, ibin;

    // Populate bins    

    for(int i=0; i<n_ghost; i++)
    {
        bin_x_idx = floor( (sys_x[i] - minpos[0] ) / search_dist );
        bin_y_idx = floor( (sys_y[i] - minpos[1] ) / search_dist );
        bin_z_idx = floor( (sys_z[i] - minpos[2] ) / search_dist );
        
        if ( bin_x_idx < 0 || bin_y_idx < 0 || bin_z_idx < 0 ) 
        {
            cout << "ERROR: Negative bin computed" << endl; // Check .xyz box lengths  and atom coords
            exit(0);
        }
        
        // Calculate bin BIN_IDX of the atom.
        int ibin = bin_x_idx + bin_y_idx * nbins_x + bin_z_idx * nbins_x* nbins_y;
        
        if ( ibin >= total_bins ) 
        {
            cout << "Error: ibin out of range\n";
            cout << ibin << " " << total_bins << endl;
            exit(1);
        }

        // Push the atom into the bin 
        bin[ibin].push_back(i);
        
    }

    // Generate neighbor lists on basis of bins
    
    for(int ai=0; ai<n_atoms; ai++)
    {
        bin_x_idx = floor( (sys_x[ai] - minpos[0]) / search_dist ) ;
        bin_y_idx = floor( (sys_y[ai] - minpos[1]) / search_dist ) ;
        bin_z_idx = floor( (sys_z[ai] - minpos[2]) / search_dist ) ;

        if ( bin_x_idx < 1 || bin_y_idx < 1 || bin_z_idx < 1 ) 
        {
            cout << "ERROR: Bad bin computed" << endl; // Check .xyz box lengths  and atom coords
            cout << bin_x_idx << " " << bin_y_idx << " " << bin_z_idx << endl;
            cout << sys_x[ai] << " " << sys_y[ai] << " " << sys_z[ai] << endl;
            exit(0);
        }

        // Loop over relevant bins only, not all atoms.
        
        int ibin, aj, ajj, ajend;
        
        for (int i=bin_x_idx-1; i<= bin_x_idx+1; i++)    //BIN_IDX_a1.X 
        {
            for (int j=bin_y_idx-1; j<=bin_y_idx+1; j++ ) //BIN_IDX_a1.Y
            {
                for (int k=bin_z_idx-1; k<=bin_z_idx+1; k++ ) // BIN_IDX_a1.Z
                {
                    ibin = i + j * nbins_x + k * nbins_x * nbins_y;

                    if (ibin >= total_bins) 
                    {
                        cout << "Error: binning BIN_IDX out of range\n";
                        cout << "BIN_IDX.X = " << i << "BIN_IDX.Y = " << j << "BIN_IDX.Z = " << k << endl;
                        exit(1);
                    }

                    ajend = bin[ibin].size();
                    
                    for (int aj=0; aj<ajend; aj++) 
                    {
                        ajj = bin[ibin][aj];

                        if ( ajj == ai ) 
                            continue;

                        if ( ai <= sys_parent[ajj]) 
                            if (get_dist(ai,ajj) < search_dist )
                                neighlist_2b[ai].push_back(ajj);        
                    }
                }
            }
        }
    }    

    if ((poly_orders[1] == 0)&&(poly_orders[2]==0))
        return;    
    
    // Make the 3- and 4-b neighbor lists

    bool valid_3mer;
    bool valid_4mer;
    
    vector<int> tmp_3mer(3);
    vector<int> tmp_4mer(4);
    
    int jj, kk, ll;
    
    double dist;
    double dist_ij, dist_ik, dist_il, dist_jk, dist_jl, dist_kl;
    
    for(int i=0; i<n_atoms; i++)
    {
        for(int j=0; j<neighlist_2b[i].size(); j++) // Neighbors of i
        {
            valid_3mer = true;
            valid_4mer = true;
            
            jj = neighlist_2b[i][j];
            
            dist_ij = get_dist(i,jj);

            if (dist_ij >= max_3b_cut)
                valid_3mer = false;
            if (dist_ij >= max_4b_cut)
                valid_4mer = false;
                
            if (!valid_3mer && !valid_4mer)
                continue;
            
            for(int k=0; k<neighlist_2b[i].size(); k++)
            {                                
                kk = neighlist_2b[i][k];
                
                if (jj == kk)
                    continue;
                if (sys_parent[jj] > sys_parent[kk])
                    continue;
                
                dist_ik = get_dist(i,kk); // Check i/k distance
            
                if (dist_ik >= max_3b_cut)
                    valid_3mer = false;
                if (dist_ik >= max_4b_cut)
                    valid_4mer = false;
                    
                dist_jk = get_dist(jj,kk); // Check j/k distance
        
                if (dist_jk >= max_3b_cut)
                    valid_3mer = false;
                if (dist_jk >= max_4b_cut)
                    valid_4mer = false;                                        
                
                if (!valid_3mer && !valid_4mer)
                {
                    if(dist_ij<max_3b_cut)
                        valid_3mer = true;
                    
                    if(dist_ij<max_4b_cut)
                        valid_4mer = true;
                
                    continue;
                }
                
                // If we're here then we have a valid 3-mer ... add it to the 3b neighbor list    
                
                tmp_3mer[0] = i;
                tmp_3mer[1] = jj;
                tmp_3mer[2] = kk;
                
                if (valid_3mer)
                    neighlist_3b.push_back(tmp_3mer);
                
                // Continue on to 4-body list
                
                if (poly_orders[2] == 0)
                    continue;
                
                if (!valid_4mer)
                {
                    if(dist_ij<max_4b_cut)
                        valid_4mer = true;
                    continue;
                }
                
                for(int l=0; l<neighlist_2b[i].size(); l++) 
                {                                
                    ll = neighlist_2b[i][l];
                
                    if (jj == ll)
                        continue;
                    if (kk == ll)
                        continue;
                    if (sys_parent[jj] > sys_parent[ll])
                        continue;                
                    if (sys_parent[kk] > sys_parent[ll])
                        continue;                

                    if (get_dist(i ,ll) >= max_4b_cut) // Check i/l distance
                        continue;                

                    if (get_dist(jj,ll) >= max_4b_cut) // Check j/l distance
                        continue;    

                    if (get_dist(kk,ll) >= max_4b_cut) // Check k/l distance
                        continue;        
                    
                    // If we're here then we have a valid 4-mer ... add it to the 4b neighbor list    
                
                    tmp_4mer[0] = i;
                    tmp_4mer[1] = jj;
                    tmp_4mer[2] = kk;
                    tmp_4mer[3] = ll;
                
                    if (valid_4mer)
                        neighlist_4b.push_back(tmp_4mer);        
                }                                                                        
            }
        }
    }
    /*
    cout << "2B neighbor list is of length:" << neighlist_2b.size() << endl;
    cout << "3B neighbor list is of length:" << neighlist_3b.size() << endl;
    cout << "4B neighbor list is of length:" << neighlist_4b.size() << endl;
    */

}

void simulation_system::run_checks(const vector<double>& max_cuts, vector<int>&poly_orders)
{
    // Sanity check 1: Are the cell vectors long enough?
    
    for(int i=0;i<max_cuts.size(); i++)
    {
        if ( 
            (max_cuts[i] > 2*latcon_a * (2*n_layers + 1)) || 
            (max_cuts[i] > 2*latcon_b * (2*n_layers + 1)) ||
            (max_cuts[i] > 2*latcon_c * (2*n_layers + 1))
            )
        {
            cout << "ERROR: Layered system is smaller than 2x the model " << i+2 <<"-body maximum outer cutoff." << endl;
            cout << "Please report this error to the developers." << endl;
            cout << "Model maximum cutoff: " << max_cuts[i] << endl;
            cout << "Layered system lattice constant (a): " << latcon_a * (2*n_layers + 1) << endl;
            cout << "Layered system lattice constant (b): " << latcon_b * (2*n_layers + 1) << endl;
            cout << "Layered system lattice constant (c): " << latcon_c * (2*n_layers + 1) << endl;
            exit(0);
            
        }
    }

    
    // Sanity check 2: Does the system have enough atoms?
    
    int bodiedness = 2;
    if (poly_orders[1] > 0)
        bodiedness++;
    if (poly_orders[2] > 0)
        bodiedness++;
    
    if (bodiedness > sys_x.size())
    {
        cout << "ERROR: Layered system contains too few atoms." << endl;
        cout << "   Model bodiedness:            " << bodiedness << endl;
        cout << "   No. atoms in layered system: " << sys_x.size() << endl;
        exit(0);
    }
}
    
// serial_chimes_interface member functions

serial_chimes_interface::serial_chimes_interface(bool small)
{
    // For small systems, allow explicit replication prior to ghost atom construction
    // This should ONLY be done for perfectly crystalline systems
    
    allow_replication = small;
    
    // Initialize Pointers, etc for chimes calculator interfacing (2-body only for now)
    // To set up for many body calculations, see the LAMMPS implementation

    dist_3b.resize(3);
    dist_4b.resize(6);
    
    dr   .resize(3);
    dr_3b.resize(3*3) ;
    dr_4b.resize(6*3);
    
    force_ptr_2b.resize(2,std::vector<double*>(3));
    
    typ_idxs_2b.resize(2);
    typ_idxs_3b.resize(3);
    typ_idxs_4b.resize(4);
    
    max_2b_cut = 0.0;
    max_3b_cut = 0.0;
    max_4b_cut = 0.0;
    
}
serial_chimes_interface::~serial_chimes_interface()
{}

void serial_chimes_interface::init_chimesFF(string chimesFF_paramfile, int rank)
{
    // Initialize the chimesFF object, read parameters
    
    init(rank);
    read_parameters(chimesFF_paramfile);
    set_atomtypes(type_list);
    build_pair_int_trip_map() ; 
    build_pair_int_quad_map() ;
}

void serial_chimes_interface::build_neigh_lists(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in)
{
    neigh.init(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in, max_cutoff_2B(true), allow_replication);
    neigh.reorient();
    neigh.build_layered_system(atmtyps, poly_orders, max_cutoff_2B(true), max_cutoff_3B(true), max_cutoff_4B(true));
    neigh.set_atomtyp_indices(type_list);
    neigh.build_neigh_lists(poly_orders, neighlist_2b, neighlist_3b, neighlist_4b, max_cutoff_2B(true), max_cutoff_3B(true), max_cutoff_4B(true));
}

void serial_chimes_interface::calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<vector<double> > & force, vector<double> & stress)
{   
    // Read system, set up lattice constants/hmats

    // Determine the max outer cutoff (MUST be 2-body, based on ChIMES logic)

    // Initialize private members (LEF)
    max_2b_cut = max_cutoff_2B(true) ;
    max_3b_cut = max_cutoff_3B(true) ;
    max_4b_cut = max_cutoff_4B(true) ;

    vector<double> stress_chimes(6,0.0) ; // Switch Chimes to a packed stressed tensor.
    
    sys.init(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in, max_2b_cut, allow_replication);   
    
    sys.build_layered_system(atmtyps,poly_orders, max_2b_cut, max_3b_cut, max_4b_cut);

    sys.set_atomtyp_indices(type_list);
    
    sys.run_checks({max_2b_cut,max_3b_cut,max_4b_cut},poly_orders);

    build_neigh_lists(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in);

    
    // Setup vars
    
    int ii, jj, kk, ll;
    
    vector<double> force_4b(4*CHDIM) ;
    vector<double> force_3b(3*CHDIM) ;
    vector<double> force_2b(2*CHDIM) ;
    chimes2BTmp chimes_2btmp(poly_orders[0]) ;
    chimes3BTmp chimes_3btmp(poly_orders[1]) ;
    chimes4BTmp chimes_4btmp(poly_orders[2]) ;      
    
    ////////////////////////
    // interate over 1- and 2b's 
    ////////////////////////

    for(int i=0; i<sys.n_atoms; i++)
    {
        compute_1B(sys.sys_atmtyp_indices[i], energy);

        for(int j=0; j<neighlist_2b[i].size(); j++) // Neighbors of i
        {
            jj = neighlist_2b[i][j];

            dist = sys.get_dist(i,jj,dr); // Populates dr, which is passed by ref (overloaded)
            
            typ_idxs_2b[0] = sys.sys_atmtyp_indices[i ];
            typ_idxs_2b[1] = sys.sys_atmtyp_indices[jj];
            
            for (int idx=0; idx<2*CHDIM; idx++) {
                force_2b[idx] = 0.0 ;
            }
            
            compute_2B(dist, dr, typ_idxs_2b, force_2b, stress_chimes, energy, chimes_2btmp);

            for (int idx=0; idx<3; idx++)
            {
                force[sys.sys_rep_parent[i]][idx]                  += force_2b[0*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[jj]]][idx] += force_2b[1*CHDIM+idx] ;     
            }
        }
    }
    
    ////////////////////////
    // interate over 3b's 
    ////////////////////////
    
    if (poly_orders[1] > 0 )
    {
        for(int i=0; i<neighlist_3b.size(); i++)
        {
            ii = neighlist_3b[i][0];
            jj = neighlist_3b[i][1];
            kk = neighlist_3b[i][2];
        
            dist_3b[0] = sys.get_dist(ii,jj,&dr_3b[0]); 
            dist_3b[1] = sys.get_dist(ii,kk,&dr_3b[3]); 
            dist_3b[2] = sys.get_dist(jj,kk,&dr_3b[6]); 
        
            typ_idxs_3b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_3b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_3b[2] = sys.sys_atmtyp_indices[kk];
        
            for (int idx=0; idx<3*CHDIM; idx++)
            {
                force_3b[idx] = 0.0 ;               
            }    
        
            compute_3B(dist_3b, dr_3b, typ_idxs_3b, force_3b, stress_chimes, energy, chimes_3btmp);

            for (int idx=0; idx<3; idx++) {
                force[sys.sys_rep_parent[sys.sys_parent[ii]]][idx] += force_3b[0*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[jj]]][idx] += force_3b[1*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[kk]]][idx] += force_3b[2*CHDIM+idx] ;
            }
        }
    }

    ////////////////////////
    // interate over 4b's 
    ////////////////////////

    if (poly_orders[2] > 0 )
    {
        for(int i=0; i<neighlist_4b.size(); i++)
        {
            ii = neighlist_4b[i][0];
            jj = neighlist_4b[i][1];
            kk = neighlist_4b[i][2];
            ll = neighlist_4b[i][3];
        
            dist_4b[0] = sys.get_dist(ii,jj,&dr_4b[0*CHDIM]); 
            dist_4b[1] = sys.get_dist(ii,kk,&dr_4b[1*CHDIM]); 
            dist_4b[2] = sys.get_dist(ii,ll,&dr_4b[2*CHDIM]); 
            dist_4b[3] = sys.get_dist(jj,kk,&dr_4b[3*CHDIM]); 
            dist_4b[4] = sys.get_dist(jj,ll,&dr_4b[4*CHDIM]); 
            dist_4b[5] = sys.get_dist(kk,ll,&dr_4b[5*CHDIM]);         

            typ_idxs_4b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_4b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_4b[2] = sys.sys_atmtyp_indices[kk];
            typ_idxs_4b[3] = sys.sys_atmtyp_indices[ll];        
        
            for (int idx=0; idx<4*CHDIM; idx++)
            {
                force_4b[idx] = 0.0 ;
            }    
        
            compute_4B(dist_4b, dr_4b, typ_idxs_4b, force_4b, stress_chimes, energy, chimes_4btmp);

            for (int idx=0; idx<3; idx++)
            {
                force[sys.sys_rep_parent[sys.sys_parent[ii]]][idx] += force_4b[0*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[jj]]][idx] += force_4b[1*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[kk]]][idx] += force_4b[2*CHDIM+idx] ;
                force[sys.sys_rep_parent[sys.sys_parent[ll]]][idx] += force_4b[3*CHDIM+idx] ;
            }    
        }    
    }

    // Correct for use of replicates, if applicable
    
    energy /= pow(sys.n_replicates+1.0,3.0);
   
    ////////////////////////
    // Finish pressure calculation
    ////////////////////////

    // Chimes calculates only unique elements of stress tensor.
    stress[0] = stress_chimes[0] ; // xx
    stress[1] = stress_chimes[1] ; // xy
    stress[2] = stress_chimes[2] ; // xz
    stress[3] = stress_chimes[1] ; // yx
    stress[4] = stress_chimes[3] ; // yy
    stress[5] = stress_chimes[4] ; // yz
    stress[6] = stress_chimes[2] ; // zx
    stress[7] = stress_chimes[4] ; // zy
    stress[8] = stress_chimes[5] ; // zz
    
    for (int idx=0; idx<9; idx++)
        stress[idx] /= sys.vol;  
}









