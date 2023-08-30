// program to generate a lammps input file
// for a simple polymer with proteins

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

#define PI 3.14159265358979
//#define DNATYPE 1
//#define PROTTYPE 2

using namespace std;

struct xyz {
  // structure for position coordinates
  int id,
    mol,
    type;
  double x,y,z;
};

struct bnd {
  // structure for bonds
  int one,two,
    type;
};

struct ang {
  // structure for angles
  int one,two,three,
    type;
};

int main() {

  int N,         // number of DNA beads
    Nprot,       // number of proteins
    Ndumb,       // number of dumbell proteins
    prot_flag,   // a flag for if there are proteins
    dumb_flag,   // a flag for if there are dumbell proteins
    lx,ly,lz,    // size of box
    seed,        // seed for random numbers
    id,          // id for atoms
    mol;         // molecule id

  int DNATYPE=1, //  atom types for DNA beads
    PROTTYPE=2,  //  atom types for proteins and dubells
    DUMBTYPEA=3,
    DUMBTYPEB=4;
  
  double theta,  // an angle
    phi;         // an angle

  vector<xyz> atoms;  // positions of atoms

  vector<bnd> bonds;  // list of bonds
  vector<ang> angles; // lust of angles

  double dx,dy,dz;
  xyz new_atom,
    last_atom,
    lastlast_atom;
  bnd new_bond;
  ang new_angle;

  char fn[50];   // file name
  ofstream ouf;  // output stream


  // get parameters
  cout<<"Size of box, x,y,z. Enter three values:"<<endl;
  cin>>lx>>ly>>lz;
  cout<<"Number of DNA beads:"<<endl;
  cin>>N;
  cout<<"Number of sphere proteins (or enter 0):"<<endl;
  cin>>Nprot;
  cout<<"Number of dumbell proteins (or enter 0):"<<endl;
  cin>>Ndumb;
  cout<<"Enter seed for random numbers."<<endl;
  cin>>seed;

  if (Nprot>0) {prot_flag=1;}
  else {
    prot_flag=0;
    DUMBTYPEA--; // if there are no sphere protein, set dumbells as type 2 and 3
    DUMBTYPEB--;
  }

  if (Ndumb>0) {dumb_flag=1;}
  else {dumb_flag=0;}

  // seed randoms
  srand(seed);

  // generate DNA beads
  id=1;
  new_atom.x=0.0; new_atom.y=0.0; new_atom.z=0.0; new_atom.id=id;
  new_atom.id=id; new_atom.type=DNATYPE; new_atom.mol=1;
  atoms.push_back( new_atom );
  last_atom=atoms.back();
  for (int i=1;i<N;i++) {

    do { // generate a random position, but check it is inside the box
      theta=double(rand())/double(RAND_MAX)*PI;
      phi=double(rand())/double(RAND_MAX)*2.0*PI;
      dx=atoms[i-1].x+sin(theta)*cos(phi);
      dy=atoms[i-1].y+sin(theta)*sin(phi);
      dz=atoms[i-1].z+cos(theta);
    } while (abs(dx)>lx||abs(dy)>ly||abs(dz)>lz);

    id++; // add a new atom to the vector
    new_atom.id=id; new_atom.type=1; new_atom.mol=1;
    new_atom.x=dx; new_atom.y=dy; new_atom.z=dz;
    atoms.push_back( new_atom );

    // add a new bond
    new_bond.one=last_atom.id;
    new_bond.two=atoms.back().id;
    new_bond.type=1;
    bonds.push_back( new_bond );

    if (i>1) {
      // add a new angle
      new_angle.one=lastlast_atom.id;
      new_angle.two=last_atom.id;
      new_angle.three=atoms.back().id;
      new_angle.type=1;
      angles.push_back( new_angle );
    }

    // set up for next atom
    lastlast_atom=last_atom;
    last_atom=atoms.back();

  }


  // Now do sphere proteins
  if (prot_flag) {
    mol=1;
    for (int i=0;i<Nprot;i++) {
      id++; mol++;
      new_atom.id=id; new_atom.type=PROTTYPE; new_atom.mol=mol;
      new_atom.x=lx*double(rand())/double(RAND_MAX)-lx*0.5;
      new_atom.y=ly*double(rand())/double(RAND_MAX)-ly*0.5;
      new_atom.z=lz*double(rand())/double(RAND_MAX)-lz*0.5;
      atoms.push_back( new_atom );
    }

  }

  // Now do dumbells
  if (dumb_flag) {
    double randX, randY, randZ;
    // id numbering carrys on
    // mol numbering carrys on, unless there were no proteins, in which case we have to set it
    if (prot_flag==0) { mol=1; }
    for (int i=0;i<Ndumb;i++) {
      mol++; // increment molecule number
      randX=(lx-3)*( double(rand())/double(RAND_MAX)-0.5 );   // don't want to intiate a dumbell at the boundary, so subtract 3 from lx
      randY=(ly-3)*( double(rand())/double(RAND_MAX)-0.5 );
      randZ=(lz-3)*( double(rand())/double(RAND_MAX)-0.5 );

      // add first atom of dumbell - type A
      id++;
      new_atom.id=id; new_atom.type=DUMBTYPEA; new_atom.mol=mol;
      new_atom.x=randX;
      new_atom.y=randY;
      new_atom.z=randZ;
      atoms.push_back( new_atom );

      // add second atom of dumbell - type B
      id++;
      new_atom.id=id; new_atom.type=DUMBTYPEB; new_atom.mol=mol;
      new_atom.x=randX+1.0;  // shift by 1 along x
      new_atom.y=randY;
      new_atom.z=randZ;
      atoms.push_back( new_atom );      
    }
  }

  // generate lammps file
  sprintf(fn,"lammps.input");
  ouf.open(fn);
  ouf<<" LAMMPS data file for gas of spheres"<<endl;
  ouf<<endl;
  ouf<<" "<<atoms.size()<<" atoms"<<endl;
  ouf<<" "<<bonds.size()<<" bonds"<<endl;
  ouf<<" "<<angles.size()<<" angles"<<endl;
  ouf<<endl;
  ouf<<" "<<1+prot_flag+2*dumb_flag<<" atom types"<<endl;
  ouf<<" 1 bond types"<<endl;
  ouf<<" 1 angle types"<<endl;
  ouf<<endl;
  ouf<<" "<<-0.5*lx<<" "<<0.5*lx<<" xlo xhi"<<endl;
  ouf<<" "<<-0.5*ly<<" "<<0.5*ly<<" ylo yhi"<<endl;
  ouf<<" "<<-0.5*lz<<" "<<0.5*lz<<" zlo zhi"<<endl;
  ouf<<endl;
  ouf<<endl<<" Masses"<<endl<<endl;
  ouf<<" "<<"1 1"<<endl;
  if (prot_flag) {
    ouf<<" "<<PROTTYPE<<" 1"<<endl;
  }
  if (dumb_flag) {
    ouf<<" "<<DUMBTYPEA<<" 1"<<endl;
    ouf<<" "<<DUMBTYPEB<<" 1"<<endl;
  }
  ouf<<endl<<" Atoms"<<endl<<endl;
  for (int i=0;i<atoms.size();i++) {
    ouf<<atoms[i].id<<" "<<atoms[i].mol<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" 0 0 0"<<endl;
  }
  ouf<<endl<<" Velocities"<<endl<<endl;
  for (int i=0;i<atoms.size();i++) {
    ouf<<" "<<i+1<<" 0 0 0"<<endl;
  }
  ouf<<endl<<" Bonds"<<endl<<endl;
  for (int i=0;i<bonds.size();i++) {
    ouf<<" "<<i+1<<" "<<bonds[i].type<<" "<<bonds[i].one<<" "<<bonds[i].two<<endl;
  }
  ouf<<endl<<" Angles"<<endl<<endl;
  for (int i=0;i<angles.size();i++) {
    ouf<<" "<<i+1<<" "<<angles[i].type<<" "<<angles[i].one<<" "<<angles[i].two<<" "<<angles[i].three<<endl;
  }

  // tidy up
  ouf.close();
  
}
