#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <limits>

using namespace std;

const int SIZE = 399999;
const std::string filenameA = "capsid.xyz";
const std::string filenameB = "P2.xyz";

const double pi = 3.14159265358979;
const int NTAB = 32;
const int NDIV = 1 + (2147483563 - 1) / NTAB;
const double EPS = 1.2e-21;
const double RNMX = 1.0 - EPS;
int idum = -873; 

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

struct Atom {
    char type;
    double coords[3];
};

// Function to read data from a file into an array of atoms
void readData(const string& filename, Atom atoms[SIZE], int& numatoms) {
    ifstream inFile(filename); // Open the file for reading
    string line;

    // Check if the file opened successfully
    if (!inFile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    // Check if the first line contains the number of atoms
    if (getline(inFile, line)) {
        istringstream iss(line);
        iss >> numatoms;
        cout << "Number of atoms: " << numatoms << endl;
    } else {
        cerr << "Error: Could not read the number of atoms from the file." << endl;
        inFile.close();
        return;
    }

    // Skip the second line (e.g., "generated by VMD")
    if (!getline(inFile, line)) {
        cerr << "Error: File ended unexpectedly while skipping the second line." << endl;
        inFile.close();
        return;
    }

    // Loop through each atom and read its details from subsequent lines
    for (int i = 0; i < numatoms; ++i) {
        if (getline(inFile, line)) {
            istringstream iss(line);
            iss >> atoms[i].type >> atoms[i].coords[0] >> atoms[i].coords[1] >> atoms[i].coords[2];
        } else {
            cerr << "Error: File ended unexpectedly while reading atom data." << endl;
            break;
        }
    }

    inFile.close(); // Close the file after reading
}

// Function to store initial coordinates to an output file
void InitialCoordinates(const string& filename, Atom atoms[SIZE], int numatoms) {
    ofstream outFile(filename);  // Open the output file for writing
    outFile << numatoms << "\nThis is an xyz file\n";  // Write number of atoms and a header line

    // Loop through each atom and write its details to the file
    for (int i = 0; i < numatoms; ++i) {
        // Write atom type and coordinates to the file
        outFile << atoms[i].type << "      " << atoms[i].coords[0] << "     " << atoms[i].coords[1] << "        " << atoms[i].coords[2] << "\n";
    }

    outFile.close();  // Close the output file
}



// Function to generate random numbers using the ran3 algorithm
float ran3(int *idum) {
    // Define static variables to retain state between function calls
    static int inext, inextp;
    static long ma[56];
    static int iff=0;
    long mj, mk;
    int i, ii, k;

    // Initialize the random number generator if required
    if (*idum < 0 || iff == 0) {
        iff=1;
        mj=MSEED-(*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        // Initialize the state vector
        for (i=1; i<54; i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        // Shuffle the state vector
        for (k=1; k<=4; k++)
            for (i=1; i<=55; i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        // Reset pointers and seed
        inext=0;
        inextp=31;
        *idum=1;
    }

    // Update pointers for the next random number
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;

    // Return the random number scaled by a factor
    return mj*FAC;
}

void MoveRandomRotateXYZMoveBack(const string& filename, Atom atoms[SIZE], int numatoms) {
    ifstream inFile(filename); // Open input file
    string junk;

    // Read the original number of atoms from the file
    int originalNumAtoms;
    inFile >> originalNumAtoms;
    getline(inFile, junk); // Consume rest of the line

    // Calculate the translation values to move atoms to origin
    double movex = -atoms[0].coords[0];
    double movey = -atoms[0].coords[1];
    double movez = -atoms[0].coords[2];

    // Translate all atoms to the origin
    for (int i = 0; i < numatoms; ++i) {
        atoms[i].coords[0] += movex;
        atoms[i].coords[1] += movey;
        atoms[i].coords[2] += movez;
    }

    // Generate random angles for rotation
    double alpha = ran3(&idum) * 2.0 * pi;
    double beta = ran3(&idum) * 2.0 * pi;
    double gamma = ran3(&idum) * 2.0 * pi;
    cout << alpha << "  " << beta << "  " << gamma;

    // Calculate the rotation matrix
    double matrix[3][3] = {
        {cos(beta) * cos(gamma), cos(beta) * sin(gamma), -sin(beta)},
        {sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma),
         sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma),
         sin(alpha) * cos(beta)},
        {cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma),
         cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma),
         cos(alpha) * cos(beta)}
    };

    // Apply rotation to all atoms
    for (int i = 0; i < numatoms; ++i) {
        double x = atoms[i].coords[0];
        double y = atoms[i].coords[1];
        double z = atoms[i].coords[2];

        atoms[i].coords[0] = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
        atoms[i].coords[1] = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
        atoms[i].coords[2] = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;
    }

    // Write rotation matrix to a file
    ofstream matrixFile("rotation_matrix.txt");
    if (matrixFile.is_open()) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i < 3 && j < 3) {
                    matrixFile << matrix[i][j] << " ";
                } else if (i == 3 && j == 3) {
                    matrixFile << "1";
                } else {
                    matrixFile << "0";
                }
                if (j < 3) {
                    matrixFile << " ";
                }
            }
            matrixFile << endl;
        }
    }

    // Translate all atoms back to their original positions
    for (int i = 0; i < numatoms; ++i) {
        atoms[i].coords[0] -= movex;
        atoms[i].coords[1] -= movey;
        atoms[i].coords[2] -= movez;
    }

    // Write updated atom positions to the output file
    ofstream outFile(filename);
    outFile << numatoms << "\nThis is an xyz file\n";
    for (int i = 0; i < numatoms; ++i) {
        outFile << atoms[i].type << "   " << atoms[i].coords[0] << "    " << atoms[i].coords[1] << "    " << atoms[i].coords[2] << "\n";
    }

    // Close input and output files
    inFile.close();
    outFile.close();
}


double calculateMinimumDistance(const Atom atomsA[], int numatomsA, const Atom atomsB[], int numatomsB) {
    //double mindist = numeric_limits<double>::infinity(); // Initialize minimum distance to infinity
    double mindist = numeric_limits<double>::infinity();

    double sqrx, sqry, sqrz; // Variables to store squared differences in coordinates

    // Iterate over atoms in set A
    for (int i = 0; i < numatomsA; ++i) {
        // Iterate over atoms in set B
        for (int j = 0; j < numatomsB; ++j) {
            // Calculate squared differences in coordinates
            sqrx = (atomsB[j].coords[0] - atomsA[i].coords[0]) * (atomsB[j].coords[0] - atomsA[i].coords[0]);
            sqry = (atomsB[j].coords[1] - atomsA[i].coords[1]) * (atomsB[j].coords[1] - atomsA[i].coords[1]);
            sqrz = (atomsB[j].coords[2] - atomsA[i].coords[2]) * (atomsB[j].coords[2] - atomsA[i].coords[2]);
            double sqrd = (sqrx + sqry + sqrz); // Sum of squared differences

            // Calculate distance between atoms
            double dist = sqrt(sqrd);

            // Update minimum distance if the calculated distance is smaller
            if (dist < mindist) {
                mindist = dist;
            }
        }
    }

    return mindist; // Return the minimum distance
}


int main () {
    const double mindist_threshold = 0.50; // Threshold for minimum distance
    static Atom atomsA[SIZE];
    static Atom atomsB[SIZE];
    static Atom initialAtomsB[SIZE];
    int numatomsA, numatomsB;

    // Read data from files into arrays of atoms
    readData(filenameA, atomsA, numatomsA);
    readData(filenameB, atomsB, numatomsB);  
    
    // Copy atomsB to initialAtomsB
    for (int i = 0; i < numatomsB; ++i) {
        initialAtomsB[i] = atomsB[i];
    }

    // Write initial coordinates to a file
    InitialCoordinates("initial_coordinates.xyz", initialAtomsB, numatomsB);

    // Loop until minimum distance exceeds threshold
    while (true) {
        // Calculate minimum distance between atomsA and atomsB
        double mindist = calculateMinimumDistance(atomsA, numatomsA, atomsB, numatomsB);
        cout << "Minimum distance is " << mindist << endl;

        // Check if minimum distance is within threshold
        if (mindist <= mindist_threshold) {
            // Write initial coordinates to a file
            InitialCoordinates("P2.xyz", initialAtomsB, numatomsB);
            
            // Move, rotate, and move back atomsB
            MoveRandomRotateXYZMoveBack(filenameB, atomsB, numatomsB);
            cout << "Moved and rotated XYZ" << endl;
        } else {
            // Minimum distance exceeds threshold, stop the loop
            cout << "Minimum distance (" << mindist << ") exceeds threshold. Stopping." << endl;
            break;
        }
    }

    return 0;
}