#include <iostream>       // Input-output stream library
#include <fstream>        // File stream library
#include <sstream>        // String stream library
#include <string>         // String library
#include <vector>         // Vector library
#include <cmath>          // Math functions library
#include <algorithm>      // Algorithm library
#include <iterator>       // Iterator library

using namespace std;    // Standard namespace for C++ libraries

const int SIZE = 1000;   // Maximum size for arrays

// Struct to represent an atom
struct Atom {
    string type;         // Atom type
    double coords[3];    // Atom coordinates
};

// Struct to store topology information of a residue
struct ResidueInfo {
    string line;         // Line containing residue information
};

// Function to read topology information from a file
vector<ResidueInfo> readTopologyInfo(const string& filename) {
    vector<ResidueInfo> topologyInfo;   // Vector to store topology information
    ifstream file(filename.c_str());    // Open file stream

    string line;    // String to store each line of the file
    bool readTopology = false;  // Flag to indicate whether topology information is being read

    // Read lines from the file until the end is reached
    while (getline(file, line)) {
        // Check if "RESI" is found to start reading topology information
        if (line.find("RESI") != string::npos) {
            readTopology = true;
        }
        // Check if "END" is found to stop reading topology information
        if (readTopology && line.find("END") != string::npos) {
            break;
        }
        // Store topology information lines in the vector
        if (readTopology) {
            topologyInfo.push_back({line});
        }
    }

    file.close();  // Close the file
    return topologyInfo;  // Return the vector containing topology information
}

// Function to read atom data from a PDB file
int readData(const string& filename, Atom atoms[SIZE]) {
    ifstream inFile(filename);   // Open file stream
    string line;    // String to store each line of the file
    int numatoms = 0;   // Counter for the number of atoms read

    // Read lines from the file until the end is reached
    while (getline(inFile, line)) {
        // Check if the line represents an atom ("HETATM" prefix)
        if (line.substr(0,6) == "HETATM" || line.substr(0, 4) == "ATOM") {
            // Parse atom type and coordinates from the line
            string atomType = line.substr(12, 4);
            double x, y, z;
            istringstream iss(line.substr(30));
            iss >> x >> y >> z;

            // Store atom data in the Atom array
            atoms[numatoms].type = atomType;
            atoms[numatoms].coords[0] = x;
            atoms[numatoms].coords[1] = y;
            atoms[numatoms].coords[2] = z;

            numatoms++;  // Increment the atom counter

            // Check if maximum number of atoms is reached
            if (numatoms >= SIZE) {
                cerr << "Maximum number of atoms reached." << endl;
                break;
            }
        } 
    }
    inFile.close();  // Close the file
    return numatoms;  // Return the number of atoms read
}

// Function to calculate distance between two atoms
double Distance(const Atom& atom1, const Atom& atom2) {
    // Calculate the differences in coordinates
    double dx = atom2.coords[0] - atom1.coords[0];
    double dy = atom2.coords[1] - atom1.coords[1];
    double dz = atom2.coords[2] - atom1.coords[2];
    // Calculate and return the Euclidean distance
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate angle formed by three atoms
double Angle(const Atom& atom1, const Atom& atom2, const Atom& atom3) {
    // Calculate the vectors between atoms
    double vec1[3] = {atom1.coords[0] - atom2.coords[0],atom1.coords[1] - atom2.coords[1], atom1.coords[2] - atom2.coords[2]};
    double vec2[3] = {atom3.coords[0] - atom2.coords[0],atom3.coords[1] - atom2.coords[1], atom3.coords[2] - atom2.coords[2]};
    
    // Calculate dot product and magnitudes of vectors
    double dotProduct = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    double magnitude1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2]);
    double magnitude2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2]);
    
    // Calculate the angle in radians and convert to degrees
    double angleRad = acos(dotProduct / (magnitude1 * magnitude2));
    double angleDeg = angleRad * 180.0 / M_PI;

    return angleDeg;  // Return the angle in degrees
}

#include <cmath> // Include the cmath library for trigonometric functions

// Function to calculate dihedral angle formed by four atoms
double dihedral(const Atom& atom1, const Atom& atom2, const Atom& atom3, const Atom& atom4) {
    // Calculate vectors between atoms
    double vec1[3] = {atom1.coords[0] - atom2.coords[0], atom1.coords[1] - atom2.coords[1], atom1.coords[2] - atom2.coords[2]};
    double vec2[3] = {atom3.coords[0] - atom2.coords[0], atom3.coords[1] - atom2.coords[1], atom3.coords[2] - atom2.coords[2]};
    double vec3[3] = {atom2.coords[0] - atom3.coords[0], atom2.coords[1] - atom3.coords[1], atom2.coords[2] - atom3.coords[2]};
    double vec4[3] = {atom4.coords[0] - atom3.coords[0], atom4.coords[1] - atom3.coords[1], atom4.coords[2] - atom3.coords[2]};
   
    // Calculate normals to the planes defined by vectors
    double normal1[3] = {
        vec1[1] * vec2[2] - vec1[2] * vec2[1],
        vec1[2] * vec2[0] - vec1[0] * vec2[2],
        vec1[0] * vec2[1] - vec1[1] * vec2[0]
    };
    double normal2[3] = {
        vec3[1] * vec4[2] - vec3[2] * vec4[1],
        vec3[2] * vec4[0] - vec3[0] * vec4[2],
        vec3[0] * vec4[1] - vec3[1] * vec4[0]
    };

    // Calculate dot product and magnitudes of normals
    double ndotProduct = normal1[0] * normal2[0] + normal1[1] * normal2[1] + normal1[2] * normal2[2];
    double nmagnitude1 = sqrt(normal1[0] * normal1[0] + normal1[1] * normal1[1] + normal1[2] * normal1[2]);
    double nmagnitude2 = sqrt(normal2[0] * normal2[0] + normal2[1] * normal2[1] + normal2[2] * normal2[2]);

    // Calculate cross product of normal vectors
    double crossProduct[3] = {
        normal1[1] * normal2[2] - normal1[2] * normal2[1],
        normal1[2] * normal2[0] - normal1[0] * normal2[2],
        normal1[0] * normal2[1] - normal1[1] * normal2[0]
    };

    // Calculate the magnitude of the cross product
    double crossMagnitude = sqrt(crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1] + crossProduct[2] * crossProduct[2]);

    // Check for division by zero
    if (nmagnitude1 == 0 || nmagnitude2 == 0 || crossMagnitude == 0) {
        return 0.0;  // Return 0 if division by zero occurs
    }

    // Calculate the cosine of the angle
    double cosAngle = ndotProduct / (nmagnitude1 * nmagnitude2);

    // Calculate the angle in radians using asin to ensure correct sign
    double phiRad = acos(cosAngle);

    // Determine the sign of phi based on the sign of signDot
    double signDot = normal1[0] * vec3[0] + normal1[1] * vec3[1] + normal1[2] * vec3[2];
    double phiDeg;
    if (signDot <= 0) {
        phiDeg = phiRad * 180.0 / M_PI;
    } else {
        phiDeg = -phiRad * 180.0 / M_PI;
    }

    return phiDeg;  // Return the dihedral angle in degrees
}



// Function to write internal coordinate table to an output stream
void writeInternalCoordinateTable(const Atom atoms[], int numAtoms, std::ostream& out) {
    // Loop through combinations of four atoms to calculate internal coordinates
    for (int i = 0; i < numAtoms - 3; ++i) {
        for (int j = i + 1; j < numAtoms - 2; ++j) {
            for (int k = j + 1; k < numAtoms - 1; ++k) {
                for (int l = k + 1; l < numAtoms; ++l) {
                    // Calculate distances and angles
                    double distAB = Distance(atoms[i], atoms[j]);
                    double angleABC = Angle(atoms[i], atoms[j], atoms[k]);
                    double dihedralABCD = dihedral(atoms[i], atoms[j], atoms[k], atoms[l]);
                    double angleBCD = Angle(atoms[j], atoms[k], atoms[l]);
                    double distCD = Distance(atoms[k], atoms[l]);

                    // Debug prints
                    std::cout << "distAB: " << distAB << " angleABC: " << angleABC << " dihedralABCD: " << dihedralABCD << " angleBCD: " << angleBCD << " distCD: " << distCD << std::endl;


                    // Check if distances are less than or equal to 3
                    if (distAB <= 2.0 && distCD <= 2.0) {
                        // Write internal coordinate table entry to output stream
                        out << "IC " << atoms[i].type << " " << atoms[j].type << " " << atoms[k].type << " " << atoms[l].type << "  " << distAB << " " << angleABC << "  " << dihedralABCD << " " << angleBCD << "  " << distCD << std::endl;
                    }
                }
            }
        }
    }
    // Add an additional line at the end
    out << "END" << std::endl;
}


// Function to insert topology information into a PDB file
void insertTopologyInfo(const string& pdbFilename, const vector<ResidueInfo>& topologyInfo, const Atom atoms[], int numAtoms) {
    ifstream pdbFile(pdbFilename);    // Open original PDB file for reading
    ofstream tempFile("temp.pdb");    // Create temporary file for writing

    string line;  // String to store each line of the file

    // Copy lines from original PDB file until "END" is found
    while (getline(pdbFile, line)) {
        if (line.find("END") != string::npos) {
            break;  // Stop copying when "END" is found
        }
        tempFile << line << endl;  // Write line to temporary file
    }
    // Insert topology information into the temporary file
    for (const auto& info : topologyInfo) {
        tempFile << info.line << endl;  // Write topology information line by line
    }

    // Write the internal coordinate table to the temporary file
    writeInternalCoordinateTable(atoms, numAtoms, tempFile);

    // Write the remaining lines from the original PDB file
    while (getline(pdbFile, line)) {
        tempFile << line << endl;  // Write remaining lines to temporary file
    }

    pdbFile.close();    // Close the original PDB file
    tempFile.close();    // Close the temporary file

    // Rename the temporary file to replace the original PDB file
    if (rename("temp.pdb", pdbFilename.c_str()) != 0) {
        cerr << "Error renaming temporary file." << endl;  // Display error message if renaming fails
    }
}

// Main function
int main() {
    std::string residueFilename = "*****.str";  // File containing topology information
    string pdbFilename = "topology_IC.inp";  // Original PDB file name

    Atom atoms[SIZE];  // Array to store atom data

    // Read atom data from PDB file
    int numAtoms = readData("*****.pdb", atoms);
    std::cout << "Number of atoms read: " << numAtoms << std::endl;

    // Read topology information from file
    auto topologyInfo = readTopologyInfo(residueFilename);
    std::cout << "Topology information read" << std::endl;

    // Insert topology information into PDB file
    insertTopologyInfo(pdbFilename, topologyInfo, atoms, numAtoms);
    std::cout << "Topology information inserted" << std::endl;

    // Open PDB file for appending
    std::ofstream pdbFileAppend(pdbFilename, std::ios_base::app);
    if (!pdbFileAppend.is_open()) {
        cerr << "Error opening pdb file for appending." << endl;  // Display error message if file opening fails
        return 1;
        
    }
}
