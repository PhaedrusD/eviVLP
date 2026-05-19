#include "FileProcessor.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

void FileProcessor::ProcessPythonFile(const std::string& filePath, const std::string& baseName) {
    std::ifstream inFile(filePath);
    if (!inFile) {
        std::cerr << "Error: Cannot open file " << filePath << std::endl;
        return;
    }
    
    std::ostringstream content;
    content << inFile.rdbuf();
    inFile.close();

    std::string fileContent = content.str();
    size_t pos = fileContent.find("****");
    while (pos != std::string::npos) {
        fileContent.replace(pos, 4, baseName);
        pos = fileContent.find("****", pos + baseName.length());
    }

    std::ofstream outFile(filePath);
    if (!outFile) {
        std::cerr << "Error: Cannot write to file " << filePath << std::endl;
        return;
    }
    outFile << fileContent;
    outFile.close();
}

void FileProcessor::ResetPythonFile(const std::string& filePath) {
    std::string templatePath;
    if (filePath.find("TMV") != std::string::npos) {
        templatePath = "Patch_orient_TMV_template.py";
    } else {
        templatePath = "Patch_orient_QB_template.py";
    }
    
    std::ifstream templateFile(templatePath);
    if (!templateFile) {
        std::cerr << "Error: Cannot open template file " << templatePath << std::endl;
        return;
    }
    
    std::ofstream outFile(filePath);
    if (!outFile) {
        std::cerr << "Error: Cannot write to file " << filePath << std::endl;
        return;
    }
    
    outFile << templateFile.rdbuf();
    templateFile.close();
    outFile.close();
}

std::string FileProcessor::GetResidueNumberFromPDB(const std::string& pdbFilePath) {
    std::ifstream pdbFile(pdbFilePath);
    if (!pdbFile) {
        std::cerr << "Error: Cannot open PDB file " << pdbFilePath << std::endl;
        return "";
    }
    
    std::string line;
    
    while (std::getline(pdbFile, line)) {
        if (line.substr(0, 4) == "ATOM") {
            if (line.length() >= 26) {
                std::string residueNum = line.substr(22, 4);
                
                // Trim whitespace
                size_t start = residueNum.find_first_not_of(" ");
                if (start != std::string::npos) {
                    size_t end = residueNum.find_last_not_of(" ");
                    residueNum = residueNum.substr(start, end - start + 1);
                }
                
                pdbFile.close();
                return residueNum;
            }
        }
    }
    
    pdbFile.close();
    return "";
}

void FileProcessor::ProcessPythonFileWithResidue(const std::string& filePath, const std::string& baseName, const std::string& pdbFilePath) {
    ProcessPythonFile(filePath, baseName);
    
    std::string residueNumber = GetResidueNumberFromPDB(pdbFilePath);
    
    if (residueNumber.empty()) {
        std::cerr << "Warning: No residue number found in PDB file" << std::endl;
        return;
    }
    
    std::ifstream inFile(filePath);
    if (!inFile) {
        std::cerr << "Error: Cannot open file " << filePath << std::endl;
        return;
    }
    
    std::ostringstream content;
    content << inFile.rdbuf();
    inFile.close();

    std::string fileContent = content.str();
    
    size_t pos = fileContent.find("r*s");
    while (pos != std::string::npos) {
        fileContent.replace(pos, 3, residueNumber);
        pos = fileContent.find("r*s", pos + residueNumber.length());
    }

    std::ofstream outFile(filePath);
    if (!outFile) {
        std::cerr << "Error: Cannot write to file " << filePath << std::endl;
        return;
    }
    outFile << fileContent;
    outFile.close();
}

void FileProcessor::ProcessQbetaFile(const std::string& pdbFilePath, const std::string& selectedNumber) {
    std::string residueNumber = GetResidueNumberFromPDB(pdbFilePath);
    
    if (residueNumber.empty()) {
        std::cerr << "Warning: No residue number found in PDB file for Qbeta" << std::endl;
        return;
    }
    
    // Map selected number to replacement value
    std::string replacementValue;
    bool isAll = false;
    if (selectedNumber == "6") {
        replacementValue = "10";
    } else if (selectedNumber == "12") {
        replacementValue = "5";
    } else if (selectedNumber == "20") {
        replacementValue = "3";
    } else if (selectedNumber == "All") {
        replacementValue = "1";
        isAll = true;
    } else {
        std::cerr << "Warning: Invalid number selection" << std::endl;
        return;
    }

    std::ifstream inFile("Qbeta_pgn_writer.py");
    if (!inFile) {
        std::cerr << "Error: Cannot open Qbeta_pgn_writer.py" << std::endl;
        return;
    }
    
    std::ostringstream content;
    content << inFile.rdbuf();
    inFile.close();

    std::string fileContent = content.str();
    
    size_t pos = fileContent.find("r*s");
    while (pos != std::string::npos) {
        fileContent.replace(pos, 3, residueNumber);
        pos = fileContent.find("r*s", pos + residueNumber.length());
    }

    // Replace ** with the mapped value
    pos = fileContent.find("**");
    while (pos != std::string::npos) {
        fileContent.replace(pos, 2, replacementValue);
        pos = fileContent.find("**", pos + replacementValue.length());
    }

    // If "All" is selected, uncomment the specified lines
    if (isAll) {
        std::istringstream iss(fileContent);
        std::ostringstream oss;
        std::string line;
        int lineNum = 0;
        
        while (std::getline(iss, line)) {
            lineNum++;
            // Uncomment lines 58-66 and 76-84
            if ((lineNum >= 58 && lineNum <= 66) || (lineNum >= 76 && lineNum <= 84)) {
                // Remove leading "# " if present
                if (line.length() >= 2 && line[0] == '#' && line[1] == ' ') {
                    line = line.substr(2);
                }
            }
            oss << line << '\n';
        }
        fileContent = oss.str();
    }

    std::ofstream outFile("Qbeta_pgn_writer.py");
    if (!outFile) {
        std::cerr << "Error: Cannot write to Qbeta_pgn_writer.py" << std::endl;
        return;
    }
    outFile << fileContent;
    outFile.close();
}

void FileProcessor::ResetQbetaFile() {
    std::ifstream templateFile("Qbeta_pgn_writer_template.py");
    if (!templateFile) {
        std::cerr << "Error: Cannot open Qbeta_pgn_writer_template.py" << std::endl;
        return;
    }
    
    std::ofstream outFile("Qbeta_pgn_writer.py");
    if (!outFile) {
        std::cerr << "Error: Cannot write to Qbeta_pgn_writer.py" << std::endl;
        return;
    }
    
    outFile << templateFile.rdbuf();
    templateFile.close();
    outFile.close();
}

void FileProcessor::ProcessTMVFile(const std::string& pdbFilePath) {
    std::cout << "Processing TMV file with PDB: " << pdbFilePath << std::endl;
    std::string residueNumber = GetResidueNumberFromPDB(pdbFilePath);
    std::cout << "Residue number extracted: '" << residueNumber << "'" << std::endl;
    
    if (residueNumber.empty()) {
        std::cerr << "Warning: No residue number found in PDB file for TMV" << std::endl;
        return;
    }
    
    std::ifstream inFile("TMV_pgn_writer.py");
    if (!inFile) {
        std::cerr << "Error: Cannot open TMV_pgn_writer.py" << std::endl;
        return;
    }
    
    std::ostringstream content;
    content << inFile.rdbuf();
    inFile.close();

    std::string fileContent = content.str();
    
    size_t pos = fileContent.find("r*s");
    while (pos != std::string::npos) {
        fileContent.replace(pos, 3, residueNumber);
        pos = fileContent.find("r*s", pos + residueNumber.length());
    }

    std::ofstream outFile("TMV_pgn_writer.py");
    if (!outFile) {
        std::cerr << "Error: Cannot write to TMV_pgn_writer.py" << std::endl;
        return;
    }
    outFile << fileContent;
    outFile.close();
}

void FileProcessor::ResetTMVFile() {
    std::ifstream templateFile("TMV_pgn_writer_template.py");
    if (!templateFile) {
        std::cerr << "Error: Cannot open TMV_pgn_writer_template.py" << std::endl;
        return;
    }
    
    std::ofstream outFile("TMV_pgn_writer.py");
    if (!outFile) {
        std::cerr << "Error: Cannot write to TMV_pgn_writer.py" << std::endl;
        return;
    }
    
    outFile << templateFile.rdbuf();
    templateFile.close();
    outFile.close();
}