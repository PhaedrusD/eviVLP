// FileProcessor.cpp
#include "FileProcessor.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

FileProcessor::FileProcessor(const std::string& templateFilePath)
    : templateFilePath(templateFilePath) {}

bool FileProcessor::ProcessAminoAcidFile(const std::string& aminoAcid, 
                                        const std::string& firstBondedAtom, 
                                        const std::string& thirdBondedAtom) {
    // Get the correct filename based on amino acid selection
    std::string fileName = GetAminoAcidFileName(aminoAcid, false);
    
    if (fileName.empty()) {
        std::cerr << "Invalid amino acid selection: " << aminoAcid << "\n";
        return false;
    }
    
    std::cout << "Processing file: " << fileName << "\n";
    std::cout << "First bonded atom: " << firstBondedAtom << "\n";
    std::cout << "Third bonded atom: " << thirdBondedAtom << "\n";
    
    // Replace *T1 with first bonded atom
    if (!firstBondedAtom.empty()) {
        if (!ReplaceInFile(fileName, "*T1", firstBondedAtom)) {
            std::cerr << "Failed to replace TEST_ATOM1 in " << fileName << "\n";
            return false;
        }
    }
    
    // Replace *T3 with third bonded atom
   if (!thirdBondedAtom.empty()) {
        if (!ReplaceInFile(fileName, "*T3", thirdBondedAtom)) {
            std::cerr << "Failed to replace TEST_ATOM3 in " << fileName << "\n";
            return false;
        }
    }
    
    std::cout << "Successfully processed " << fileName << "\n";
    return true;
}

void FileProcessor::ResetAminoAcidFile(const std::string& aminoAcid) {
    std::string templateFileName = GetAminoAcidFileName(aminoAcid, true);
    std::string targetFileName = GetAminoAcidFileName(aminoAcid, false);
    
    if (templateFileName.empty() || targetFileName.empty()) {
        std::cerr << "Invalid amino acid selection for reset: " << aminoAcid << "\n";
        return;
    }
    
    std::ifstream templateFile(templateFileName);
    if (!templateFile.is_open()) {
        std::cerr << "Failed to open template file: " << templateFileName << "\n";
        return;
    }
    
    std::ofstream targetFile(targetFileName);
    if (!targetFile.is_open()) {
        std::cerr << "Failed to open target file: " << targetFileName << "\n";
        return;
    }
    
    targetFile << templateFile.rdbuf();
    std::cout << "Reset " << targetFileName << " from " << templateFileName << "\n";
}

bool FileProcessor::ProcessPatchMakerFile(const std::string& aminoAcid,
                                         const std::vector<std::string>& bondedAtoms,
                                         const std::vector<std::string>& deletedAtoms) {
    std::string fileName = GetPatchMakerFileName(aminoAcid, false);
    
    if (fileName.empty()) {
        std::cerr << "Invalid amino acid selection: " << aminoAcid << "\n";
        return false;
    }
    
    std::cout << "Processing patch_maker file: " << fileName << "\n";
    
    // Create atom_name array (bonded + deleted atoms)
    std::vector<std::string> atomNameArray = CreateAtomNameArray(bondedAtoms, deletedAtoms);
    
    // Create atom_info array (bonded atoms only)
    std::vector<std::string> atomInfoArray = CreateAtomInfoArray(bondedAtoms);
    
    // Read the file
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << fileName << "\n";
        return false;
    }
    
    std::ostringstream buffer;
    buffer << file.rdbuf();
    std::string fileContent = buffer.str();
    file.close();
    
    // Replace atom_name array
    std::ostringstream atomNameStr;
    atomNameStr << "atom_name = [";
    for (size_t i = 0; i < atomNameArray.size(); ++i) {
        atomNameStr << "'" << atomNameArray[i] << "'";
        if (i < atomNameArray.size() - 1) {
            atomNameStr << ", ";
        }
    }
    atomNameStr << "]";
    
    // Replace atom_info array
    std::ostringstream atomInfoStr;
    atomInfoStr << "atom_info = [";
    for (size_t i = 0; i < atomInfoArray.size(); ++i) {
        atomInfoStr << "'" << atomInfoArray[i] << "'";
        if (i < atomInfoArray.size() - 1) {
            atomInfoStr << ", ";
        }
    }
    atomInfoStr << "]";
    
    // Find and replace the atom_name line
    size_t atomNamePos = fileContent.find("atom_name = [");
    if (atomNamePos != std::string::npos) {
        size_t endPos = fileContent.find("]", atomNamePos);
        if (endPos != std::string::npos) {
            fileContent.replace(atomNamePos, endPos - atomNamePos + 1, atomNameStr.str());
            std::cout << "Replaced atom_name array\n";
        }
    }
    
    // Find and replace the atom_info line
    size_t atomInfoPos = fileContent.find("atom_info = [");
    if (atomInfoPos != std::string::npos) {
        size_t endPos = fileContent.find("]", atomInfoPos);
        if (endPos != std::string::npos) {
            fileContent.replace(atomInfoPos, endPos - atomInfoPos + 1, atomInfoStr.str());
            std::cout << "Replaced atom_info array\n";
        }
    }
    
    // Replace individual TEST_ATOM placeholders
    // Replace *T1 with first bonded atom if available
    if (!bondedAtoms.empty() && !bondedAtoms[0].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
    while ((pos = fileContent.find("*T1", pos)) != std::string::npos) {
        fileContent.replace(pos, 3, bondedAtoms[0]); // 3 is length of "*T1"
            pos += bondedAtoms[0].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of TEST_ATOM1 with " << bondedAtoms[0] << "\n";
        }
    }
    
    // Replace TEST_ATOM2 with second bonded atom if available
    if (bondedAtoms.size() > 1 && !bondedAtoms[1].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
        while ((pos = fileContent.find("*T2", pos)) != std::string::npos) {
            fileContent.replace(pos, 3, bondedAtoms[1]); // 3 is length of "*T2"
            pos += bondedAtoms[1].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of *T2 with " << bondedAtoms[1] << "\n";
        }
    }
    
    // Replace TEST_ATOM3 with third bonded atom if available
    if (bondedAtoms.size() > 2 && !bondedAtoms[2].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
        while ((pos = fileContent.find("*T3", pos)) != std::string::npos) {
            fileContent.replace(pos, 3, bondedAtoms[2]); // 3 is length of "*T3"
            pos += bondedAtoms[2].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of TEST_ATOM3 with " << bondedAtoms[2] << "\n";
        }
    }
    // Replace TEST_ATOM4 with third bonded atom if available
    if (bondedAtoms.size() > 2 && !bondedAtoms[2].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
        while ((pos = fileContent.find("*T4", pos)) != std::string::npos) {
            fileContent.replace(pos, 3, bondedAtoms[2]); // 3 is length of "*T3"
            pos += bondedAtoms[2].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of TEST_ATOM3 with " << bondedAtoms[2] << "\n";
        }
    }

    // Replace TEST_ATOM5 with third bonded atom if available
    if (bondedAtoms.size() > 2 && !bondedAtoms[2].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
        while ((pos = fileContent.find("*T5", pos)) != std::string::npos) {
            fileContent.replace(pos, 3, bondedAtoms[2]); // 3 is length of "*T3"
            pos += bondedAtoms[2].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of TEST_ATOM3 with " << bondedAtoms[2] << "\n";
        }
    }

    // Replace TEST_ATOM6 with third bonded atom if available
    if (bondedAtoms.size() > 2 && !bondedAtoms[2].empty()) {
        size_t pos = 0;
        int replacementCount = 0;
        while ((pos = fileContent.find("*T6", pos)) != std::string::npos) {
            fileContent.replace(pos, 3, bondedAtoms[2]); // 3 is length of "*T3"
            pos += bondedAtoms[2].length();
            replacementCount++;
        }
        if (replacementCount > 0) {
            std::cout << "Replaced " << replacementCount << " occurrences of TEST_ATOM3 with " << bondedAtoms[2] << "\n";
        }
    }
    // Write the updated content back to the file
    std::ofstream outFile(fileName);
    if (!outFile.is_open()) {
        std::cerr << "Failed to write to " << fileName << "\n";
        return false;
    }
    
    outFile << fileContent;
    std::cout << "Successfully processed " << fileName << "\n";
    return true;
}

void FileProcessor::ResetPatchMakerFile(const std::string& aminoAcid) {
    std::string templateFileName = GetPatchMakerFileName(aminoAcid, true);
    std::string targetFileName = GetPatchMakerFileName(aminoAcid, false);
    
    if (templateFileName.empty() || targetFileName.empty()) {
        std::cerr << "Invalid amino acid selection for patch_maker reset: " << aminoAcid << "\n";
        return;
    }
    
    std::ifstream templateFile(templateFileName);
    if (!templateFile.is_open()) {
        std::cerr << "Failed to open template file: " << templateFileName << "\n";
        return;
    }
    
    std::ofstream targetFile(targetFileName);
    if (!targetFile.is_open()) {
        std::cerr << "Failed to open target file: " << targetFileName << "\n";
        return;
    }
    
    targetFile << templateFile.rdbuf();
    std::cout << "Reset " << targetFileName << " from " << templateFileName << "\n";
}

bool FileProcessor::ReplaceInFile(const std::string& filePath, const std::string& target, const std::string& replacement) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filePath << "\n";
        return false;
    }

    std::ostringstream buffer;
    buffer << file.rdbuf();
    std::string fileContent = buffer.str();
    file.close();

    // Count replacements for verification
    int replacementCount = 0;
    size_t pos = 0;
    while ((pos = fileContent.find(target, pos)) != std::string::npos) {
        fileContent.replace(pos, target.length(), replacement);
        pos += replacement.length();
        replacementCount++;
    }

    if (replacementCount == 0) {
        std::cout << "Warning: No occurrences of '" << target << "' found in " << filePath << "\n";
    } else {
        std::cout << "Replaced " << replacementCount << " occurrences of '" << target << "' with '" << replacement << "' in " << filePath << "\n";
    }

    // Write the updated content back to the file
    std::ofstream outFile(filePath);
    if (!outFile.is_open()) {
        std::cerr << "Failed to write to " << filePath << "\n";
        return false;
    }

    outFile << fileContent;
    return true;
}

std::string FileProcessor::GetAminoAcidFileName(const std::string& aminoAcid, bool isTemplate) {
    // Convert amino acid to lowercase for comparison
    std::string lowerAminoAcid = aminoAcid;
    std::transform(lowerAminoAcid.begin(), lowerAminoAcid.end(), lowerAminoAcid.begin(), ::tolower);
    
    std::string baseName;
    if (lowerAminoAcid == "lysine") {
        baseName = "orient_place_drug_lys";
    } else if (lowerAminoAcid == "cysteine") {
        baseName = "orient_place_drug_cys";
    } else if (lowerAminoAcid == "tyrosine") {
        baseName = "orient_place_drug_tyr";
    } else if (lowerAminoAcid == "N-Terminus") {
        baseName = "orient_place_drug_ntrm";
    } else {
        return ""; // Invalid amino acid
    }
    
    if (isTemplate) {
        return baseName + "_template.py";
    } else {
        return baseName + ".py";
    }
}

std::string FileProcessor::GetPatchMakerFileName(const std::string& aminoAcid, bool isTemplate) {
    // Convert amino acid to lowercase for comparison
    std::string lowerAminoAcid = aminoAcid;
    std::transform(lowerAminoAcid.begin(), lowerAminoAcid.end(), lowerAminoAcid.begin(), ::tolower);
    
    std::string baseName;
    if (lowerAminoAcid == "lysine") {
        baseName = "patch_maker_lys";
    } else if (lowerAminoAcid == "cysteine") {
        baseName = "patch_maker_cys";
    } else if (lowerAminoAcid == "tyrosine") {
        baseName = "patch_maker_tyr";
    } else if (lowerAminoAcid == "N-Terminus") {
        baseName = "patch_maker_ntrm";
    } else {
        return ""; // Invalid amino acid
    }
    
    if (isTemplate) {
        return baseName + "_template.py";
    } else {
        return baseName + ".py";
    }
}

std::vector<std::string> FileProcessor::CreateAtomNameArray(
    const std::vector<std::string>& bondedAtoms,
    const std::vector<std::string>& deletedAtoms) 
{
    std::vector<std::string> atomNameArray;
    
    // First add all bonded atoms (these will be replaced later)
    for (const auto& atom : bondedAtoms) {
        if (!atom.empty()) {
            atomNameArray.push_back("2" + atom);
        }
    }
    
    // Then add all deleted atoms 
    for (const auto& atom : deletedAtoms) {
        if (!atom.empty()) {
            atomNameArray.push_back("2" + atom);
        }
    }
    
    // NO PADDING - return only the real atoms
    return atomNameArray;
}

std::vector<std::string> FileProcessor::CreateAtomInfoArray(
    const std::vector<std::string>& bondedAtoms) 
{
    std::vector<std::string> atomInfoArray;
    
    // Add exactly the bonded atoms (will be replaced later)
    for (const auto& atom : bondedAtoms) {
        if (!atom.empty()) {
            atomInfoArray.push_back("2" + atom);
        }
    }
    
    // Pad to exactly 3 elements if needed
    while (atomInfoArray.size() < 3) {
        atomInfoArray.push_back(""); // Empty string as placeholder
    }
    
    return atomInfoArray;
}

