// FileProcessor.h
#ifndef FILEPROCESSOR_H
#define FILEPROCESSOR_H

#include <string>
#include <vector>

class FileProcessor {
public:
    FileProcessor(const std::string& templateFilePath);
    
    // Existing methods
    bool ProcessAminoAcidFile(const std::string& aminoAcid, 
                             const std::string& firstBondedAtom, 
                             const std::string& thirdBondedAtom);
    void ResetAminoAcidFile(const std::string& aminoAcid);
    
    // New method for patch_maker processing
    bool ProcessPatchMakerFile(const std::string& aminoAcid,
                              const std::vector<std::string>& bondedAtoms,
                              const std::vector<std::string>& deletedAtoms);
    void ResetPatchMakerFile(const std::string& aminoAcid);
    
private:
    std::string templateFilePath;
    
    // Existing helper methods
    bool ReplaceInFile(const std::string& filePath, const std::string& target, const std::string& replacement);
    std::string GetAminoAcidFileName(const std::string& aminoAcid, bool isTemplate);
    
    // New helper methods for patch_maker
    std::string GetPatchMakerFileName(const std::string& aminoAcid, bool isTemplate);
    std::vector<std::string> CreateAtomNameArray(const std::vector<std::string>& bondedAtoms, 
                                                const std::vector<std::string>& deletedAtoms);
    std::vector<std::string> CreateAtomInfoArray(const std::vector<std::string>& bondedAtoms);
};

#endif // FILEPROCESSOR_H
