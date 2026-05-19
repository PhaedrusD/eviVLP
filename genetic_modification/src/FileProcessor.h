#ifndef FILEPROCESSOR_H
#define FILEPROCESSOR_H

#include <string>

class FileProcessor {
public:
    void ProcessPythonFile(const std::string& filePath, const std::string& baseName);
    void ResetPythonFile(const std::string& filePath);
    std::string GetResidueNumberFromPDB(const std::string& pdbFilePath);
    void ProcessPythonFileWithResidue(const std::string& filePath, const std::string& baseName, const std::string& pdbFilePath);
    void ProcessQbetaFile(const std::string& pdbFilePath, const std::string& selectedNumber);  // Add second parameter
    void ResetQbetaFile();
    void ProcessTMVFile(const std::string& pdbFilePath);
    void ResetTMVFile();
};

#endif // FILEPROCESSOR_H