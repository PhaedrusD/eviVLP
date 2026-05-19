#include <wx/wx.h>
#include <wx/progdlg.h>
#include <wx/notebook.h>
#include "FileProcessor.h" // Include your custom file processing class
#include <cstdlib> // for system()

class MyApp : public wxApp {
public:
    virtual bool OnInit();
};

class MyFrame : public wxFrame {
private:
    wxString selectedPDB;
    wxString selectedVLP;
    wxString selectedNumber;
    wxNotebook* notebook;
    
    // Panel 1 controls
    wxComboBox* vlpComboBox;
    wxButton* selectProteinBtn;
    wxButton* processPDBBtn;
    // Panel 2 controls
    wxComboBox* selectNumberCombo;
    wxButton* selectNumberBtn;
    wxButton* displayPDBBtn;
    wxButton* displayZIPBtn;

    void OnSelectVLP(wxCommandEvent& event);
    void OnSelectProtein(wxCommandEvent& event);
    void OnProcessPDB(wxCommandEvent& event);
    void OnSelectNumber(wxCommandEvent& event);
    void OnDisplayPDB(wxCommandEvent& event);
    void OnZipFiles(wxCommandEvent& event);

public:
    MyFrame(const wxString& title);
};

wxIMPLEMENT_APP(MyApp);

bool MyApp::OnInit() {
    MyFrame* frame = new MyFrame("VLP Genetic Modification");
    frame->Show(true);
    return true;
}

MyFrame::MyFrame(const wxString& title) 
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(400, 400)) {

    // Create notebook for panels
    notebook = new wxNotebook(this, wxID_ANY);

    // ===== PANEL 1: VLP and Protein Selection =====
    wxPanel* panel1 = new wxPanel(notebook, wxID_ANY);
    
    // VLP Selection Label and Dropdown
    wxStaticText* vlpLabel = new wxStaticText(panel1, wxID_ANY, "Select VLP Type:", 
                                              wxPoint(50, 20));
    wxArrayString vlpChoices;
    vlpChoices.Add("TMV");
    vlpChoices.Add("QBeta");
    vlpComboBox = new wxComboBox(panel1, wxID_ANY, "Select VLP", 
                                 wxPoint(50, 45), wxSize(200, 30), 
                                 vlpChoices, wxCB_READONLY);
    
    // Select Protein Button
    selectProteinBtn = new wxButton(panel1, wxID_ANY, "Select Protein File", 
                                    wxPoint(50, 95), wxSize(200, 30));
    
        // Process PDB button
    processPDBBtn = new wxButton(panel1, wxID_ANY, "Process PDB", 
                                 wxPoint(50, 135), wxSize(200, 30));
    // Add panel to notebook
    notebook->AddPage(panel1, "VLP & Protein Selection");

    // ===== PANEL 2: Processing Options =====
    wxPanel* panel2 = new wxPanel(notebook, wxID_ANY);
    
    // Select Number with dropdown
    wxArrayString numChoices;
    numChoices.Add("6");
    numChoices.Add("12");
    numChoices.Add("20");
    numChoices.Add("All");
    selectNumberCombo = new wxComboBox(panel2, wxID_ANY, "Select Number", 
                                       wxPoint(50, 50), wxSize(200, 30), 
                                       numChoices, wxCB_READONLY);
    selectNumberBtn = new wxButton(panel2, wxID_ANY, "Run Selected Option", 
                                   wxPoint(50, 90), wxSize(200, 30));
    
    // Display PDB button
    displayPDBBtn = new wxButton(panel2, wxID_ANY, "Display PDB", 
                                 wxPoint(50, 130), wxSize(200, 30));
    
        // Process PDB button
    displayZIPBtn = new wxButton(panel2, wxID_ANY, "Get Files", 
                                 wxPoint(50, 170), wxSize(200, 30));

    // Add panel to notebook
    notebook->AddPage(panel2, "Processing");

    // Bind events
    vlpComboBox->Bind(wxEVT_COMBOBOX, &MyFrame::OnSelectVLP, this);
    selectProteinBtn->Bind(wxEVT_BUTTON, &MyFrame::OnSelectProtein, this);
    processPDBBtn->Bind(wxEVT_BUTTON, &MyFrame::OnProcessPDB, this);
    selectNumberBtn->Bind(wxEVT_BUTTON, &MyFrame::OnSelectNumber, this);
    displayPDBBtn->Bind(wxEVT_BUTTON, &MyFrame::OnDisplayPDB, this);
    displayZIPBtn->Bind(wxEVT_BUTTON, &MyFrame::OnZipFiles, this);
}

void MyFrame::OnSelectVLP(wxCommandEvent& event) {
    selectedVLP = vlpComboBox->GetStringSelection();
    wxMessageBox("Selected VLP: " + selectedVLP);
}

void MyFrame::OnSelectProtein(wxCommandEvent& event) {
    wxFileDialog openFileDialog(this, _("Select PDB file"), "", "", 
                                "PDB files (*.pdb)|*.pdb", wxFD_OPEN | wxFD_FILE_MUST_EXIST);

    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;

    selectedPDB = openFileDialog.GetPath();
    wxMessageBox("Selected Protein File: " + selectedPDB);
}

void MyFrame::OnProcessPDB(wxCommandEvent& event) {
    if (selectedPDB.IsEmpty()) {
        wxMessageBox("Please select a PDB file first.");
        return;
    }
    
    if (selectedVLP.IsEmpty()) {
        wxMessageBox("Please select a VLP type first.");
        return;
    }

    // Create FileProcessor instance
    FileProcessor processor;
    
    // Extract base name from the selected PDB file
    wxString baseName = selectedPDB.BeforeLast('.');
    
    if (selectedVLP == "QBeta") {
        // QBeta processing
        processor.ProcessPythonFileWithResidue("Patch_orient_QB.py", baseName.ToStdString(), selectedPDB.ToStdString());
        system("Python3 Patch_orient_QB.py");
        processor.ResetPythonFile("Patch_orient_QB.py");
    } else if (selectedVLP == "TMV") {
        // TMV processing
        processor.ProcessPythonFileWithResidue("Patch_orient_TMV.py", baseName.ToStdString(), selectedPDB.ToStdString());
        system("Python3 Patch_orient_TMV.py");
        processor.ResetPythonFile("Patch_orient_TMV.py");
    }
    
    // Common steps for both
    system("./rotate");
    system("Python3 rotate_protein.py");
    
    wxMessageBox("PDB processing completed successfully.");
}

void MyFrame::OnSelectNumber(wxCommandEvent& event) {
    selectedNumber = selectNumberCombo->GetStringSelection();

    if (selectedVLP.IsEmpty()) {
        wxMessageBox("Please select a VLP type first.");
        return;
    }
    
    // Create Capsid
    FileProcessor processor;
    
    if (selectedVLP == "QBeta") {
        system("Python3 Qbeta_capsid_maker.py");
        processor.ProcessQbetaFile(selectedPDB.ToStdString(), selectedNumber.ToStdString());
        system("Python3 Qbeta_pgn_writer.py");
        processor.ResetQbetaFile();
    } else if (selectedVLP == "TMV") {
        system("Python3 TMV_capsid_maker.py");
        processor.ProcessTMVFile(selectedPDB.ToStdString());
        system("Python3 TMV_pgn_writer.py");
        processor.ResetTMVFile();
    }
    
    wxMessageBox("Selected number operations completed.");
}

void MyFrame::OnDisplayPDB(wxCommandEvent& event) {
    if (selectedVLP.IsEmpty()) {
        wxMessageBox("Please select a VLP type first.");
        return;
    }
    
    wxString pgnFile, tclFile;
    
    if (selectedVLP == "QBeta") {
        pgnFile = "QB_capsid.pgn";
        tclFile = "load_capsid.tcl";
    } else if (selectedVLP == "TMV") {
        pgnFile = "TMV_capsid.pgn";
        tclFile = "load_rod.tcl"; // Change if TMV uses different TCL file
    }
    
    // First command: Run the PGN script
    wxString cmd1 = wxString::Format("/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e %s", pgnFile);
    int result1 = system(cmd1.mb_str());
    
    if (result1 != 0) {
        wxMessageBox("Error: Failed to run " + pgnFile + " script.");
        return;
    }
    
    // Second command: Load and display the capsid
    wxString cmd2 = wxString::Format("/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -e %s", tclFile);
    int result2 = system(cmd2.mb_str());
    
    if (result2 != 0) {
        wxMessageBox("Error: Failed to load capsid visualization.");
        return;
    }
}

void MyFrame::OnZipFiles(wxCommandEvent& event)
{
    wxMessageBox("Starting file compression...\n");
    
    // Define the files to zip
    wxArrayString filesToZip;
    filesToZip.Add("rod.psf");
    filesToZip.Add("rod.pdb");
    
    // Check if all files exist
    bool allFilesExist = true;
    for (size_t i = 0; i < filesToZip.GetCount(); i++)
    {
        if (!wxFileExists(filesToZip[i]))
        {
            wxMessageBox("Error: File not found: " + filesToZip[i] + "\n");
            allFilesExist = false;
        }
    }
    
    if (!allFilesExist)
    {
        wxMessageBox("Some required files are missing. Please check the status log.", 
                     "Missing Files", wxOK | wxICON_ERROR);
        return;
    }
    
    // Create timestamp for unique filename
    wxDateTime now = wxDateTime::Now();
    wxString timestamp = now.Format("%Y%m%d_%H%M%S");
    wxString zipFileName = "capsid_files_" + timestamp + ".zip";
    
    // Build the zip command
    // Using system zip command (available on macOS and Linux)
    wxString command = "zip " + zipFileName;
    for (size_t i = 0; i < filesToZip.GetCount(); i++)
    {
        command += " " + filesToZip[i];
    }
    
    wxMessageBox("Running command: " + command + "\n");
    
    int result = system(command.mb_str());
    
    if (result == 0 && wxFileExists(zipFileName))
    {
        wxMessageBox("Files successfully compressed to: " + zipFileName + "\n");
        wxMessageBox("Archive contains:\n");
        for (size_t i = 0; i < filesToZip.GetCount(); i++)
        {
            wxMessageBox("  - " + filesToZip[i] + "\n");
        }
        
        wxMessageBox("Files successfully compressed to:\n" + zipFileName, 
                     "Success", wxOK | wxICON_INFORMATION);
        
        // Open the directory containing the zip file
        #ifdef __WXMSW__
            wxExecute("explorer /select," + wxGetCwd() + "\\" + zipFileName);
        #elif __WXMAC__
            wxExecute("open -R " + wxGetCwd() + "/" + zipFileName);
        #elif __WXGTK__
            wxExecute("xdg-open " + wxGetCwd());
        #endif
    }
    else
    {
        wxMessageBox("Error: Failed to create zip file.\n");
        wxMessageBox("Failed to create zip file. Make sure the zip command is available on your system.", 
                     "Error", wxOK | wxICON_ERROR);
    }
}