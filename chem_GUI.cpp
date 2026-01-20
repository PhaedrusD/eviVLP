#include <wx/wx.h>
#include <wx/filedlg.h>
#include <wx/choice.h>
#include <wx/notebook.h>
#include "FileProcessor.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

class MyFrame : public wxFrame
{
public:
    MyFrame(const wxString& title);

private:
    wxString pdbFilePath;
    wxString pyFilePath;
    wxString incomingPdbPath;
    wxString mol2FilePath;
    wxString strFilePath;
    wxTextCtrl* statusText;
    wxChoice* aminoAcidChoice;
    wxChoice* vlpTypeChoice;  // New: VLP type selector
    
    // Text controls for bonded and deleted atoms
    wxTextCtrl* bondedAtom1Text;
    wxTextCtrl* bondedAtom2Text;
    wxTextCtrl* bondedAtom3Text;
    wxTextCtrl* deletedAtom1Text;
    wxTextCtrl* deletedAtom2Text;
    wxTextCtrl* deletedAtom3Text;
    
    FileProcessor processor;

    // Panel creation methods
    wxPanel* CreateSetupPanel(wxNotebook* notebook);
    wxPanel* CreateAtomConfigPanel(wxNotebook* notebook);
    wxPanel* CreateTopologyPanel(wxNotebook* notebook);
    wxPanel* CreateCapsidPanel(wxNotebook* notebook);

    // Existing methods
    void SelectFile(wxString& target, const wxString& message, const wxString& wildcard);
    void OnSelectAminoAcidAndPDB(wxCommandEvent& event);
    void OnDownloadOneMol2(wxCommandEvent& event);
    void OnUploadOneStr(wxCommandEvent& event);
    void OnRunRenaming(wxCommandEvent& event);
    void OnSelectBondedAndDeleted(wxCommandEvent& event);
    void OnProcessFile(wxCommandEvent& event);
    void OnProcessPatchMaker(wxCommandEvent& event);
    void OnCreateTopology(wxCommandEvent& event);
    void OnDownloadNAAMol2(wxCommandEvent& event);
    void OnUploadNAAStr(wxCommandEvent& event);
    void OnCreateViewCapsid(wxCommandEvent& event);
    void OnMinimizeCapsid(wxCommandEvent& event);
    void CreateViewQbeta();  // New helper method
    void CreateViewTMV();    // New helper method
    void OnZipFiles(wxCommandEvent& event);

    bool CopyFile(const wxString& source, const wxString& destination);
    wxString GetAminoAcidAbbreviation(const wxString& fullName);

    wxDECLARE_EVENT_TABLE();
};

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_BUTTON(1001, MyFrame::OnSelectAminoAcidAndPDB)
    EVT_BUTTON(1002, MyFrame::OnRunRenaming)
    EVT_BUTTON(1003, MyFrame::OnDownloadOneMol2)
    EVT_BUTTON(1004, MyFrame::OnUploadOneStr)
    EVT_BUTTON(1005, MyFrame::OnSelectBondedAndDeleted)
    EVT_BUTTON(1006, MyFrame::OnProcessFile)
    EVT_BUTTON(1007, MyFrame::OnProcessPatchMaker)
    EVT_BUTTON(1008, MyFrame::OnCreateTopology)
    EVT_BUTTON(1009, MyFrame::OnDownloadNAAMol2)
    EVT_BUTTON(1010, MyFrame::OnUploadNAAStr)
    EVT_BUTTON(1011, MyFrame::OnCreateViewCapsid)
    EVT_BUTTON(1012, MyFrame::OnMinimizeCapsid)
    EVT_BUTTON(1013, MyFrame::OnZipFiles)  // New button
wxEND_EVENT_TABLE()

MyFrame::MyFrame(const wxString& title)
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(800, 700)),
      processor("IC_table_noCGENFF.cpp")
{
    // Create main sizer
    wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
    
    // Create notebook (tabbed interface)
    wxNotebook* notebook = new wxNotebook(this, wxID_ANY);
    
    // Create all tabs
    wxPanel* setupPanel = CreateSetupPanel(notebook);
    wxPanel* atomConfigPanel = CreateAtomConfigPanel(notebook);
    wxPanel* topologyPanel = CreateTopologyPanel(notebook);
    wxPanel* capsidPanel = CreateCapsidPanel(notebook);
    
    // Add tabs to notebook
    notebook->AddPage(setupPanel, "1. Initial Setup", true);
    notebook->AddPage(atomConfigPanel, "2. Atom Configuration", false);
    notebook->AddPage(topologyPanel, "3. Topology & Patch", false);
    notebook->AddPage(capsidPanel, "4. Capsid Operations", false);
    
    // Create status text control (shared across all tabs)
    statusText = new wxTextCtrl(this, wxID_ANY, "Status Log:\n", 
                               wxDefaultPosition, wxSize(-1, 200), 
                               wxTE_MULTILINE | wxTE_READONLY);
    
    // Add notebook and status to main sizer
    mainSizer->Add(notebook, 1, wxEXPAND | wxALL, 5);
    mainSizer->Add(new wxStaticText(this, wxID_ANY, "Status Log:"), 
                   0, wxLEFT | wxRIGHT | wxTOP, 5);
    mainSizer->Add(statusText, 0, wxEXPAND | wxALL, 5);
    
    SetSizer(mainSizer);
}

wxPanel* MyFrame::CreateSetupPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    // Add some spacing at top
    sizer->AddSpacer(10);
    
    // Amino acid selection section
    wxStaticBoxSizer* aminoAcidBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Amino Acid Selection");
    wxString aminoAcids[] = {"Lysine", "Cysteine", "Tyrosine", "N-Terminus"};
    aminoAcidChoice = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxDefaultSize, 4, aminoAcids);
    aminoAcidChoice->SetSelection(0);
    aminoAcidBox->Add(aminoAcidChoice, 0, wxALL | wxEXPAND, 5);
    sizer->Add(aminoAcidBox, 0, wxALL | wxEXPAND, 10);
    
    // File selection section
    wxStaticBoxSizer* fileBox = new wxStaticBoxSizer(wxVERTICAL, panel, "File Operations");
    wxButton* selectDrugButton = new wxButton(panel, 1001, "Select Incoming Drug");
    fileBox->Add(selectDrugButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(fileBox, 0, wxALL | wxEXPAND, 10);
    
    // MOL2/STR workflow section
    wxStaticBoxSizer* mol2StrBox = new wxStaticBoxSizer(wxVERTICAL, panel, "CGenFF Workflow");
    wxStaticText* instructionText = new wxStaticText(panel, wxID_ANY, 
        "1. Download one.mol2\n"
        "2. Upload to CGenFF website: https://cgenff.silcsbio.com/\n"
        "3. Download the one.str file\n"
        "4. Upload one.str back here\n"
        "5. Run Renaming");
    instructionText->Wrap(400);
    mol2StrBox->Add(instructionText, 0, wxALL, 5);
    
    wxButton* downloadOneMol2Button = new wxButton(panel, 1003, "Download one.mol2");  // ID 1003
    wxButton* uploadOneStrButton = new wxButton(panel, 1004, "Upload one.str");        // ID 1004
    wxButton* runRenamingButton = new wxButton(panel, 1002, "Run Renaming");           // ID 1002
    mol2StrBox->Add(downloadOneMol2Button, 0, wxALL | wxEXPAND, 5);
    mol2StrBox->Add(uploadOneStrButton, 0, wxALL | wxEXPAND, 5);
    mol2StrBox->Add(runRenamingButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(mol2StrBox, 0, wxALL | wxEXPAND, 10);
    
    sizer->AddStretchSpacer();
    panel->SetSizer(sizer);
    return panel;
}

wxPanel* MyFrame::CreateAtomConfigPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    sizer->AddSpacer(10);
    
    // Bonded atoms section
    wxStaticBoxSizer* bondedBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Bonded Atoms");
    wxStaticText* bondedInstr = new wxStaticText(panel, wxID_ANY,
        "Select atom to be directly bonded to the amino acid,\n"
        "then a second and third atoms consecutively bonded to that atom.");
    bondedInstr->Wrap(400);
    bondedBox->Add(bondedInstr, 0, wxALL, 5);
    
    wxFlexGridSizer* bondedGrid = new wxFlexGridSizer(3, 2, 5, 10);
    bondedGrid->AddGrowableCol(1, 1);
    
    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "1st (directly bonded):"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom1Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    bondedGrid->Add(bondedAtom1Text, 1, wxEXPAND);
    
    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "2nd (linear):"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom2Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    bondedGrid->Add(bondedAtom2Text, 1, wxEXPAND);
    
    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "3rd (linear):"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom3Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    bondedGrid->Add(bondedAtom3Text, 1, wxEXPAND);
    
    bondedBox->Add(bondedGrid, 0, wxALL | wxEXPAND, 5);
    sizer->Add(bondedBox, 0, wxALL | wxEXPAND, 10);
    
    // Deleted atoms section
    wxStaticBoxSizer* deletedBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Deleted Atoms");
    wxStaticText* deletedInstr = new wxStaticText(panel, wxID_ANY,
        "Select all atoms to be deleted (up to three):");
    deletedBox->Add(deletedInstr, 0, wxALL, 5);
    
    wxFlexGridSizer* deletedGrid = new wxFlexGridSizer(3, 2, 5, 10);
    deletedGrid->AddGrowableCol(1, 1);
    
    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "1st:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom1Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    deletedGrid->Add(deletedAtom1Text, 1, wxEXPAND);
    
    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "2nd:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom2Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    deletedGrid->Add(deletedAtom2Text, 1, wxEXPAND);
    
    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "3rd:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom3Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(150, -1));
    deletedGrid->Add(deletedAtom3Text, 1, wxEXPAND);
    
    deletedBox->Add(deletedGrid, 0, wxALL | wxEXPAND, 5);
    sizer->Add(deletedBox, 0, wxALL | wxEXPAND, 10);
    
    // Action buttons
    wxButton* selectBondedDeletedButton = new wxButton(panel, 1005, "Confirm Atom Selection");
    wxButton* processFileButton = new wxButton(panel, 1006, "Process and Orient Drug");
    sizer->Add(selectBondedDeletedButton, 0, wxALL | wxEXPAND, 10);
    sizer->Add(processFileButton, 0, wxALL | wxEXPAND, 10);
    
    sizer->AddStretchSpacer();
    panel->SetSizer(sizer);
    return panel;
}

wxPanel* MyFrame::CreateTopologyPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    sizer->AddSpacer(10);
    
    // Patch maker section
    wxStaticBoxSizer* patchBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Patch Generation");
    wxButton* processPatchMakerButton = new wxButton(panel, 1007, "Process Patch Maker");
    patchBox->Add(processPatchMakerButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(patchBox, 0, wxALL | wxEXPAND, 10);
    
    // NAA MOL2/STR workflow section
    wxStaticBoxSizer* naaMol2StrBox = new wxStaticBoxSizer(wxVERTICAL, panel, "NAA CGenFF Workflow");
    wxStaticText* naaInstructionText = new wxStaticText(panel, wxID_ANY, 
        "1. Download NAA.mol2\n"
        "2. Upload to CGenFF website: https://cgenff.silcsbio.com/\n"
        "3. Download the NAA.str file\n"
        "4. Upload NAA.str back here");
    naaInstructionText->Wrap(400);
    naaMol2StrBox->Add(naaInstructionText, 0, wxALL, 5);
    
    wxButton* downloadNAAMol2Button = new wxButton(panel, 1009, "Download NAA.mol2");
    wxButton* uploadNAAStrButton = new wxButton(panel, 1010, "Upload NAA.str");
    naaMol2StrBox->Add(downloadNAAMol2Button, 0, wxALL | wxEXPAND, 5);
    naaMol2StrBox->Add(uploadNAAStrButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(naaMol2StrBox, 0, wxALL | wxEXPAND, 10);
    
    // Topology creation section
    wxStaticBoxSizer* topologyBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Novel Amino Acid Creation");
    wxButton* createTopologyButton = new wxButton(panel, 1008, "Create Novel Amino Acid");
    topologyBox->Add(createTopologyButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(topologyBox, 0, wxALL | wxEXPAND, 10);
    sizer->AddStretchSpacer();
    panel->SetSizer(sizer);
    return panel;
}

wxPanel* MyFrame::CreateCapsidPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    sizer->AddSpacer(10);
    
    // VLP type selection
    wxStaticBoxSizer* vlpBox = new wxStaticBoxSizer(wxVERTICAL, panel, "VLP Type Selection");
    wxString vlpTypes[] = {"Qbeta", "TMV"};  // Added TMV
    vlpTypeChoice = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxDefaultSize, 2, vlpTypes);  // Changed from 1 to 2
    vlpTypeChoice->SetSelection(0);
    vlpBox->Add(vlpTypeChoice, 0, wxALL | wxEXPAND, 5);
    sizer->Add(vlpBox, 0, wxALL | wxEXPAND, 10);
    
    // Capsid operations section
        wxStaticBoxSizer* capsidBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Capsid Operations");
    wxButton* createViewCapsidButton = new wxButton(panel, 1011, "Create And View Capsid");
    wxButton* minimizeCapsidButton = new wxButton(panel, 1012, "Minimize Capsid");
    wxButton* zipFilesButton = new wxButton(panel, 1013, "Package Files (Zip)");  // New button
    capsidBox->Add(createViewCapsidButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(minimizeCapsidButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(zipFilesButton, 0, wxALL | wxEXPAND, 5);  // Add the new button
    sizer->Add(capsidBox, 0, wxALL | wxEXPAND, 10);
    
    sizer->AddStretchSpacer();
    panel->SetSizer(sizer);
    return panel;
}

// Helper function to get amino acid abbreviation
wxString MyFrame::GetAminoAcidAbbreviation(const wxString& fullName)
{
    if (fullName == "Lysine") return "lys";
    if (fullName == "Cysteine") return "cys";
    if (fullName == "Tyrosine") return "tyr";
    if (fullName == "N-Terminus") return "ntrm";
    return fullName.Lower();
}

// SelectFile function definition
void MyFrame::SelectFile(wxString& target, const wxString& message, const wxString& wildcard)
{
    wxFileDialog openFileDialog(this, message, "", "", wildcard, wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (openFileDialog.ShowModal() == wxID_CANCEL) return;
    target = openFileDialog.GetPath();
    statusText->AppendText("Selected file: " + target + "\n");
}

// Helper function to copy file
bool MyFrame::CopyFile(const wxString& source, const wxString& destination)
{
    std::ifstream src(source.ToStdString(), std::ios::binary);
    std::ofstream dst(destination.ToStdString(), std::ios::binary);
    
    if (!src.is_open() || !dst.is_open()) {
        return false;
    }
    
    dst << src.rdbuf();
    return true;
}

void MyFrame::OnDownloadOneMol2(wxCommandEvent& event)
{
    if (!wxFileExists("one.pdb"))
    {
        wxMessageBox("one.pdb file not found. Please run the renaming step first.", "Error", wxOK | wxICON_ERROR);
        return;
    }
    
    wxString command = "obabel one.pdb -O one.mol2";
    statusText->AppendText("Running command: " + command + "\n");
    
    int result = wxExecute(command, wxEXEC_SYNC);
    
    if (result == 0 && wxFileExists("one.mol2"))
    {
        statusText->AppendText("one.mol2 file generated successfully!\n");
        statusText->AppendText("\n*** NEXT STEP ***\n");
        statusText->AppendText("1. Locate one.mol2 in your working directory\n");
        statusText->AppendText("2. Upload it to the CGenFF website: https://cgenff.silcsbio.com/\n");
        statusText->AppendText("3. Download the one.str file from CGenFF\n");
        statusText->AppendText("4. Click 'Upload one.str' button to select the downloaded file\n\n");
        
        #ifdef __WXMSW__
            wxExecute("explorer /select," + wxGetCwd() + "\\one.mol2");
        #elif __WXMAC__
            wxExecute("open -R " + wxGetCwd() + "/one.mol2");
        #elif __WXGTK__
            wxExecute("xdg-open " + wxGetCwd());
        #endif
    }
    else
    {
        statusText->AppendText("Error: one.mol2 generation failed.\n");
        wxMessageBox("Failed to generate one.mol2 file. Make sure obabel is installed.", "Error", wxOK | wxICON_ERROR);
    }
}

void MyFrame::OnUploadOneStr(wxCommandEvent& event)
{
    wxFileDialog openFileDialog(this, "Select one.str file", "", "",
                               "STR files (*.str)|*.str", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    
    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;
    
    wxString selectedStrPath = openFileDialog.GetPath();
    wxString targetPath = wxGetCwd() + "/one.str";
    
    if (wxCopyFile(selectedStrPath, targetPath, true))
    {
        statusText->AppendText("one.str file uploaded successfully: " + selectedStrPath + "\n");
        statusText->AppendText("File copied to: " + targetPath + "\n");
        statusText->AppendText("You can now proceed to the next step.\n\n");
    }
    else
    {
        statusText->AppendText("Error: Failed to copy one.str file.\n");
        wxMessageBox("Failed to upload one.str file.", "Error", wxOK | wxICON_ERROR);
    }
}

void MyFrame::OnDownloadNAAMol2(wxCommandEvent& event)
{
    if (!wxFileExists("patched.pdb"))
    {
        wxMessageBox("patched.pdb file not found. Please complete the topology creation step first.", "Error", wxOK | wxICON_ERROR);
        return;
    }
    
    wxString command = "obabel patched.pdb -O naa.mol2";
    statusText->AppendText("Running command: " + command + "\n");
    
    int result = wxExecute(command, wxEXEC_SYNC);
    
    if (result == 0 && wxFileExists("naa.mol2"))
    {
        statusText->AppendText("naa.mol2 file generated successfully!\n");
        statusText->AppendText("\n*** NEXT STEP ***\n");
        statusText->AppendText("1. Locate naa.mol2 in your working directory\n");
        statusText->AppendText("2. Upload it to the CGenFF website: https://cgenff.silcsbio.com/\n");
        statusText->AppendText("3. Download the naa.str file from CGenFF\n");
        statusText->AppendText("4. Click 'Upload naa.str' button to select the downloaded file\n\n");
        
        #ifdef __WXMSW__
            wxExecute("explorer /select," + wxGetCwd() + "\\naa.mol2");
        #elif __WXMAC__
            wxExecute("open -R " + wxGetCwd() + "/naa.mol2");
        #elif __WXGTK__
            wxExecute("xdg-open " + wxGetCwd());
        #endif
    }
    else
    {
        statusText->AppendText("Error: naa.mol2 generation failed.\n");
        wxMessageBox("Failed to generate naa.mol2 file. Make sure obabel is installed.", "Error", wxOK | wxICON_ERROR);
    }
}

void MyFrame::OnUploadNAAStr(wxCommandEvent& event)
{
    wxFileDialog openFileDialog(this, "Select naa.str file", "", "",
                               "STR files (*.str)|*.str", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    
    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;
    
    wxString selectedStrPath = openFileDialog.GetPath();
    wxString targetPath = wxGetCwd() + "/naa.str";
    
    if (wxCopyFile(selectedStrPath, targetPath, true))
    {
        statusText->AppendText("naa.str file uploaded successfully: " + selectedStrPath + "\n");
        statusText->AppendText("File copied to: " + targetPath + "\n");
        statusText->AppendText("You can now proceed to create and view the capsid.\n\n");
    }
    else
    {
        statusText->AppendText("Error: Failed to copy naa.str file.\n");
        wxMessageBox("Failed to upload naa.str file.", "Error", wxOK | wxICON_ERROR);
    }
}

void MyFrame::OnSelectAminoAcidAndPDB(wxCommandEvent& event)
{
    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);
    statusText->AppendText("Selected amino acid: " + selectedAminoAcid + "\n");

    SelectFile(incomingPdbPath, _("Select incoming PDB file"), "PDB files (*.pdb)|*.pdb");
    
    if (!incomingPdbPath.empty()) {
        wxString onePdbPath = "one.pdb";
        if (CopyFile(incomingPdbPath, onePdbPath)) {
            statusText->AppendText("Successfully copied " + incomingPdbPath + " to one.pdb\n");
            
            wxString command = "obabel one.pdb -O one.mol2";
            statusText->AppendText("Running command: " + command + "\n");
            
            int result = system(command.mb_str());
            
            if (result == 0) {
                statusText->AppendText("Successfully converted one.pdb to one.mol2\n");
            } else {
                statusText->AppendText("Error: obabel conversion failed. Make sure obabel is installed and in your PATH.\n");
            }
        } else {
            statusText->AppendText("Error: Failed to copy file to one.pdb\n");
        }
    } else {
        statusText->AppendText("No incoming PDB file selected.\n");
    }
}

void MyFrame::OnRunRenaming(wxCommandEvent& event)
{
    statusText->AppendText("Starting atom renaming process...\n");
    
    statusText->AppendText("Resetting top_all36_cgenff_CBD.rtf from template...\n");
    
    std::ifstream templateFile("top_all36_cgenff_CBD_template.rtf");
    if (!templateFile.is_open()) {
        statusText->AppendText("Error: Failed to open template file top_all36_cgenff_CBD_template.rtf\n");
        return;
    }
    
    std::ofstream targetFile("top_all36_cgenff_CBD.rtf");
    if (!targetFile.is_open()) {
        statusText->AppendText("Error: Failed to create target file top_all36_cgenff_CBD.rtf\n");
        templateFile.close();
        return;
    }
    
    targetFile << templateFile.rdbuf();
    templateFile.close();
    targetFile.close();
    
    statusText->AppendText("Successfully reset top_all36_cgenff_CBD.rtf from template.\n");
    
    wxString command = "python3 fix_atom_naming_pdb.py";
    statusText->AppendText("Running command: " + command + "\n");
    
    int result = system(command.mb_str());
    
    if (result == 0) {
        statusText->AppendText("Atom renaming script executed successfully.\n");
    } else {
        statusText->AppendText("Error: Atom renaming script failed. Make sure fix_atom_naming_pdb.py exists and python3 is installed.\n");
        return;
    }
    
    statusText->AppendText("Running ./IC_table_one command...\n");
    
    wxString icTableCommand = "./IC_table";
    int icResult = system(icTableCommand.mb_str());
    
    if (icResult == 0) {
        statusText->AppendText("IC_table command executed successfully.\n");
        statusText->AppendText("Atom renaming process completed successfully!\n");
    } else {
        statusText->AppendText("Error: IC_table command failed. Make sure IC_table executable exists and has proper permissions.\n");
    }
}

void MyFrame::OnSelectBondedAndDeleted(wxCommandEvent& event)
{
    wxString bonded1 = bondedAtom1Text->GetValue().Trim();
    wxString bonded2 = bondedAtom2Text->GetValue().Trim();
    wxString bonded3 = bondedAtom3Text->GetValue().Trim();
    wxString deleted1 = deletedAtom1Text->GetValue().Trim();
    wxString deleted2 = deletedAtom2Text->GetValue().Trim();
    wxString deleted3 = deletedAtom3Text->GetValue().Trim();
    
    if (bonded1.IsEmpty()) {
        statusText->AppendText("Error: At least the first bonded atom must be specified.\n");
        return;
    }
    
    statusText->AppendText("=== Bonded and Deleted Atoms Selected ===\n");
    statusText->AppendText("Bonded Atoms:\n");
    statusText->AppendText("  1st (directly bonded): " + bonded1 + "\n");
    if (!bonded2.IsEmpty()) {
        statusText->AppendText("  2nd (linear): " + bonded2 + "\n");
    }
    if (!bonded3.IsEmpty()) {
        statusText->AppendText("  3rd (linear): " + bonded3 + "\n");
    }
    
    statusText->AppendText("Deleted Atoms:\n");
    if (!deleted1.IsEmpty()) {
        statusText->AppendText("  1st: " + deleted1 + "\n");
    }
    if (!deleted2.IsEmpty()) {
        statusText->AppendText("  2nd: " + deleted2 + "\n");
    }
    if (!deleted3.IsEmpty()) {
        statusText->AppendText("  3rd: " + deleted3 + "\n");
    }
    
    statusText->AppendText("Atom selection completed successfully.\n\n");
}

void MyFrame::OnProcessFile(wxCommandEvent& event)
{
    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);
    
    wxString firstBondedAtom = bondedAtom1Text->GetValue().Trim();
    wxString thirdBondedAtom = bondedAtom3Text->GetValue().Trim();
    
    if (firstBondedAtom.IsEmpty()) {
        statusText->AppendText("Error: First bonded atom must be specified for processing.\n");
        return;
    }
    
    statusText->AppendText("Processing amino acid file...\n");
    statusText->AppendText("Amino acid: " + selectedAminoAcid + "\n");
    statusText->AppendText("First bonded atom: " + firstBondedAtom + "\n");
    
    if (!thirdBondedAtom.IsEmpty()) {
        statusText->AppendText("Third bonded atom: " + thirdBondedAtom + "\n");
    } else {
        statusText->AppendText("Third bonded atom: Not specified\n");
    }
    
    bool success = processor.ProcessAminoAcidFile(std::string(selectedAminoAcid.mb_str()),
                                                  std::string(firstBondedAtom.mb_str()),
                                                  std::string(thirdBondedAtom.mb_str()));
    
    if (success) {
        statusText->AppendText("Amino acid file processed successfully!\n");
        statusText->AppendText("TEST_ATOM1 replaced with: " + firstBondedAtom + "\n");
        if (!thirdBondedAtom.IsEmpty()) {
            statusText->AppendText("TEST_ATOM3 replaced with: " + thirdBondedAtom + "\n");
        }
        
        wxString aminoAcidAbbrev = GetAminoAcidAbbreviation(selectedAminoAcid);
        wxString pythonCommand = "python3 orient_place_drug_" + aminoAcidAbbrev + ".py";
        statusText->AppendText("Running command: " + pythonCommand + "\n");
        
        int result = system(pythonCommand.mb_str());
        
        if (result == 0) {
            statusText->AppendText("Orient place drug script executed successfully.\n");
            
            statusText->AppendText("Resetting amino acid file...\n");
            processor.ResetAminoAcidFile(std::string(selectedAminoAcid.mb_str()));
            statusText->AppendText("Amino acid file reset completed.\n");
        } else {
            statusText->AppendText("Error: Orient place drug script failed. Make sure orient_place_drug_" + aminoAcidAbbrev + ".py exists and python3 is installed.\n");
            statusText->AppendText("Skipping file reset due to script failure.\n");
        }
    } else {
        statusText->AppendText("Error: Failed to process amino acid file.\n");
    }
}

void MyFrame::OnProcessPatchMaker(wxCommandEvent& event)
{
    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);
    
    wxString bonded1 = bondedAtom1Text->GetValue().Trim();
    wxString bonded2 = bondedAtom2Text->GetValue().Trim();
    wxString bonded3 = bondedAtom3Text->GetValue().Trim();
    
    wxString deleted1 = deletedAtom1Text->GetValue().Trim();
    wxString deleted2 = deletedAtom2Text->GetValue().Trim();
    wxString deleted3 = deletedAtom3Text->GetValue().Trim();
    
    if (bonded1.IsEmpty()) {
        statusText->AppendText("Error: At least the first bonded atom must be specified for patch maker processing.\n");
        return;
    }
    
    std::vector<std::string> bondedAtoms;
    std::vector<std::string> deletedAtoms;
    
    if (!bonded1.IsEmpty()) {
        bondedAtoms.push_back(std::string(bonded1.mb_str()));
    }
    if (!bonded2.IsEmpty()) {
        bondedAtoms.push_back(std::string(bonded2.mb_str()));
    }
    if (!bonded3.IsEmpty()) {
        bondedAtoms.push_back(std::string(bonded3.mb_str()));
    }
    
    if (!deleted1.IsEmpty()) {
        deletedAtoms.push_back(std::string(deleted1.mb_str()));
    }
    if (!deleted2.IsEmpty()) {
        deletedAtoms.push_back(std::string(deleted2.mb_str()));
    }
    if (!deleted3.IsEmpty()) {
        deletedAtoms.push_back(std::string(deleted3.mb_str()));
    }
    
    statusText->AppendText("Processing patch maker file...\n");
    statusText->AppendText("Amino acid: " + selectedAminoAcid + "\n");
    statusText->AppendText("Bonded atoms: ");
    for (const auto& atom : bondedAtoms) {
        statusText->AppendText(wxString(atom) + " ");
    }
    statusText->AppendText("\n");
    
    if (!deletedAtoms.empty()) {
        statusText->AppendText("Deleted atoms: ");
        for (const auto& atom : deletedAtoms) {
            statusText->AppendText(wxString(atom) + " ");
        }
        statusText->AppendText("\n");
    }
    
    bool success = processor.ProcessPatchMakerFile(std::string(selectedAminoAcid.mb_str()),
                                                   bondedAtoms,
                                                   deletedAtoms);
    
    if (success) {
        statusText->AppendText("Patch maker file processed successfully!\n");
        
        wxString aminoAcidAbbrev = GetAminoAcidAbbreviation(selectedAminoAcid);
        wxString pythonCommand = "python3 patch_maker_" + aminoAcidAbbrev + ".py";
        statusText->AppendText("Running command: " + pythonCommand + "\n");
        
        int result = system(pythonCommand.mb_str());
        
        if (result == 0) {
            statusText->AppendText("Patch maker script executed successfully.\n");
            
            statusText->AppendText("Resetting patch maker file...\n");
            processor.ResetPatchMakerFile(std::string(selectedAminoAcid.mb_str()));
            statusText->AppendText("Patch maker file reset completed.\n");
            
            wxString vmdCommand = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e patched_" + aminoAcidAbbrev + ".pgn";
            statusText->AppendText("Running VMD command: " + vmdCommand + "\n");
            
            int vmdResult = system(vmdCommand.mb_str());
            
            if (vmdResult == 0) {
                statusText->AppendText("VMD command executed successfully.\n");
            } else {
                statusText->AppendText("Error: VMD command failed.\n");
            }
            
            wxString pythonCommand = "python3 element_column_fix.py patched.pdb naa.pdb";
            statusText->AppendText("Running command: " + pythonCommand + "\n");
            int result = wxExecute(pythonCommand, wxEXEC_SYNC);
            if (result == 0)
            {
                statusText->AppendText("Python script executed successfully.\n");
            }
            else
            {
                statusText->AppendText("Error: Python script execution failed.\n");
            }
            wxString obabelCommand = "obabel naa.pdb -O naa.mol2";
            statusText->AppendText("Running obabel command: " + obabelCommand + "\n");
            
            int obabelResult = system(obabelCommand.mb_str());
            
            if (obabelResult == 0) {
                statusText->AppendText("Obabel command executed successfully.\n");
            } else {
                statusText->AppendText("Error: Obabel command failed.\n");
            }
            
        } else {
            statusText->AppendText("Error: Patch maker script failed. Make sure patch_maker_" + aminoAcidAbbrev + ".py exists and python3 is installed.\n");
        }
    } else {
        statusText->AppendText("Error: Failed to process patch maker file.\n");
    }
}

void MyFrame::OnCreateTopology(wxCommandEvent& event)
{
    statusText->AppendText("Starting final processing...\n");
    
    wxString command1 = "python3 noX.py";
    statusText->AppendText("Running command: " + command1 + "\n");
    
    int result1 = system(command1.mb_str());
    
    if (result1 == 0) {
        statusText->AppendText("noX.py executed successfully.\n");
        
        wxString command2 = "python3 noX_psf.py";
        statusText->AppendText("Running command: " + command2 + "\n");
        
        int result2 = system(command2.mb_str());
        
        if (result2 == 0) {
            statusText->AppendText("noX_psf.py executed successfully.\n");
            
            wxString command0 = "python3 fix_atom_naming_pdb_naa.py";
            statusText->AppendText("Running command: " + command0 + "\n");
            
            int result0 = system(command0.mb_str());
            
            if (result0 == 0) {
                statusText->AppendText("fix_atom_naming_pdb_naa.py executed successfully.\n");
                
                wxString command3 = "./IC_naa";
                statusText->AppendText("Running command: " + command3 + "\n");
                
                int result3 = system(command3.mb_str());
                
                if (result3 == 0) {
                    statusText->AppendText("IC_naa executed successfully.\n");
                    statusText->AppendText("Final processing completed successfully!\n");
                } else {
                    statusText->AppendText("Error: IC_naa failed. Make sure IC_naa executable exists and has proper permissions.\n");
                }
            } else {
                statusText->AppendText("Error: fix_atom_naming_pdb_naa.py failed. Make sure fix_atom_naming_pdb_naa.py exists and python3 is installed.\n");
                statusText->AppendText("Skipping IC_naa due to script failure.\n");
            }
        } else {
            statusText->AppendText("Error: noX_psf.py failed. Make sure noX_psf.py exists and python3 is installed.\n");
            statusText->AppendText("Skipping remaining steps due to script failure.\n");
        }
    } else {
        statusText->AppendText("Error: noX.py failed. Make sure noX.py exists and python3 is installed.\n");
        statusText->AppendText("Skipping remaining steps due to script failure.\n");
    }
}

void MyFrame::OnCreateViewCapsid(wxCommandEvent& event)
{
    // Check which VLP type is selected
    int vlpSelection = vlpTypeChoice->GetSelection();
    
    if (vlpSelection == 0) {
        // Qbeta workflow
        CreateViewQbeta();
    } else if (vlpSelection == 1) {
        // TMV workflow
        CreateViewTMV();
    } else {
        statusText->AppendText("Error: Unknown VLP type selected.\n");
    }

}
void MyFrame::CreateViewQbeta()
{
    // All your existing OnCreateViewCapsid code goes here unchanged
    wxString command1 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e fix_a.pgn";
    statusText->AppendText("Running command: " + command1 + "\n");
    
    int result1 = system(command1.mb_str());
    
    if (result1 == 0) {
        statusText->AppendText("fix_a.pgn executed successfully.\n");
        
        wxString command2 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e fix_b.pgn";
        statusText->AppendText("Running command: " + command2 + "\n");
        
        int result2 = system(command2.mb_str());
        
        if (result2 == 0) {
            statusText->AppendText("fix_b.pgn executed successfully.\n");
            
            wxString command3 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e fix_c.pgn";
            statusText->AppendText("Running command: " + command3 + "\n");
            
            int result3 = system(command3.mb_str());
            
            if (result3 == 0) {
                statusText->AppendText("fix_c.pgn executed successfully.\n");
                
                wxString command4 = "python3 Qbeta_capsid_maker.py";
                statusText->AppendText("Running command: " + command4 + "\n");
                
                int result4 = system(command4.mb_str());
                
                if (result4 == 0) {
                    statusText->AppendText("Qbeta_capsid_maker.py executed successfully.\n");
                    
                    wxString command5 = "python3 Qbeta_pgn_writer.py";
                    statusText->AppendText("Running command: " + command5 + "\n");
                    
                    int result5 = system(command5.mb_str());
                    
                    if (result5 == 0) {
                        wxString command6 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e QB_capsid.pgn";
                        wxString command7 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -e capsid_zoom_naa.tcl";
                        statusText->AppendText("Running command: " + command6 + "\n");

                        int result6 = system(command6.mb_str());

                        if (result6 == 0) {
                            statusText->AppendText("QB_capsid.pgn executed successfully.\n");
                            
                            statusText->AppendText("Running command: " + command7 + "\n");
                            int result7 = system(command7.mb_str());
                            
                            if (result7 == 0) {
                                statusText->AppendText("capsid_zoom_naa.tcl executed successfully.\n");
                                statusText->AppendText("Capsid creation and viewing process completed successfully!\n");
                            } else {
                                statusText->AppendText("Error: capsid_zoom_naa.tcl failed. Make sure the file exists and VMD is installed correctly.\n");
                            }
                        } else {
                            statusText->AppendText("Error: QB_capsid.pgn failed. Make sure the file exists and VMD is installed correctly.\n");
                        }
                    } else {
                        statusText->AppendText("Error: Qbeta_pgn_writer.py failed. Make sure the script exists and python3 is installed.\n");
                    }
                } else {
                    statusText->AppendText("Error: Qbeta_capsid_maker.py failed. Make sure the script exists and python3 is installed.\n");
                }
            } else {
                statusText->AppendText("Error: fix_c.pgn failed. Make sure the file exists and VMD is installed correctly.\n");
            }
        } else {
            statusText->AppendText("Error: fix_b.pgn failed. Make sure the file exists and VMD is installed correctly.\n");
        }
    } else {
        statusText->AppendText("Error: fix_a.pgn failed. Make sure the file exists and VMD is installed correctly.\n");
    }
}

void MyFrame::CreateViewTMV()
{
    statusText->AppendText("Starting TMV capsid creation and viewing process...\n");
    
    wxString command1 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e TMV_tyr.pgn";
    statusText->AppendText("Running command: " + command1 + "\n");
    int result1 = system(command1.mb_str());
    
    if (result1 == 0) {
        statusText->AppendText("TMV_tyr.pgn executed successfully.\n");  // Fixed message
        
        wxString command2 = "python3 TMV_capsid_maker.py";
        statusText->AppendText("Running command: " + command2 + "\n");
        int result2 = system(command2.mb_str());
        
        if (result2 == 0) {  // Properly nested
            statusText->AppendText("TMV_capsid_maker.py executed successfully.\n");
            
            wxString command3 = "python3 TMV_pgn_writer.py";
            statusText->AppendText("Running command: " + command3 + "\n");
            
            int result3 = system(command3.mb_str());  // Fixed: was result5, should be result3
            
            if (result3 == 0) {  // Properly nested
                wxString command4 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e TMV_all_out.pgn";
                wxString command5 = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -e capsid_zoom_naa_TMV.tcl";
                statusText->AppendText("Running command: " + command4 + "\n");

                int result4 = system(command4.mb_str());  // Fixed: was result6, should be result4

                if (result4 == 0) {  // Fixed: was result6
                    statusText->AppendText("TMV_all_out.pgn executed successfully.\n");  // Fixed message
                    
                    statusText->AppendText("Running command: " + command5 + "\n");
                    int result5 = system(command5.mb_str());  // Fixed: was result7, should be result5
                    
                    if (result5 == 0) {  // Fixed: was result7
                        statusText->AppendText("capsid_zoom_naa.tcl executed successfully.\n");
                        statusText->AppendText("TMV capsid creation and viewing process completed successfully!\n");
                    } else {
                        statusText->AppendText("Error: capsid_zoom_naa.tcl failed.\n");
                    }
                } else {
                    statusText->AppendText("Error: TMV_all_out.pgn failed.\n");
                }
            } else {
                statusText->AppendText("Error: TMV_pgn_writer.py failed.\n");
            }
        } else {
            statusText->AppendText("Error: TMV_capsid_maker.py failed.\n");
        }
    } else {
        statusText->AppendText("Error: TMV_tyr.pgn failed.\n");
    }
}

void MyFrame::OnMinimizeCapsid(wxCommandEvent& event)
{
    wxString command9 = "~/Downloads/NAMD_3.0b5_MacOS-universal-multicore/namd3 +p8 minimize.conf output.log";
    
    statusText->AppendText("Minimize Capsid button clicked!\n");
    statusText->AppendText("Running command: " + command9 + "\n");
    
    int result = system(command9.mb_str());
    
    if (result == 0) {
        statusText->AppendText("NAMD minimize command executed successfully.\n");
    } else {
        statusText->AppendText("Error: NAMD minimize command failed. Make sure NAMD is installed and the path is correct.\n");
    }
}

void MyFrame::OnZipFiles(wxCommandEvent& event)
{
    statusText->AppendText("Starting file compression...\n");
    
    // Define the files to zip
    wxArrayString filesToZip;
    filesToZip.Add("minimize.conf");
    filesToZip.Add("QB_naa_capsid_wb.psf");
    filesToZip.Add("QB_naa_capsid_wb.pdb");
    filesToZip.Add("par_GUI.prm");
    filesToZip.Add("top_all36_cgenff_CBD.rtf");
    
    // Check if all files exist
    bool allFilesExist = true;
    for (size_t i = 0; i < filesToZip.GetCount(); i++)
    {
        if (!wxFileExists(filesToZip[i]))
        {
            statusText->AppendText("Error: File not found: " + filesToZip[i] + "\n");
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
    
    statusText->AppendText("Running command: " + command + "\n");
    
    int result = system(command.mb_str());
    
    if (result == 0 && wxFileExists(zipFileName))
    {
        statusText->AppendText("Files successfully compressed to: " + zipFileName + "\n");
        statusText->AppendText("Archive contains:\n");
        for (size_t i = 0; i < filesToZip.GetCount(); i++)
        {
            statusText->AppendText("  - " + filesToZip[i] + "\n");
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
        statusText->AppendText("Error: Failed to create zip file.\n");
        wxMessageBox("Failed to create zip file. Make sure the zip command is available on your system.", 
                     "Error", wxOK | wxICON_ERROR);
    }
}

class MyApp : public wxApp
{
public:
    virtual bool OnInit() {
        MyFrame* frame = new MyFrame("Create Novel Amino Acid");
        frame->Show(true);
        return true;
    }
};

wxIMPLEMENT_APP(MyApp);
