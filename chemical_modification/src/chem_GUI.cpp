#include <wx/wx.h>
#include <wx/filedlg.h>
#include <wx/choice.h>
#include <wx/notebook.h>
#include <wx/splitter.h>
#include <wx/tglbtn.h>
#include "FileProcessor.h"
#include "MoleculeViewer.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <wx/zipstrm.h>
#include <wx/wfstream.h>
#include <wx/stdpaths.h>
#include <wx/filename.h>
#include <wx/textfile.h>
#include <wx/file.h>

class MyFrame : public wxFrame
{
public:
    MyFrame(const wxString& title);

private:
    wxString pdbFilePath;
    wxString pyFilePath;
    wxString incomingPdbPath;
    wxString mol2FilePath;
    wxString chargeOverrideFile;   // optional hand-curated charges; see OnSelectChargeOverride
    wxStaticText* chargeOverrideLabel = nullptr;
    wxString strFilePath;
    wxTextCtrl* statusText;
    wxChoice* aminoAcidChoice;
    wxChoice* vlpTypeChoice;
    
    // Text controls for bonded and deleted atoms
    wxTextCtrl* bondedAtom1Text;
    wxTextCtrl* bondedAtom2Text;
    wxTextCtrl* bondedAtom3Text;
    wxTextCtrl* deletedAtom1Text;
    wxTextCtrl* deletedAtom2Text;
    wxTextCtrl* deletedAtom3Text;

    // 3D molecule viewer
    MoleculeViewer* molViewer;
    wxToggleButton* bondedModeBtn;
    wxToggleButton* deletedModeBtn;
    wxStaticText* viewerStatusLabel;
    
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
    bool RunOrientStep();   // fill *T1/*T3 in the orient script, run it, reset it
    int  RunVmd(const wxString& command, const wxString& expectedOutput);
    wxString NtrmVariant();          // "ntrm" (Ala, Qbeta) or "ntrm_ser" (Ser, TMV)
    bool StageNtrmReference();       // copy the right model residue to ntrm_ref.pdb
    void OnCreateTopology(wxCommandEvent& event);
    void OnDownloadNAAMol2(wxCommandEvent& event);
    void OnUploadNAAStr(wxCommandEvent& event);
    void OnDownloadNtrmMol2(wxCommandEvent& event);
    void OnSelectChargeOverride(wxCommandEvent& event);
    void OnUploadNtrmStr(wxCommandEvent& event);
    bool MakeMol2(const wxString& pdb, const wxString& psf, const wxString& mol2);
    void SaveMol2Copy(const wxString& mol2);
    void OnCreateViewCapsid(wxCommandEvent& event);
    void OnMinimizeCapsid(wxCommandEvent& event);
    void OnMinimizeMonomer(wxCommandEvent& event);
    bool MinimizeOneMonomer(const wxString& stem);
    void OnBuildMonomer(wxCommandEvent& event);
    wxArrayString MonomerStems();    // fix_a/b/c (Qbeta) or TMV_naa (TMV)
    bool BuildMonomers();            // fill the pgn templates and run psfgen
    void ShowLogTail(const wxString& logFile, int lines = 18);
    void CreateViewQbeta();
    void CreateViewTMV();
    void OnZipFiles(wxCommandEvent& event);

    // New viewer-related methods
    void OnLoadMolecule(wxCommandEvent& event);
    void OnBondedModeToggle(wxCommandEvent& event);
    void OnDeletedModeToggle(wxCommandEvent& event);
    void OnClearBondedSelection(wxCommandEvent& event);
    void OnClearDeletedSelection(wxCommandEvent& event);
    void OnSelectionChanged(SelectionMode mode, const std::vector<std::string>& names);
    void UpdateViewerStatus();

    bool CopyFile(const wxString& source, const wxString& destination);
    wxString GetAminoAcidAbbreviation(const wxString& fullName);

    wxDECLARE_EVENT_TABLE();
};

enum {
    ID_SELECT_DRUG = 1001,
    ID_RUN_RENAMING = 1002,
    ID_DOWNLOAD_MOL2 = 1003,
    ID_UPLOAD_STR = 1004,
    ID_CONFIRM_ATOMS = 1005,
    ID_PROCESS_ORIENT = 1006,
    ID_PATCH_MAKER = 1007,
    ID_CREATE_TOPOLOGY = 1008,
    ID_DOWNLOAD_NAA_MOL2 = 1009,
    ID_UPLOAD_NAA_STR = 1010,
    ID_CREATE_VIEW_CAPSID = 1011,
    ID_MINIMIZE_CAPSID = 1012,
    ID_ZIP_FILES = 1013,
    ID_MINIMIZE_MONOMER = 1032,
    ID_BUILD_MONOMER = 1033,
    ID_LOAD_MOLECULE = 1014,
    ID_BONDED_MODE = 1015,
    ID_DELETED_MODE = 1016,
    ID_CLEAR_BONDED = 1017,
    ID_DOWNLOAD_NTRM_MOL2 = 1030,
    ID_UPLOAD_NTRM_STR = 1031,
    ID_CLEAR_DELETED = 1018,
    ID_APPLY_CHARGE_OVERRIDES = 1034
};

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_BUTTON(ID_SELECT_DRUG, MyFrame::OnSelectAminoAcidAndPDB)
    EVT_BUTTON(ID_RUN_RENAMING, MyFrame::OnRunRenaming)
    EVT_BUTTON(ID_DOWNLOAD_MOL2, MyFrame::OnDownloadOneMol2)
    EVT_BUTTON(ID_UPLOAD_STR, MyFrame::OnUploadOneStr)
    EVT_BUTTON(ID_CONFIRM_ATOMS, MyFrame::OnSelectBondedAndDeleted)
    EVT_BUTTON(ID_PROCESS_ORIENT, MyFrame::OnProcessFile)
    EVT_BUTTON(ID_PATCH_MAKER, MyFrame::OnProcessPatchMaker)
    EVT_BUTTON(ID_CREATE_TOPOLOGY, MyFrame::OnCreateTopology)
    EVT_BUTTON(ID_DOWNLOAD_NAA_MOL2, MyFrame::OnDownloadNAAMol2)
    EVT_BUTTON(ID_UPLOAD_NAA_STR, MyFrame::OnUploadNAAStr)
    EVT_BUTTON(ID_DOWNLOAD_NTRM_MOL2, MyFrame::OnDownloadNtrmMol2)
    EVT_BUTTON(ID_APPLY_CHARGE_OVERRIDES, MyFrame::OnSelectChargeOverride)
    EVT_BUTTON(ID_UPLOAD_NTRM_STR, MyFrame::OnUploadNtrmStr)
    EVT_BUTTON(ID_CREATE_VIEW_CAPSID, MyFrame::OnCreateViewCapsid)
    EVT_BUTTON(ID_MINIMIZE_CAPSID, MyFrame::OnMinimizeCapsid)
    EVT_BUTTON(ID_BUILD_MONOMER, MyFrame::OnBuildMonomer)
    EVT_BUTTON(ID_MINIMIZE_MONOMER, MyFrame::OnMinimizeMonomer)
    EVT_BUTTON(ID_ZIP_FILES, MyFrame::OnZipFiles)
    EVT_BUTTON(ID_LOAD_MOLECULE, MyFrame::OnLoadMolecule)
    EVT_TOGGLEBUTTON(ID_BONDED_MODE, MyFrame::OnBondedModeToggle)
    EVT_TOGGLEBUTTON(ID_DELETED_MODE, MyFrame::OnDeletedModeToggle)
    EVT_BUTTON(ID_CLEAR_BONDED, MyFrame::OnClearBondedSelection)
    EVT_BUTTON(ID_CLEAR_DELETED, MyFrame::OnClearDeletedSelection)
wxEND_EVENT_TABLE()

MyFrame::MyFrame(const wxString& title)
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(900, 800)),
      processor("IC_table_noCGENFF.cpp")
{
    wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
    
    wxNotebook* notebook = new wxNotebook(this, wxID_ANY);
    
    wxPanel* setupPanel = CreateSetupPanel(notebook);
    wxPanel* atomConfigPanel = CreateAtomConfigPanel(notebook);
    wxPanel* topologyPanel = CreateTopologyPanel(notebook);
    wxPanel* capsidPanel = CreateCapsidPanel(notebook);
    
    notebook->AddPage(setupPanel, "1. Initial Setup", true);
    notebook->AddPage(atomConfigPanel, "2. Atom Configuration", false);
    notebook->AddPage(topologyPanel, "3. Topology & Patch", false);
    notebook->AddPage(capsidPanel, "4. Capsid Operations", false);
    
    statusText = new wxTextCtrl(this, wxID_ANY, "Status Log:\n", 
                               wxDefaultPosition, wxSize(-1, 150), 
                               wxTE_MULTILINE | wxTE_READONLY);
    
    mainSizer->Add(notebook, 1, wxEXPAND | wxALL, 5);
    mainSizer->Add(new wxStaticText(this, wxID_ANY, "Status Log:"), 
                   0, wxLEFT | wxRIGHT | wxTOP, 5);
    mainSizer->Add(statusText, 0, wxEXPAND | wxALL, 5);
    
    SetSizer(mainSizer);

    // --- environment check: which python3 / obabel will the pipeline use? ---
    {
        wxArrayString out, err;
        statusText->AppendText("Working directory: " + wxGetCwd() + "\n");
        if (wxExecute("python3 -c \"import sys, numpy; print(sys.executable, 'numpy', numpy.__version__)\"",
                      out, err, wxEXEC_SYNC) == 0 && !out.IsEmpty())
            statusText->AppendText("python3 OK: " + out[0] + "\n");
        else
            statusText->AppendText("WARNING: 'python3 -c \"import numpy\"' failed. Every pipeline script needs NumPy. "
                                   "Fix the active Python environment (e.g. pip install --force-reinstall numpy) "
                                   "or launch the GUI from an environment where it works.\n");
        out.Clear(); err.Clear();
        if (wxExecute("obabel -V", out, err, wxEXEC_SYNC) == 0 && !out.IsEmpty())
            statusText->AppendText("obabel OK: " + out[0] + "\n\n");
        else
            statusText->AppendText("WARNING: obabel not found on PATH.\n\n");
    }
}

wxPanel* MyFrame::CreateSetupPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    sizer->AddSpacer(10);
    
    // Amino acid selection section
    wxStaticBoxSizer* aminoAcidBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Amino Acid Selection");
    wxString aminoAcids[] = {"Lysine", "Cysteine", "Tyrosine"};
    aminoAcidChoice = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxDefaultSize, 3, aminoAcids);
    aminoAcidBox->Add(aminoAcidChoice, 0, wxALL | wxEXPAND, 5);
    sizer->Add(aminoAcidBox, 0, wxALL | wxEXPAND, 10);
    
    // File selection section
    wxStaticBoxSizer* fileBox = new wxStaticBoxSizer(wxVERTICAL, panel, "File Operations");
    wxButton* selectDrugButton = new wxButton(panel, ID_SELECT_DRUG, "Select Incoming Drug");
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
    
    wxButton* downloadOneMol2Button = new wxButton(panel, ID_DOWNLOAD_MOL2, "Download one.mol2");
    wxButton* uploadOneStrButton = new wxButton(panel, ID_UPLOAD_STR, "Import CGenFF result (zip/str)");
    wxButton* runRenamingButton = new wxButton(panel, ID_RUN_RENAMING, "Run Renaming");
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
    wxBoxSizer* mainSizer = new wxBoxSizer(wxHORIZONTAL);

    // ── Left side: 3D molecule viewer ────────────────────────────────────────
    wxBoxSizer* viewerSizer = new wxBoxSizer(wxVERTICAL);

    wxStaticText* viewerTitle = new wxStaticText(panel, wxID_ANY, "Molecule Viewer");
    wxFont titleFont = viewerTitle->GetFont();
    titleFont.SetWeight(wxFONTWEIGHT_BOLD);
    titleFont.SetPointSize(titleFont.GetPointSize() + 1);
    viewerTitle->SetFont(titleFont);
    viewerSizer->Add(viewerTitle, 0, wxALL, 5);

    // Viewer controls row
    wxBoxSizer* viewerCtrlSizer = new wxBoxSizer(wxHORIZONTAL);
    wxButton* loadMolBtn = new wxButton(panel, ID_LOAD_MOLECULE, "Load PDB",
                                        wxDefaultPosition, wxSize(90, -1));
    viewerCtrlSizer->Add(loadMolBtn, 0, wxRIGHT, 5);

    viewerStatusLabel = new wxStaticText(panel, wxID_ANY, "No molecule loaded");
    viewerCtrlSizer->Add(viewerStatusLabel, 1, wxALIGN_CENTER_VERTICAL | wxLEFT, 5);
    viewerSizer->Add(viewerCtrlSizer, 0, wxEXPAND | wxLEFT | wxRIGHT, 5);

    // The GL canvas
    molViewer = new MoleculeViewer(panel, wxID_ANY, wxDefaultPosition, wxSize(400, 350));
    viewerSizer->Add(molViewer, 1, wxEXPAND | wxALL, 5);

    // Instructions below viewer
    wxStaticText* viewerHelp = new wxStaticText(panel, wxID_ANY,
        "Left-drag: rotate | Right-drag: pan | Scroll: zoom\n"
        "Click atom to select | Click again to deselect");
    viewerHelp->SetForegroundColour(wxColour(100, 100, 100));
    wxFont smallFont = viewerHelp->GetFont();
    smallFont.SetPointSize(smallFont.GetPointSize() - 1);
    viewerHelp->SetFont(smallFont);
    viewerSizer->Add(viewerHelp, 0, wxLEFT | wxRIGHT | wxBOTTOM, 5);

    mainSizer->Add(viewerSizer, 1, wxEXPAND | wxALL, 5);

    // ── Right side: atom selection controls ──────────────────────────────────
    wxBoxSizer* rightSizer = new wxBoxSizer(wxVERTICAL);
    rightSizer->AddSpacer(5);

    // ── Selection mode toggles ───────────────────────────────────────────────
    wxStaticBoxSizer* modeBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Selection Mode");
    wxBoxSizer* modeBtnSizer = new wxBoxSizer(wxHORIZONTAL);

    bondedModeBtn = new wxToggleButton(panel, ID_BONDED_MODE, "Select Bonded Atoms");
    deletedModeBtn = new wxToggleButton(panel, ID_DELETED_MODE, "Select Deleted Atoms");
    bondedModeBtn->SetValue(true);  // Start in bonded mode

    modeBtnSizer->Add(bondedModeBtn, 1, wxRIGHT, 3);
    modeBtnSizer->Add(deletedModeBtn, 1, wxLEFT, 3);
    modeBox->Add(modeBtnSizer, 0, wxEXPAND | wxALL, 5);

    wxStaticText* modeHelp = new wxStaticText(panel, wxID_ANY,
        "Toggle mode, then click atoms in the 3D viewer.\n"
        "You can also type atom names manually below.");
    modeHelp->SetForegroundColour(wxColour(100, 100, 100));
    modeHelp->SetFont(smallFont);
    modeBox->Add(modeHelp, 0, wxALL, 5);
    rightSizer->Add(modeBox, 0, wxEXPAND | wxALL, 5);

    // ── Bonded atoms section ─────────────────────────────────────────────────
    wxStaticBoxSizer* bondedBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Bonded Atoms");
    wxStaticText* bondedInstr = new wxStaticText(panel, wxID_ANY,
        "1st: directly bonded to amino acid\n"
        "2nd, 3rd: consecutively bonded (for vector)");
    bondedInstr->SetFont(smallFont);
    bondedInstr->SetForegroundColour(wxColour(100, 100, 100));
    bondedBox->Add(bondedInstr, 0, wxALL, 3);

    wxFlexGridSizer* bondedGrid = new wxFlexGridSizer(3, 3, 3, 5);
    bondedGrid->AddGrowableCol(1, 1);

    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "1st:"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom1Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    bondedGrid->Add(bondedAtom1Text, 1, wxEXPAND);
    // Color indicator
    wxPanel* bondedColor1 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    bondedColor1->SetBackgroundColour(wxColour(25, 255, 50));
    bondedGrid->Add(bondedColor1, 0, wxALIGN_CENTER_VERTICAL);

    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "2nd:"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom2Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    bondedGrid->Add(bondedAtom2Text, 1, wxEXPAND);
    wxPanel* bondedColor2 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    bondedColor2->SetBackgroundColour(wxColour(25, 255, 50));
    bondedGrid->Add(bondedColor2, 0, wxALIGN_CENTER_VERTICAL);

    bondedGrid->Add(new wxStaticText(panel, wxID_ANY, "3rd:"), 0, wxALIGN_CENTER_VERTICAL);
    bondedAtom3Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    bondedGrid->Add(bondedAtom3Text, 1, wxEXPAND);
    wxPanel* bondedColor3 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    bondedColor3->SetBackgroundColour(wxColour(25, 255, 50));
    bondedGrid->Add(bondedColor3, 0, wxALIGN_CENTER_VERTICAL);

    bondedBox->Add(bondedGrid, 0, wxALL | wxEXPAND, 5);

    wxButton* clearBondedBtn = new wxButton(panel, ID_CLEAR_BONDED, "Clear Bonded",
                                            wxDefaultPosition, wxSize(-1, -1));
    bondedBox->Add(clearBondedBtn, 0, wxALL | wxALIGN_RIGHT, 3);

    rightSizer->Add(bondedBox, 0, wxEXPAND | wxALL, 5);

    // ── Deleted atoms section ────────────────────────────────────────────────
    wxStaticBoxSizer* deletedBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Deleted Atoms");
    wxStaticText* deletedInstr = new wxStaticText(panel, wxID_ANY,
        "Select atoms to be removed (up to 3):");
    deletedInstr->SetFont(smallFont);
    deletedInstr->SetForegroundColour(wxColour(100, 100, 100));
    deletedBox->Add(deletedInstr, 0, wxALL, 3);

    wxFlexGridSizer* deletedGrid = new wxFlexGridSizer(3, 3, 3, 5);
    deletedGrid->AddGrowableCol(1, 1);

    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "1st:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom1Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    deletedGrid->Add(deletedAtom1Text, 1, wxEXPAND);
    wxPanel* deletedColor1 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    deletedColor1->SetBackgroundColour(wxColour(255, 50, 25));
    deletedGrid->Add(deletedColor1, 0, wxALIGN_CENTER_VERTICAL);

    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "2nd:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom2Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    deletedGrid->Add(deletedAtom2Text, 1, wxEXPAND);
    wxPanel* deletedColor2 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    deletedColor2->SetBackgroundColour(wxColour(255, 50, 25));
    deletedGrid->Add(deletedColor2, 0, wxALIGN_CENTER_VERTICAL);

    deletedGrid->Add(new wxStaticText(panel, wxID_ANY, "3rd:"), 0, wxALIGN_CENTER_VERTICAL);
    deletedAtom3Text = new wxTextCtrl(panel, wxID_ANY, "", wxDefaultPosition, wxSize(100, -1));
    deletedGrid->Add(deletedAtom3Text, 1, wxEXPAND);
    wxPanel* deletedColor3 = new wxPanel(panel, wxID_ANY, wxDefaultPosition, wxSize(14, 14));
    deletedColor3->SetBackgroundColour(wxColour(255, 50, 25));
    deletedGrid->Add(deletedColor3, 0, wxALIGN_CENTER_VERTICAL);

    deletedBox->Add(deletedGrid, 0, wxALL | wxEXPAND, 5);

    wxButton* clearDeletedBtn = new wxButton(panel, ID_CLEAR_DELETED, "Clear Deleted",
                                             wxDefaultPosition, wxSize(-1, -1));
    deletedBox->Add(clearDeletedBtn, 0, wxALL | wxALIGN_RIGHT, 3);

    rightSizer->Add(deletedBox, 0, wxEXPAND | wxALL, 5);

    // ── Action buttons ───────────────────────────────────────────────────────
    wxButton* confirmAtomsButton = new wxButton(panel, ID_CONFIRM_ATOMS, "Confirm Atom Selection");
    wxButton* processFileButton = new wxButton(panel, ID_PROCESS_ORIENT, "Process and Orient Drug");
    rightSizer->Add(confirmAtomsButton, 0, wxEXPAND | wxALL, 5);
    rightSizer->Add(processFileButton, 0, wxEXPAND | wxALL, 5);

    rightSizer->AddStretchSpacer();
    mainSizer->Add(rightSizer, 0, wxEXPAND | wxALL, 5);

    panel->SetSizer(mainSizer);

    // ── Set up selection callback ────────────────────────────────────────────
    molViewer->SetSelectionCallback(
        [this](SelectionMode mode, const std::vector<std::string>& names) {
            OnSelectionChanged(mode, names);
        }
    );

    return panel;
}

wxPanel* MyFrame::CreateTopologyPanel(wxNotebook* notebook)
{
    wxPanel* panel = new wxPanel(notebook);
    wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
    
    sizer->AddSpacer(10);
    
    wxStaticBoxSizer* patchBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Patch Generation");
    wxButton* processPatchMakerButton = new wxButton(panel, ID_PATCH_MAKER, "Process Patch Maker");
    patchBox->Add(processPatchMakerButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(patchBox, 0, wxALL | wxEXPAND, 10);
    
    wxStaticBoxSizer* naaMol2StrBox = new wxStaticBoxSizer(wxVERTICAL, panel, "NAA CGenFF Workflow");
    wxStaticText* naaInstructionText = new wxStaticText(panel, wxID_ANY, 
        "1. Download naa.mol2 and upload it to https://cgenff.silcsbio.com/\n"
        "   *** Tick 'include parameters that are already in CGenFF' ***\n"
        "   Without it the .str holds only the new, analogy-derived terms, and\n"
        "   every other term must come from your local par_all36_cgenff.prm -\n"
        "   so a newer service library than your local release breaks the run.\n"
        "2. Import the returned zip here (naa.str is extracted)\n"
        "LYSINE ONLY: a lysine modification also modifies the N-terminus.\n"
        "   Repeat steps 1-2 for ntrm.mol2 -> ntrm.str. Both are required.");
    naaInstructionText->Wrap(400);
    naaMol2StrBox->Add(naaInstructionText, 0, wxALL, 5);
    
    wxButton* downloadNAAMol2Button = new wxButton(panel, ID_DOWNLOAD_NAA_MOL2, "Download NAA.mol2");
    wxButton* uploadNAAStrButton = new wxButton(panel, ID_UPLOAD_NAA_STR, "Import CGenFF result (zip/str)");
    naaMol2StrBox->Add(downloadNAAMol2Button, 0, wxALL | wxEXPAND, 5);
    naaMol2StrBox->Add(uploadNAAStrButton, 0, wxALL | wxEXPAND, 5);
    wxButton* downloadNtrmMol2Button = new wxButton(panel, ID_DOWNLOAD_NTRM_MOL2, "Download ntrm.mol2 (Lysine runs)");
    wxButton* uploadNtrmStrButton = new wxButton(panel, ID_UPLOAD_NTRM_STR, "Import N-terminus CGenFF result (zip/str)");
    naaMol2StrBox->Add(downloadNtrmMol2Button, 0, wxALL | wxEXPAND, 5);
    naaMol2StrBox->Add(uploadNtrmStrButton, 0, wxALL | wxEXPAND, 5);

    // Optional step, so it sits to the side rather than in the main column:
    // only molecules with chemistry CGenFF cannot type (a persistent radical is
    // the usual case) need hand-curated charges. A run without an override file
    // is unaffected.
    wxBoxSizer* overrideRow = new wxBoxSizer(wxHORIZONTAL);
    chargeOverrideLabel = new wxStaticText(panel, wxID_ANY,
        "Optional: no charge override file selected.\n"
        "Pick one only if CGenFF cannot type part of your\n"
        "molecule (a radical, say). Applied automatically\n"
        "during 'Create Novel Amino Acid'.");
    chargeOverrideLabel->Wrap(260);
    wxButton* applyOverridesButton = new wxButton(panel, ID_APPLY_CHARGE_OVERRIDES,
        "Charge\noverride file...", wxDefaultPosition, wxSize(120, 40));
    overrideRow->Add(chargeOverrideLabel, 1, wxALIGN_CENTER_VERTICAL | wxLEFT | wxRIGHT, 5);
    overrideRow->Add(applyOverridesButton, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);
    naaMol2StrBox->Add(overrideRow, 0, wxEXPAND | wxTOP, 5);

    sizer->Add(naaMol2StrBox, 0, wxALL | wxEXPAND, 10);
    
    wxStaticBoxSizer* topologyBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Novel Amino Acid Creation");
    wxButton* createTopologyButton = new wxButton(panel, ID_CREATE_TOPOLOGY, "Create Novel Amino Acid");
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
    
    wxStaticBoxSizer* vlpBox = new wxStaticBoxSizer(wxVERTICAL, panel, "VLP Type Selection");
    wxString vlpTypes[] = {"Qbeta", "TMV"};
    vlpTypeChoice = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxDefaultSize, 2, vlpTypes);
    vlpTypeChoice->SetSelection(0);
    vlpBox->Add(vlpTypeChoice, 0, wxALL | wxEXPAND, 5);
    sizer->Add(vlpBox, 0, wxALL | wxEXPAND, 10);
    
    wxStaticBoxSizer* capsidBox = new wxStaticBoxSizer(wxVERTICAL, panel, "Capsid Operations");
    wxButton* buildMonomerButton = new wxButton(panel, ID_BUILD_MONOMER, "1. Build Modified Monomer(s)");
    wxButton* minimizeMonomerButton = new wxButton(panel, ID_MINIMIZE_MONOMER, "2. Minimize Monomer (protein fixed)");
    wxButton* createViewCapsidButton = new wxButton(panel, ID_CREATE_VIEW_CAPSID, "3. Assemble And View Capsid");
    wxButton* minimizeCapsidButton = new wxButton(panel, ID_MINIMIZE_CAPSID, "Minimize Capsid");
    wxButton* zipFilesButton = new wxButton(panel, ID_ZIP_FILES, "Package Files (Zip)");
    capsidBox->Add(buildMonomerButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(minimizeMonomerButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(createViewCapsidButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(minimizeCapsidButton, 0, wxALL | wxEXPAND, 5);
    capsidBox->Add(zipFilesButton, 0, wxALL | wxEXPAND, 5);
    sizer->Add(capsidBox, 0, wxALL | wxEXPAND, 10);
    
    sizer->AddStretchSpacer();
    panel->SetSizer(sizer);
    return panel;
}

// ── Viewer-related event handlers ────────────────────────────────────────────

void MyFrame::OnLoadMolecule(wxCommandEvent& event)
{
    // Try to auto-load one.pdb if it exists, otherwise ask user
    wxString pdbPath = "one.pdb";
    
    if (!wxFileExists(pdbPath)) {
        wxFileDialog openFileDialog(this, "Select PDB file to view", "", "",
                                    "PDB files (*.pdb)|*.pdb", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
        if (openFileDialog.ShowModal() == wxID_CANCEL)
            return;
        pdbPath = openFileDialog.GetPath();
    }

    if (molViewer->LoadPDB(pdbPath.ToStdString())) {
        statusText->AppendText("Loaded molecule from: " + pdbPath + "\n");
        UpdateViewerStatus();
    } else {
        statusText->AppendText("Error: Failed to load molecule from: " + pdbPath + "\n");
        wxMessageBox("Failed to load PDB file.", "Error", wxOK | wxICON_ERROR);
    }
}

void MyFrame::OnBondedModeToggle(wxCommandEvent& event)
{
    bondedModeBtn->SetValue(true);
    deletedModeBtn->SetValue(false);
    molViewer->SetSelectionMode(SelectionMode::BONDED);
    statusText->AppendText("Selection mode: Bonded Atoms (green)\n");
}

void MyFrame::OnDeletedModeToggle(wxCommandEvent& event)
{
    deletedModeBtn->SetValue(true);
    bondedModeBtn->SetValue(false);
    molViewer->SetSelectionMode(SelectionMode::DELETED);
    statusText->AppendText("Selection mode: Deleted Atoms (red)\n");
}

void MyFrame::OnClearBondedSelection(wxCommandEvent& event)
{
    molViewer->ClearBondedSelections();
    bondedAtom1Text->SetValue("");
    bondedAtom2Text->SetValue("");
    bondedAtom3Text->SetValue("");
    statusText->AppendText("Cleared bonded atom selections.\n");
}

void MyFrame::OnClearDeletedSelection(wxCommandEvent& event)
{
    molViewer->ClearDeletedSelections();
    deletedAtom1Text->SetValue("");
    deletedAtom2Text->SetValue("");
    deletedAtom3Text->SetValue("");
    statusText->AppendText("Cleared deleted atom selections.\n");
}

void MyFrame::OnSelectionChanged(SelectionMode mode, const std::vector<std::string>& names)
{
    std::cout << "OnSelectionChanged called with " << names.size() << " names, mode="
              << (mode == SelectionMode::BONDED ? "BONDED" : "DELETED") << std::endl;

    if (mode == SelectionMode::BONDED) {
        bondedAtom1Text->SetValue(names.size() > 0 ? wxString(names[0]) : "");
        bondedAtom2Text->SetValue(names.size() > 1 ? wxString(names[1]) : "");
        bondedAtom3Text->SetValue(names.size() > 2 ? wxString(names[2]) : "");

        statusText->AppendText("Bonded atoms updated: ");
        for (const auto& n : names) statusText->AppendText(wxString(n) + " ");
        statusText->AppendText("\n");
    } else {
        deletedAtom1Text->SetValue(names.size() > 0 ? wxString(names[0]) : "");
        deletedAtom2Text->SetValue(names.size() > 1 ? wxString(names[1]) : "");
        deletedAtom3Text->SetValue(names.size() > 2 ? wxString(names[2]) : "");

        statusText->AppendText("Deleted atoms updated: ");
        for (const auto& n : names) statusText->AppendText(wxString(n) + " ");
        statusText->AppendText("\n");
    }

    std::cout << "OnSelectionChanged completed" << std::endl;
}

void MyFrame::UpdateViewerStatus()
{
    if (molViewer->HasMolecule()) {
        wxString status = wxString::Format("Loaded: %d atoms, %d bonds",
                                           molViewer->GetAtomCount(),
                                           molViewer->GetBondCount());
        viewerStatusLabel->SetLabel(status);
    } else {
        viewerStatusLabel->SetLabel("No molecule loaded");
    }
}

// ── Helper functions (unchanged) ─────────────────────────────────────────────

wxString MyFrame::GetAminoAcidAbbreviation(const wxString& fullName)
{
    if (fullName == "Lysine") return "lys";
    if (fullName == "Cysteine") return "cys";
    if (fullName == "Tyrosine") return "tyr";
    if (fullName == "N-Terminus") return "ntrm";
    return fullName.Lower();
}

void MyFrame::SelectFile(wxString& target, const wxString& message, const wxString& wildcard)
{
    wxFileDialog openFileDialog(this, message, "", "", wildcard, wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (openFileDialog.ShowModal() == wxID_CANCEL) return;
    target = openFileDialog.GetPath();
    statusText->AppendText("Selected file: " + target + "\n");
}

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

// ── All existing event handlers below (unchanged) ────────────────────────────

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
        statusText->AppendText("1. Choose where to save one.mol2 (dialog follows)\n");
        statusText->AppendText("2. Upload it to the CGenFF website: https://cgenff.silcsbio.com/\n");
        statusText->AppendText("3. Download the result zip from CGenFF (no need to unzip it)\n");
        statusText->AppendText("4. Click 'Import CGenFF result' and pick that zip; one.str is extracted and the zip is trashed\n\n");
        
        // Let the user save a copy where the CGenFF upload can find it
        wxFileDialog saveDialog(this, "Save one.mol2 for CGenFF upload", "", "one.mol2",
                                "MOL2 files (*.mol2)|*.mol2", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
        if (saveDialog.ShowModal() == wxID_OK)
        {
            wxString dest = saveDialog.GetPath();
            if (wxCopyFile("one.mol2", dest, true))
                statusText->AppendText("one.mol2 saved to: " + dest + "\n\n");
            else
                statusText->AppendText("Error: could not write " + dest + "\n");
        }
        else
        {
            statusText->AppendText("Save cancelled; one.mol2 remains in the working directory: " + wxGetCwd() + "\n\n");
        }
    }
    else
    {
        statusText->AppendText("Error: one.mol2 generation failed.\n");
        wxMessageBox("Failed to generate one.mol2 file. Make sure obabel is installed.", "Error", wxOK | wxICON_ERROR);
    }
}


// ── CGenFF result import ──────────────────────────────────────────────────────
// Move a file to the Trash (recoverable) instead of deleting it outright.
static bool MoveToTrash(const wxString& path)
{
    if (!wxFileExists(path)) return false;
#ifdef __WXMAC__
    wxString script = "tell application \"Finder\" to delete POSIX file \"" + path + "\"";
    wxArrayString out, err;
    return wxExecute("osascript -e '" + script + "'", out, err, wxEXEC_SYNC) == 0;
#else
    return wxRemoveFile(path);
#endif
}

// Accept either a bare .str or the .zip returned by cgenff.silcsbio.com.
// Writes <targetName> in the working directory. On success the zip (and a
// copy of <mol2Name> sitting next to it, if any) is moved to the Trash.
static bool ImportCGenFFResult(wxWindow* parent, wxTextCtrl* log,
                               const wxString& targetName, const wxString& mol2Name)
{
    wxString downloads = wxStandardPaths::Get().GetUserDir(wxStandardPaths::Dir_Downloads);
    wxFileDialog dlg(parent, "Select CGenFF result (" + targetName + " or the downloaded .zip)",
                     downloads, "", "CGenFF result (*.zip;*.str)|*.zip;*.str|All files (*.*)|*.*",
                     wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (dlg.ShowModal() == wxID_CANCEL) return false;

    wxString src = dlg.GetPath();
    wxString target = wxGetCwd() + "/" + targetName;
    wxString staging = target + ".incoming";
    wxFileName srcFn(src);

    // Verify the stream file really describes mol2Name before accepting it.
    // The naa and ntrm zips come back from CGenFF with near-identical names.
    auto verifyAndInstall = [&](void) -> bool {
        wxString other = (mol2Name == "naa.mol2") ? "ntrm.mol2" : "naa.mol2";
        wxString cmd = "python3 check_str_matches_mol2.py \"" + staging + "\" " + mol2Name + " " + other;
        wxArrayString out, err;
        long rc = wxExecute(cmd, out, err, wxEXEC_SYNC);
        for (const auto& l : out) log->AppendText("  " + l + "\n");
        if (rc != 0) {
            wxRemoveFile(staging);
            log->AppendText("Import REJECTED: this is not the CGenFF result for " + mol2Name + ".\n\n");
            wxMessageBox("That stream file does not match " + mol2Name + ".\n\nSee the status log - it probably belongs "
                         "to the other molecule. Nothing was changed.", "Wrong CGenFF result", wxOK | wxICON_WARNING);
            return false;
        }
        if (wxFileExists(target)) wxRemoveFile(target);
        if (!wxRenameFile(staging, target)) { log->AppendText("Error: could not write " + target + "\n"); return false; }
        return true;
    };

    if (srcFn.GetExt().Lower() != "zip")
    {
        if (!wxCopyFile(src, staging, true)) { log->AppendText("Error: could not copy " + src + "\n"); return false; }
        if (!verifyAndInstall()) return false;
        log->AppendText(targetName + " copied from " + src + "\n\n");
        return true;
    }

    // --- zip: find the first *.str entry ---
    wxFileInputStream fin(src);
    if (!fin.IsOk()) { log->AppendText("Error: cannot open " + src + "\n"); return false; }
    wxZipInputStream zip(fin);
    bool found = false;
    wxString entryName;
    for (wxZipEntry* e = zip.GetNextEntry(); e; e = zip.GetNextEntry())
    {
        wxString name = e->GetName();
        bool isStr = name.Lower().EndsWith(".str") && !name.StartsWith("__MACOSX");
        if (isStr && zip.OpenEntry(*e))
        {
            wxFileOutputStream fout(staging);
            if (!fout.IsOk()) { delete e; break; }
            zip.Read(fout);
            found = true; entryName = name;
            delete e;
            break;
        }
        delete e;
    }
    if (!found) { log->AppendText("Error: no .str file inside " + src + "\n"); return false; }
    log->AppendText("Extracted " + entryName + " from " + srcFn.GetFullName() + "\n");
    if (!verifyAndInstall()) return false;
    log->AppendText("Installed as " + target + "\n");

    // --- tidy up: trash the zip and the mol2 copy that was uploaded, if next to it ---
    if (MoveToTrash(src)) log->AppendText("Moved " + srcFn.GetFullName() + " to Trash\n");
    wxString mol2Copy = srcFn.GetPath() + "/" + mol2Name;
    if (wxFileExists(mol2Copy) && MoveToTrash(mol2Copy))
        log->AppendText("Moved " + mol2Name + " (upload copy) to Trash\n");
    log->AppendText("You can now proceed to the next step.\n\n");
    return true;
}

void MyFrame::OnUploadOneStr(wxCommandEvent& event)
{
    if (!ImportCGenFFResult(this, statusText, "one.str", "one.mol2"))
        wxMessageBox("No one.str was imported.", "CGenFF import", wxOK | wxICON_INFORMATION);
}


// obabel -> <stem>_raw.mol2, then reconcile atoms/bonds against the PSF
// (psfgen's topology is the truth; obabel's distance-based perception can
// bond the junction N to a nearby O instead of the chosen carbon).
bool MyFrame::MakeMol2(const wxString& pdb, const wxString& psf, const wxString& mol2)
{
    if (!wxFileExists(pdb) || !wxFileExists(psf)) {
        statusText->AppendText("Error: need " + pdb + " and " + psf + " to build " + mol2 + "\n");
        return false;
    }
    wxString raw = mol2.BeforeLast('.') + "_raw.mol2";
    wxString cmd1 = "obabel " + pdb + " -O " + raw;
    statusText->AppendText("Running command: " + cmd1 + "\n");
    if (system(cmd1.mb_str()) != 0 || !wxFileExists(raw)) {
        statusText->AppendText("Error: obabel failed for " + pdb + "\n");
        return false;
    }
    wxString cmd2 = "python3 reconcile_mol2_with_psf.py " + raw + " " + psf + " " + mol2;
    statusText->AppendText("Running command: " + cmd2 + "\n");
    if (system(cmd2.mb_str()) != 0 || !wxFileExists(mol2)) {
        statusText->AppendText("Error: reconcile_mol2_with_psf.py failed for " + mol2 + "\n");
        return false;
    }
    statusText->AppendText(mol2 + " written (bonds taken from " + psf + ").\n");
    return true;
}

void MyFrame::SaveMol2Copy(const wxString& mol2)
{
    wxFileDialog saveDialog(this, "Save " + mol2 + " for CGenFF upload", "", mol2,
                            "MOL2 files (*.mol2)|*.mol2", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
    if (saveDialog.ShowModal() == wxID_OK) {
        wxString dest = saveDialog.GetPath();
        if (wxCopyFile(mol2, dest, true)) statusText->AppendText(mol2 + " saved to: " + dest + "\n\n");
        else statusText->AppendText("Error: could not write " + dest + "\n");
    } else {
        statusText->AppendText("Save cancelled; " + mol2 + " remains in " + wxGetCwd() + "\n\n");
    }
}

// Hand-curated atom types and partial charges, for groups CGenFF has no
// chemistry for. CGenFF types the nearest closed-shell analogue instead - a
// persistent radical's oxygen comes back as an alkoxide, which also saddles the
// residue with a formal charge it should not have. The values must come from the
// literature, so the GUI supplies only the mechanism; the numbers live in an
// override .prm the user selects here.
//
// This button only RECORDS the file. The override cannot be applied on demand,
// because nad.str is created and consumed inside a single OnCreateTopology run:
// noX_psf.py writes it, the rename retypes it, and ./IC_naa then reads it to
// build the RESI NAA block in top_all36_cgenff_CBD.rtf - which is where psfgen
// gets the charges that end up in the PSF. Applying it before that run would be
// overwritten; after, the topology is already written. So OnCreateTopology
// applies it in the one correct position: after the rename (so the protein-side
// types are settled) and before the parameter steps (so they see the corrected
// atom types). With no file selected, nothing changes and the run behaves as
// it always did.
void MyFrame::OnSelectChargeOverride(wxCommandEvent& event)
{
    wxFileDialog openDialog(this, "Select the charge override file", "", "",
                            "Parameter files (*.prm)|*.prm|All files (*.*)|*.*",
                            wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if (openDialog.ShowModal() == wxID_CANCEL) return;

    chargeOverrideFile = openDialog.GetPath();
    if (chargeOverrideLabel) {
        chargeOverrideLabel->SetLabel("Charge overrides WILL be applied from:\n" +
                                      wxFileName(chargeOverrideFile).GetFullName() +
                                      "\nduring 'Create Novel Amino Acid'.");
        chargeOverrideLabel->Wrap(260);
        chargeOverrideLabel->GetParent()->Layout();
    }
    statusText->AppendText("Charge override file selected: " + chargeOverrideFile +
                           "\nIt will be applied during 'Create Novel Amino Acid', "
                           "after the atom-type renaming and before the parameter steps.\n\n");
}

void MyFrame::OnDownloadNtrmMol2(wxCommandEvent& event)
{
    if (!wxFileExists("ntrm_fixed.pdb") || !wxFileExists("ntrm.psf")) {
        wxMessageBox("No N-terminus structure yet. Run the Patch Maker step on a Lysine job first.",
                     "ntrm.mol2", wxOK | wxICON_INFORMATION);
        return;
    }
    if (MakeMol2("ntrm_fixed.pdb", "ntrm.psf", "ntrm.mol2")) SaveMol2Copy("ntrm.mol2");
}

void MyFrame::OnUploadNtrmStr(wxCommandEvent& event)
{
    if (!ImportCGenFFResult(this, statusText, "ntrm.str", "ntrm.mol2"))
        wxMessageBox("No ntrm.str was imported.", "CGenFF import", wxOK | wxICON_INFORMATION);
}

void MyFrame::OnDownloadNAAMol2(wxCommandEvent& event)
{
    if (!wxFileExists("naa.pdb") || !wxFileExists("patched.psf")) {
        wxMessageBox("naa.pdb / patched.psf not found. Run the Patch Maker step first.", "naa.mol2", wxOK | wxICON_ERROR);
        return;
    }
    if (MakeMol2("naa.pdb", "patched.psf", "naa.mol2")) SaveMol2Copy("naa.mol2");
}

void MyFrame::OnUploadNAAStr(wxCommandEvent& event)
{
    if (!ImportCGenFFResult(this, statusText, "naa.str", "naa.mol2"))
        wxMessageBox("No naa.str was imported.", "CGenFF import", wxOK | wxICON_INFORMATION);
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

    // par_GUI.prm accumulates cross-FF blocks from every job; start it fresh too
    {
        std::ifstream prmTemplate("par_GUI_template.prm");
        if (!prmTemplate.is_open()) {
            statusText->AppendText("Error: Failed to open template file par_GUI_template.prm\n");
            return;
        }
        std::ofstream prmTarget("par_GUI.prm");
        prmTarget << prmTemplate.rdbuf();
        statusText->AppendText("Successfully reset par_GUI.prm from template.\n");
    }
    
    wxString command = "python3 fix_atom_naming_pdb.py";
    statusText->AppendText("Running command: " + command + "\n");
    
    int result = system(command.mb_str());
    
    if (result == 0) {
        statusText->AppendText("Atom renaming script executed successfully.\n");
    } else {
        statusText->AppendText("Error: Atom renaming script failed. Make sure fix_atom_naming_pdb.py exists and python3 is installed.\n");
        return;
    }
    wxString bbsCommand = "python3 fix_atom_naming_pdb_BBS.py";
    statusText->AppendText("Running command: " + bbsCommand + "\n");
    
    int bbsResult = system(bbsCommand.mb_str());
    
    if (bbsResult == 0) {
        statusText->AppendText("BBS atom renaming script executed successfully.\n");
    } else {
        statusText->AppendText("Error: BBS atom renaming script failed. Make sure fix_atom_naming_pdb_BBS.py exists.\n");
        return;
    }
    statusText->AppendText("Running ./IC_table_one command...\n");
    
    wxString icTableCommand = "./IC_table_one";   // built from IC_table_one.cpp
    int icResult = system(icTableCommand.mb_str());
    
    if (icResult == 0) {
        statusText->AppendText("IC_table command executed successfully.\n");
        statusText->AppendText("Atom renaming process completed successfully!\n");

        // Auto-load the renamed molecule into the viewer
        if (wxFileExists("one.pdb") && molViewer) {
            if (molViewer->LoadPDB("one.pdb")) {
                statusText->AppendText("Molecule automatically loaded into 3D viewer.\n");
                UpdateViewerStatus();
            }
        }
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


// The N-terminus model residue depends on the VLP: Qbeta residue 1 is Ala
// (RESI BBQ, BBS_QB.pdb), TMV residue 1 is Ser (RESI BBT, BBS_TMV.pdb). The
// scripts always read ntrm_ref.pdb, so only this staging step differs.
wxString MyFrame::NtrmVariant()
{
    return (vlpTypeChoice->GetSelection() == 1) ? "ntrm_ser" : "ntrm";
}

bool MyFrame::StageNtrmReference()
{
    bool tmv = (vlpTypeChoice->GetSelection() == 1);
    wxString src = tmv ? "BBS_TMV.pdb" : "BBS_QB.pdb";
    wxString resi = tmv ? "BBT (Ser)" : "BBQ (Ala)";
    if (!wxFileExists(src)) {
        statusText->AppendText("Error: N-terminus reference " + src + " not found.\n");
        return false;
    }
    if (tmv) {
        // RESI BBT must be installed in the topology file or psfgen cannot build it
        wxTextFile rtf;
        bool haveBBT = false;
        if (rtf.Open("top_all36_cgenff_CBD.rtf")) {
            for (wxString l = rtf.GetFirstLine(); !rtf.Eof(); l = rtf.GetNextLine())
                if (l.StartsWith("RESI BBT")) { haveBBT = true; break; }
            rtf.Close();
        }
        if (!haveBBT) {
            statusText->AppendText("Error: RESI BBT (TMV Ser N-terminus) is not in top_all36_cgenff_CBD.rtf.\n"
                                   "One-time setup: send BBS_TMV.mol2 to CGenFF, then run\n"
                                   "  python3 install_nterm_reference.py BBT.str\n"
                                   "and click 'Run Renaming' to refresh the working topology file.\n");
            wxMessageBox("The TMV serine N-terminus model (RESI BBT) is not installed yet.\n\n"
                         "Send BBS_TMV.mol2 to CGenFF, then run:\n"
                         "  python3 install_nterm_reference.py BBT.str\n\n"
                         "See the status log.", "N-terminus model missing", wxOK | wxICON_ERROR);
            return false;
        }
    }
    if (!wxCopyFile(src, "ntrm_ref.pdb", true)) {
        statusText->AppendText("Error: could not stage " + src + " as ntrm_ref.pdb\n");
        return false;
    }
    statusText->AppendText("N-terminus model: " + src + " -> ntrm_ref.pdb, RESI " + resi + "\n");
    statusText->AppendText("  Assumes a FREE (unacetylated) N-terminal alpha-amine, as expected for\n"
                           "  recombinant expression in E. coli, which lacks the eukaryotic NatA/NatB/NatC\n"
                           "  N-terminal acetyltransferases. A capsid from a native plant/eukaryotic source\n"
                           "  may be N-terminally acetylated and therefore unreactive to amine chemistry -\n"
                           "  in that case the N-terminus should not be modified at all.\n");
    return true;
}

// VMD returns 0 even when psfgen aborts inside the script, so success is
// judged by whether the file the script was meant to write actually appeared.
int MyFrame::RunVmd(const wxString& command, const wxString& expectedOutput)
{
    if (!expectedOutput.IsEmpty() && wxFileExists(expectedOutput)) wxRemoveFile(expectedOutput);
    int rc = system(command.mb_str());
    if (rc != 0) return rc;
    if (!expectedOutput.IsEmpty() && !wxFileExists(expectedOutput)) {
        statusText->AppendText("Error: VMD finished but did not write " + expectedOutput +
                               " - psfgen failed inside the script (look for 'add atom failed' / 'unknown atom type' in the terminal output).\n");
        return 1;
    }
    return 0;
}

void MyFrame::OnProcessFile(wxCommandEvent& event)
{
    RunOrientStep();
}

// Fill the bonded-atom placeholders in orient_place_drug_<aa>.py, run it
// (and the N-terminus variant for Lysine), then reset the script from its
// template. Returns true only if one_aligned.pdb was produced.
bool MyFrame::RunOrientStep()
{
    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);
    wxString firstBondedAtom = bondedAtom1Text->GetValue().Trim();
    wxString thirdBondedAtom = bondedAtom3Text->GetValue().Trim();

    if (firstBondedAtom.IsEmpty()) {
        statusText->AppendText("Error: First bonded atom must be specified for orientation.\n");
        return false;
    }
    if (!wxFileExists("one.pdb")) {
        statusText->AppendText("Error: one.pdb not found. Load the incoming molecule first.\n");
        return false;
    }

    statusText->AppendText("=== Orienting " + selectedAminoAcid + " / drug (*T1=" + firstBondedAtom +
                           (thirdBondedAtom.IsEmpty() ? wxString("") : ", *T3=" + thirdBondedAtom) + ") ===\n");

    // Always start from the template so a previous run's names cannot linger
    processor.ResetAminoAcidFile(std::string(selectedAminoAcid.mb_str()));
    bool success = processor.ProcessAminoAcidFile(std::string(selectedAminoAcid.mb_str()),
                                                  std::string(firstBondedAtom.mb_str()),
                                                  std::string(thirdBondedAtom.mb_str()));
    if (!success) {
        statusText->AppendText("Error: could not fill the orient script placeholders.\n");
        return false;
    }

    wxString aminoAcidAbbrev = GetAminoAcidAbbreviation(selectedAminoAcid);
    wxString pythonCommand = "python3 orient_place_drug_" + aminoAcidAbbrev + ".py";
    statusText->AppendText("Running command: " + pythonCommand + "\n");
    int result = system(pythonCommand.mb_str());
    bool ok = (result == 0) && wxFileExists("one_aligned.pdb");

    if (ok && selectedAminoAcid == "Lysine") {
        // N-terminus copy: same placeholders, separate script
        std::string ntrmScript = "orient_place_drug_ntrm.py";
        processor.ResetAminoAcidFile("N-Terminus");
        processor.ProcessAminoAcidFile("N-Terminus",
                                       std::string(firstBondedAtom.mb_str()),
                                       std::string(thirdBondedAtom.mb_str()));
        if (!StageNtrmReference()) return false;
        wxString ntrmCommand = "python3 orient_place_drug_ntrm.py";
        statusText->AppendText("Running command: " + ntrmCommand + "\n");
        if (system(ntrmCommand.mb_str()) != 0) {
            statusText->AppendText("Error: orient_place_drug_ntrm.py failed.\n");
            ok = false;
        }
        processor.ResetAminoAcidFile("N-Terminus");
    }

    processor.ResetAminoAcidFile(std::string(selectedAminoAcid.mb_str()));
    if (ok)
        statusText->AppendText("Orientation complete: one_aligned.pdb written.\n\n");
    else
        statusText->AppendText("Error: orient script failed (see output above); one_aligned.pdb not updated.\n\n");
    return ok;
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

    // Re-runnable: drop patches/residues a previous Patch Maker run appended
    {
        wxString strip = "python3 strip_rtf_residues.py ZAP ZPP NAA NPC ntrm_fix naa ntrm";
        statusText->AppendText("Running command: " + strip + "\n");
        system(strip.mb_str());
    }

    // one_aligned.pdb must come from the current one.pdb; (re)orient if not
    {
        wxFileName alignedFn("one_aligned.pdb"), oneFn("one.pdb");
        bool stale = !alignedFn.FileExists() ||
                     (oneFn.FileExists() && oneFn.GetModificationTime() > alignedFn.GetModificationTime());
        if (stale) {
            statusText->AppendText("one_aligned.pdb is missing or older than one.pdb - running orientation first.\n");
            if (!RunOrientStep()) return;
        }
    }
    
    std::vector<std::string> bondedAtoms;
    std::vector<std::string> deletedAtoms;
    
    if (!bonded1.IsEmpty()) bondedAtoms.push_back(std::string(bonded1.mb_str()));
    if (!bonded2.IsEmpty()) bondedAtoms.push_back(std::string(bonded2.mb_str()));
    if (!bonded3.IsEmpty()) bondedAtoms.push_back(std::string(bonded3.mb_str()));
    
    if (!deleted1.IsEmpty()) deletedAtoms.push_back(std::string(deleted1.mb_str()));
    if (!deleted2.IsEmpty()) deletedAtoms.push_back(std::string(deleted2.mb_str()));
    if (!deleted3.IsEmpty()) deletedAtoms.push_back(std::string(deleted3.mb_str()));
    
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
            
            int vmdResult = RunVmd(vmdCommand, "patched.pdb");
            
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
            if (!MakeMol2("naa.pdb", "patched.psf", "naa.mol2")) {
                statusText->AppendText("Error: naa.mol2 not produced.\n");
            }
// --- N-terminus parallel processing (for Lysine only) ---
            if (selectedAminoAcid == "Lysine") {
                if (!StageNtrmReference()) return;
                // Process ntrm patch maker template with the same bonded/deleted atoms
                bool ntrmSuccess = processor.ProcessPatchMakerFile("N-Terminus",
                                                                    bondedAtoms,
                                                                    deletedAtoms);
                if (!ntrmSuccess) {
                    statusText->AppendText("Error: Failed to process N-terminus patch maker file.\n");
                    return;
                }
                wxString ntrmPatchCommand = "python3 patch_maker_ntrm.py";
                statusText->AppendText("Running command: " + ntrmPatchCommand + "\n");
                int ntrmPatchResult = system(ntrmPatchCommand.mb_str());
                
                if (ntrmPatchResult != 0) {
                    statusText->AppendText("Error: patch_maker_ntrm.py failed.\n");
                    return;
                }
                statusText->AppendText("N-terminus patch maker executed successfully.\n");
                
                // Reset the ntrm patch maker template
                processor.ResetPatchMakerFile("N-Terminus");

                wxString ntrmVmdCommand = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e patched_ntrm.pgn";
                statusText->AppendText("Running N-terminus VMD command: " + ntrmVmdCommand + "\n");
                
                int ntrmVmdResult = RunVmd(ntrmVmdCommand, "ntrm.pdb");
                
                if (ntrmVmdResult == 0) {
                    statusText->AppendText("N-terminus VMD command executed successfully.\n");
                    
                    wxString ntrmElementFix = "python3 element_column_fix_NPC.py ntrm.pdb ntrm_fixed.pdb";
                    statusText->AppendText("Running command: " + ntrmElementFix + "\n");
                    int ntrmFixResult = wxExecute(ntrmElementFix, wxEXEC_SYNC);
                    
                    if (ntrmFixResult == 0) {
                        statusText->AppendText("N-terminus element fix executed successfully.\n");
                        
                        if (MakeMol2("ntrm_fixed.pdb", "ntrm.psf", "ntrm.mol2")) {
                            statusText->AppendText("\n*** LYSINE RUN: TWO CGenFF JOBS NEEDED ***\n"
                                                   "A lysine modification also modifies the N-terminus.\n"
                                                   "Upload BOTH naa.mol2 and ntrm.mol2 to CGenFF and import both results\n"
                                                   "(naa.str and ntrm.str) before 'Create Novel Amino Acid'.\n\n");
                            wxMessageBox("Lysine runs modify the N-terminus as well.\n\n"
                                         "Send BOTH naa.mol2 and ntrm.mol2 to CGenFF and import both results "
                                         "(naa.str and ntrm.str) before creating the novel amino acid.",
                                         "Two CGenFF jobs required", wxOK | wxICON_INFORMATION);
                        } else {
                            statusText->AppendText("Error: ntrm.mol2 not produced.\n");
                        }
                    } else {
                        statusText->AppendText("Error: N-terminus element fix failed.\n");
                    }
                } else {
                    statusText->AppendText("Error: N-terminus VMD command failed.\n");
                }
            }
            
        } else {
            statusText->AppendText("Error: Patch maker script failed...\n");
        }
    } else {
        statusText->AppendText("Error: Failed to process patch maker file.\n");
    }
}

void MyFrame::OnCreateTopology(wxCommandEvent& event)
{
    statusText->AppendText("Starting final processing...\n");

    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);
    wxString aminoAcidAbbrev = GetAminoAcidAbbreviation(selectedAminoAcid);

    // --- preflight: the CGenFF results must exist and be newer than the mol2 they came from ---
    auto strReady = [this](const wxString& mol2, const wxString& str) -> bool {
        if (!wxFileExists(str)) {
            statusText->AppendText("Error: " + str + " not found. Send " + mol2 + " to CGenFF and import the result first.\n");
            return false;
        }
        wxFileName m(mol2), t(str);
        if (m.FileExists() && m.GetModificationTime() > t.GetModificationTime()) {
            statusText->AppendText("Error: " + str + " is older than " + mol2 + " - it belongs to a previous molecule. "
                                   "Re-run CGenFF on the current " + mol2 + " and import the new result.\n");
            return false;
        }
        return true;
    };
    bool ready = strReady("naa.mol2", "naa.str");
    if (selectedAminoAcid == "Lysine") ready = strReady("ntrm.mol2", "ntrm.str") && ready;
    if (!ready) {
        wxMessageBox("CGenFF results are missing or stale - see the status log.", "Cannot create novel amino acid", wxOK | wxICON_ERROR);
        return;
    }

    // Make this step re-runnable: drop residues a previous attempt appended
    {
        wxString strip = "python3 strip_rtf_residues.py NAA NPC ntrm_fix naa ntrm";
        statusText->AppendText("Running command: " + strip + "\n");
        system(strip.mb_str());
    }

    // --- NAA path (existing) ---
    wxString command1 = "python3 noX.py";
    statusText->AppendText("Running command: " + command1 + "\n");
    int result1 = system(command1.mb_str());
    
    if (result1 != 0) {
        statusText->AppendText("Error: noX.py failed.\n");
        return;
    }
    statusText->AppendText("noX.py executed successfully.\n");
    
    wxString command2 = "python3 noX_psf.py";
    statusText->AppendText("Running command: " + command2 + "\n");
    int result2 = system(command2.mb_str());
    
    if (result2 != 0) {
        statusText->AppendText("Error: noX_psf.py failed.\n");
        return;
    }
    statusText->AppendText("noX_psf.py executed successfully.\n");
    
    wxString renameTypes = "python3 rename_atom_types_naa.py " + aminoAcidAbbrev + " nad.str";
    statusText->AppendText("Running command: " + renameTypes + "\n");
    int renameResult = system(renameTypes.mb_str());
    
    if (renameResult != 0) {
        statusText->AppendText("Error: rename_atom_types_naa.py failed.\n");
        return;
    }
    statusText->AppendText("Atom type renaming executed successfully.\n");

    // Optional hand-curated charges (see OnSelectChargeOverride). Must land here:
    // after the rename, so the protein-side types are final, and before the
    // parameter steps and ./IC_naa, which read the types and charges from this
    // same nad.str.
    if (!chargeOverrideFile.IsEmpty()) {
        wxString ovr = "python3 apply_charge_overrides.py nad.str \"" +
                       chargeOverrideFile + "\" > charge_overrides.log 2>&1";
        statusText->AppendText("Running command: " + ovr + "\n");
        int ovrResult = system(ovr.mb_str());
        ShowLogTail("charge_overrides.log", 40);
        if (ovrResult != 0) {
            statusText->AppendText("Error: apply_charge_overrides.py refused or failed. "
                                   "nad.str is unchanged and the run is stopping, because "
                                   "continuing would build the residue with CGenFF's "
                                   "uncorrected charges.\n");
            wxMessageBox("Charge overrides were not applied - see the log in the status "
                         "pane. The run has stopped so the residue is not built with the "
                         "wrong charges.", "Overrides not applied", wxOK | wxICON_ERROR);
            return;
        }
        statusText->AppendText("Charge overrides applied to nad.str.\n");
    }
    
    wxString command0 = "python3 fix_atom_naming_pdb_naa.py";
    statusText->AppendText("Running command: " + command0 + "\n");
    int result0 = system(command0.mb_str());
    
    if (result0 != 0) {
        statusText->AppendText("Error: fix_atom_naming_pdb_naa.py failed.\n");
        return;
    }
    statusText->AppendText("fix_atom_naming_pdb_naa.py executed successfully.\n");

    // Novel CGenFF + cross-force-field parameters for the modified side chain.
    // par_GUI.prm is the third parameter file a NAMD run loads, alongside
    // par_all36m_prot.prm and par_all36_cgenff.prm, so it must carry every term
    // the released sets lack. Without this the junction dihedrals stayed in
    // nad.str and never reached the simulation.
    {
        wxString addNaa = "python3 add_naa_params_to_prm.py nad.str " + aminoAcidAbbrev + " NAA";
        statusText->AppendText("Running command: " + addNaa + "\n");
        if (system(addNaa.mb_str()) != 0) {
            statusText->AppendText("Error: add_naa_params_to_prm.py failed.\n");
            return;
        }
        statusText->AppendText("NAA parameters added to par_GUI.prm.\n");

        // Renaming protein-side types orphans any term CGenFF satisfied from its
        // own release library: it is absent from the .str, so nothing re-emits it
        // under the new type. This copies those values across.
        wxString xff = "python3 complete_cross_ff.py naa.str nad.str NAA";
        statusText->AppendText("Running command: " + xff + "\n");
        if (system(xff.mb_str()) != 0) {
            statusText->AppendText("Error: complete_cross_ff.py failed for NAA.\n");
            return;
        }

        // Catch a malformed parameter file here rather than in NAMD, where it
        // surfaces as "FATAL ERROR: ABNORMAL EOF FOUND - buffer=*END*".
        wxString checkPrm = "python3 validate_prm.py par_GUI.prm";
        statusText->AppendText("Running command: " + checkPrm + "\n");
        if (system(checkPrm.mb_str()) != 0) {
            statusText->AppendText("Error: par_GUI.prm failed validation (see terminal). "
                                   "NAMD would abort while reading it.\n");
            wxMessageBox("par_GUI.prm is malformed - see the terminal output. "
                         "Click 'Run Renaming' to rebuild it from par_GUI_template.prm.",
                         "Parameter file invalid", wxOK | wxICON_ERROR);
            return;
        }
        statusText->AppendText("par_GUI.prm validated.\n");
    }

    wxString command3 = "./IC_naa";
    statusText->AppendText("Running command: " + command3 + "\n");
    int result3 = system(command3.mb_str());
    
    if (result3 != 0) {
        statusText->AppendText("Error: IC_naa failed.\n");
        return;
    }
    statusText->AppendText("IC_naa executed successfully.\n");

    {
        wxString v = "python3 validate_rtf_ic.py NAA";
        statusText->AppendText("Running command: " + v + "\n");
        if (system(v.mb_str()) != 0) {
            statusText->AppendText("Error: RESI NAA failed validation (see terminal). "
                                   "psfgen would place atoms from degenerate IC frames.\n");
            wxMessageBox("RESI NAA failed its internal-coordinate check - see the terminal output. "
                         "Do not build a capsid from this topology.", "Topology validation failed",
                         wxOK | wxICON_ERROR);
            return;
        }
        statusText->AppendText("RESI NAA internal-coordinate check passed.\n");
    }
    
    // --- NPC path (N-terminus, for Lysine) ---
    
    if (selectedAminoAcid == "Lysine") {
        wxString ntrmNoX = "python3 noX_ntrm.py";
        statusText->AppendText("Running command: " + ntrmNoX + "\n");
        int ntrmNoXResult = system(ntrmNoX.mb_str());
        
        if (ntrmNoXResult != 0) {
            statusText->AppendText("Error: noX_ntrm.py failed.\n");
            return;
        }
        statusText->AppendText("noX_ntrm.py executed successfully.\n");
        
        wxString ntrmNoXPsf = "python3 noX_psf_ntrm.py";
        statusText->AppendText("Running command: " + ntrmNoXPsf + "\n");
        int ntrmNoXPsfResult = system(ntrmNoXPsf.mb_str());
        
        if (ntrmNoXPsfResult != 0) {
            statusText->AppendText("Error: noX_psf_ntrm.py failed.\n");
            return;
        }
        statusText->AppendText("noX_psf_ntrm.py executed successfully.\n");
        
        wxString chargeFix = "python3 charge_fix_ntrm.py";
        statusText->AppendText("Running command: " + chargeFix + "\n");
        int chargeFixResult = system(chargeFix.mb_str());
        
        if (chargeFixResult != 0) {
            statusText->AppendText("Error: charge_fix_ntrm.py failed.\n");
            return;
        }
    statusText->AppendText("charge_fix_ntrm.py executed successfully.\n");

        wxString renameTypes = "python3 rename_atom_types_naa.py " + NtrmVariant() + " ntrm_clean.str";
        statusText->AppendText("Running command: " + renameTypes + "\n");
        int renameResult = system(renameTypes.mb_str());
        
        if (renameResult != 0) {
            statusText->AppendText("Error: rename_atom_types_naa.py failed.\n");
            return;
        }
        statusText->AppendText("rename_atom_types_naa.py executed successfully.\n");

        wxString ntrmFixNaming = "python3 fix_atom_naming_pdb_ntrm.py";
        statusText->AppendText("Running command: " + ntrmFixNaming + "\n");
        int ntrmFixResult = system(ntrmFixNaming.mb_str());
        
        if (ntrmFixResult != 0) {
            statusText->AppendText("Error: fix_atom_naming_pdb_ntrm.py failed.\n");
            return;
        }
        statusText->AppendText("fix_atom_naming_pdb_ntrm.py executed successfully.\n");
        
        wxString addParams = "python3 add_npc_params_to_prm.py " + NtrmVariant();
        statusText->AppendText("Running command: " + addParams + "\n");
        int addParamsResult = system(addParams.mb_str());

        if (addParamsResult != 0) {
            statusText->AppendText("Error: add_npc_params_to_prm.py failed.\n");
            return;
        }
statusText->AppendText("NPC cross-FF parameters added to par_GUI.prm.\n");

        wxString xffN = "python3 complete_cross_ff.py ntrm.str ntrm_clean.str NPC";
        statusText->AppendText("Running command: " + xffN + "\n");
        if (system(xffN.mb_str()) != 0) {
            statusText->AppendText("Error: complete_cross_ff.py failed for NPC.\n");
            return;
        }

        wxString ntrmIC = "./IC_ntrm_NPC";
        statusText->AppendText("Running command: " + ntrmIC + "\n");
        int ntrmICResult = system(ntrmIC.mb_str());
        
        if (ntrmICResult != 0) {
            statusText->AppendText("Error: IC_ntrm_NPC failed.\n");
            return;
        }
        statusText->AppendText("IC_ntrm_NPC executed successfully.\n");

        {
            wxString v = "python3 validate_rtf_ic.py NPC";
            statusText->AppendText("Running command: " + v + "\n");
            if (system(v.mb_str()) != 0) {
                statusText->AppendText("Error: RESI NPC failed validation (see terminal).\n");
                wxMessageBox("RESI NPC failed its internal-coordinate check - see the terminal output.",
                             "Topology validation failed", wxOK | wxICON_ERROR);
                return;
            }
            statusText->AppendText("RESI NPC internal-coordinate check passed.\n");
        }
    }
    
    statusText->AppendText("Final processing completed successfully!\n");
}


// Print the end of a NAMD/VMD log into the status box. A missing parameter or a
// bad config aborts NAMD during setup, and the reason is only in its output.
void MyFrame::ShowLogTail(const wxString& logFile, int lines)
{
    wxTextFile log;
    if (!wxFileExists(logFile) || !log.Open(logFile)) {
        statusText->AppendText("  (no log at " + logFile + " - NAMD may not have started; "
                               "check the binary path)\n");
        return;
    }
    int n = log.GetLineCount();
    statusText->AppendText("  --- last lines of " + logFile + " ---\n");
    for (int i = (n > lines ? n - lines : 0); i < n; i++)
        statusText->AppendText("  " + log.GetLine(i) + "\n");
    statusText->AppendText("  --- end of log ---\n");
    log.Close();
}

// ── stage 1: build the modified monomer(s) ──────────────────────────────────
// The coat-protein templates only LABEL the target residues NAA/NPC; the real
// residue, with its drug adduct, is constructed here by psfgen from RESI NAA /
// RESI NPC and guesscoord. Splitting this out from assembly is what lets the
// monomer be minimised before the symmetry matrices replicate it.
wxArrayString MyFrame::MonomerStems()
{
    wxArrayString stems;
    if (vlpTypeChoice->GetSelection() == 1) {
        stems.Add("TMV_naa");
    } else {
        stems.Add("fix_a");
        stems.Add("fix_b");
        stems.Add("fix_c");
    }
    return stems;
}

bool MyFrame::BuildMonomers()
{
    int selection = aminoAcidChoice->GetSelection();
    wxString selectedAminoAcid = aminoAcidChoice->GetString(selection);

    if (!processor.ProcessPgnFile(std::string(selectedAminoAcid.mb_str()))) {
        statusText->AppendText("Error: could not generate the monomer pgn file(s) from templates.\n");
        return false;
    }

    bool tmv = (vlpTypeChoice->GetSelection() == 1);
    if (tmv && !wxFileExists("TMV_naa.pgn")) {
        statusText->AppendText("Error: no TMV coat-protein template for " + selectedAminoAcid +
                               ". TMV supports Tyrosine (tmv_Y139.pdb) and Lysine "
                               "(tmv_K53_K68.pdb) only.\n");
        wxMessageBox("TMV has no coat-protein template for " + selectedAminoAcid + ".",
                     "TMV template missing", wxOK | wxICON_ERROR);
        return false;
    }

    wxArrayString stems = MonomerStems();
    for (size_t i = 0; i < stems.GetCount(); i++) {
        wxString pgn = tmv ? "TMV_naa.pgn" : (stems[i] + ".pgn");
        wxString cmd = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/"
                       "startup.command.csh -dispdev text -e " + pgn;
        statusText->AppendText("Running command: " + cmd + "\n");
        if (RunVmd(cmd, stems[i] + ".pdb") != 0) {
            statusText->AppendText("Error: " + pgn + " failed.\n");
            return false;
        }
        statusText->AppendText(pgn + " built " + stems[i] + ".pdb / .psf\n");
        // a fresh build invalidates any earlier minimisation of this monomer
        if (wxFileExists(stems[i] + "_premin.pdb")) wxRemoveFile(stems[i] + "_premin.pdb");
        if (wxFileExists(stems[i] + "_min.coor")) wxRemoveFile(stems[i] + "_min.coor");
    }
    return true;
}

void MyFrame::OnBuildMonomer(wxCommandEvent& event)
{
    statusText->AppendText("=== Building the modified monomer(s) ===\n");
    if (!BuildMonomers()) {
        wxMessageBox("Monomer build failed - see the status log.", "Build failed",
                     wxOK | wxICON_ERROR);
        return;
    }
    statusText->AppendText("Monomer(s) built. Next: 'Minimize Monomer' to relieve the "
                           "geometry guesscoord built, then 'Assemble And View Capsid'.\n\n");
}

void MyFrame::OnCreateViewCapsid(wxCommandEvent& event)
{
    int vlpSelection = vlpTypeChoice->GetSelection();
    if (vlpSelection != 0 && vlpSelection != 1) {
        statusText->AppendText("Error: Unknown VLP type selected.\n");
        return;
    }

    // Build the monomers only if they are missing. A monomer that is already
    // present may have been minimised, and rebuilding would silently discard
    // that work, so it is left alone.
    wxArrayString stems = MonomerStems();
    wxArrayString missing;
    for (size_t i = 0; i < stems.GetCount(); i++)
        if (!wxFileExists(stems[i] + ".pdb")) missing.Add(stems[i]);

    if (!missing.IsEmpty()) {
        statusText->AppendText("Monomer(s) not built yet (" +
                               wxJoin(missing, ',') + ") - building them now.\n");
        if (!BuildMonomers()) return;
    } else {
        int minimised = 0;
        for (size_t i = 0; i < stems.GetCount(); i++)
            if (wxFileExists(stems[i] + "_premin.pdb")) minimised++;
        statusText->AppendText(wxString::Format(
            "Using the existing monomer(s); %d of %d have been minimised.\n",
            minimised, (int)stems.GetCount()));
    }

    if (vlpSelection == 0) CreateViewQbeta();
    else                   CreateViewTMV();
}

void MyFrame::CreateViewQbeta()
{
    // Monomers must already exist; BuildMonomers() is called by the caller only
    // when they are absent, so a minimised monomer is never silently rebuilt.
    statusText->AppendText("Assembling the Qbeta capsid from fix_a/b/c.pdb...\n");

    wxString maker = "python3 Qbeta_capsid_maker.py";
    statusText->AppendText("Running command: " + maker + "\n");
    if (system(maker.mb_str()) != 0) {
        statusText->AppendText("Error: Qbeta_capsid_maker.py failed.\n");
        return;
    }

    wxString writer = "python3 Qbeta_pgn_writer.py";
    statusText->AppendText("Running command: " + writer + "\n");
    if (system(writer.mb_str()) != 0) {
        statusText->AppendText("Error: Qbeta_pgn_writer.py failed.\n");
        return;
    }

    wxString assemble = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e " "QB_capsid.pgn";
    statusText->AppendText("Running command: " + assemble + "\n");
    if (RunVmd(assemble, "capsid_naa.pdb") != 0) {
        statusText->AppendText("Error: QB_capsid.pgn failed.\n");
        return;
    }
    statusText->AppendText("QB_capsid.pgn executed successfully.\n");

    wxString view = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -e capsid_zoom_naa.tcl";
    statusText->AppendText("Running command: " + view + "\n");
    if (RunVmd(view, "") != 0) {
        statusText->AppendText("Error: capsid_zoom_naa.tcl failed.\n");
        return;
    }
    statusText->AppendText("Qbeta capsid assembled and displayed.\n\n");
}

void MyFrame::CreateViewTMV()
{
    statusText->AppendText("Assembling the TMV rod from TMV_naa.pdb...\n");

    wxString maker = "python3 TMV_capsid_maker.py";
    statusText->AppendText("Running command: " + maker + "\n");
    if (system(maker.mb_str()) != 0) {
        statusText->AppendText("Error: TMV_capsid_maker.py failed.\n");
        return;
    }

    wxString writer = "python3 TMV_pgn_writer.py";
    statusText->AppendText("Running command: " + writer + "\n");
    if (system(writer.mb_str()) != 0) {
        statusText->AppendText("Error: TMV_pgn_writer.py failed.\n");
        return;
    }

    wxString assemble = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -dispdev text -e " "TMV_all_out.pgn";
    statusText->AppendText("Running command: " + assemble + "\n");
    if (RunVmd(assemble, "TMV_naa_capsid.pdb") != 0) {
        statusText->AppendText("Error: TMV_all_out.pgn failed.\n");
        return;
    }
    statusText->AppendText("TMV_all_out.pgn executed successfully.\n");

    wxString view = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh -e capsid_zoom_naa_TMV.tcl";
    statusText->AppendText("Running command: " + view + "\n");
    if (RunVmd(view, "") != 0) {
        statusText->AppendText("Error: capsid_zoom_naa_TMV.tcl failed.\n");
        return;
    }
    statusText->AppendText("TMV rod assembled and displayed.\n\n");
}


// ── monomer minimisation, before the capsid is assembled ────────────────────
// Relieves the geometry guesscoord built for the modified residues while every
// crystallographic heavy atom stays fixed, so the symmetry matrices that
// assemble the capsid remain valid. Clashes BETWEEN neighbouring subunits'
// adducts cannot be seen here - that is what "Minimize Capsid" is for.
bool MyFrame::MinimizeOneMonomer(const wxString& stem)
{
    if (!wxFileExists(stem + ".pdb") || !wxFileExists(stem + ".psf")) {
        statusText->AppendText("Error: " + stem + ".pdb / .psf not found. Build the monomer first.\n");
        return false;
    }

    // 1. flag which atoms may move (psfgen marks guessed coordinates occupancy 0)
    wxString flags = "python3 write_fixed_atoms.py " + stem + ".pdb " + stem +
                     "_fix.pdb --only-residues NAA NPC";
    statusText->AppendText("Running command: " + flags + "\n");
    if (system(flags.mb_str()) != 0) {
        statusText->AppendText("Error: write_fixed_atoms.py failed for " + stem + ".\n");
        return false;
    }

    // 2. fill the NAMD template for this monomer
    {
        wxTextFile tpl;
        if (!tpl.Open("minimize_monomer_template.conf")) {
            statusText->AppendText("Error: minimize_monomer_template.conf not found.\n");
            return false;
        }
        wxString conf;
        for (wxString l = tpl.GetFirstLine(); ; l = tpl.GetNextLine()) {
            l.Replace("{{STEM}}", stem, true);
            conf << l << "\n";
            if (tpl.Eof()) break;
        }
        tpl.Close();
        wxFile out("minimize_monomer.conf", wxFile::write);
        out.Write(conf);
    }

    // 3. minimise
    // NAMD takes ONE config file and writes to stdout, so the log needs a shell
    // redirect; passing it as a second argument gives NAMD two config files and
    // it writes nothing at all.
    wxString logFile = stem + "_min.log";
    wxString namd = "~/Downloads/NAMD_3.0b5_MacOS-universal-multicore/namd3 +p8 "
                    "minimize_monomer.conf > " + logFile + " 2>&1";
    statusText->AppendText("Running command: " + namd + "\n");
    int rc = system(namd.mb_str());
    if (rc != 0 || !wxFileExists(stem + "_min.coor")) {
        statusText->AppendText("Error: NAMD minimisation failed for " + stem +
                               " (exit " + wxString::Format("%d", rc) + ").\n");
        ShowLogTail(logFile);
        return false;
    }

    // 4. write the minimised coordinates back over the monomer PDB, keeping the original
    if (!wxFileExists(stem + "_premin.pdb"))
        wxCopyFile(stem + ".pdb", stem + "_premin.pdb", false);
    {
        wxTextFile tpl;
        if (!tpl.Open("coor_to_pdb_template.tcl")) {
            statusText->AppendText("Error: coor_to_pdb_template.tcl not found.\n");
            return false;
        }
        wxString tcl;
        for (wxString l = tpl.GetFirstLine(); ; l = tpl.GetNextLine()) {
            l.Replace("{{STEM}}", stem, true);
            tcl << l << "\n";
            if (tpl.Eof()) break;
        }
        tpl.Close();
        wxFile out("coor_to_pdb.tcl", wxFile::write);
        out.Write(tcl);
    }
    wxString vmd = "/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/MacOS/startup.command.csh "
                   "-dispdev text -e coor_to_pdb.tcl";
    statusText->AppendText("Running command: " + vmd + "\n");
    system(vmd.mb_str());

    statusText->AppendText(stem + ".pdb updated from the minimisation (original kept as " +
                           stem + "_premin.pdb).\n\n");
    return true;
}

void MyFrame::OnMinimizeMonomer(wxCommandEvent& event)
{
    statusText->AppendText("=== Minimising the modified residues, protein fixed ===\n");

    wxArrayString stems = MonomerStems();
    int done = 0;
    for (size_t i = 0; i < stems.GetCount(); i++) {
        if (!wxFileExists(stems[i] + ".pdb")) {
            statusText->AppendText("Skipping " + stems[i] + " (not built).\n");
            continue;
        }
        if (!MinimizeOneMonomer(stems[i])) {
            wxMessageBox("Monomer minimisation failed for " + stems[i] +
                         " - see the status log.", "Minimisation failed", wxOK | wxICON_ERROR);
            return;
        }
        done++;
    }
    if (!done) {
        statusText->AppendText("No monomers found to minimise. Run the capsid step's first "
                               "stage (or Create Novel Amino Acid) first.\n");
        return;
    }
    statusText->AppendText("Monomer minimisation complete for " + wxString::Format("%d", done) +
                           " subunit(s). Now run 'Create And View Capsid' to replicate and "
                           "assemble, then 'Minimize Capsid' to resolve contacts between "
                           "neighbouring subunits.\n\n");
}

void MyFrame::OnMinimizeCapsid(wxCommandEvent& event)
{
    // Each VLP has its own minimize conf: box geometry, structure name, and cell
    // differ between Qbeta (cubic 320 Å) and TMV (rectangular rod 209x209x121 Å).
    bool tmv = (vlpTypeChoice->GetSelection() == 1);
    wxString confFile  = tmv ? "TMV_minimize.conf"  : "minimize.conf";
    wxString wbPdb     = tmv ? "TMV_naa_capsid_wb.pdb" : "QB_naa_capsid_wb.pdb";
    wxString logFile   = "minimize_capsid.log";

    statusText->AppendText("Minimize Capsid: using " + confFile + "\n");

    if (!wxFileExists(wbPdb)) {
        statusText->AppendText("Error: " + wbPdb + " not found. Run Create/View Capsid "
                               "first (it solvates and ionizes the capsid).\n");
        return;
    }
    // NAMD does not read CRYST1, so derive the periodic cell from the solvated
    // PDB itself (box from CRYST1, origin from the geometric centre). The conf
    // sources the <stem>_cell.namd this writes and falls back to its literal
    // values if the file is missing.
    wxString cellCmd = "python3 capsid_cell.py " + wbPdb;
    statusText->AppendText("Running command: " + cellCmd + "\n");
    if (system(cellCmd.mb_str()) != 0) {
        statusText->AppendText("Warning: capsid_cell.py failed; " + confFile +
                               " will use its built-in cell values.\n");
    }

    wxString command9 = "~/Downloads/NAMD_3.0b5_MacOS-universal-multicore/namd3 +p8 " +
                        confFile + " > " + logFile + " 2>&1";
    statusText->AppendText("Running command: " + command9 + "\n");

    int result = system(command9.mb_str());

    if (result == 0) {
        statusText->AppendText("NAMD minimize command executed successfully "
                               "(log: " + logFile + ").\n");
    } else {
        statusText->AppendText("Error: NAMD minimize command failed (exit " +
                               wxString::Format("%d", result) + ").\n");
        ShowLogTail(logFile);
    }
}

void MyFrame::OnZipFiles(wxCommandEvent& event)
{
    statusText->AppendText("Starting file compression...\n");
    
    // The solvated capsid is named after the VLP, so the zip must follow the
    // selection rather than always naming the Qbeta files.
    wxString capsidStem = (vlpTypeChoice->GetSelection() == 1)
                              ? "TMV_naa_capsid_wb" : "QB_naa_capsid_wb";
    statusText->AppendText("Packaging " + capsidStem + " (" +
                           vlpTypeChoice->GetString(vlpTypeChoice->GetSelection()) + ")\n");

    wxArrayString filesToZip;
    filesToZip.Add(capsidStem == "TMV_naa_capsid_wb" ? "TMV_minimize.conf" : "minimize.conf");
    if (wxFileExists(capsidStem + "_cell.namd"))
        filesToZip.Add(capsidStem + "_cell.namd");   // periodic cell the conf sources
    filesToZip.Add(capsidStem + ".psf");
    filesToZip.Add(capsidStem + ".pdb");
    filesToZip.Add("par_GUI.prm");          // novel CGenFF + cross-FF parameters
    filesToZip.Add("par_all36_cgenff.prm"); // the other three parameter files the conf reads
    filesToZip.Add("par_all36m_prot.prm");
    filesToZip.Add("par_water_ions.inp");
    filesToZip.Add("top_all36_cgenff_CBD.rtf");  // RESI NAA / RESI NPC topology
    
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
    
    wxDateTime now = wxDateTime::Now();
    wxString timestamp = now.Format("%Y%m%d_%H%M%S");
    wxString zipFileName = "capsid_files_" + timestamp + ".zip";
    
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
        
        wxMessageBox("Files successfully compressed to:\n" + zipFileName, 
                     "Success", wxOK | wxICON_INFORMATION);
        
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
        wxMessageBox("Failed to create zip file.", "Error", wxOK | wxICON_ERROR);
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