#include <wx/wx.h>
#include "FileProcessor.h" // Include your custom file processing class
#include <cstdlib> // for system()

class MyApp : public wxApp {
public:
    virtual bool OnInit();
};

class MyFrame : public wxFrame {
private:
    wxString selectedPDB;

    void OnSelectFile(wxCommandEvent& event);
    void OnProcessFile(wxCommandEvent& event);
    void OnRunPython(wxCommandEvent& event);
    void GetMatrix(wxCommandEvent& event);
    void RotateProtein(wxCommandEvent& event);
    void CreateCapsid(wxCommandEvent& event);
    void CreatePGN(wxCommandEvent& event);
    void RunPGN(wxCommandEvent& event);

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
    : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(300, 400)) {

    wxPanel* panel = new wxPanel(this, wxID_ANY);

    wxButton* selectFileBtn = new wxButton(panel, wxID_ANY, "Select PDB File", wxPoint(50, 10));
    wxButton* processFileBtn = new wxButton(panel, wxID_ANY, "Process Python File", wxPoint(50, 50));
    wxButton* runPythonBtn = new wxButton(panel, wxID_ANY, "Run Python Program", wxPoint(50, 90));
    wxButton* GetMatrixBtn = new wxButton(panel, wxID_ANY, "Get Rotation Matrix", wxPoint(50, 130));
    wxButton* RotateProteinBtn = new wxButton(panel, wxID_ANY, "Implement Rotation", wxPoint(50, 170));
    wxButton* CreateCapsidBtn = new wxButton(panel, wxID_ANY, "Create Capsid", wxPoint(50, 210));
    wxButton* CreatePGNBtn = new wxButton(panel, wxID_ANY, "Create PGN File", wxPoint(50, 250));
    wxButton* RunPGNBtn = new wxButton(panel, wxID_ANY, "Run PGN File", wxPoint(50, 290));

    selectFileBtn->Bind(wxEVT_BUTTON, &MyFrame::OnSelectFile, this);
    processFileBtn->Bind(wxEVT_BUTTON, &MyFrame::OnProcessFile, this);
    runPythonBtn->Bind(wxEVT_BUTTON, &MyFrame::OnRunPython, this);
    GetMatrixBtn->Bind(wxEVT_BUTTON, &MyFrame::GetMatrix, this);
    RotateProteinBtn->Bind(wxEVT_BUTTON, &MyFrame::RotateProtein, this);
    CreateCapsidBtn->Bind(wxEVT_BUTTON, &MyFrame::CreateCapsid, this);
    CreatePGNBtn->Bind(wxEVT_BUTTON, &MyFrame::CreatePGN, this);
    RunPGNBtn->Bind(wxEVT_BUTTON, &MyFrame::CreatePGN, this);
}

void MyFrame::OnSelectFile(wxCommandEvent& event) {
    wxFileDialog openFileDialog(this, _("Select PDB file"), "", "", 
                                "PDB files (*.pdb)|*.pdb", wxFD_OPEN | wxFD_FILE_MUST_EXIST);

    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;

    selectedPDB = openFileDialog.GetPath();
    wxMessageBox("Selected File: " + selectedPDB);
}

void MyFrame::OnProcessFile(wxCommandEvent& event) {
    if (selectedPDB.IsEmpty()) {
        wxMessageBox("Please select a PDB file first.");
        return;
    }

    wxString baseName = selectedPDB.BeforeLast('.'); // Get base name
    FileProcessor::ProcessPythonFile("Patch_orient_QB.py", baseName.ToStdString());
    wxMessageBox("File processed successfully.");
}

void MyFrame::OnRunPython(wxCommandEvent& event) {
    if (selectedPDB.IsEmpty()) {
        wxMessageBox("Please select and process a PDB file first.");
        return;
    }

    system("Python3 Patch_orient_QB.py"); // Run the Python program
    FileProcessor::ResetPythonFile("Patch_orient_QB.py");
    wxMessageBox("Python program executed and file reset.");
}

void  MyFrame::GetMatrix(wxCommandEvent& event) {
    system("./rot_mat");
}

void  MyFrame::RotateProtein(wxCommandEvent& event) {
    system("Python3 rotate_protein.py");
}

void MyFrame::CreateCapsid(wxCommandEvent& event) {
    system("Python3 Qbeta_capsid_maker.py");
}

void MyFrame::CreatePGN(wxCommandEvent& event) {
    system("Python3 Qbeta_pgn_writer.py");
}

void MyFrame::RunPGN(wxCommandEvent& event) {
    system("vmd -dispdev text -e QB_capsid.pgn");
}
