#ifndef MOLECULEVIEWER_H
#define MOLECULEVIEWER_H

#include <wx/wx.h>
#include <wx/glcanvas.h>
#include <vector>
#include <string>
#include <cmath>
#include <functional>

// ── Data structures ──────────────────────────────────────────────────────────

struct Atom {
    std::string name;      // e.g. "C1", "N2", "O3"
    std::string element;   // e.g. "C", "N", "O"
    float x, y, z;         // 3-D coordinates from PDB
    int serial;            // atom serial number
};

struct Bond {
    int atomIndex1;
    int atomIndex2;
};

// ── Selection mode ───────────────────────────────────────────────────────────

enum class SelectionMode {
    BONDED,   // selecting the 3 bonded atoms
    DELETED   // selecting the 3 deleted atoms
};

// ── Molecule viewer widget ───────────────────────────────────────────────────

class MoleculeViewer : public wxGLCanvas
{
public:
    MoleculeViewer(wxWindow* parent, wxWindowID id = wxID_ANY,
                   const wxPoint& pos = wxDefaultPosition,
                   const wxSize& size = wxDefaultSize);
    ~MoleculeViewer();

    // Load molecule from a PDB file; returns false on failure
    bool LoadPDB(const std::string& filePath);

    // Clear everything
    void ClearMolecule();

    // Selection
    void SetSelectionMode(SelectionMode mode);
    SelectionMode GetSelectionMode() const { return selectionMode; }
    void ClearSelections();
    void ClearBondedSelections();
    void ClearDeletedSelections();

    // Retrieve selected atom names
    std::vector<std::string> GetBondedAtomNames() const;
    std::vector<std::string> GetDeletedAtomNames() const;

    // Callback fired whenever selection changes
    // signature: void callback(SelectionMode mode, vector<string> names)
    using SelectionCallback = std::function<void(SelectionMode, const std::vector<std::string>&)>;
    void SetSelectionCallback(SelectionCallback cb) { selectionCallback = cb; }

    // Query
    bool HasMolecule() const { return !atoms.empty(); }
    int  GetAtomCount()  const { return (int)atoms.size(); }
    int  GetBondCount()  const { return (int)bonds.size(); }

private:
    // ── OpenGL helpers ───────────────────────────────────────────────────────
    void InitGL();
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnMouse(wxMouseEvent& event);

    void Render();
    void DrawAtom(const Atom& atom, bool selected, bool hovered, int selectionOrder);
    void DrawBond(const Bond& bond);

    // Sphere / cylinder approximations using GLU-style subdivision
    void DrawSphere(float cx, float cy, float cz, float radius, int slices, int stacks);
    void DrawCylinder(float x1, float y1, float z1,
                      float x2, float y2, float z2, float radius, int slices);

    // Hit-testing: which atom is under (screenX, screenY)?  -1 if none
    int  PickAtom(int screenX, int screenY);

    // Atom colour by element
    void GetElementColor(const std::string& element, float& r, float& g, float& b);
    float GetElementRadius(const std::string& element);

    // Bond detection (distance-based)
    void DetectBonds();

    // Centre / scale the view around the loaded molecule
    void ComputeBounds();

    // ── Data ─────────────────────────────────────────────────────────────────
    wxGLContext* glContext;
    bool         glInitialized;

    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    // Camera
    float rotX, rotY;          // rotation angles (degrees)
    float zoom;                // distance multiplier
    float panX, panY;          // panning offset
    float centerX, centerY, centerZ; // molecule center
    float molRadius;           // bounding-sphere radius

    // Mouse tracking
    bool  dragging;
    bool  panning;
    int   lastMouseX, lastMouseY;
    int   hoveredAtom;         // index or -1

    // Selection state
    SelectionMode selectionMode;
    std::vector<int> bondedSelected;   // up to 3 indices
    std::vector<int> deletedSelected;  // up to 3 indices

    SelectionCallback selectionCallback;

    wxDECLARE_EVENT_TABLE();
};

#endif // MOLECULEVIEWER_H