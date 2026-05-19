// MoleculeViewer.cpp
#include "MoleculeViewer.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ── GL canvas attributes ─────────────────────────────────────────────────────

static const int GL_ATTRIBS[] = {
    WX_GL_RGBA,
    WX_GL_DOUBLEBUFFER,
    WX_GL_DEPTH_SIZE, 16,
    0   // removed multisampling to avoid issues on some Macs
};

// ── Event table ──────────────────────────────────────────────────────────────

wxBEGIN_EVENT_TABLE(MoleculeViewer, wxGLCanvas)
    EVT_PAINT(MoleculeViewer::OnPaint)
    EVT_SIZE(MoleculeViewer::OnSize)
    EVT_LEFT_DOWN(MoleculeViewer::OnMouse)
    EVT_LEFT_UP(MoleculeViewer::OnMouse)
    EVT_RIGHT_DOWN(MoleculeViewer::OnMouse)
    EVT_RIGHT_UP(MoleculeViewer::OnMouse)
    EVT_MOUSEWHEEL(MoleculeViewer::OnMouse)
    EVT_MOTION(MoleculeViewer::OnMouse)
wxEND_EVENT_TABLE()

// ── Constructor / destructor ─────────────────────────────────────────────────

MoleculeViewer::MoleculeViewer(wxWindow* parent, wxWindowID id,
                               const wxPoint& pos, const wxSize& size)
    : wxGLCanvas(parent, id, GL_ATTRIBS, pos, size,
                 wxFULL_REPAINT_ON_RESIZE)
    , glContext(nullptr)
    , glInitialized(false)
    , rotX(20.0f), rotY(-30.0f)
    , zoom(1.0f)
    , panX(0.0f), panY(0.0f)
    , centerX(0), centerY(0), centerZ(0)
    , molRadius(10.0f)
    , dragging(false), panning(false)
    , lastMouseX(0), lastMouseY(0)
    , hoveredAtom(-1)
    , selectionMode(SelectionMode::BONDED)
{
    glContext = new wxGLContext(this);
}

MoleculeViewer::~MoleculeViewer()
{
    delete glContext;
}

// ── OpenGL initialization ────────────────────────────────────────────────────

void MoleculeViewer::InitGL()
{
    SetCurrent(*glContext);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    float lightPos[]     = { 1.0f, 1.0f, 1.0f, 0.0f };
    float lightAmbient[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    float lightDiffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    float lightSpecular[]= { 0.4f, 0.4f, 0.4f, 1.0f };

    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

    glClearColor(0.15f, 0.15f, 0.2f, 1.0f);

    glInitialized = true;
}

// ── PDB loading ──────────────────────────────────────────────────────────────

bool MoleculeViewer::LoadPDB(const std::string& filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "MoleculeViewer: cannot open " << filePath << "\n";
        return false;
    }

    atoms.clear();
    bonds.clear();
    bondedSelected.clear();
    deletedSelected.clear();
    hoveredAtom = -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line.size() < 54) continue;
        std::string recType = line.substr(0, 6);
        recType.erase(recType.find_last_not_of(' ') + 1);

        if (recType != "ATOM" && recType != "HETATM") continue;

        Atom a;
        a.serial = std::atoi(line.substr(6, 5).c_str());

        a.name = line.substr(12, 4);
        a.name.erase(0, a.name.find_first_not_of(' '));
        a.name.erase(a.name.find_last_not_of(' ') + 1);

        a.x = std::atof(line.substr(30, 8).c_str());
        a.y = std::atof(line.substr(38, 8).c_str());
        a.z = std::atof(line.substr(46, 8).c_str());

        if (line.size() >= 78) {
            a.element = line.substr(76, 2);
            a.element.erase(0, a.element.find_first_not_of(' '));
            a.element.erase(a.element.find_last_not_of(' ') + 1);
        }
        if (a.element.empty()) {
            for (char c : a.name) {
                if (std::isalpha(c)) {
                    a.element = std::string(1, std::toupper(c));
                    break;
                }
            }
        }

        atoms.push_back(a);
    }

    if (atoms.empty()) {
        std::cerr << "MoleculeViewer: no atoms found in " << filePath << "\n";
        return false;
    }

    std::cout << "MoleculeViewer: loaded " << atoms.size() << " atoms from " << filePath << "\n";

    DetectBonds();
    ComputeBounds();

    rotX = 20.0f;
    rotY = -30.0f;
    zoom = 1.0f;
    panX = panY = 0.0f;

    Refresh(false);
    return true;
}

void MoleculeViewer::ClearMolecule()
{
    atoms.clear();
    bonds.clear();
    bondedSelected.clear();
    deletedSelected.clear();
    hoveredAtom = -1;
    Refresh(false);
}

// ── Bond detection ───────────────────────────────────────────────────────────

void MoleculeViewer::DetectBonds()
{
    bonds.clear();

    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            float dx = atoms[i].x - atoms[j].x;
            float dy = atoms[i].y - atoms[j].y;
            float dz = atoms[i].z - atoms[j].z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);

            float maxDist = 1.9f;
            if (atoms[i].element == "S" || atoms[j].element == "S")
                maxDist = 2.3f;
            if (atoms[i].element == "P" || atoms[j].element == "P")
                maxDist = 2.2f;

            if (dist < 0.4f) continue;

            if (dist <= maxDist) {
                bonds.push_back({(int)i, (int)j});
            }
        }
    }

    std::cout << "MoleculeViewer: detected " << bonds.size() << " bonds\n";
}

// ── Bounding box ─────────────────────────────────────────────────────────────

void MoleculeViewer::ComputeBounds()
{
    if (atoms.empty()) return;

    float minX = atoms[0].x, maxX = atoms[0].x;
    float minY = atoms[0].y, maxY = atoms[0].y;
    float minZ = atoms[0].z, maxZ = atoms[0].z;

    for (const auto& a : atoms) {
        if (a.x < minX) minX = a.x; if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y; if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z; if (a.z > maxZ) maxZ = a.z;
    }

    centerX = (minX + maxX) * 0.5f;
    centerY = (minY + maxY) * 0.5f;
    centerZ = (minZ + maxZ) * 0.5f;

    float dx = maxX - minX;
    float dy = maxY - minY;
    float dz = maxZ - minZ;
    molRadius = std::sqrt(dx*dx + dy*dy + dz*dz) * 0.5f;
    if (molRadius < 1.0f) molRadius = 1.0f;
}

// ── Selection ────────────────────────────────────────────────────────────────

void MoleculeViewer::SetSelectionMode(SelectionMode mode)
{
    selectionMode = mode;
    Refresh(false);
}

void MoleculeViewer::ClearSelections()
{
    bondedSelected.clear();
    deletedSelected.clear();
    Refresh(false);
}

void MoleculeViewer::ClearBondedSelections()
{
    bondedSelected.clear();
    Refresh(false);
}

void MoleculeViewer::ClearDeletedSelections()
{
    deletedSelected.clear();
    Refresh(false);
}

std::vector<std::string> MoleculeViewer::GetBondedAtomNames() const
{
    std::vector<std::string> names;
    for (int idx : bondedSelected) {
        if (idx >= 0 && idx < (int)atoms.size())
            names.push_back(atoms[idx].name);
    }
    return names;
}

std::vector<std::string> MoleculeViewer::GetDeletedAtomNames() const
{
    std::vector<std::string> names;
    for (int idx : deletedSelected) {
        if (idx >= 0 && idx < (int)atoms.size())
            names.push_back(atoms[idx].name);
    }
    return names;
}

// ── Rendering ────────────────────────────────────────────────────────────────

void MoleculeViewer::OnPaint(wxPaintEvent& WXUNUSED(event))
{
    wxPaintDC dc(this);  // must create even if unused, to validate the window

    if (!glContext) return;
    if (!glInitialized) InitGL();
    SetCurrent(*glContext);

    Render();
    SwapBuffers();
}

void MoleculeViewer::OnSize(wxSizeEvent& event)
{
    if (!IsShownOnScreen()) return;
    if (!glContext) return;

    SetCurrent(*glContext);
    if (!glInitialized) InitGL();

    // Use the actual pixel size for GL viewport
    wxSize sz = GetClientSize();
    double scale = GetContentScaleFactor();
    glViewport(0, 0, (int)(sz.GetWidth() * scale), (int)(sz.GetHeight() * scale));

    Refresh(false);  // false = don't erase background, reduces flicker
}

void MoleculeViewer::Render()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (atoms.empty()) return;

    wxSize sz = GetClientSize();
    float aspect = (float)sz.GetWidth() / std::max((float)sz.GetHeight(), 1.0f);

    // Projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float viewDist = molRadius * 3.0f / zoom;
    gluPerspective(45.0, aspect, 0.1, viewDist * 4.0);

    // Modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(panX, panY, -viewDist);
    glRotatef(rotX, 1, 0, 0);
    glRotatef(rotY, 0, 1, 0);
    glTranslatef(-centerX, -centerY, -centerZ);

    // Draw bonds
    for (const auto& bond : bonds) {
        DrawBond(bond);
    }

    // Draw atoms
    for (size_t i = 0; i < atoms.size(); ++i) {
        bool isBondedSel = std::find(bondedSelected.begin(), bondedSelected.end(), (int)i)
                           != bondedSelected.end();
        bool isDeletedSel = std::find(deletedSelected.begin(), deletedSelected.end(), (int)i)
                            != deletedSelected.end();
        bool selected = isBondedSel || isDeletedSel;
        bool hovered  = ((int)i == hoveredAtom);

        int order = -1;
        if (isBondedSel) {
            for (size_t k = 0; k < bondedSelected.size(); ++k)
                if (bondedSelected[k] == (int)i) { order = (int)k + 1; break; }
        } else if (isDeletedSel) {
            for (size_t k = 0; k < deletedSelected.size(); ++k)
                if (deletedSelected[k] == (int)i) { order = (int)k + 1; break; }
        }

        DrawAtom(atoms[i], selected, hovered, order);
    }
}

void MoleculeViewer::DrawAtom(const Atom& atom, bool selected, bool hovered, int selectionOrder)
{
    float r, g, b;
    GetElementColor(atom.element, r, g, b);
    float radius = GetElementRadius(atom.element) * 0.35f;

    if (selected) {
        bool isBondedSel = false;
        for (int idx : bondedSelected)
            if (&atoms[idx] == &atom) { isBondedSel = true; break; }

        if (isBondedSel) {
            r = 0.1f; g = 1.0f; b = 0.2f;  // green
        } else {
            r = 1.0f; g = 0.2f; b = 0.1f;  // red
        }
        radius *= 1.3f;
    } else if (hovered) {
        r = 1.0f; g = 1.0f; b = 0.3f;  // bright yellow on hover
        radius *= 1.25f;
    }

    glColor3f(r, g, b);
    DrawSphere(atom.x, atom.y, atom.z, radius, 10, 8);  // reduced poly count
}

void MoleculeViewer::DrawBond(const Bond& bond)
{
    const Atom& a1 = atoms[bond.atomIndex1];
    const Atom& a2 = atoms[bond.atomIndex2];

    glColor3f(0.6f, 0.6f, 0.6f);
    DrawCylinder(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z, 0.08f, 6);  // reduced slices
}

// ── Primitive drawing ────────────────────────────────────────────────────────

void MoleculeViewer::DrawSphere(float cx, float cy, float cz, float radius,
                                int slices, int stacks)
{
    glPushMatrix();
    glTranslatef(cx, cy, cz);

    for (int i = 0; i < stacks; ++i) {
        float lat0 = (float)M_PI * (-0.5f + (float)i / stacks);
        float lat1 = (float)M_PI * (-0.5f + (float)(i + 1) / stacks);
        float y0 = std::sin(lat0), yr0 = std::cos(lat0);
        float y1 = std::sin(lat1), yr1 = std::cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= slices; ++j) {
            float lng = 2.0f * (float)M_PI * (float)j / slices;
            float x = std::cos(lng), z = std::sin(lng);

            glNormal3f(x * yr0, y0, z * yr0);
            glVertex3f(radius * x * yr0, radius * y0, radius * z * yr0);

            glNormal3f(x * yr1, y1, z * yr1);
            glVertex3f(radius * x * yr1, radius * y1, radius * z * yr1);
        }
        glEnd();
    }

    glPopMatrix();
}

void MoleculeViewer::DrawCylinder(float x1, float y1, float z1,
                                  float x2, float y2, float z2,
                                  float radius, int slices)
{
    float dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
    float length = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (length < 0.001f) return;

    float ax = dx / length, ay = dy / length, az = dz / length;

    float px, py, pz;
    if (std::fabs(ax) < 0.9f) {
        px = 0; py = -az; pz = ay;
    } else {
        px = -az; py = 0; pz = ax;
    }
    float pl = std::sqrt(px*px + py*py + pz*pz);
    px /= pl; py /= pl; pz /= pl;

    float qx = ay * pz - az * py;
    float qy = az * px - ax * pz;
    float qz = ax * py - ay * px;

    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i <= slices; ++i) {
        float angle = 2.0f * (float)M_PI * i / slices;
        float c = std::cos(angle), s = std::sin(angle);

        float nx = c * px + s * qx;
        float ny = c * py + s * qy;
        float nz = c * pz + s * qz;

        glNormal3f(nx, ny, nz);
        glVertex3f(x1 + radius * nx, y1 + radius * ny, z1 + radius * nz);
        glVertex3f(x2 + radius * nx, y2 + radius * ny, z2 + radius * nz);
    }
    glEnd();
}

// ── Mouse interaction ────────────────────────────────────────────────────────

void MoleculeViewer::OnMouse(wxMouseEvent& event)
{
    int mx = event.GetX();
    int my = event.GetY();

    // Scroll → zoom
    if (event.GetWheelRotation() != 0) {
        float delta = event.GetWheelRotation() > 0 ? 1.1f : 0.9f;
        zoom *= delta;
        zoom = std::max(0.1f, std::min(zoom, 20.0f));
        Refresh(false);
        return;
    }

    // Left click → pick atom or start drag
    if (event.LeftDown()) {
        int picked = PickAtom(mx, my);

        if (picked >= 0) {
            std::cout << "Picked atom: " << atoms[picked].name
                      << " (index " << picked << ")" << std::endl;

            std::vector<int>& sel = (selectionMode == SelectionMode::BONDED)
                                    ? bondedSelected : deletedSelected;

            auto it = std::find(sel.begin(), sel.end(), picked);
            if (it != sel.end()) {
                sel.erase(it);
            } else if (sel.size() < 3) {
                sel.push_back(picked);
            }

            if (selectionCallback) {
                auto names = (selectionMode == SelectionMode::BONDED)
                             ? GetBondedAtomNames() : GetDeletedAtomNames();
                std::cout << "Firing callback with " << names.size() << " names" << std::endl;
                selectionCallback(selectionMode, names);
                std::cout << "Callback returned" << std::endl;
            } else {
                std::cout << "WARNING: No callback set!" << std::endl;
            }

            Refresh(false);
        } else {
            dragging = true;
            lastMouseX = mx;
            lastMouseY = my;
        }
        return;
    }

    // Right button → pan
    if (event.RightDown()) {
        panning = true;
        lastMouseX = mx;
        lastMouseY = my;
        return;
    }

    // Dragging → rotate
    if (event.Dragging() && dragging) {
        float ddx = (float)(mx - lastMouseX);
        float ddy = (float)(my - lastMouseY);
        rotY += ddx * 0.5f;
        rotX += ddy * 0.5f;
        lastMouseX = mx;
        lastMouseY = my;
        Refresh(false);
        return;
    }

    // Panning
    if (event.Dragging() && panning) {
        float ddx = (float)(mx - lastMouseX);
        float ddy = (float)(my - lastMouseY);
        float scale = molRadius * 0.005f / zoom;
        panX += ddx * scale;
        panY -= ddy * scale;
        lastMouseX = mx;
        lastMouseY = my;
        Refresh(false);
        return;
    }

    // Release
    if (event.LeftUp()) { dragging = false; return; }
    if (event.RightUp()) { panning = false; return; }

    // All other events (motion, enter, leave, etc.) — skip them
    event.Skip();
}

// ── Atom picking ─────────────────────────────────────────────────────────────
//
// The challenge on Retina Macs: wxWidgets mouse coords are in "points",
// glViewport is in physical pixels, and gluProject uses the viewport.
// Rather than trying to convert between spaces, we override the viewport
// to match point space, so gluProject output matches mouse coords directly.

int MoleculeViewer::PickAtom(int screenX, int screenY)
{
    if (atoms.empty()) return -1;
    if (!glInitialized) return -1;

    SetCurrent(*glContext);

    double scaleFactor = GetContentScaleFactor();
    wxSize sz = GetClientSize();  // in logical points

    // Reconstruct the same matrices that Render() uses
    float aspect = (float)sz.GetWidth() / std::max((float)sz.GetHeight(), 1.0f);
    float viewDist = molRadius * 3.0f / zoom;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, aspect, 0.1, viewDist * 4.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(panX, panY, -viewDist);
    glRotatef(rotX, 1, 0, 0);
    glRotatef(rotY, 0, 1, 0);
    glTranslatef(-centerX, -centerY, -centerZ);

    GLdouble modelview[16], projection[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);

    // Use a FAKE viewport in logical points so gluProject output matches mouse coords
    GLint pointViewport[4] = { 0, 0, sz.GetWidth(), sz.GetHeight() };

    // Mouse Y flipped from top-left to bottom-left, in point space
    double glMouseY = (double)sz.GetHeight() - (double)screenY;

    // Debug output
    static int debugCount = 0;
    if (debugCount < 10) {
        std::cout << "PickAtom: mouse=(" << screenX << "," << screenY << ")"
                  << " pointVP=(" << pointViewport[2] << "x" << pointViewport[3] << ")"
                  << " glMouseY=" << glMouseY
                  << " scale=" << scaleFactor << "\n";

        if (!atoms.empty()) {
            GLdouble px, py, pz;
            gluProject(atoms[0].x, atoms[0].y, atoms[0].z,
                       modelview, projection, pointViewport, &px, &py, &pz);
            std::cout << "  atom[0] '" << atoms[0].name
                      << "' projects to (" << px << "," << py << ")\n";
        }
        debugCount++;
    }

    int bestIdx = -1;
    float bestDistSq = 99999.0f;
    float bestZ = 1e9f;
    float pickThresholdSq = 25.0f * 25.0f;  // 25 point radius

    for (size_t i = 0; i < atoms.size(); ++i) {
        GLdouble projX, projY, projZ;
        gluProject(atoms[i].x, atoms[i].y, atoms[i].z,
                   modelview, projection, pointViewport,
                   &projX, &projY, &projZ);

        float dx = (float)(projX - (double)screenX);
        float dy = (float)(projY - glMouseY);
        float distSq = dx*dx + dy*dy;

        if (distSq < pickThresholdSq) {
            if (distSq < bestDistSq ||
                (distSq < bestDistSq + 100.0f && projZ < bestZ)) {
                bestDistSq = distSq;
                bestZ = (float)projZ;
                bestIdx = (int)i;
            }
        }
    }

    return bestIdx;
}

// ── Element properties ───────────────────────────────────────────────────────

void MoleculeViewer::GetElementColor(const std::string& element, float& r, float& g, float& b)
{
    if (element == "C")       { r = 0.5f;  g = 0.5f;  b = 0.5f;  }
    else if (element == "N")  { r = 0.2f;  g = 0.3f;  b = 1.0f;  }
    else if (element == "O")  { r = 1.0f;  g = 0.2f;  b = 0.2f;  }
    else if (element == "S")  { r = 1.0f;  g = 0.85f; b = 0.2f;  }
    else if (element == "P")  { r = 1.0f;  g = 0.5f;  b = 0.0f;  }
    else if (element == "H")  { r = 0.9f;  g = 0.9f;  b = 0.9f;  }
    else if (element == "F")  { r = 0.0f;  g = 0.8f;  b = 0.0f;  }
    else if (element == "CL") { r = 0.0f;  g = 0.9f;  b = 0.0f;  }
    else if (element == "BR") { r = 0.6f;  g = 0.1f;  b = 0.1f;  }
    else if (element == "I")  { r = 0.4f;  g = 0.0f;  b = 0.7f;  }
    else if (element == "FE") { r = 0.8f;  g = 0.5f;  b = 0.0f;  }
    else                      { r = 0.7f;  g = 0.5f;  b = 0.7f;  }
}

float MoleculeViewer::GetElementRadius(const std::string& element)
{
    if (element == "H")  return 1.2f;
    if (element == "C")  return 1.7f;
    if (element == "N")  return 1.55f;
    if (element == "O")  return 1.52f;
    if (element == "S")  return 1.8f;
    if (element == "P")  return 1.8f;
    if (element == "F")  return 1.47f;
    if (element == "CL") return 1.75f;
    if (element == "BR") return 1.85f;
    if (element == "I")  return 1.98f;
    return 1.7f;
}