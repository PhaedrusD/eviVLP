# eviVLP
Programs designed to prepare and modify pdb files of virus-like particles (VLP) for molecular dynamics simulations
VLPs are used in immunology and imaging and consist of many repeated proteins forming a large structure. The proteins are called coat proteins. The modifications can be genetic, which will add another protein to some or all coat proteins, or chemical, which add small drug like molecules or imaging agents to amino acids.
The GUI is written in C++, but also calls python and a few other programs. 
For functionality, many of the C++ codes must be compiled before running the GUI. 
To automate the outputing of a mol2 file, obabel must be installed.
