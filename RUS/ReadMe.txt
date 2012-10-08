========================================================================
    CONSOLE APPLICATION : RUS Project Overview
========================================================================

AppWizard has created this RUS application for you.

This file contains a summary of what you will find in each of the files that
make up your RUS application.


RUS.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

RUS.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

RUS.cpp
    This is the main application source file.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named RUS.pch and a precompiled types file named StdAfx.obj.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////

Resonant ultrasound calculations. 

Timing: 
with a population of 40, and taking advantage of the blocks (and with chrome and foobar open...), basis size of 10 takes 35ms per eigenvalue-sovling step, and 3472ms per generation per thread.

Old way, same parameters, but calculating without the block-advantage, takes 13440ms per generation per thread. 

Same as first (40, blocks, basis 10, etc..) but NOW i calculate potential energy only within the blocks. BAM. 1443ms per gen per thread.

calcualte the gradient integrals only once at the start... now it's 1032ms per gen per thread