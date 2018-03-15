# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=OSM - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to OSM - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OSM - Win32 Release" && "$(CFG)" != "OSM - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "OSM.mak" CFG="OSM - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OSM - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "OSM - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "OSM - Win32 Debug"
RSC=rc.exe
F90=fl32.exe

!IF  "$(CFG)" == "OSM - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\OSM.exe"

CLEAN : 
	-@erase ".\Release\OSM.exe"
	-@erase ".\Release\ReflectION.obj"
	-@erase ".\Release\Radial.obj"
	-@erase ".\Release\param.obj"
	-@erase ".\Release\Pedestal.obj"
	-@erase ".\Release\Qdiv.obj"
	-@erase ".\Release\Denlim.obj"
	-@erase ".\Release\Distr.obj"
	-@erase ".\Release\divert3.obj"
	-@erase ".\Release\Atomic.obj"
	-@erase ".\Release\ATCOOL2.OBJ"
	-@erase ".\Release\edgecalc.obj"
	-@erase ".\Release\Limit2.obj"
	-@erase ".\Release\Editdiv.obj"
	-@erase ".\Release\Expint.obj"
	-@erase ".\Release\ncefits.obj"
	-@erase ".\Release\Neutdist.obj"
	-@erase ".\Release\Divstab.obj"
	-@erase ".\Release\Svion.obj"
	-@erase ".\Release\transm.obj"
	-@erase ".\Release\sgesv_jm.obj"
	-@erase ".\Release\Neutfrac.obj"
	-@erase ".\Release\poloidal.obj"
	-@erase ".\Release\Soldiv.obj"
	-@erase ".\Release\Dsigmav.obj"
	-@erase ".\Release\Divrad.obj"
	-@erase ".\Release\DivSol.obj"
	-@erase ".\Release\Interp.obj"
	-@erase ".\Release\Neutxpt.obj"
	-@erase ".\Release\rotate.obj"
	-@erase ".\Release\cxrcefits.obj"
	-@erase ".\Release\Radintg2.obj"
	-@erase ".\Release\Geometry.obj"
	-@erase ".\Release\soldata.obj"
	-@erase ".\Release\Disrupt.obj"
	-@erase ".\Release\Molecule.obj"
	-@erase ".\Release\TOROTATE.OBJ"
	-@erase ".\Release\CEFITS.OBJ"
	-@erase ".\Release\EDGEROTRAN.OBJ"
	-@erase ".\Release\width.obj"
	-@erase ".\Release\Besj0.obj"
	-@erase ".\Release\Corad.obj"
	-@erase ".\Release\SIGMAV.OBJ"
	-@erase ".\Release\Neutral3.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OSM.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/OSM.pdb" /machine:I386 /out:"$(OUTDIR)/OSM.exe" 
LINK32_OBJS= \
	"$(INTDIR)/ReflectION.obj" \
	"$(INTDIR)/Radial.obj" \
	"$(INTDIR)/param.obj" \
	"$(INTDIR)/Pedestal.obj" \
	"$(INTDIR)/Qdiv.obj" \
	"$(INTDIR)/Denlim.obj" \
	"$(INTDIR)/Distr.obj" \
	"$(INTDIR)/divert3.obj" \
	"$(INTDIR)/Atomic.obj" \
	"$(INTDIR)/ATCOOL2.OBJ" \
	"$(INTDIR)/edgecalc.obj" \
	"$(INTDIR)/Limit2.obj" \
	"$(INTDIR)/Editdiv.obj" \
	"$(INTDIR)/Expint.obj" \
	"$(INTDIR)/ncefits.obj" \
	"$(INTDIR)/Neutdist.obj" \
	"$(INTDIR)/Divstab.obj" \
	"$(INTDIR)/Svion.obj" \
	"$(INTDIR)/transm.obj" \
	"$(INTDIR)/sgesv_jm.obj" \
	"$(INTDIR)/Neutfrac.obj" \
	"$(INTDIR)/poloidal.obj" \
	"$(INTDIR)/Soldiv.obj" \
	"$(INTDIR)/Dsigmav.obj" \
	"$(INTDIR)/Divrad.obj" \
	"$(INTDIR)/DivSol.obj" \
	"$(INTDIR)/Interp.obj" \
	"$(INTDIR)/Neutxpt.obj" \
	"$(INTDIR)/rotate.obj" \
	"$(INTDIR)/cxrcefits.obj" \
	"$(INTDIR)/Radintg2.obj" \
	"$(INTDIR)/Geometry.obj" \
	"$(INTDIR)/soldata.obj" \
	"$(INTDIR)/Disrupt.obj" \
	"$(INTDIR)/Molecule.obj" \
	"$(INTDIR)/TOROTATE.OBJ" \
	"$(INTDIR)/CEFITS.OBJ" \
	"$(INTDIR)/EDGEROTRAN.OBJ" \
	"$(INTDIR)/width.obj" \
	"$(INTDIR)/Besj0.obj" \
	"$(INTDIR)/Corad.obj" \
	"$(INTDIR)/SIGMAV.OBJ" \
	"$(INTDIR)/Neutral3.obj"

"$(OUTDIR)\OSM.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OSM - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\OSM.exe"

CLEAN : 
	-@erase ".\Debug\OSM.exe"
	-@erase ".\Debug\Neutral3.obj"
	-@erase ".\Debug\Neutdist.obj"
	-@erase ".\Debug\Limit2.obj"
	-@erase ".\Debug\Qdiv.obj"
	-@erase ".\Debug\Expint.obj"
	-@erase ".\Debug\SIGMAV.OBJ"
	-@erase ".\Debug\Pedestal.obj"
	-@erase ".\Debug\poloidal.obj"
	-@erase ".\Debug\EDGEROTRAN.OBJ"
	-@erase ".\Debug\Denlim.obj"
	-@erase ".\Debug\ATCOOL2.OBJ"
	-@erase ".\Debug\Soldiv.obj"
	-@erase ".\Debug\width.obj"
	-@erase ".\Debug\DivSol.obj"
	-@erase ".\Debug\edgecalc.obj"
	-@erase ".\Debug\Svion.obj"
	-@erase ".\Debug\ReflectION.obj"
	-@erase ".\Debug\Atomic.obj"
	-@erase ".\Debug\Editdiv.obj"
	-@erase ".\Debug\ncefits.obj"
	-@erase ".\Debug\Divstab.obj"
	-@erase ".\Debug\cxrcefits.obj"
	-@erase ".\Debug\Distr.obj"
	-@erase ".\Debug\sgesv_jm.obj"
	-@erase ".\Debug\Neutfrac.obj"
	-@erase ".\Debug\Dsigmav.obj"
	-@erase ".\Debug\transm.obj"
	-@erase ".\Debug\Neutxpt.obj"
	-@erase ".\Debug\Besj0.obj"
	-@erase ".\Debug\Divrad.obj"
	-@erase ".\Debug\Radial.obj"
	-@erase ".\Debug\divert3.obj"
	-@erase ".\Debug\Interp.obj"
	-@erase ".\Debug\Radintg2.obj"
	-@erase ".\Debug\soldata.obj"
	-@erase ".\Debug\Corad.obj"
	-@erase ".\Debug\Geometry.obj"
	-@erase ".\Debug\Disrupt.obj"
	-@erase ".\Debug\rotate.obj"
	-@erase ".\Debug\Molecule.obj"
	-@erase ".\Debug\TOROTATE.OBJ"
	-@erase ".\Debug\param.obj"
	-@erase ".\Debug\CEFITS.OBJ"
	-@erase ".\Debug\OSM.ilk"
	-@erase ".\Debug\OSM.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/OSM.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OSM.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/OSM.pdb" /debug /machine:I386 /out:"$(OUTDIR)/OSM.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Neutral3.obj" \
	"$(INTDIR)/Neutdist.obj" \
	"$(INTDIR)/Limit2.obj" \
	"$(INTDIR)/Qdiv.obj" \
	"$(INTDIR)/Expint.obj" \
	"$(INTDIR)/SIGMAV.OBJ" \
	"$(INTDIR)/Pedestal.obj" \
	"$(INTDIR)/poloidal.obj" \
	"$(INTDIR)/EDGEROTRAN.OBJ" \
	"$(INTDIR)/Denlim.obj" \
	"$(INTDIR)/ATCOOL2.OBJ" \
	"$(INTDIR)/Soldiv.obj" \
	"$(INTDIR)/width.obj" \
	"$(INTDIR)/DivSol.obj" \
	"$(INTDIR)/edgecalc.obj" \
	"$(INTDIR)/Svion.obj" \
	"$(INTDIR)/ReflectION.obj" \
	"$(INTDIR)/Atomic.obj" \
	"$(INTDIR)/Editdiv.obj" \
	"$(INTDIR)/ncefits.obj" \
	"$(INTDIR)/Divstab.obj" \
	"$(INTDIR)/cxrcefits.obj" \
	"$(INTDIR)/Distr.obj" \
	"$(INTDIR)/sgesv_jm.obj" \
	"$(INTDIR)/Neutfrac.obj" \
	"$(INTDIR)/Dsigmav.obj" \
	"$(INTDIR)/transm.obj" \
	"$(INTDIR)/Neutxpt.obj" \
	"$(INTDIR)/Besj0.obj" \
	"$(INTDIR)/Divrad.obj" \
	"$(INTDIR)/Radial.obj" \
	"$(INTDIR)/divert3.obj" \
	"$(INTDIR)/Interp.obj" \
	"$(INTDIR)/Radintg2.obj" \
	"$(INTDIR)/soldata.obj" \
	"$(INTDIR)/Corad.obj" \
	"$(INTDIR)/Geometry.obj" \
	"$(INTDIR)/Disrupt.obj" \
	"$(INTDIR)/rotate.obj" \
	"$(INTDIR)/Molecule.obj" \
	"$(INTDIR)/TOROTATE.OBJ" \
	"$(INTDIR)/param.obj" \
	"$(INTDIR)/CEFITS.OBJ"

"$(OUTDIR)\OSM.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "OSM - Win32 Release"
# Name "OSM - Win32 Debug"

!IF  "$(CFG)" == "OSM - Win32 Release"

!ELSEIF  "$(CFG)" == "OSM - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\width.for
DEP_F90_WIDTH=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\width.obj" : $(SOURCE) $(DEP_F90_WIDTH) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\transm.for
DEP_F90_TRANS=\
	".\geometry.fi"\
	

"$(INTDIR)\transm.obj" : $(SOURCE) $(DEP_F90_TRANS) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\TOROTATE.FOR
DEP_F90_TOROT=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\TOROTATE.OBJ" : $(SOURCE) $(DEP_F90_TOROT) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Svion.for
DEP_F90_SVION=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Svion.obj" : $(SOURCE) $(DEP_F90_SVION) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Soldiv.for
DEP_F90_SOLDI=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Soldiv.obj" : $(SOURCE) $(DEP_F90_SOLDI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\soldata.for
DEP_F90_SOLDA=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\soldata.obj" : $(SOURCE) $(DEP_F90_SOLDA) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sgesv_jm.for

"$(INTDIR)\sgesv_jm.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\rotate.for
DEP_F90_ROTAT=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\rotate.obj" : $(SOURCE) $(DEP_F90_ROTAT) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ReflectION.for
DEP_F90_REFLE=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\ReflectION.obj" : $(SOURCE) $(DEP_F90_REFLE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Radintg2.for
DEP_F90_RADIN=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Radintg2.obj" : $(SOURCE) $(DEP_F90_RADIN) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Radial.for
DEP_F90_RADIA=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Radial.obj" : $(SOURCE) $(DEP_F90_RADIA) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Qdiv.for
DEP_F90_QDIV_=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Qdiv.obj" : $(SOURCE) $(DEP_F90_QDIV_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\poloidal.for
DEP_F90_POLOI=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\poloidal.obj" : $(SOURCE) $(DEP_F90_POLOI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Pedestal.for
DEP_F90_PEDES=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Pedestal.obj" : $(SOURCE) $(DEP_F90_PEDES) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\param.for
DEP_F90_PARAM=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\param.obj" : $(SOURCE) $(DEP_F90_PARAM) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Neutxpt.FOR
DEP_F90_NEUTX=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Neutxpt.obj" : $(SOURCE) $(DEP_F90_NEUTX) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Neutral3.for
DEP_F90_NEUTR=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Neutral3.obj" : $(SOURCE) $(DEP_F90_NEUTR) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Neutfrac.for
DEP_F90_NEUTF=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Neutfrac.obj" : $(SOURCE) $(DEP_F90_NEUTF) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Neutdist.for
DEP_F90_NEUTD=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Neutdist.obj" : $(SOURCE) $(DEP_F90_NEUTD) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ncefits.for

"$(INTDIR)\ncefits.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Molecule.for
DEP_F90_MOLEC=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Molecule.obj" : $(SOURCE) $(DEP_F90_MOLEC) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Limit2.for
DEP_F90_LIMIT=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Limit2.obj" : $(SOURCE) $(DEP_F90_LIMIT) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Interp.for
DEP_F90_INTER=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Interp.obj" : $(SOURCE) $(DEP_F90_INTER) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Geometry.for
DEP_F90_GEOME=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Geometry.obj" : $(SOURCE) $(DEP_F90_GEOME) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Expint.for

"$(INTDIR)\Expint.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Editdiv.for
DEP_F90_EDITD=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Editdiv.obj" : $(SOURCE) $(DEP_F90_EDITD) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\EDGEROTRAN.FOR
DEP_F90_EDGER=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\EDGEROTRAN.OBJ" : $(SOURCE) $(DEP_F90_EDGER) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\edgecalc.for
DEP_F90_EDGEC=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\edgecalc.obj" : $(SOURCE) $(DEP_F90_EDGEC) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Dsigmav.f

"$(INTDIR)\Dsigmav.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Divstab.for
DEP_F90_DIVST=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Divstab.obj" : $(SOURCE) $(DEP_F90_DIVST) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Divrad.for
DEP_F90_DIVRA=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Divrad.obj" : $(SOURCE) $(DEP_F90_DIVRA) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\divert3.for
DEP_F90_DIVER=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\divert3.obj" : $(SOURCE) $(DEP_F90_DIVER) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Distr.for
DEP_F90_DISTR=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Distr.obj" : $(SOURCE) $(DEP_F90_DISTR) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Disrupt.for
DEP_F90_DISRU=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Disrupt.obj" : $(SOURCE) $(DEP_F90_DISRU) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Denlim.for
DEP_F90_DENLI=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Denlim.obj" : $(SOURCE) $(DEP_F90_DENLI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\cxrcefits.for

"$(INTDIR)\cxrcefits.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Corad.for
DEP_F90_CORAD=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Corad.obj" : $(SOURCE) $(DEP_F90_CORAD) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\CEFITS.FOR

"$(INTDIR)\CEFITS.OBJ" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Besj0.for

"$(INTDIR)\Besj0.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Atomic.for
DEP_F90_ATOMI=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\Atomic.obj" : $(SOURCE) $(DEP_F90_ATOMI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\DivSol.for
DEP_F90_DIVSO=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\DivSol.obj" : $(SOURCE) $(DEP_F90_DIVSO) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ATCOOL2.FOR
DEP_F90_ATCOO=\
	".\SOLDIV.FI"\
	

"$(INTDIR)\ATCOOL2.OBJ" : $(SOURCE) $(DEP_F90_ATCOO) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\SIGMAV.F

"$(INTDIR)\SIGMAV.OBJ" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
