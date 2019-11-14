cat("Generate nlmixr installer")
library(nlmixr)
nsi.lauch.stub <- "
CRCCheck On
RequestExecutionLevel user
; Best Compression
SetCompress Auto
SetCompressor /SOLID lzma
SetCompressorDictSize 32
SetDatablockOptimize On
;SetCompress off
Name \"<%=name%>\"
Icon \"<%=icon%>\"
OutFile \"<%=name%>.exe\"
AutoCloseWindow true
Caption \"Loading nlmixr\"
Subcaption 3 \" \"
ChangeUI all \"${NSISDIR}\\Contrib\\UIs\\LoadingBar_Icon.exe\"
XPStyle on
Section Main sec_main
WriteRegStr HKCU \"Software\\nlmixr<%=archext%>\" \"\" \"$EXEDIR\"
WriteRegStr HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\" \"InstallPath\" \"$EXEDIR\\R\"
Exec \"$EXEDIR\\R\\bin\\<%=Rdir%>\\Rgui.exe\"
SectionEnd"

shiny.lauch.stub <- "
CRCCheck On
RequestExecutionLevel user
; Best Compression
SetCompress Auto
SetCompressor /SOLID lzma
SetCompressorDictSize 32
SetDatablockOptimize On
;SetCompress off
Name \"<%=shiny.name%>\"
Icon \"<%=icon%>\"
OutFile \"<%=shiny.name%>.exe\"
AutoCloseWindow true
Caption \"Starting ShinyMixR\"
Subcaption 3 \" \"
ChangeUI all \"${NSISDIR}\\Contrib\\UIs\\LoadingBar_Icon.exe\"
XPStyle on
Section Main sec_main
WriteRegStr HKCU \"Software\\nlmixr<%=archext%>\" \"\" \"$EXEDIR\"
WriteRegStr HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\" \"InstallPath\" \"$EXEDIR\\R\"
Exec '$EXEDIR\\R\\bin\\<%=Rdir%>\\R.exe -e options(keep.source=TRUE);library(shinyMixR);nlmixr:::.setRoot();run_shinymixr(launch.browser=TRUE)'
SectionEnd"


update.lauch.stub <- "
CRCCheck On
RequestExecutionLevel user
; Best Compression
SetCompress Auto
SetCompressor /SOLID lzma
SetCompressorDictSize 32
SetDatablockOptimize On
!include \"MUI2.nsh\"
!include \"MUI_EXTRAPAGES.nsh\"
!include \"update.nsdinc\"
!define MUI_HEADERIMAGE_BITMAP \"nlmixr-header.bmp\"
!define MUI_HEADERIMAGE_BITMAP_NOSTRETCH
!define MUI_HEADERIMAGE_UNBITMAP \"nlmixr-header.bmp\" ; 150x57 pixels
!define MUI_PAGE_HEADER_TEXT \"nlmixr\"
!define MUI_PAGE_HEADER_SUBTEXT \"Nonlinear Mixed Effects Models in R\"
BrandingText \"nlmixr - Nonlinear Mixed Effects Models in R\"

;SetCompress off
Name \"Update\"
Icon \"Oxygen-Icons.org-Oxygen-Apps-system-software-update.ico\"
!define MUI_ICON \"Oxygen-Icons.org-Oxygen-Apps-system-software-update.ico\"
OutFile \"update.exe\"
AutoCloseWindow true
Caption \"Updating RxODE/nlmixr\"
Function fnc_update_Validate
  ${NSD_GetState} $hCtl_update_CheckBox1 $R0
  ${If} $R0 == ${BST_CHECKED}
    System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"binOpt\", \"true\").r0'
  ${Else}
    System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"binOpt\", \"false\").r0'
  ${EndIf}
  ${NSD_GetText} $hCtl_update_nlmixr $R0
  System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"nlmixrRef\", \"$R0\").r0'
  ${NSD_GetText} $hCtl_update_RxODE $R0
  System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"rxodeRef\", \"$R0\").r0'
  ${NSD_GetState} $hCtl_update_CheckBox2 $R0
  ${If} $R0 == ${BST_CHECKED}
    System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"useCRAN\", \"true\").r0'
  ${Else}
    System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"useCRAN\", \"false\").r0'
  ${EndIf}
  WriteRegStr HKCU \"Software\\nlmixr<%=archext%>\" \"\" \"$EXEDIR\"
  WriteRegStr HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\" \"InstallPath\" \"$EXEDIR\\R\"
  System::Call 'Kernel32::SetEnvironmentVariableA(t, t) i(\"HOME\", \"$TEMP\").r0'
  Exec '$EXEDIR\\R\\bin\\Rscript.exe \"$EXEDIR\\R\\update.R\"'
FunctionEnd
Page custom fnc_update_Show fnc_update_Validate
Section Main sec_main
SectionEnd
"

nsi.stub <- "
CRCCheck On
RequestExecutionLevel user
; Best Compression
SetCompress Auto
SetCompressor /SOLID lzma
SetCompressorDictSize 32
SetDatablockOptimize On
;SetCompress off
!include \"MUI2.nsh\"
!include \"MUI_EXTRAPAGES.nsh\"
## For some reason, this doesn't work... :(
## !include LogicLib.nsh
## !include StrContains.nsh
## Function .onInit
##     ${StrContains} $0 \"64\" $SysDir
##     StrCmp $0 \"\" bit32
##         Goto bit64
##     bit32:
##         MessageBox MB_OK \"This installer is for 64 bit systems only $SysDir\"
##         Abort
##     bit64:
##         Goto done
##     done:
##FunctionEnd

Name \"<%=name%>\"
!define MUI_ICON \"<%=icon%>\"
OutFile \"<%=name%>_<%=nlmixr.ver%>_<%=arch%>_install.exe\"
InstallDir \"c:\\R\\nlmixr<%=archext%>_<%=nlmixr.ver%>\"
InstallDirRegKey HKCU \"Software\\nlmixr<%=archext%>\\<%=nlmixr.ver%>\" \"\"
!define MUI_HEADERIMAGE

!define MUI_HEADERIMAGE_BITMAP \"<%=header%>\"
!define MUI_HEADERIMAGE_BITMAP_NOSTRETCH
!define MUI_HEADERIMAGE_UNBITMAP \"<%=header%>\" ; 150x57 pixels

!define MUI_WELCOMEFINISHPAGE_BITMAP \"<%=welcome%>\" ;164x314 pixels
!define MUI_UNWELCOMEFINISHPAGE_BITMAP \"<%=welcome%>\" ;164x314 pixels

!define MUI_ABORTWARNING
!define MUI_UNABORTWARNING

!define MUI_PAGE_HEADER_TEXT \"nlmixr\"
!define MUI_PAGE_HEADER_SUBTEXT \"Nonlinear Mixed Effects Models in R\"
BrandingText \"<%=name%> - Nonlinear Mixed Effects Models in R\"

!define MUI_COMPONENTSPAGE_SMALLDESC
!define MUI_HEADERIMAGE_RIGHT

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE \"<%=lic%>\"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_SHOWREADME $INSTDIR\\installation-notes.rtf
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH
!insertmacro MUI_LANGUAGE \"English\"
Section Main sec_main
WriteRegStr HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\" \"InstallPath\" \"$INSTDIR\\R\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"Current Version\" \"<%=rver%>nlmixr<%=archext%>\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"Current Version\" \"<%=rver%>nlmixr<%=archext%>\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"InstallPath\" \"$INSTDIR\\R\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\" \"Current Version\" \"<%=rtoolsver%>\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\" \"InstallPath\" \"$INSTDIR\\rtools\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"FullVersion\" \"<%=fullr%>\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"InstallPath\" \"$INSTDIR\\R\\Rtools\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"MinRVersion\" \"<%=minr%>\"
SetOutPath \"$INSTDIR\"
File \"nlmixr.exe\"
File \"shinyMixR.exe\"
File \"update.exe\"
File \"installation-notes.rtf\"
SetOutPath \"$INSTDIR\\R\"
File \"update.R\"
File /r <%=R%>\\*

##CreateDirectory \"c:\\R\\nlmixr<%=arch%>-<%=nlmixr.ver%>\"
<%=shortcuts%>

;Store installation folder
WriteRegStr HKCU \"Software\\nlmixr<%=archext%>\" \"\" $INSTDIR
;Create uninstaller
WriteUninstaller \"$INSTDIR\\Uninstall.exe\"
SectionEnd

Section \"Uninstall\"
RmDir /r \"$INSTDIR\\rtools\"
RmDir /r \"$INSTDIR\\python\"
RmDir /r \"$INSTDIR\\R\"
Delete \"$SMPROGRAMS\\nlmixr\\*.lnk\"
RmDir /r \"$SMPROGRAMS\\nlmixr\"

Delete \"$INSTDIR\\Uninstall.exe\"
Delete \"$INSTDIR\\Uninstall.exe\"
RmDir /r \"$INSTDIR\"
Delete \"$DESKTOP\\nlmixr R *.lnk\"
DeleteRegKey HKCU \"Software\\nlmixr<%=archext%>\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\FullVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\MinRVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\"
SectionEnd"

buildInstaller <- function(name="nlmixr"){
    rtools <- gsub("/", "\\", RxODE:::.rxRtoolsBaseWin(), fixed=TRUE);
    python <- gsub("/", "\\", RxODE:::.rxPythonBaseWin(), fixed=TRUE);
    R <- gsub("/", "\\", Sys.getenv("R_HOME"), fixed=TRUE);
    lic <- gsub("/", "\\", devtools::package_file("LICENSE"), fixed=TRUE);
    readme <- gsub("/", "\\", devtools::package_file("build/installation-notes.rtf"), fixed=TRUE);
    header <- gsub("/", "\\", devtools::package_file("build/nlmixr-header.bmp"), fixed=TRUE)
    welcome <- gsub("/", "\\", devtools::package_file("build/nlmixr-welcome.bmp"), fixed=TRUE)
    icon <- gsub("/", "\\", devtools::package_file("build/icon_red.ico"), fixed=TRUE)
    rver <- paste(R.version$major,R.version$minor,sep=".");
    ## rtools.curr <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HLM", view = "32-bit", maxdepth = 2);
    full.ver <- gsub("Rtools version ", "", readLines(file.path(RxODE:::.rxRtoolsBaseWin(), "VERSION.txt")))
    min.rver <- gsub("([0-9]+[.][0-9]+).*", "\\1", full.ver);
    rtools.curr <- min.rver;
    rtools.cur.ver <- rtools.curr
    arch <- R.version$arch;
    nlmixr.ver <- sessionInfo()$otherPkgs$nlmixr$Version;
    archext <- ifelse(.Platform$r_arch == "i386", "32", "")
    if (archext == "32"){
        shortcut <- sprintf("CreateShortCut \"$DESKTOP\\nlmixr R (32 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R %s (32 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"", nlmixr.ver);
        Rdir <- "i386"
    } else {
        shortcut <- sprintf("CreateShortCut \"$DESKTOP\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"", nlmixr.ver);
        Rdir <- "x64"
    }
    rtoolsver <- rtools.cur.ver;
    minr <- min.rver;
    fullr <- full.ver;
    shortcuts <- shortcut;
    dr <- gsub("/", "\\", devtools::package_file("build"), fixed=TRUE)
    dir <- dr;
    brew::brew(text=update.lauch.stub, output=file.path(dr,"update.nsi"));
    system(sprintf("makensis %s", file.path(dr, "update.nsi")));
    exe <- file.path(dr, "nlmixr.nsi");
    brew::brew(text=nsi.lauch.stub, output=file.path(dr, "nlmixr.nsi"));
    system(sprintf("makensis %s", file.path(dr, "nlmixr.nsi")));
    icon <- gsub("/", "\\", devtools::package_file("build/shinyMixR.ico"), fixed=TRUE)
    shiny.name <- "shinyMixR"
    brew::brew(text = shiny.lauch.stub, output=file.path(dr, "shinyMixR.nsi"))
    system(sprintf("makensis %s", file.path(dr, "shinyMixR.nsi")));
    icon <- gsub("/", "\\", devtools::package_file("build/icon_red.ico"), fixed=TRUE)
    ## unlink(file.path(dr, "nlmixr.nsi"))
    dr <- normalizePath(file.path(dr, sprintf("%s%s.nsi", name, archext)))
    brew::brew(text=nsi.stub, output=dr)
    system(sprintf("makensis %s", dr));
    ## unlink(dr)
    ## unlink(exe)
}

buildInstaller()
