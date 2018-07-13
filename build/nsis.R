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
InstallDir \"$LOCALAPPDATA\\nlmixr<%=archext%>\"
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
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"InstallPath\" \"$INSTDIR\\rtools\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"MinRVersion\" \"<%=minr%>\"
SetOutPath \"$INSTDIR\"
File \"nlmixr.exe\"
SetOutPath \"$INSTDIR\\python\"
File /r <%=python%>\\*
SetOutPath \"$INSTDIR\\rtools\"
File /r <%=rtools%>\\*
SetOutPath \"$INSTDIR\\R\"
File /r <%=R%>\\*

CreateDirectory \"$SMPROGRAMS\\nlmixr\"
<%=shortcuts%>

;Store installation folder
WriteRegStr HKCU \"Software\\nlmixr<%=archext%>\" \"\" $INSTDIR
;Create uninstaller
WriteUninstaller \"$INSTDIR\\etc\\Uninstall.exe\"
SectionEnd

Section \"Uninstall\"
RmDir /r \"$INSTDIR\\rtools\"
RmDir /r \"$INSTDIR\\python\"
RmDir /r \"$INSTDIR\\R\"
Delete \"$SMPROGRAMS\\nlmixr\\*.lnk\"
RmDir /r \"$SMPROGRAMS\\nlmixr\"

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
    header <- gsub("/", "\\", devtools::package_file("build/nlmixr-header.bmp"), fixed=TRUE)
    welcome <- gsub("/", "\\", devtools::package_file("build/nlmixr-welcome.bmp"), fixed=TRUE)
    icon <- gsub("/", "\\", devtools::package_file("build/icon_red.ico"), fixed=TRUE)
    rver <- paste(R.version$major,R.version$minor,sep=".");
    rtools.curr <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HLM", view = "32-bit", maxdepth = 2);
    rtools.cur.ver <- rtools.curr$`Current Version`;
    full.ver <- rtools.curr[[rtools.cur.ver]][["FullVersion"]];
    min.rver <- rtools.curr[[rtools.cur.ver]][["MinRVersion"]];
    arch <- R.version$arch;
    nlmixr.ver <- sessionInfo()$otherPkgs$nlmixr$Version;
    archext <- ifelse(.Platform$r_arch == "i386", "32", "")
    if (archext == "32"){
        shortcut <- "CreateShortCut \"$DESKTOP\\nlmixr R (32 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R (32 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"";
        Rdir <- "i386"
    } else {
        shortcut <- "CreateShortCut \"$DESKTOP\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\nlmixr.exe\"";
        Rdir <- "x64"
    }
    rtoolsver <- rtools.cur.ver;
    minr <- min.rver;
    fullr <- full.ver;
    shortcuts <- shortcut;
    dr <- gsub("/", "\\", devtools::package_file("build"), fixed=TRUE)
    dir <- dr;
    exe <- file.path(dr, "nlmixr.nsi");
    brew::brew(text=nsi.lauch.stub, output=file.path(dr, "nlmixr.nsi"));
    system(sprintf("makensis %s", file.path(dr, "nlmixr.nsi")));
    ## unlink(file.path(dr, "nlmixr.nsi"))
    dr <- normalizePath(file.path(dr, sprintf("%s%s.nsi", name, archext)))
    brew::brew(text=nsi.stub, output=dr)
    system(sprintf("makensis %s", dr));
    ## unlink(dr)
    ## unlink(exe)
}

buildInstaller()
