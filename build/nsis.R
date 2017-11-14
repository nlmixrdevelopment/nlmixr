cat("Generate nlmixr installer")

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
File \"<%=icon%>\"
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
WriteUninstaller \"$INSTDIR\\Uninstall.exe\"
SectionEnd

Section \"Uninstall\"
RmDir /r \"$INSTDIR\\rtools\"
RmDir /r \"$INSTDIR\\python\"
RmDir /r \"$INSTDIR\\R\"
Delete \"$SMPROGRAMS\\nlmixr\\*.lnk\"
RmDir /r \"$SMPROGRAMS\\nlmixr\"

Delete \"$INSTDIR\\Uninstall.exe\"
RmDir /r \"$INSTDIR\"

DeleteRegKey HKCU \"Software\\nlmixr<%=archext%>\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr<%=archext%>\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\FullVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\MinRVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\"
SectionEnd"

buildInstaller <- function(name="nlmixr"){
    rtools <- gsub("/", "\\", RxODE:::rxRtoolsBaseWin(), fixed=TRUE);
    python <- gsub("/", "\\", RxODE:::rxPythonBaseWin(), fixed=TRUE);
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
        shortcut <- "CreateShortCut \"$DESKTOP\\nlmixr R (32 bit).lnk\" \"$INSTDIR\\R\\bin\\i386\\Rgui.exe\" \"\" \"$INSTDIR\\icon_red.ico\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R (32 bit).lnk\" \"$INSTDIR\\R\\bin\\i386\\Rgui.exe\" \"\" \"$INSTDIR\\icon_red.ico\"";
    } else {
        shortcut <- "CreateShortCut \"$DESKTOP\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\R\\bin\\x64\\Rgui.exe\" \"\" \"$INSTDIR\\icon_red.ico\"\nCreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R (64 bit).lnk\" \"$INSTDIR\\R\\bin\\x64\\Rgui.exe\" \"\" \"$INSTDIR\\icon_red.ico\""
    }
    nsis <- gsub("<%=icon%>",icon,
                 gsub("<%=welcome%>", welcome,
                      gsub("<%=header%>",header,
                           gsub("<%=name%>", name,
                                gsub("<%=lic%>", lic,
                                     gsub("<%=R%>",R,
                                          gsub("<%=python%>", python,
                                               gsub("<%=rtools%>", rtools, nsi.stub, fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE)
    nsis <- gsub("<%=rver%>", rver, nsis, fixed=TRUE)
    nsis <- gsub("<%=rtoolsver%>", rtools.cur.ver, nsis, fixed=TRUE)
    nsis <- gsub("<%=minr%>", min.rver, nsis, fixed=TRUE)
    nsis <- gsub("<%=fullr%>", full.ver, nsis, fixed=TRUE)
    nsis <- gsub("<%=arch%>", arch, nsis, fixed=TRUE)
    nsis <- gsub("<%=nlmixr.ver%>", nlmixr.ver, nsis, fixed=TRUE);
    nsis <- gsub("<%=shortcuts%>", shortcut, nsis, fixed=TRUE);
    nsis <- gsub("<%=archext%>", archext, nsis, fixed=TRUE);
    dr <- gsub("/", "\\", devtools::package_file("build"), fixed=TRUE)
    nsis <- gsub("<%=dir%>", dr, nsis, fixed=TRUE);
    sink(file.path(devtools::package_file("build"), sprintf("%s%s.nsi", name, archex)));
    cat(nsis)
    sink()
    nsis <- sprintf("%s\\%s%s.nsi", dr, name, archex);
    system(sprintf("makensis %s", nsis));
    unlink(nsis)
}

buildInstaller()
