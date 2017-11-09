cat("Generate nlmixr installer")

nsi.stub <- "CRCCheck On
RequestExecutionLevel user
; Best Compression
SetCompress Auto
SetCompressor /SOLID lzma
SetCompressorDictSize 32
SetDatablockOptimize On
;SetCompress off
!include \"MUI2.nsh\"
Name \"<%=name%>\"
!define MUI_ICON \"<%=icon%>\"
OutFile \"<%=name%>-install.exe\"
InstallDir \"$LOCALAPPDATA\\nlmixr\"
InstallDirRegKey HKCU \"Software\\nlmixr\" \"\"
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
WriteRegStr HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr\" \"InstallPath\" \"$INSTDIR\\R\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"Current Version\" \"<%=rver%>nlmixr\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"Current Version\" \"<%=rver%>nlmixr\"
WriteRegStr HKCU \"Software\\R-core\\R\" \"InstallPath\" \"$INSTDIR\\R\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\" \"Current Version\" \"<%=rtoolsver%>\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\" \"InstallPath\" \"$INSTDIR\\rtools\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"FullVersion\" \"<%=fullr%>\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"InstallPath\" \"$INSTDIR\\rtools\"
WriteRegStr HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\" \"MinRVersion\" \"<%=minr%>\"

SetOutPath \"$INSTDIR\\python\"
File /r <%=python%>\\*
SetOutPath \"$INSTDIR\\rtools\"
File /r <%=rtools%>\\*
SetOutPath \"$INSTDIR\\R\"
File /r <%=R%>\\*

CreateDirectory \"$SMPROGRAMS\\nlmixr\"
CreateShortCut \"$SMPROGRAMS\\nlmixr\\nlmixr R.lnk\" \"$INSTDIR\\R\\bin\\x64\\Rgui.exe\"

;Store installation folder
WriteRegStr HKCU \"Software\\nlmixr\" \"\" $INSTDIR
;Create uninstaller
WriteUninstaller \"$INSTDIR\\Uninstall.exe\"
SectionEnd

Section \"Uninstall\"
RmDir /r \"$INSTDIR\\rtools\"
RmDir /r \"$INSTDIR\\python\"
RmDir /r \"$INSTDIR\\R\"
Delete \"$SMPROGRAMS\\nlmixr\\nlmixr R.lnk\"
RmDir /r \"$SMPROGRAMS\\nlmixr\"

Delete \"$INSTDIR\\Uninstall.exe\"
RMDir \"$INSTDIR\"

DeleteRegKey HKCU \"Software\\nlmixr\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\R\\<%=rver%>nlmixr\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\FullVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\InstallPath\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\\MinRVersion\"
DeleteRegKey HKCU \"Software\\R-core\\Rtools\\<%=rtoolsver%>\"
SectionEnd"

buildInstaller <- function(name="nlmixr"){
    rtools <- gsub("/", "\\", RxODE:::rxRtoolsBaseWin(), fixed=TRUE);
    python <- gsub("/", "\\", RxODE:::rxPythonBaseWin(), fixed=TRUE)
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
    sink(sprintf("%s.nsi", name));
    cat(nsis)
    sink()
    system(sprintf("makensis %s.nsi", name));
}

buildInstaller()
