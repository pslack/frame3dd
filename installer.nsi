; NSIS script to create a binary installer for the Windows platform.
; by John Pye, 2007.
;
; Originally based on example2.nsi from the NSIS distribution.

;--------------------------------

; The name of the installer
Name "FRAME3DD"

!include LogicLib.nsh
!include nsis\registerExtension.nsh

; The file to write
!ifdef OUTFILE
OutFile "${OUTFILE}"
!else
OutFile frame3dd-setup.exe
!endif


;SetCompressor /FINAL zlib
SetCompressor /SOLID lzma

; The default installation directory
InstallDir $PROGRAMFILES\FRAME3DD

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\FRAME3DD" "Install_Dir"

;--------------------------------

; Pages

Page components
Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------

; The stuff to install
Section "FRAME3DD (required)"
	SectionIn RO

	DetailPrint "--- FRAME3DD application ---"

	; Set output path to the installation directory.
	SetOutPath $INSTDIR
	File "build\frame3dd.exe"
	File "README.txt"
	File "README-win32.txt"
	File "TODO.txt"
	File "LICENSE.txt"
	File "ChangeLog.txt"
	
	SetOutPath "$INSTDIR\examples"
	File "examples\exA.frm"
	File "examples\exB.frm"
	File "examples\exC.frm"
	File "examples\exD.frm"
	File "examples\exE.frm"
	File "examples\exF.frm"
	File "examples\exG.frm"
	File "examples\exH.frm"
	File "examples\exI.frm"
	
	SetOutPath "$INSTDIR\doc"
	File "doc\user-manual.html"
	File "doc\version.html"
	SetOutPath "$INSTDIR\doc\img"
	File "doc\img\tors.jpg"
	File "doc\img\psoe.jpg"
	File "doc\img\CSV1.png"
	File "doc\img\CSV2.png"
	File "doc\img\CSV3.png"
	File "doc\img\CSV4.png"
		
	; Write the installation path into the registry
	WriteRegStr HKLM SOFTWARE\FRAME3DD "Install_Dir" "$INSTDIR"

	; Write the uninstall keys for Windows
	WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRAME3DD" "DisplayName" "FRAME3DD"
	WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRAME3DD" "UninstallString" '"$INSTDIR\uninstall.exe"'
	WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRAME3DD" "NoModify" 1
	WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRAME3DD" "NoRepair" 1
	WriteUninstaller "uninstall.exe"

	
SectionEnd

Section "Microstran Viewer"
	DetailPrint "--- Microstran Viewer ---"
	SetOutPath $INSTDIR
	File "build\microstranparser.dll"
	File "build\arc2iv.exe"
	File "build\forcebalance.exe"
	File "src\microstran\properties.txt"
	
	${registerExtension} "$INSTDIR\arc2iv.exe" ".arc" "Microstran Archive File"
	
	; Record the fact that we've got the Microstran components installed
	WriteRegDWORD HKLM SOFTWARE\FRAME3DD "MSTRANP_INSTALLED" 1

SectionEnd



;---------------------------------

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"
  
  WriteRegDWORD HKLM "SOFTWARE\FRAME3DD" "StartMenu" 1

  CreateDirectory "$SMPROGRAMS\FRAME3DD"
  CreateShortCut "$SMPROGRAMS\FRAME3DD\FRAME3DD README.lnk" "$INSTDIR\README.txt" "" "$INSTDIR\README.txt" 0
  CreateShortCut "$SMPROGRAMS\FRAME3DD\README for Windows users.lnk" "$INSTDIR\README-win32.txt" "" "$INSTDIR\README-win32.txt" 0
  CreateShortCut "$SMPROGRAMS\FRAME3DD\FRAME3DD Documentation.lnk" "$INSTDIR\doc\user-manual.html" "" "$INSTDIR\doc\user-manual.html" 0
  CreateShortCut "$SMPROGRAMS\FRAME3DD\Examples.lnk" "$INSTDIR\examples" "" "$INSTDIR\examples" 0
  CreateShortCut "$SMPROGRAMS\FRAME3DD\ChangeLog.lnk" "$INSTDIR\ChangeLog.txt" "" "$INSTDIR\ChangeLog.txt" 0
  CreateShortCut "$SMPROGRAMS\FRAME3DD\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  
SectionEnd

;------------------------------------------------------------------
; UNINSTALLER

Section "Uninstall"
  
;--- start menu ---

	ReadRegDWORD $1 HKLM "SOFTWARE\FRAME3DD" "StartMenu"
	${If} $1 == 1
		; Remove shortcuts, if any
		DetailPrint "--- REMOVING START MENU SHORTCUTS ---"
		RmDir /r "$SMPROGRAMS\FRAME3DD"
	${EndIf}

;--- MSTRANP components ---

	ReadRegDWORD $1 HKLM "SOFTWARE\FRAME3DD" "MSTRANP_INSTALLED"	
	${If} $1 == 1
		; Remove Microstran parser components
		DetailPrint "--- REMOVING MICROSTRAN PARSER COMPONENTS ---"
		${unregisterExtension} ".arc" "Microstran Archive File"
		Delete "$INSTDIR\microstranparser.dll"
		Delete "$INSTDIR\arc2iv.exe"
		Delete "$INSTDIR\forcebalance.exe"
		Delete "$INSTDIR\properties.txt"

	${EndIf}

;--- common components ---

	DetailPrint "--- REMOVING FRAME3DD application ---"
	; Remove registry keys

	DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRAME3DD"
	DeleteRegKey HKLM "SOFTWARE\FRAME3DD"

	; Remove files and uninstaller

	Delete "$INSTDIR\frame3dd.exe"
	Delete "$INSTDIR\microstranparser.dll"
	Delete "$INSTDIR\LICENSE.txt"
	Delete "$INSTDIR\TODO.txt"
	Delete "$INSTDIR\README.txt"
	Delete "$INSTDIR\README-win32.txt"
	Delete "$INSTDIR\ChangeLog.txt"

	Delete "$INSTDIR\examples\exA.frm"
	Delete "$INSTDIR\examples\exB.frm"
	Delete "$INSTDIR\examples\exC.frm"
	Delete "$INSTDIR\examples\exD.frm"
	Delete "$INSTDIR\examples\exE.frm"
	Delete "$INSTDIR\examples\exF.frm"
	Delete "$INSTDIR\examples\exG.frm"
	Delete "$INSTDIR\examples\exH.frm"
	Delete "$INSTDIR\examples\exI.frm"	
	RMDir "$INSTDIR\examples"
	
	Delete "$INSTDIR\doc\user-manual.html"
	Delete "$INSTDIR\doc\img\tors.jpg"
	Delete "$INSTDIR\doc\img\psoe.jpg"
	Delete "$INSTDIR\doc\img\CSV1.png"
	Delete "$INSTDIR\doc\img\CSV2.png"
	Delete "$INSTDIR\doc\img\CSV3.png"
	Delete "$INSTDIR\doc\img\CSV4.png"
	RmDir "$INSTDIR\doc\img"
	RmDir "$INSTDIR\doc"
	
	Delete "$INSTDIR\uninstall.exe"
	RMDir "$INSTDIR"

SectionEnd
