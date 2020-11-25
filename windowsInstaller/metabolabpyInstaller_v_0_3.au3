#AutoIt3Wrapper_icon="C:\metabolabpy\metabolabpy\icon\icon.ico"

#include <ButtonConstants.au3>
#include <GUIConstantsEx.au3>
#include <MsgBoxConstants.au3>
#include <WinAPIFiles.au3>
#include <InetConstants.au3>
#include <WindowsConstants.au3>


ConsoleWrite($cmdLine[0] & @CRLF)
$gui1 = GUICreate("MetaboLabPy Installer v0.3", 500, 300, 100, 100) ; ==> draw GUI
$nextButton = GUICtrlCreateButton("Next",440, 270, 50, 20)
$cancelButton = GUICtrlCreateButton("Cancel", 380, 270, 50, 20)
$installMiniconda = GUICtrlCreateCheckbox("Install Miniconda3", 10, 10,250, 20)
GUICtrlSetState($installMiniconda, $GUI_CHECKED)
GUICtrlCreateLabel("Miniconda3 path", 10, 40, 120, 20)
$minicondaPath = GUICtrlCreateEdit("", 140, 40, 270, 20)
$setMinicondaPath = GUICtrlCreateButton("Select", 420, 40, 70, 20)
$makeMinicondaStandard = GUICtrlCreateCheckbox("Make Miniconda3 your standard python distribution", 40, 100, 400, 20)
GUICtrlSetState($makeMinicondaStandard, $GUI_CHECKED)
$addMinicondaToPath = GUICtrlCreateCheckbox("Add Miniconda3 to your PATH variable", 40, 130, 400, 20)
GUICtrlSetState($addMinicondaToPath, $GUI_CHECKED)
GUISetState() ; show GUI


while 1
   $msg = GUIGetMsg() ; check for GUI activity
   Select
	  case $msg = -3 ; -3 = close by clicking x (top right)
		 Exit
	  case $msg = $setMinicondaPath
		 _select_folder()
	  case $msg = $cancelButton
		 Exit
	  case $msg = $nextButton
		 If GUICtrlRead($minicondaPath) = "" Then
			$isChecked = GUICtrlRead($installMiniconda) = $GUI_CHECKED
			warningMiniconda3PathEmpty($isChecked)
		 Else
			ExitLoop
		 EndIf
   EndSelect
WEnd

$installDir = GUICtrlRead($minicondaPath)
If GUICtrlRead($addMinicondaToPath) = $GUI_CHECKED Then
   $addToPath = 1
Else
   $addToPath = 0
EndIf

If GUICtrlRead($makeMinicondaStandard) = $GUI_CHECKED Then
   $registerPython = 1
Else
   $registerPython = 0
EndIf

If GUICtrlRead($installMiniconda) = $GUI_CHECKED Then
   $installMinicondaVar = True
Else
   $installMinicondaVar = False
EndIf

GUIDelete($gui1)
$gui2 = GUICreate("MetaboLabPy Installer v0.3", 500, 300, 100, 100) ; ==> draw GUI
GUICtrlCreateLabel("MetaboLabPy Installation in progress, please wait...", 10, 10, 450, 20)
$messageBoard = GUICtrlCreateEdit("", 10, 40, 450, 220)
$finishButton = GUICtrlCreateButton("Finish",440, 270, 50, 20)
GUICtrlSetState($finishButton, $GUI_DISABLE)
WinSetOnTop($gui2, "", 1)
GUISetState() ; show GUI
InstallMetaboLabPy($installDir, $installMinicondaVar, $addToPath, $registerPython)

GUICtrlSetState($finishButton, $GUI_ENABLE)

while 1
   $msg = GUIGetMsg() ; check for GUI activity
   Select
	  case $msg = -3 ; -3 = close by clicking x (top right)
		 Exit
	  case $msg = $finishButton
		 ;SelfDelete()
		 Exit
   EndSelect
WEnd


FUNC _select_folder()
    $sDir = FileSelectFolder("Select Miniconda3 folder", "", 0)
    If @error <> 1 Then
        GUICtrlSetData($minicondaPath,$sDir)
    Else
        ToolTip("user abort",0,0)
        Sleep(500)
        ToolTip("",0,0)
    EndIf
 EndFunc

Func InstallMetaboLabPy($installDir, $installMinicondaVar, $addToPath, $registerPython)
   If $installMinicondaVar Then
	  $msgbData = GUICtrlRead($messageBoard) & "Downloading Miniconda3..." & @CRLF
	  GUICtrlSetData($messageBoard, $msgbData)

	  ; Save the downloaded file to the temporary folder.
	  Local $sFilePath = @TempDir & "\Miniconda3.exe"

	  ; Download the file by waiting for it to complete. The option of 'get the file from the local cache' has been selected.
	  Local $iBytesSize = InetGet("https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe", $sFilePath, $INET_FORCERELOAD)

	  $msgbData = GUICtrlRead($messageBoard) & "Installing Miniconda3..." & @CRLF
	  GUICtrlSetData($messageBoard, $msgbData)

	  ; Run installer
	  $installString = $sFilePath & " /InstallationType=JustMe /AddToPath=" & $addToPath & " /RegisterPython=" & $registerPython & " /S /D=" & $installDir
	  RunWait($installString)

; Delete the file.
	  FileDelete($sFilePath)
	  $msgbData = GUICtrlRead($messageBoard) & "Miniconda3 Installation finished" & @CRLF
	  GUICtrlSetData($messageBoard, $msgbData)
   Else
	  ConsoleWrite("No install!" & @CRLF)
   EndIf
   $path = EnvGet("path")
   EnvSet("path",$installDir & ";" & $installDir & "\Library\mingw-w64\bin;" & $installDir & "\Library\usr\bin;" & $installDir & "\Library\bin;" & $installDir & "\Scripts;" & $path)
   EnvUpdate()
   $msgbData = GUICtrlRead($messageBoard) & "Creating metabolabpy environment..." & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
   $condaExe = $installDir & "\Scripts\conda.exe"
   $condaBat = $installDir & "\condabin\conda.bat"
   $cmdStr = "cmd.exe" & ' /c "' & $condaExe & ' create --name metabolabpy -y"'
   RunWait($cmdStr)
   $path = EnvGet("path")
   EnvSet("path",$installDir & "\envs\metabolabpy;" & $path)
   EnvUpdate()
   $cmdStr = "cmd.exe" & ' /c "' & $condaBat & ' init cmd.exe"'
   RunWait($cmdStr)
   $msgbData = GUICtrlRead($messageBoard) & "Activating metabolabpy environment..." & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
   $cmdStr = "cmd.exe" & ' /c "' & $condaExe & ' init cmd.exe&&'& $condaExe & ' activate metabolabpy"'
   RunWait($cmdStr)
   $msgbData = GUICtrlRead($messageBoard) & "Installing Python 3.7.4, NumPy with Intel acceleration and MetaboLabPy..." & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
   $cmdStr = "cmd.exe" & ' /c " ' & $condaBat & ' activate metabolabpy && ' & $condaExe & ' install python==3.7.4 -y && ' & $condaExe & ' install numpy -y && pip install metabolabpy"'
   $msgbData = GUICtrlRead($messageBoard) & $cmdStr & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
   RunWait($cmdStr)
   $msgbData = GUICtrlRead($messageBoard) & "Creating Starter Script and Desktop lnk..." & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
   $batFileHandle = FileOpen($installDir & "\ml.bat", 2)
   FileWrite($batFileHandle, "start /min " & $installDir & "\ml_exec.bat")
   $sFilePath = @DesktopDir & "\MetaboLabPy.lnk"
   $bat2FileHandle = FileOpen($installDir & "\ml_exec.bat", 2)
   FileWrite($bat2FileHandle, "conda activate metabolabpy && metabolabpy && exit")
   $sFilePath = @DesktopDir & "\MetaboLabPy.lnk"
   FileCreateShortcut($installDir & "\ml.bat", $sFilePath, @DesktopDir, "", "", $installDir & "\envs\metabolabpy\Lib\site-packages\metabolabpy\icon\icon.ico", @SW_SHOWMINIMIZED)
   $msgbData = GUICtrlRead($messageBoard) & "Installation finshed" & @CRLF
   GUICtrlSetData($messageBoard, $msgbData)
EndFunc   ;==>Example

Func warningMiniconda3PathEmpty($isChecked)
   If $isChecked Then
	  MsgBox($MB_ICONWARNING, "MetaboLabPy Install Error", "No miniconda3 install path provided, please try again")
   Else
	  MsgBox($MB_ICONWARNING, "MetaboLabPy Install Error", "Please provide the location of your miniconda3 installation, please try again")
   EndIf
EndFunc

Func warningMetabolabpyPathEmpty()
   MsgBox($MB_ICONWARNING, "MetaboLabPy Install Error", "No metabolabpy git repository path provided, please try again")
EndFunc

Func SelfDelete()
    Local $fp
    Local $Buffer
    Local $lzBatchFile

    ;Batch file.
    $lzBatchFile = "delme.bat"

    $Buffer = ":Repeat" & @CRLF
    $Buffer &= "attrib " & @ScriptName & " -r -s" & @CRLF
    $Buffer &= "del " & @ScriptName & @CRLF
    $Buffer &= "if exist " & @ScriptName & " goto Repeat" & @CRLF
    $Buffer &= "del " & $lzBatchFile & @CRLF
    $Buffer &= "exit" & @CRLF

    ;Write contents of batch to file.
    $fp = FileOpen($lzBatchFile, 2)
    FileWrite($fp, $Buffer)
    FileClose($fp)

    ;Clear up
    $Buffer = ""
    ;Run the script.
    Run($lzBatchFile, @ScriptDir, 0)

EndFunc   ;==>SelfDelete
