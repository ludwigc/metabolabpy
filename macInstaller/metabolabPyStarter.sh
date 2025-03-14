#!/usr/bin/osascript

on is_running(appName)
    tell application "System Events" to (name of processes) contains appName
end is_running

set termRunning to is_running("Terminal")

tell Application "Terminal"
    set newTab to do script "conda activate metabolabpy; metabolabpy"
    set newWindow to id of front window

    tell window id newWindow
        set index to 1
        set visible to false
    end tell
    repeat
        delay 0.1
        if not busy of newTab then exit repeat
    end repeat
    if termRunning then
        tell Application "Terminal" to close
    else
        close (every window whose id ≠ newWindow)
        tell Application "Terminal" to quit
    end if
end tell
