@echo off
title CritterGripper NiceGUI
cd /d "%~dp0"

echo.
echo  CritterGripper — starting NiceGUI (not Streamlit)
echo  Open:  http://localhost:8765  (landing page — use Enter workflow)
echo  (Keep this window open while you use the app.)
echo.

python main.py
if errorlevel 1 pause
