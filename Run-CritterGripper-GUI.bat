@echo off
title CritterGripper GUI
cd /d "%~dp0"

echo.
echo  CritterGripper — starting Streamlit
echo  Bookmark:  http://localhost:8501
echo  (Use this exact URL; 127.0.0.1 is a different origin for some browsers.)
echo.

streamlit run bender_streamlit_gui.py
if errorlevel 1 pause
