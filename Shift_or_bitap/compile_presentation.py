#!/usr/bin/env python3
"""
Compile Shift-Or Presentation (LaTeX/Beamer) to PDF

Usage:
    python compile_presentation.py

Requirements:
    - pdflatex installed (TeX distribution with beamer)
    - shift_or_presentation.tex in current directory
"""

import subprocess
import os
import sys

def check_pdflatex():
    """Check if pdflatex is installed"""
    try:
        subprocess.run(['pdflatex', '--version'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def compile_presentation():
    """Compile Beamer presentation to PDF"""
    
    tex_file = 'shift_or_presentation.tex'
    
    if not os.path.exists(tex_file):
        print(f"Error: {tex_file} not found!")
        print("Make sure shift_or_presentation.tex is in the current directory.")
        return False
    
    if not check_pdflatex():
        print("Error: pdflatex not found!")
        print("\nInstall TeX distribution:")
        print("  Ubuntu/Debian: sudo apt-get install texlive-latex-extra texlive-fonts-recommended")
        print("  macOS: brew install mactex")
        print("  Windows: Download from https://miktex.org/")
        return False
    
    print("=" * 70)
    print("COMPILING SHIFT-OR PRESENTATION (BEAMER)")
    print("=" * 70)
    
    # Compilation pass (Beamer typically needs only one pass)
    print("\n[1/2] Compiling presentation...")
    result = subprocess.run(
        ['pdflatex', '-interaction=nonstopmode', tex_file],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print("✗ Compilation failed!")
        print(result.stdout[-1000:])  # Show last 1000 chars
        return False
    
    print("✓ Compilation completed")
    
    # Clean up auxiliary files
    print("[2/2] Cleaning up auxiliary files...")
    aux_extensions = ['.aux', '.log', '.out', '.nav', '.snm', '.toc']
    for ext in aux_extensions:
        aux_file = tex_file.replace('.tex', ext)
        if os.path.exists(aux_file):
            os.remove(aux_file)
    
    # Check if PDF was created
    pdf_file = tex_file.replace('.tex', '.pdf')
    if os.path.exists(pdf_file):
        pdf_size = os.path.getsize(pdf_file) / (1024 * 1024)  # Size in MB
        print(f"\n✓ PRESENTATION COMPILED SUCCESSFULLY!")
        print(f"  File: {pdf_file}")
        print(f"  Size: {pdf_size:.2f} MB")
        print(f"  Slides: 11 (Title + 10 content slides)")
        print(f"\nYou can now:")
        print(f"  - View slides: xdg-open {pdf_file}")
        print(f"  - Present with PDF viewer (F5 for fullscreen)")
        print(f"  - Share with colleagues/advisor")
        return True
    else:
        print("✗ PDF file not created!")
        return False

if __name__ == '__main__':
    success = compile_presentation()
    sys.exit(0 if success else 1)
