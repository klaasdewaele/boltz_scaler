#!/usr/bin/env python3
"""
Simple launcher script for the Boltz GUI with dependency checking
"""

import sys
import subprocess
import os

def check_and_install_dependencies():
    """Check for required packages and offer to install if missing"""
    required_packages = [
        'PySide6',
        'Bio',
        'pandas', 
        'numpy'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            if package == 'Bio':
                import Bio
            else:
                __import__(package)
            print(f"✓ {package} is available")
        except ImportError:
            print(f"✗ {package} is missing")
            missing_packages.append(package)
    
    if missing_packages:
        print(f"\nMissing packages: {', '.join(missing_packages)}")
        print("Install with: pip install -r gui_requirements.txt")
        
        response = input("\nWould you like to install the missing packages now? (y/n): ")
        if response.lower() in ['y', 'yes']:
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', 'gui_requirements.txt'])
                print("✓ Dependencies installed successfully!")
                return True
            except subprocess.CalledProcessError as e:
                print(f"✗ Failed to install dependencies: {e}")
                return False
        else:
            print("Please install the dependencies manually before running the GUI.")
            return False
    
    return True

def main():
    """Main launcher function"""
    print("Boltz Protein Prediction GUI Launcher")
    print("=" * 40)
    
    # Check if GUI file exists
    gui_file = "boltz_gui.py"
    if not os.path.exists(gui_file):
        print(f"✗ GUI file '{gui_file}' not found in current directory")
        return 1
    
    # Check dependencies
    if not check_and_install_dependencies():
        return 1
    
    print("\nLaunching Boltz GUI...")
    try:
        # Import and run the GUI
        import boltz_gui
        boltz_gui.main()
    except Exception as e:
        print(f"✗ Failed to launch GUI: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())