#!/usr/bin/env python3
"""
Boltz Protein Prediction GUI

A PySide6-based graphical user interface for the Boltz protein structure prediction pipeline.
Provides an intuitive interface for loading FASTA files, visualizing proteins, and configuring
predictions without requiring command line usage.

Usage:
    python boltz_gui.py

Requirements:
    PySide6, biopython, pandas, numpy
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import string

try:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget,
        QPushButton, QFileDialog, QTableWidget, QTableWidgetItem,
        QLabel, QGroupBox, QSplitter, QTextEdit, QHeaderView,
        QFrame, QMessageBox, QProgressBar, QTabWidget,
        QSpinBox, QCheckBox, QLineEdit, QFormLayout
    )
    from PySide6.QtCore import Qt, QMimeData, QTimer, QThread, Signal
    from PySide6.QtGui import QDragEnterEvent, QDropEvent, QColor, QFont
except ImportError:
    print("PySide6 is required. Please install it with: pip install PySide6")
    sys.exit(1)

try:
    from Bio import SeqIO
    import pandas as pd
    import numpy as np
    from protein_parser import parse_protein_files, ProteinRecord
except ImportError:
    print("Required packages missing. Please install: pip install biopython pandas numpy")
    print("Also ensure protein_parser.py is in the same directory as this GUI.")
    sys.exit(1)


class ProteinData:
    """Data class to hold protein information"""
    def __init__(self, identifier: str, name: str, sequence: str, source_file: str, file_index: int):
        self.identifier = identifier
        self.name = name
        self.sequence = sequence
        self.source_file = source_file
        self.file_index = file_index
        self.length = len(sequence)



class ProteinTableWidget(QTableWidget):
    """Custom table widget for displaying protein information with color coding"""
    
    def __init__(self):
        super().__init__()
        self.setup_table()
        self.file_colors = []  # Store colors for each file
        
    def setup_table(self):
        """Initialize the table structure"""
        self.setColumnCount(5)
        headers = ["ID", "Protein Name", "Length", "Source File", "Sequence"]
        self.setHorizontalHeaderLabels(headers)
        
        # Set column widths
        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Fixed)  # ID
        header.setSectionResizeMode(1, QHeaderView.Stretch)  # Name
        header.setSectionResizeMode(2, QHeaderView.Fixed)  # Length
        header.setSectionResizeMode(3, QHeaderView.Stretch)  # Source
        header.setSectionResizeMode(4, QHeaderView.Stretch)  # Sequence
        
        self.setColumnWidth(0, 50)   # ID
        self.setColumnWidth(2, 80)   # Length
        
        # Table styling
        self.setAlternatingRowColors(True)
        self.setSelectionBehavior(QTableWidget.SelectRows)
        self.setStyleSheet("""
            QTableWidget {
                gridline-color: #d0d0d0;
                background-color: white;
            }
            QTableWidget::item {
                padding: 8px;
            }
        """)
    
    def add_proteins(self, proteins: List[ProteinData], file_colors: List[QColor]):
        """Add proteins to the table with appropriate color coding"""
        self.file_colors = file_colors
        
        # Track protein count per file for hue variation
        file_protein_counts = {}
        
        for protein in proteins:
            row_pos = self.rowCount()
            self.insertRow(row_pos)
            
            # Get base color for this file
            base_color = file_colors[protein.file_index]
            
            # Count proteins per file for hue variation
            if protein.file_index not in file_protein_counts:
                file_protein_counts[protein.file_index] = 0
            file_protein_counts[protein.file_index] += 1
            
            # Create varied hue for each protein within the same file
            bg_color = self.create_protein_color(base_color, file_protein_counts[protein.file_index])
            
            # Create table items
            id_item = QTableWidgetItem(protein.identifier)
            name_item = QTableWidgetItem(protein.name)
            length_item = QTableWidgetItem(str(protein.length))
            source_item = QTableWidgetItem(os.path.basename(protein.source_file))
            sequence_item = QTableWidgetItem(self.truncate_sequence(protein.sequence))
            
            # Apply color coding with black text
            for item in [id_item, name_item, length_item, source_item, sequence_item]:
                item.setBackground(bg_color)
                item.setForeground(QColor(0, 0, 0))  # Black text
                # Store full protein data in the first column for easy access
                if item == id_item:
                    item.setData(Qt.UserRole, protein)
            
            # Set items in table
            self.setItem(row_pos, 0, id_item)
            self.setItem(row_pos, 1, name_item)
            self.setItem(row_pos, 2, length_item)
            self.setItem(row_pos, 3, source_item)
            self.setItem(row_pos, 4, sequence_item)
    
    def create_protein_color(self, base_color: QColor, protein_index: int) -> QColor:
        """Create a unique hue variation for each protein within a file"""
        # Get HSV values from base color
        hue = base_color.hue()
        saturation = base_color.saturation()
        
        # Vary the lightness for each protein (keeping it light for readability)
        # Start at 90% lightness and decrease by steps
        lightness = max(240 - (protein_index * 15), 200)  # Keep it light
        
        # Create new color with varied lightness
        color = QColor()
        color.setHsv(hue, saturation // 3, lightness)  # Reduce saturation for lighter appearance
        return color
    
    def truncate_sequence(self, sequence: str, max_length: int = 50) -> str:
        """Truncate sequence for display purposes"""
        if len(sequence) <= max_length:
            return sequence
        return f"{sequence[:max_length]}..."
    
    def clear_proteins(self):
        """Clear all proteins from the table"""
        self.setRowCount(0)
        self.file_colors = []


class FileInfoWidget(QWidget):
    """Widget to display loaded file information"""
    
    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.file_colors = []
    
    def setup_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        
        self.title_label = QLabel("Files:")
        self.title_label.setFont(QFont("Arial", 10, QFont.Bold))
        layout.addWidget(self.title_label)
        
        self.files_widget = QWidget()
        self.files_layout = QVBoxLayout(self.files_widget)
        self.files_layout.setContentsMargins(0, 0, 0, 0)
        self.files_layout.setSpacing(3)
        layout.addWidget(self.files_widget)
        
        layout.addStretch()  # Push content to top
        self.setLayout(layout)
    
    def update_files(self, file_paths: List[str], file_colors: List[QColor]):
        """Update the display with loaded files and their colors"""
        self.file_colors = file_colors
        
        # Clear existing widgets
        for i in reversed(range(self.files_layout.count())):
            child = self.files_layout.itemAt(i).widget()
            if child:
                child.deleteLater()
        
        # Add simple file entries
        for i, (file_path, color) in enumerate(zip(file_paths, file_colors)):
            file_layout = QHBoxLayout()
            file_layout.setContentsMargins(0, 0, 0, 0)
            
            # Color indicator
            color_label = QLabel("●")
            color_label.setStyleSheet(f"color: {color.name()}; font-size: 16px;")
            color_label.setFixedWidth(20)
            file_layout.addWidget(color_label)
            
            # File name
            file_name = QLabel(os.path.basename(file_path))
            file_name.setStyleSheet("font-size: 11px;")
            file_layout.addWidget(file_name)
            
            file_layout.addStretch()
            
            # Container widget
            container = QWidget()
            container.setLayout(file_layout)
            self.files_layout.addWidget(container)


class BoltzGUI(QMainWindow):
    """Main GUI application for Boltz protein prediction pipeline"""
    
    def __init__(self):
        super().__init__()
        self.proteins: List[ProteinData] = []
        self.loaded_files: List[str] = []
        self.file_colors: List[QColor] = []
        self.setup_ui()
        self.setup_colors()
        
        # Enable drag and drop on the main window
        self.setAcceptDrops(True)
        
    def setup_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("Boltz Protein Structure Prediction GUI")
        self.setGeometry(100, 100, 1200, 800)
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Title
        title = QLabel("Boltz Protein Structure Prediction Pipeline")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("padding: 10px;")
        main_layout.addWidget(title)
        
        # Create tabs
        tabs = QTabWidget()
        main_layout.addWidget(tabs)
        
        # File Loading Tab
        file_tab = QWidget()
        file_layout = QVBoxLayout(file_tab)
        
        # Top section with loading controls - keep it minimal
        top_section = QWidget()
        top_layout = QHBoxLayout(top_section)
        top_layout.setContentsMargins(10, 5, 10, 5)
        
        # User-friendly instruction
        instruction_label = QLabel("Click here to load FASTA files or text files for proteins to be predicted")
        instruction_label.setStyleSheet("color: #666; font-size: 12px;")
        top_layout.addWidget(instruction_label)
        
        top_layout.addStretch()
        
        # Load button
        load_button = QPushButton("Load Files")
        load_button.clicked.connect(self.load_files)
        load_button.setFixedSize(120, 30)
        top_layout.addWidget(load_button)
        
        file_layout.addWidget(top_section)
        
        # Main content area - this will expand to fill remaining space
        content_splitter = QSplitter(Qt.Horizontal)
        
        # File info widget (left, smaller)
        self.file_info_widget = FileInfoWidget()
        self.file_info_widget.setMaximumWidth(200)
        self.file_info_widget.setMinimumWidth(150)
        content_splitter.addWidget(self.file_info_widget)
        
        # Protein table (right, larger - emphasis on this)
        table_group = QGroupBox("Loaded Proteins")
        table_layout = QVBoxLayout(table_group)
        table_layout.setContentsMargins(5, 15, 5, 5)
        
        self.protein_table = ProteinTableWidget()
        table_layout.addWidget(self.protein_table)
        
        content_splitter.addWidget(table_group)
        
        # Set splitter proportions - emphasize the table
        content_splitter.setStretchFactor(0, 0)  # File info fixed
        content_splitter.setStretchFactor(1, 1)  # Protein table gets all remaining space
        
        # Add the splitter and make it expand to fill remaining vertical space
        file_layout.addWidget(content_splitter, 1)  # stretch factor = 1
        
        tabs.addTab(file_tab, "File Loading & Protein Visualization")
        
        # Prediction Configuration Tab (placeholder for future expansion)
        config_tab = QWidget()
        config_layout = QVBoxLayout(config_tab)
        
        config_group = QGroupBox("Prediction Configuration")
        config_form = QFormLayout(config_group)
        
        self.n_mer_spin = QSpinBox()
        self.n_mer_spin.setRange(1, 10)
        self.n_mer_spin.setValue(2)
        config_form.addRow("N-mer size (r):", self.n_mer_spin)
        
        self.output_dir_edit = QLineEdit("./boltz_gui_output")
        config_form.addRow("Output directory:", self.output_dir_edit)
        
        self.recycling_steps_spin = QSpinBox()
        self.recycling_steps_spin.setRange(1, 10)
        self.recycling_steps_spin.setValue(3)
        config_form.addRow("Recycling steps:", self.recycling_steps_spin)
        
        self.gpu_count_spin = QSpinBox()
        self.gpu_count_spin.setRange(1, 8)
        self.gpu_count_spin.setValue(1)
        config_form.addRow("GPU count:", self.gpu_count_spin)
        
        self.skip_msa_check = QCheckBox()
        config_form.addRow("Skip MSA generation:", self.skip_msa_check)
        
        config_layout.addWidget(config_group)
        
        # Run button (placeholder)
        run_button = QPushButton("Generate Prediction Pipeline")
        run_button.setMinimumHeight(40)
        config_layout.addWidget(run_button)
        
        tabs.addTab(config_tab, "Prediction Configuration")
        
        # Status bar
        self.statusBar().showMessage("Ready - Load FASTA files to begin")
    
    def setup_colors(self):
        """Initialize color palette for different files"""
        # Define a set of distinct colors for different files
        self.color_palette = [
            QColor("#e74c3c"),  # Red
            QColor("#2ecc71"),  # Green  
            QColor("#3498db"),  # Blue
            QColor("#f39c12"),  # Orange
            QColor("#9b59b6"),  # Purple
            QColor("#1abc9c"),  # Teal
            QColor("#34495e"),  # Dark gray
            QColor("#e67e22"),  # Orange
            QColor("#95a5a6"),  # Light gray
            QColor("#27ae60"),  # Dark green
        ]
    
    def load_files(self):
        """Open file dialog to load FASTA/text files"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select FASTA or Text Files",
            "",
            "Protein Files (*.fasta *.fas *.fa *.txt);;FASTA Files (*.fasta *.fas *.fa);;Text Files (*.txt);;All Files (*)"
        )
        
        if file_paths:
            self.load_files_from_paths(file_paths)
    
    def dragEnterEvent(self, event: QDragEnterEvent):
        """Handle drag enter event for the main window"""
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            valid_files = [url.toLocalFile() for url in urls 
                          if url.toLocalFile().lower().endswith(('.fasta', '.fas', '.fa', '.txt'))]
            if valid_files:
                event.acceptProposedAction()
    
    def dropEvent(self, event: QDropEvent):
        """Handle drop event for the main window"""
        urls = event.mimeData().urls()
        valid_files = [url.toLocalFile() for url in urls 
                      if url.toLocalFile().lower().endswith(('.fasta', '.fas', '.fa', '.txt'))]
        if valid_files:
            self.load_files_from_paths(valid_files)
            event.acceptProposedAction()
    
    def load_files_from_paths(self, file_paths: List[str]):
        """Load FASTA/text files from given paths using shared parser"""
        try:
            # Clear existing data
            self.proteins.clear()
            self.loaded_files.clear()
            self.file_colors.clear()
            self.protein_table.clear_proteins()
            
            # Use the shared parsing logic from protein_parser module
            protein_records = parse_protein_files(file_paths)
            
            if not protein_records:
                QMessageBox.information(self, "No Proteins Found", 
                                      "No valid protein sequences were found in the selected files.")
                return
            
            # Create file mapping for color assignment
            file_set = list(set([record.source_file for record in protein_records.values()]))
            self.loaded_files = file_set
            
            # Assign colors to files
            for i, file_path in enumerate(self.loaded_files):
                color_index = i % len(self.color_palette)
                self.file_colors.append(self.color_palette[color_index])
            
            # Convert parsed records to GUI format with identifiers
            identifier_counter = 0
            alphabet = string.ascii_uppercase
            
            for protein_id, record in protein_records.items():
                # Generate unique letter identifier for GUI
                if identifier_counter < len(alphabet):
                    gui_identifier = alphabet[identifier_counter]
                else:
                    # Handle more than 26 proteins
                    gui_identifier = f"{alphabet[(identifier_counter // 26) - 1]}{alphabet[identifier_counter % 26]}"
                
                # Find file index for color coding
                file_index = self.loaded_files.index(record.source_file)
                
                protein = ProteinData(
                    identifier=gui_identifier,
                    name=record.name,
                    sequence=record.sequence,
                    source_file=record.source_file,
                    file_index=file_index
                )
                self.proteins.append(protein)
                identifier_counter += 1
            
            # Update UI
            self.protein_table.add_proteins(self.proteins, self.file_colors)
            self.file_info_widget.update_files(self.loaded_files, self.file_colors)
            
            # Update status
            total_proteins = len(self.proteins)
            total_files = len(self.loaded_files)
            self.statusBar().showMessage(
                f"Loaded {total_proteins} proteins from {total_files} files"
            )
            
        except Exception as e:
            QMessageBox.critical(self, "Loading Error", f"Failed to load files: {str(e)}")
    
    def get_loaded_proteins(self) -> List[ProteinData]:
        """Get the currently loaded proteins for downstream processing"""
        return self.proteins.copy()
    
    def get_config(self) -> Dict:
        """Get current prediction configuration"""
        return {
            'n_mer_size': self.n_mer_spin.value(),
            'output_dir': self.output_dir_edit.text(),
            'recycling_steps': self.recycling_steps_spin.value(),
            'gpu_count': self.gpu_count_spin.value(),
            'skip_msa': self.skip_msa_check.isChecked(),
            'loaded_files': self.loaded_files.copy()
        }


def main():
    """Main application entry point"""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("Boltz Protein Prediction GUI")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("Boltz Pipeline")
    
    # Create and show the main window
    window = BoltzGUI()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()