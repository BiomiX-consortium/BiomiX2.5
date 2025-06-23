import sys
import subprocess
import os

# No need to manually define conda_python_path anymore
# We'll use sys.executable to inherit the current environment

def rest_of_the_script():
    import os
    import subprocess
    from PyQt5 import QtCore, QtGui, QtWidgets
    from PyQt5.QtCore import Qt, QTimer
    from PyQt5.QtGui import QPixmap
    from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QPushButton, QFileDialog, QSizePolicy, QHBoxLayout, QVBoxLayout, QWidget

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()

            self.setWindowTitle("BiomiX")
            icon = QtGui.QIcon()
            icon.addPixmap(
                QtGui.QPixmap("BiomiX_logo3.png"),
                QtGui.QIcon.Mode.Normal,
                QtGui.QIcon.State.Off,
            )
            self.setWindowIcon(icon)
            self.setGeometry(200, 200, 800, 600)

            central_widget = QWidget(self)
            self.setCentralWidget(central_widget)

            self.figure_label = QLabel(self)
            self.figure_label.setAlignment(Qt.AlignCenter)
            self.figure_label.setScaledContents(True)

            self.label = QLabel("Welcome to BiomiX", self)
            font = QtGui.QFont()
            font.setFamily("Arial")
            font.setPointSize(22)
            font.setBold(False)
            font.setWeight(50)
            self.label.setFont(font)
            self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

            self.label2 = QLabel("Upload the multiomics metadata file\n to start the analysis", self)
            font.setPointSize(16)
            self.label2.setFont(font)
            self.label2.setAlignment(Qt.AlignmentFlag.AlignCenter)

            self.button = QPushButton("Upload", self)
            font.setPointSize(16)
            self.button.setFont(font)
            self.button.clicked.connect(self.open_dialog)

            self.figure_path = "BiomiX_logo3.png"

            self.label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.label2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            self.figure_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

            label_layout = QVBoxLayout()
            label_layout.addWidget(self.label)
            label_layout.addWidget(self.label2)

            bottom_layout = QHBoxLayout()
            bottom_layout.addWidget(self.button)
            bottom_layout.addWidget(self.figure_label)

            main_layout = QVBoxLayout()
            main_layout.addLayout(label_layout)
            main_layout.addLayout(bottom_layout)

            central_widget.setLayout(main_layout)

            QTimer.singleShot(100, self.load_figure)

        def resizeEvent(self, event):
            self.adjust_font_size()
            self.load_figure()
            super().resizeEvent(event)

        def adjust_font_size(self):
            width = self.width()
            height = self.height()
            font_size_label = int(height / 15)
            font_size_label2 = int(height / 20)

            font = self.label.font()
            font.setPointSize(font_size_label)
            self.label.setFont(font)

            font = self.label2.font()
            font.setPointSize(font_size_label2)
            self.label2.setFont(font)

        def load_figure(self):
            pixmap = QPixmap(self.figure_path)
            pixmap = pixmap.scaled(self.figure_label.size(), Qt.AspectRatioMode.KeepAspectRatio)
            self.figure_label.setPixmap(pixmap)

        def open_dialog(self):
            dialog = QFileDialog()
            dialog.setFileMode(QFileDialog.ExistingFile)
            dialog.setNameFilter("Metadata Files (*.txt *.csv *.tsv *.xls *.xlsx)")
            if dialog.exec_():
                file_path = dialog.selectedFiles()[0]
                print(file_path)
                directory = os.path.dirname(file_path)
                self.write_directory(file_path)

                # Use current Python executable to run second script
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "MODULE_interface.py")
                python_cmd = [sys.executable, script_path]
                print(f"Launching: {python_cmd}")
                subprocess.Popen(python_cmd)

            self.close()

        def write_directory(self, directory):
            script_directory = os.path.dirname(os.path.abspath(__file__))
            output_path = os.path.join(script_directory, "directory.txt")
            with open(output_path, "w") as file:
                file.write(directory)

    if __name__ == '__main__':
        print("BiomiX interface is loading...")
        app = QApplication([])
        window = MainWindow()
        window.show()
        app.exec_()


if __name__ == "__main__":
    # Always run the main GUI script
    rest_of_the_script()
