from PyQt5.QtGui import QIntValidator
import stochastic as sde
import matplotlib.pyplot as plt

import sys
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QWidget, QLabel, QLineEdit
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtGui import QIntValidator
from PyQt5.QtCore import QSize    



class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        
        self.setMinimumSize(QSize(320, 140))    
        self.setWindowTitle("Parameter Input") 

        self.nameLabel = QLabel(self)
        self.nameLabel.setText('Start Time')
        self.line1 = QLineEdit(self)


        self.line1.move(80, 10)
        self.line1.resize(200, 32)
        self.nameLabel.move(20, 10)

        self.nameLabel = QLabel(self)
        self.nameLabel.setText('End Time')
        self.line2 = QLineEdit(self)
      
   

        self.line2.move(80, 40)
        self.line2.resize(200, 32)
        self.nameLabel.move(20, 40)  

        
        self.nameLabel = QLabel(self)
        self.nameLabel.setText('N Steps')
        self.line3 = QLineEdit(self)
    
     

        self.line3.move(80, 70)
        self.line3.resize(200, 32)
        self.nameLabel.move(20, 70)

        pybutton = QPushButton('OK', self)

        pybutton.clicked.connect(self.clickedon)
        pybutton.clicked.connect(self.close)
        pybutton.resize(100,32)
        pybutton.move(80, 105)

        pybutton = QPushButton('Cancel', self)

        pybutton.clicked.connect(self.close)
        pybutton.resize(100,32)
        pybutton.move(170, 105)


    def line_one(self):
        self.A = self.line1.text()
        return self.A

    def line_two(self):
        self.B = self.line2.text()
        return self.B

    def line_three(self):
        self.C = self.line3.text()
        return self.C

    def clickedon(self):
       
        print('Your Start Time: ' + self.line_one())
        print('Your Start Time: ' + self.line_two())
        print('Your Start Time: ' + self.line_three())
        
        
    
    def connection(self):
        self.start = int(self.A)
        self.end = int(self.B)
        self.n = int(self.C)

        gg = sde.Graph(self.start, self.end, self.n)
        gg.graph()
        plt.show()




    
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    app.exec_()
    mainWin.line_one()
    mainWin.line_two()
    mainWin.line_three()
    mainWin.connection()



