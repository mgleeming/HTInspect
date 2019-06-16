import os, sys

from PyQt4 import QtCore, QtGui
#from res.utils import *
from gui.HT_guiprocs import HT_search
from res.utils import *
from multiprocessing import Process, Queue

def main():
    runner_thread = QtCore.QThread()
    runner = Runner(start_signal = runner_thread.started)

    app = QtGui.QApplication(sys.argv)

    try:
        import qdarkstyle as qds
        app.setStyleSheet(qds.load_stylesheet(pyside = False))
        darkMode = True
    except:
        darkMode = False

    gui = HT_search(darkMode = darkMode, runner = (runner, runner_thread, Queue()))
    gui.show()
    sys.exit(app.exec_())

    return

if __name__ == '__main__':
    sys.exit(main())
