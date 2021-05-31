
__authors__ = ["S. Basu"]
__license__ = "MIT"
__date__ = "07/05/2021"

import os
import json
import pathlib
import logging

logger = logging.getLogger("sxdm")

class Abstract(object):
    def __init__(self, inData):
          self._ioDict = dict()
          self._ioDict['inData'] = json.dumps(inData, default=str)
          self._ioDict['outData'] = json.dumps(dict(), default=str)
          self._ioDict['success'] = True
          self._ioDict['WorkingDirectory'] = None
          self._logFileName = None
          return

    def get_inData(self):
        return json.loads(self._ioDict['inData'])

    def set_inData(self, inData):
        self._ioDict['inData'] = json.dumps(inData, default=str)

    jshandle = property(get_inData, set_inData)

    def get_outData(self):
        return json.loads(self._ioDict['outData'])

    def set_outData(self, results):
        self._ioDict['outData'] = json.dumps(results, default=str)

    results = property(get_outData, set_outData)

    def writeInputData(self, inData):
        # Write input data
        if self._ioDict['WorkingDirectory'] is not None:
            jsonName = "inData_" + self.__class__.__name__ + ".json"
            with open(str(self.getOutputDirectory() / jsonName), 'w') as f:
                f.write(json.dumps(inData, default=str, indent=4))
        return

    def writeOutputData(self, results):
        self.set_outData(results)
        if self._ioDict['WorkingDirectory'] is not None:
            jsonName = "outData_" + self.__class__.__name__ + ".json"
            with open(str(self.getOutputDirectory() / jsonName), 'w') as f:
                f.write(json.dumps(results, default=str, indent=4))
        return

    def setFailure(self):
        self._ioDict['success'] = False

    def is_success(self):
        return self._ioDict['success']

    def setOutputDirectory(self, path=None):
        if path is None:
            directory = self.jshandle.get('processing_directory', os.getcwd())
            self._ioDict['WorkingDirectory'] = pathlib.Path(directory)
        else:
            self._ioDict['WorkingDirectory'] = pathlib.Path(path)

    def getOutputDirectory(self):
        return self._ioDict['WorkingDirectory']

    def getLogPath(self):
        if self._logFileName is None:
            self._logFileName = self.__class__.__name__ + ".log"
        logPath = self._ioDict['WorkingDirectory'] / self._logFileName
        return logPath

    def setLogFileName(self, logFileName):
        self._logFileName = logFileName

    def getLogFileName(self):
        return self._logFileName

    def readLog(self):
        with open(str(self.getLogPath())) as f:
            log = f.read()
        return log

