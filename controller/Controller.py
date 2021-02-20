#
# Libraries and modules
import click
from model import Model

#
# Controller class
#
class Controller:
    
    ''' Attributes '''

    #
    @property
    def model(self):
        self._model


    ''' Methods '''

    #
    def __init__(self):
        self._model = Model()
