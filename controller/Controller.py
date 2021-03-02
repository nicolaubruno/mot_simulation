#
# Libraries and modules
from model import Model
from view import View

#
# Controller class
#
class Controller:
    
    ''' Attributes '''


    ''' Methods '''

    #
    def __init__(self):
        self.__model = Model()
        self.__view = View(self.__model)

    #
    def main(self):
        call_menu = True

        while call_menu:
            option_code = self.__view.initial_menu()
            call_menu = False

            #
            # Run Simulation
            if option_code == 1:
                self.__model.run_simulation()

            #
            # Show parameters of the simulation
            elif option_code == 2:
                code = self.__view.parameters()

                if code == 7:
                    call_menu = True

            #
            # Show results
            elif option_code == 3:
                pass