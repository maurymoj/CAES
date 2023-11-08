% pressure over 70 km pipeline
% pyversion c:\Anaconda3\python.exe % Sets Python 3.11 as a interpreter
CP = py.importlib.import_module('CoolProp.CoolProp');
%%
CP.PropsSI('D','P',101325,'T',298,'Air')
