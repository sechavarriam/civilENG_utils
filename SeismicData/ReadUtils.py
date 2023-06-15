
# Funtion to read accelerogram data from a .txt preprocesed file f which 
# contents a list of times an accelertions previously procesed from other formats.

def importAccelerogram(file_location):
    
    # FALTA PONER EXCEPCIONES, CHEQUEOS.... 

    f = open(file_location, 'r')

    raw = f.readlines()

    f.close()

    n = len(raw)

    time = [float(raw[t].split()[0]) for t in range(n)] # Time discretization corresponding to ag
    ag   = [float(raw[t].split()[1]) for t in range(n)]                    

    return [time, ag] 









