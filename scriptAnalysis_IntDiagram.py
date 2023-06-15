#

# Data to be used in funtion


## INPUT SECCIÓN ===============
# La sección podría ser un objeto en el futurio.

# Dimensiones sección S3
b =  1000 #[mm]
h =  500 #[mm]

d4 = 50      # [mm]
d3 = 190     # [mm]
d2 = 500-190 # [mm]
d1 = 500-50  # [mm]


# Áreas de barras =====
N6 = 284 # [mm]^2
N7 = 387 # [mm]^2
N8 = 510 # [mm]^2 
# =====================
As_1 = 6*N7
As_2 = 2*N7
As_3 = 2*N7
As_4 = 6*N7


d = [d1, d2, d3, d4]
As= [As_1, As_2, As_3, As_4]

## INPUT MATERIAL ==============
# El material podría ser un objeto en el futuro. Sin embargo
# para el análisis elástico no justiifca.


fc = 35  # [MPa] f'c
fy = 420 # [MPa]

Es = 200e3 # [MPa]

phi = 0.75 # Estribos
#phi = 0.80 # Espiral





























