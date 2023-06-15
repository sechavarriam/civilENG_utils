import matplotlib.pyplot as plt

f = open('TestData/Tsukidate-MYG004/MYG0041103111446.NS', 'r')

Date      = (f.readline().split())[-2]
Latitude  = float((f.readline().split())[-1])
Longitude = float((f.readline().split())[-1])
Depth     = float((f.readline().split())[-1])
Magnitude = float((f.readline().split())[-1])

(f.readline().split()) # Station Code
(f.readline().split()) # Station Latitude
(f.readline().split()) # Station Longitude
(f.readline().split()) # Station Height (m)
(f.readline().split()) # Record Time

SampFreq = float((((f.readline().split())[-1]).split('H'))[0])# Sampling Frequency

Duration = float((f.readline().split())[-1]) #Duration Time
Direction = (f.readline().split())[-1]

ScaleFactor = ((f.readline().split())[-1]).split('/')

Aux1 = float((ScaleFactor[0].split('('))[0])
Aux2 = float(ScaleFactor[1])
ScaleFactor = Aux1/Aux2 # [gal]

MaxAcc = float((f.readline().split())[-1])

(f.readline().split())# Last correction date.
(f.readline().split())

Acc_DATA = f.readlines()
f.close()

Aux_Acc  = [list(map(float,i.split())) for i in Acc_DATA]


Acceleration = [rawData*ScaleFactor for sublist in Aux_Acc for rawData in sublist] #[gal]
Time = [ t/SampFreq for t in range(0, len(Acceleration))]


plt.plot(Time, Acceleration)
plt.show()



















