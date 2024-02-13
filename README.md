<div align="center">
# LuminetCpp
Unofficial C++ version of Luminet project of bgmeulem //https://github.com/bgmeulem/Luminet/.
Work in progress. Do not use as a lot of  functionalities have to be implemented.
# Partial results :
<img src="https://github.com/Niohori/LuminetCpp/blob/main/Documentation/Isoradial.gif" width="800" />
<img src="https://github.com/Niohori/LuminetCpp/blob/main/Documentation/Isoradial_89.PNG" width="800" />
<img src="https://github.com/Niohori/LuminetCpp/blob/main/Documentation/Isoredshifts3.gif" width="800" />

# TODO

## Bugs
- Calculation of redshift isolines are quite fast(about 30 s for 180 different inclinations) but quality of rendering has to be improved.
- Redshift color rendering is not realistic

## Improvements
- Cleaning the code.
- Create from scratch the accretion disk 

## The aim is also to give some clarifications to the maths behind the program (in progress). 
<img src="https://github.com/Niohori/LuminetCpp/blob/main/Documentation/Math/images/Coordinates_system.PNG" width="800" />

# Dependencies
- Boost for the elliptic integrals
- Dislin for the graphical rendering //https://www.dislin.de/
- dlib for multithreading //http://www.dlib.net/api.html
