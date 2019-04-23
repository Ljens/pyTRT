# pyTRT

**Set of classes and functions for analysing Thermal Response Test (TRT) data for Shallow Geothermal applications.**

**Includes implementation of the following TRT models:**
- Infinte Line Source Model (ILSM)
    - Simulation
    - Linear Regression
    - 2 paramter optimisation
    - 3 paramter optimistaion
    - MARS segmented regression (implements R functions, see MARS folder)
- Infinite Cylindrical Source Model (ICSM)
    - 2 paramter optimisation
    - 3 paramter optimistaion
- Finite Line Source Model (FLSM)
    - 2 paramter optimisation
    - 3 paramter optimistaion
- Pile G-function Model:
    (Loveridge, Fleur, and William Powrie. "Temperature response functions (G-functions) for single pile heat exchangers." Energy 57 (2013): 554-564.)
    - 2 paramter optimisation
    - 3 paramter optimistaion

**General Implementation**
- TRT_class is the wrapper and API interface of the underlying models and fitting functions
- 



