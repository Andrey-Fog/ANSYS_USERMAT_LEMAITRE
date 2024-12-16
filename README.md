# ANSYS-USERMAT-LEMAITRE
The Lemaitre damage model with combined hardening is implemented using ANSYS user-programmable features within the USERMAT subroutine.
The Lemaitre damage model is a continuum damage mechanics approach that characterizes the progressive degradation of material stiffness due to microstructural deterioration. This model is particularly effective for predicting the initiation and evolution of damage in metals under complex loading conditions. The incorporation of a combined isotropic-kinematic hardening law allows for modeling of the cyclic plasticity and ratcheting effects observed in real materials. By accounting for damage accumulation, the model enhances the prediction of material durability and failure under fatigue and monotonic loading scenarios.

<center>

![image](https://github.com/user-attachments/assets/0239e339-26d3-4ccf-a9c3-d245d879b33a)


</center>
<b>Distribution of damage parameter depending on the displacement of the upper edge of the 3D specimen</b>
<br>
<br>
If you using this code for research or industrial purposes please cite:   

Tumanov A.V., Kosov D.A., Fedorenkov D.I. 
"Numerical and experimental methods for determining the parameters of generalized models of a damaged visco-plastic medium in durability prediction"
[DOI:10.15593/perm.mech/2024.5.10](https://doi.org/10.15593/perm.mech/2024.5.10)

 ## Extended theory  
- [Lemaitre J. (1985) "A continuous damage mechanics model for ductile fracture"](http://dx.doi.org/10.1115/1.3225775)

- [De Souza Neto E., Peric D. and Owen D. (2008) "Computational Methods for Plasticity: Theory and Applications"](http://dx.doi.org/10.1002/9780470694626)  

- [Robert Lee Gates (2018) "A Finite Element Implementation of a Ductile"](https://arxiv.org/pdf/1302.2439)



## Acknowledgment
This project was made possible thanks to the support of Russian Science Foundation (Project №  [24-29-00475](https://rscf.ru/en/project/24-29-00475/))


<br>
<br>
<br>

## Instructions for compiling and attaching USERMATLIB.DLL 

---

Install Visual Studio first and then Intel fortran compiler. When installing the compiler, select "Integrate into Visual Studio". Supported versions can be found in the ANSYS documentation in the section on User Programmable Features (UPF). Add LIB and INCLUDE variables in the system environment. Create new solution and add new fortran dll project. The name of the created library must be "USERMATLIB.DLL". Add all fortran files from Source directory to your dll project. Tune compiler according to instructions present below. After compiling connect library to ANSYS.


<br>

### Connecting to ANSYS

---

After creation the dll file you have to connect this library to ANSYS:

<br>

**1. Create environment variable named ANS_USER_PATH**

*My Computer->Properties->Advanced system settings->Advanced*  

On the tab, click on the button:

*Environment Variables->System Variables->New*

<br>

**2. In the variable value field, specify the path to the folder where library is located. Use only latin characters in the path.**

*For example:* 
>C:\Username\......\Usermatlib
   
If everything is connected correctly in the ANSYS output window at startup there will be a line 


>User link path <ANS_USER_PATH>: *path to your folder*" 

<br>

**3. After launching the ANSYS, create an user material**

*Preprocessor->Material Props->Material models*

<br>

**4. In the drop-down list of materials, select**

*Structural->Specialized Materials->User material options->user material*


And add cells. There should be 6 properties in total. Of which:

| NN  |     | Property                              | 
| --- | --- | ------------------------------------  | 
|  C1 |  -  |Young modulus                          | 
|  C2 |  -  |Puasson ratio                          | 
|  C3 |  -  |Yelding stress                         | 
|  C4 |  -  |Damage law constant                    | 
|  C5 |  -  |Damage law constant                    | 
|  C6 |  -  |Isotropic hardening constant           | 
|  C7 |  -  |Isotropic hardening constant           | 
|  C8 |  -  |Isotropic hardening constant           | 
|  C9 |  -  |Kinematic hardening constant           | 
| C10 |  -  |Kinematic hardening constant           | 
| C11 |  -  |Multiaxial function const              | 
| C12 |  -  |Multiaxial function power              | 

In command line it will be looks like present bellow
```
!*** Define parameters related to generalised model  
!* Modulus of Elasticity  
Young	= 200000   
!* Poisson ratio  
nu	= 0.3  
!* Yield Strength  
S02	= 300  
!* Damage law constants  
r 	= 3.5   
s	= 1  
!* Isotropic hardening constants  
R0	= 0.001  
Rinf 	= 3300  
gamma	= 0.4  
!* Kinematic hardening constant   
a 	= 2500   
b	= 20   
!* Multiaxial function const   
Nlconst	= 1    
Nlpower	= 1 
    
!*** add user model   
TB,USER,1,1,12,   
TBTEMP,0   
TBDATA,,Young, nu, S02, r, s, R0  
TBDATA,,Rinf, gamma, a, b,  Nlconst, Nlpower
```
**5. Add 20 state variables**  

*Preprocessor->Material Props->Material models->Structural->Specialized Materials->User material options->State Variables*

| SVAR| Storing value                                 |
| --- | ------------------------------------- |
|  1     | equivalent plastic strain at end of time increment    |
| 2-7    | Plastic strain vector                                 |
|  8     | Damage parameter                                      |
|  9     | Von-mises stress                                      |
|  10    | Equivalent elastic stress                             |
|  11-17 | Backstress vector                                     |
|  18-20 | additional variables                                  |

Or use this APDL script in preprocessor section
```
TB,STAT,1,1,20,  
TBTEMP,0  
TBDATA,,0,0,0,0,0,0  
TBDATA,,0,0,0,0,0,0  
TBDATA,,0,0,0,0,0,0   
TBDATA,,0,0  
```
**6. Access to user variable arrays**

Before starting of the solution, in the solver section (/SOL) write the line:

- to save every substeps results  
> OUTRES,SVAR,ALL

- to save only the last step  
> OUTRES,SVAR,LAST

In order for all elements of user arrays to be available, command GRAF must be used in the postprocessor (/POST) section.  

> /GRA,FULL

That's all. Further we work as with the usual scheme.

<br>
<br>
<br>

#### COMPILATOR SETTINGS (*projectname*->properties). 

| Name     |   | Value |
| ----------- | ----------- |----------- |
|Supress startup banner:| 	 - |	            Yes (/nologo) | 
|Additional include Directories:|  - |	        C:\Program Files\ANSYS Inc\v***\ansys\customize\include | 
|Optimization:| 			 - |	            Disable (/Od) |
|Preprocessor definitions:| 	 - |	        /DNOSTDCALL /DARGTRAIL /DPCWIN64_SYS /DPCWINX64_SYS /DPCWINNT_SYS /DCADOE_ANSYS /D__EFL /DFORTRAN /auto /c /Fo.\ /MD /W0  |
|Debug information Format:|		 - |            Full (/debug:full)  |
|Preprocess Source file:|		 - |            Yes (/fpp)  |
|Preprocessor Definitions to fpp only:| -|	Yes (/noD)  |
|Use Portlib Library:| 		 - |	            Yes (/4Yportlib)  |

#### LINKER SETTINGS  
| Name    |  |    Value |
| ----------- | ----------- |----------- |
|Enable incremental linking:| - |		No (/INCREMENTAL:NO)  |
|Supress startup banner:|  - |		    Yes (/nologo) | 
|Additional library Directories:|  - |	C:\Program Files\ANSYS Inc\v***\ansys\custom\lib\winx64|  
|Additional dependencies:| 	 - |	    ANSYS.LIB  |
|Generate debug info: |	 - |		    Yes (/DEBUG)  |

*** - your version of ANSYS.  
All other settings by default. Its allows me connect to ANSYS for debugging.
